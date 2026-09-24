#!/usr/bin/env perl
# Coarse storage parity and explicitly retained-variant FASTQ regressions.
#
# Reuse the ONT fixture and its controlled clusterer/mappers. Preprocessing and
# native seed/count checks use the installed SDM. A separate capture wrapper
# exercises failed-seed handling. Standard retained-quality HQ support is required.
#
# Run: LOTUS_TEST_SDM=/path/to/sdm prove -v tests/coarse_derep.t
use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/lib";
use Digest::SHA ();
use File::Basename qw(basename);
use File::Temp qw(tempdir);
use JSON::PP ();
use List::Util qw(sum0 uniq);
use Test::More;
use DerepAudit qw(audit py_int);
use LotusTest qw($ROOT $SDM case contains lacks read_text write_text append_text text_lines json_decode rc run_command after count_of every4 has_arg replace_all);

my ($version) = run_command(\%ENV, 30, $SDM, '-v');
my $HAS_COARSE = $version =~ /sdm\s+(\d+)\.(\d+)/ && ($1 <=> 3 || $2 <=> 52) >= 0;
my ($flags) = run_command(\%ENV, 30, $SDM, '-help_flags');
my $HAS_SEEDS = index($flags, '-seedSubclusters') >= 0;
my $HAS_STANDARD_HQ = $HAS_SEEDS && index($flags, 'Standard HQ output: .1.hq.fq; paired additionally .2.hq.fq') >= 0;

# Shared by the SDM stand-ins: argument lookup, last-suffix removal, JSON with
# ", " separators, and the real SDM run (killed by signal N: exit 256-N).
my $SDM_PRELUDE = <<'PERL';
#!/usr/bin/env perl
use strict; use warnings; use JSON::PP ();
$SIG{__DIE__} = sub { die @_ if $^S; print STDERR @_; exit 1 };
my @args = @ARGV;
sub env_set { return length($ENV{$_[0]} // '') }
sub has { my ($flag) = @_; return scalar grep { $_ eq $flag } @args }
sub after { my ($flag) = @_; for my $i (0 .. $#args) { return $args[$i+1] // die "$flag: missing value\n" if $args[$i] eq $flag } return undef }
sub stem { my ($dir, $name) = $_[0] =~ m{\A(.*/)?([^/]*)\z}s; my $i = rindex($name, '.');
    $name = substr($name, 0, $i) if $i > 0 && $i < length($name) - 1; return ($dir // '') . $name }
sub spew { my ($f, $t) = @_; open my $fh, '>', $f or die "$f: $!\n"; print {$fh} $t; close $fh or die "$f: $!\n" }
sub text { my ($s) = @_; utf8::decode($s); return $s }
sub dumps { my ($v) = @_; my $j = JSON::PP->new->ascii->allow_nonref;
    return ref $v eq 'HASH' ? '{' . join(', ', map { $j->encode(text($_)) . ': ' . $j->encode(text($v->{$_})) } sort keys %$v) . '}'
        : '[' . join(', ', map { $j->encode(text($_)) } @$v) . ']' }
sub real_sdm { my $sdm = $ENV{COARSE_TEST_REAL_SDM} // die "COARSE_TEST_REAL_SDM unset\n";
    system { $sdm } $sdm, @args; die "Cannot run $sdm: $!\n" if $? == -1; return $? & 127 ? 256 - ($? & 127) : $? >> 8 }
PERL

my $WRAPPER = $SDM_PRELUDE . <<'PERL';
if (@args && ($args[0] eq '-v' || $args[0] eq '-version')) {
    print(($ENV{COARSE_TEST_VERSION} // 'sdm 3.52 beta'), "\n");
}
elsif (@args == 1 && $args[0] eq '-help_flags') {
    print env_set('COARSE_TEST_NO_CAPABILITY') ? "-derepStoreQuals <0|1>\n" : "-seedSubclusters <0|1>\n";
    print "Standard HQ output: .1.hq.fq; paired additionally .2.hq.fq\n-derepCoarseClusters <0|1>\n"
        unless env_set('COARSE_TEST_LEGACY_HQ');
}
elsif (has('-optimalRead2Cluster')) {
    spew($ENV{COARSE_TEST_SEED_CALL} // die("COARSE_TEST_SEED_CALL unset\n"), dumps(\@args));
    exit 17;  # Intentional failure: never claim to emulate native seed selection.
}
else {
    my $status = real_sdm();
    if ($status == 0 && env_set('COARSE_TEST_REMOVE_MATE') && defined(my $out = after('-o_dereplicate'))) {
        my $mate = stem($out) . '.2.hq.fq';
        unlink $mate or die "$mate: $!\n";
    }
    exit $status;
}
PERL

# Snapshot SDM's preprocessing outputs and record the seed command.
my $SNAPSHOT_WRAPPER = $SDM_PRELUDE . <<'PERL';
my $status = real_sdm();
if ($status == 0 && defined(my $out = after('-o_dereplicate'))) {
    my ($dir, $name) = stem($out) =~ m{\A(.*/)?([^/]*)\z}s;
    $dir //= '';
    opendir my $dh, $dir eq '' ? '.' : $dir or die "$dir: $!\n";
    my %files;
    for my $entry (grep { index($_, "$name.") == 0 && -f "$dir$_" } readdir $dh) {
        open my $fh, '<:raw', "$dir$entry" or die "$dir$entry: $!\n";
        local $/; (my $content = <$fh> // '') =~ s/\r\n?/\n/g;
        $files{$entry} = $content;
    }
    spew($ENV{COARSE_TEST_SNAPSHOT} // die("COARSE_TEST_SNAPSHOT unset\n"), dumps(\%files));
}
spew($ENV{COARSE_TEST_SEED_CALL} // die("COARSE_TEST_SEED_CALL unset\n"), dumps(\@args)) if has('-optimalRead2Cluster');
exit $status;
PERL

# DADA2 stand-in: checks the per-run FASTQ dereplicates and emits one ASV.
my $RSCRIPT = <<'PERL';
#!/usr/bin/env perl
use strict; use warnings;
$SIG{__DIE__} = sub { die @_ if $^S; print STDERR @_; exit 1 };
sub slurp { my ($f) = @_; open my $fh, '<', $f or die "$f: $!\n"; local $/; return <$fh> // '' }
sub spew { my ($f, $t) = @_; open my $fh, '>', $f or die "$f: $!\n"; print {$fh} $t; close $fh or die "$f: $!\n" }
my @args = @ARGV;
die "expected --vanilla: @args\n" unless ($args[0] // '') eq '--vanilla';
my $out = $args[3] // die "missing output directory\n";
opendir my $dh, $out or die "$out: $!\n";
my @files = map { "$out/$_" } grep { /\Aderep\..*\.fas\z/s } readdir $dh;
die "@files\n" unless @files == 2;
die "not FASTQ: $_\n" for grep { substr(slurp($_), 0, 1) ne '@' } @files;
spew("$out/dada2.uc", '');
spew("$out/dada2_p1_errF.pdf", 'test');
spew("$out/uniqueSeqs.fna", ">OTU0\n" . ($ENV{ONT_TEST_CONSENSUS} // die "ONT_TEST_CONSENSUS unset\n") . "\n");
PERL

my $JSON_PY = JSON::PP->new->canonical->allow_nonref;
# Replace the first occurrence; a missing anchor is an error, not a no-op.
sub replace_once { my ($text, $from, $to) = @_; my $at = index($text, $from); die "anchor not found: $from\n" if $at < 0;
    substr($text, $at, length $from) = $to; return $text }
sub sha256_file { return Digest::SHA->new(256)->addfile($_[0], 'b')->hexdigest }
sub size_of { my ($header) = @_; return $header =~ /;size=(\d+);/ ? 0 + $1 : die "no ;size= in '$header'\n" }
sub vkey { return join "\0", @_ }
sub cmds { return read_text("$_[0]{out}/LotuSLogS/LotuS_cmds.log") }
# Names in one directory matching a '*' pattern (dotfiles included).
sub glob_names {
    my ($dir, $pattern) = @_;
    my $re = join '', map { $_ eq '*' ? '.*' : quotemeta } split /(\*)/, $pattern;
    opendir my $dh, $dir or return [];
    return [sort grep { !/\A\.\.?\z/ && /\A$re\z/s } readdir $dh];
}
# Typed equality on decoded JSON: '200' != 200, while false == 0.
sub json_py { return JSON::PP->new->boolean_values(0, 1)->decode(read_text($_[0])) }
sub is_py {
    my ($got, $want, $name) = @_;
    local $Test::Builder::Level = $Test::Builder::Level + 1;
    return is($JSON_PY->encode($got), $JSON_PY->encode($want), $name);
}


sub capture_wrapper {
    my ($t) = @_;
    my $wrapper = "$t->{tools}/sdm";
    chmod 0755, write_text($wrapper, $WRAPPER) or die "chmod $wrapper: $!\n";
    write_text($t->{cfg}, replace_all(read_text($t->{cfg}), $SDM, $wrapper));
    $t->{seed_call} = "$t->{root}/seed_call.json";
    @{ $t->{env} }{qw(COARSE_TEST_REAL_SDM COARSE_TEST_SEED_CALL)} = ($SDM, $t->{seed_call});
}

sub options { my ($t, $text) = @_; append_text("$t->{root}/options.txt", $text) }

sub inputs {
    my ($t, $paired, %o) = @_;
    my $counts = $o{counts} // [4, 3];
    write_text($t->{map}, "#SampleID\tfastqFile\n" . join('', map { "$_\t$_.1.fq" . ($paired ? ",$_.2.fq" : '') . "\n" } qw(s1 s2)));
    for my $s (0, 1) {
        my $sample = (qw(s1 s2))[$s];
        for my $mate (1 .. ($paired ? 2 : 1)) {
            my $fq = '';
            for my $i (0 .. $counts->[$s] - 1) {
                my $seq = $t->{seq};
                $seq = (substr($seq, 0, 1) ne 'A' ? 'A' : 'T') . substr($seq, 1) if $sample eq 's1' && $i == 0 && $mate == ($paired ? 2 : 1);
                $fq .= "\@${sample}_$i/$mate\n$seq\n+\n" . ('I' x length $seq) . "\n";
            }
            write_text("$t->{reads}/$sample.$mate.fq", $fq);
        }
    }
    my $options = write_text("$t->{root}/options.txt", "minSeqLength\t100\nmaxSeqLength\t2000\nminAvgQuality\t0\n"
        . "RejectSeqWithoutFwdPrim\tF\nRejectSeqWithoutRevPrim\tF\nTrimWindowThreshhold\t0\nmaxHomonucleotide\t100\n");
    $t->options("derepStoreQuals\t1\n") if $o{retain} // 1;
    return ['-p', 'miSeq', '-CL', 'vsearch', '-s', $options, '-derepMin', '1', '-sdmThreads', '1'];
}

sub check_handoff {
    my ($t, $paired, $identity) = @_;
    $t->capture_wrapper;
    my $result = $t->run_lotus(extra => [@{ $t->inputs($paired) }, '-coarseDerep', $identity], ok => 0);
    contains($result->{output}, 'SDM retained-variant seed extension failed');
    lacks($result->{output}, 'Fallback to');
    my $args = json_decode(read_text($t->{seed_call}));
    is(after($args, '-seedSubclusters'), '1', '-seedSubclusters 1');
    if ($paired) { is(after($args, '-merge_pairs_seed'), '1', '-merge_pairs_seed 1') }
    else { ok(!has_arg($args, '-merge_pairs_seed'), 'no -merge_pairs_seed') }
    is(after($args, '-i_qual_offset'), '33', '-i_qual_offset 33');
    my @files = split /,/, after($args, '-i_fastq') // '', -1;
    is_deeply([map { basename($_) } @files], [$paired ? ('derep.1.hq.fq', 'derep.2.hq.fq') : 'derep.1.hq.fq'], 'seed FASTQs');
    my @records = map { [text_lines(read_text($_))] } @files;
    is(int(@{ $records[0] } / 4), 2, 'two exported variants');
    is(sum0(map { size_of($_) } @{ every4($records[0], 0) }), 7, 'variant sizes sum to 7');
    if ($paired) {
        is_deeply(every4($records[0], 0), every4($records[1], 0), 'mate headers agree');
        is(scalar(uniq @{ every4($records[0], 1) }), 1, 'one R1 sequence');
        is(scalar(uniq @{ every4($records[1], 1) }), 2, 'two R2 sequences');
    }
    my $commands = $t->cmds;
    contains($commands, '-derepIdentity ' . sprintf('%.12g', $identity * 100));
    contains($commands, $_) for '-derepStoreQuals 1', '-derepStoreDiffs 0', '-derepSubclusterFasta 0', '-derepReassign 0',
        '-derepCoarseClusters 0', '-merge_pairs_derep 0', '-merge_pairs_filter 0', '-merge_pairs_demulti 0';
    is(basename(after($args, '-derep_map') // ''), 'derep.map', 'derep map');
    is_deeply(glob_names("$t->{out}/tmpFiles", '*.subclusters*.fq'), [], 'no subcluster FASTQs');
    is_deeply(glob_names("$t->{out}/tmpFiles", '*.diff'), [], 'no diff files');
}

sub metadata {
    my ($t, %o) = @_;
    my $retained = $o{retained} // 1;
    my $granularity = $o{granularity} // 'exact search keys with configured prefix consolidation';
    my $m = json_py("$t->{out}/primary/sdm_dereplication.json");
    is_py($m->{status}, 'complete', 'status');
    is_py($m->{contract}, 'standard_hq_v2', 'contract');
    is_py($m->{quality_retention}, $retained, 'quality_retention');
    is_py($m->{output_granularity}, $granularity, 'output_granularity');
    is_py($m->{hq_record_layout}, $retained ? 'exact variants' : 'representatives', 'hq_record_layout');
    is_py($m->{sdm_binary_sha256}, sha256_file($SDM), 'sdm_binary_sha256');
    is_py($m->{sdm_options_sha256}, sha256_file("$t->{root}/options.txt"), 'sdm_options_sha256');
    is_py($m->{noSearchWithMerge}, { present => 0, value => undef }, 'noSearchWithMerge');
    is_py($m->{preprocessing_merge_flags}, { map { ("merge_pairs_$_" => 0) } qw(derep filter demulti) }, 'preprocessing_merge_flags');
    if (!$m->{derep_per_sequencing_run}) {
        my $counts = audit("$t->{out}/tmpFiles/derep.fas");
        is($counts->{map_total}, $counts->{passing} + $counts->{rest}, 'map total = passing + rest');
    }
    my $manifest = read_text("$t->{out}/LotuSLogS/run_manifest.txt");
    contains($manifest, "SDM quality retention: $retained");
    contains($manifest, "SDM dereplication output: $granularity");
    return $m;
}

sub variant_counts {
    my ($t, $paired) = @_;
    my $base = "$t->{out}/tmpFiles";
    my %parents;
    for my $line (text_lines(read_text("$base/derep.map"))) {
        next if $line =~ /\A#/;
        my ($header, @samples) = split /\t/, $line, -1;
        my $count = size_of($header // '');
        is($count, sum0(map { py_int((split /:/, $_, -1)[1]) } @samples), "map row $header");
        $parents{ (split /;size=/, $header, -1)[0] } = $count;
    }
    my @r1 = text_lines(read_text("$base/derep.1.hq.fq"));
    my @r2 = $paired ? text_lines(read_text("$base/derep.2.hq.fq")) : ('') x @r1;
    is_deeply(every4(\@r1, 0), every4(\@r2, 0), 'mate headers agree') if $paired;
    my (%sums, %records);
    for (my $i = 0; $i < @r1; $i += 4) {
        my $header = substr($r1[$i], 1);
        my $key = (split /;size=/, $header, -1)[0] // '';
        my $at = rindex($key, '.sub');
        die "no .sub suffix in '$header'\n" if $at < 0;
        ok(substr($key, $at + 4) =~ /\A[0-9]+\z/, "$header: numbered variant");
        my $count = size_of($header);
        $sums{ substr($key, 0, $at) } += $count;
        my $old = $records{ vkey($r1[$i+1], $r2[$i+1]) } // [0];
        $records{ vkey($r1[$i+1], $r2[$i+1]) } = [$old->[0] + $count, $r1[$i+3], $r2[$i+3]];
    }
    is_deeply(\%sums, \%parents, 'variant sizes sum to their parents');
    is_deeply(glob_names($base, '*.subclusters*.fq'), [], 'no subcluster FASTQs');
    return \%records;
}

# Make the helpers above callable as fixture methods.
{ no strict 'refs'; *{"LotusTest::$_"} = \&{"main::$_"} for qw(capture_wrapper options inputs check_handoff metadata variant_counts cmds); }

# CoarsePreflight

case test_invalid_ani_rejected_before_output_creation_even_without_strict => sub {
    my $t = shift;
    for my $value ('0', '0.9499', '1.0001', '95', '-1', 'NaN', 'Inf', 'abc', '0.97;echo', '1e999', '') {
        subtest "value='$value'" => sub {
            my $result = $t->run_lotus(extra => ["-coarseDerep=$value", '--no-strict'], ok => 0);
            contains($result->{output}, $value ne '' ? '-coarseDerep must be a numeric ANI' : 'coarseDerep requires an argument');
            ok(!-e $t->{out}, 'no output created');
        };
    }
};

case test_inclusive_ani_bounds_and_fraction_formats => sub {
    my $t = shift;
    for my $value (qw(0.95 .975 9.9e-1 1.0)) {
        subtest "value=$value" => sub {
            contains($t->run_lotus(extra => ['-coarseDerep', $value, '-v'])->{output}, 'LotuS 3.');
        };
    }
};

case test_old_sdm_rejected_before_processing => sub {
    my $t = shift;
    $t->capture_wrapper;
    $t->{env}{COARSE_TEST_VERSION} = 'sdm 3.51 beta';
    my $result = $t->run_lotus(extra => ['-CL', 'vsearch', '-coarseDerep', '0.97'], ok => 0);
    contains($result->{output}, 'requires SDM >= 3.52');
    ok(!-e $t->{seed_call}, 'no seed call');
    ok(!-e "$t->{out}/tmpFiles", 'no tmpFiles');
};

case test_retention_without_seed_capability_rejected_before_processing => sub {
    my $t = shift;
    $t->capture_wrapper;
    $t->{env}{COARSE_TEST_NO_CAPABILITY} = '1';
    my $options = write_text("$t->{root}/retention.txt", "derepStoreQuals\t1\n");
    my $result = $t->run_lotus(extra => ['-CL', 'vsearch', '-coarseDerep', '0.97', '-s', $options], ok => 0);
    contains($result->{output}, 'does not advertise -seedSubclusters');
    ok(!-e $t->{seed_call}, 'no seed call');
    ok(!-e "$t->{out}/tmpFiles", 'no tmpFiles');
};

case test_old_retained_quality_layout_rejected_despite_seed_flag => sub {
    my $t = shift;
    $t->capture_wrapper;
    $t->{env}{COARSE_TEST_LEGACY_HQ} = '1';
    my $result = $t->run_lotus(extra => ['-CL', 'vsearch', '-coarseDerep', '0.97'], ok => 0);
    contains($result->{output}, 'does not advertise the standard retained-quality HQ layout');
    ok(!-e $t->{seed_call}, 'no seed call');
    ok(!-e "$t->{out}/tmpFiles", 'no tmpFiles');
};

case test_incompatible_modes_rejected_even_without_strict => sub {
    my $t = shift;
    for my $c (['$sdmDerepDo = 0;', 'full run with SDM dereplication'],
               ['$mergePreCluster = 1; $sdmRetainQuals = 1;', '-mergePreClusterReads 0'],
               ['$saveDemulti = 1;', 'demultiplex-only'], ['$TaxOnly = "1";', 'taxonomy-only'], ['$onlyTaxRedo = 1;', 'taxonomy-only']) {
        my ($settings, $diagnostic) = @$c;
        subtest "settings=$settings" => sub {
            my $result = $t->probe("$settings validate_option_combinations();", extra => ['-coarseDerep', '0.97', '--no-strict'], ok => 0);
            contains($result->{output}, $diagnostic);
        };
    }
};

case test_finalized_abundance_and_pair_bound_are_reported => sub {
    my $t = shift;
    my $result = $t->probe(<<'PERL');
print parse_sdm_short_report(
            "Maximum pairs with both reads accepted: 900 (90.0%) n/a\n",
            "Dereplication: 10 unique sequences (avg size 78; 780 counts, 10 merged)\n"
            . "Dereplication abundance: 1200 total in map; 780 passing; 420 in rest\n"
            . "Retention exclusions: 40 pairs with empty R2 skipped before admission\n");
PERL
    contains($result->{output}, 'Pairs accepted on both ends (upper bound): 900 (90.0%)');
    contains($result->{output}, '10 passing unique sequences from 780 reads');
    contains($result->{output}, '1,200 total in map; 780 passing (main + merged); 420 in rest');
    contains($result->{output}, 'Retention exclusions: 40 pairs with physically empty R2');
};

subtest test_count_auditor_includes_merged_and_rest_and_rejects_missing_parents => sub {
    my $root = tempdir('derep-audit-test-XXXXXXXX', TMPDIR => 1, CLEANUP => 1);
    my $main = "$root/derep.fas";
    write_text("$root/derep.map", "#SMPLS\t0:s1\t1:s2\np;size=4;\t0:3\t1:1\nm;size=2;\t0:1\t1:1\nr;size=3;\t1:3\n");
    write_text("$root/derep.fas.rest", ">r;size=3;\nACGT\n");
    for my $fastq (0, 1) {
        subtest "fastq=$fastq" => sub {
            write_text($main, $fastq ? "\@p;size=4;\nACGT\n+\nIIII\n" : ">p;size=4;\nAC\nGT\n");
            my $merged = write_text("$root/derep.merg.fas", $fastq ? "\@m;size=2;\nAGGT\n+\nIIII\n" : ">m;size=2;\nAGGT\n");
            my $counts = audit($main);
            is_deeply([@$counts{qw(map_total passing rest)}], [9, 6, 3], 'map total, passing, rest');
            is_deeply($counts->{samples}, { s1 => 4, s2 => 5 }, 'per-sample map counts');
            unlink $merged or die "$merged: $!\n";
            ok(!eval { audit($main); 1 }, 'audit fails without the merged output');
            like($@, qr/map parents absent/, 'missing parents reported');
        };
    }
};

case test_ordinary_preprocessing_merge_policy_is_preserved => sub {
    my $t = shift;
    my $result = $t->probe(<<'PERL');
$mergePreCluster = 1; $ClusterPipe = 7;
            print JSON::PP->new->canonical->encode(sdm_preprocessing_merge_options());
PERL
    contains($result->{output}, '{"merge_pairs_demulti":1,"merge_pairs_derep":1}');
};

# CoarseHandoff: every case needs SDM with standard retained-quality HQ support.
sub handoff {
    my ($name, $code) = @_;
    return case($name, $code) if $HAS_COARSE && $HAS_STANDARD_HQ;
    subtest $name => sub { plan skip_all => 'requires SDM with standard retained-quality HQ support' };
}

handoff test_single_end_reconstructed_fastq_and_percent_conversion => sub { $_[0]->check_handoff(0, '0.975') };

handoff test_paired_reconstructed_fastqs_preserve_r2_only_variants => sub { $_[0]->check_handoff(1, '0.95') };

handoff test_explicit_identity_one_still_exports_subclusters => sub { $_[0]->check_handoff(1, '1.0') };

handoff test_missing_reconstructed_mate_fails_before_seed_command => sub {
    my $t = shift;
    $t->capture_wrapper;
    $t->{env}{COARSE_TEST_REMOVE_MATE} = '1';
    my $result = $t->run_lotus(extra => [@{ $t->inputs(1) }, '-coarseDerep', '0.97'], ok => 0);
    contains($result->{output}, 'Missing or empty SDM seed FASTQ');
    ok(!-e $t->{seed_call}, 'no seed call');
    ok(!-e "$t->{out}/primary/sdm_dereplication.json", 'no dereplication metadata');
};

handoff test_default_hq_path_still_counts_without_coarse_flags => sub {
    my $t = shift;
    $t->run_lotus(extra => $t->inputs(0, retain => 0));
    $t->check_counts({ s1 => 4, s2 => 3 });
    my $commands = $t->cmds;
    contains($commands, 'derep.1.hq.fq');
    lacks($commands, $_) for qw(-seedSubclusters -derepIdentity -derepStoreQuals);
    $t->metadata(retained => 0);
};

handoff test_options_file_retention_selects_variant_reader_at_100_percent => sub {
    my $t = shift;
    my $extra = $t->inputs(1);
    $t->options("derepStoreQuals\t0\nderepStoreQuals\t1\nderepIdentity\t100\n");
    $t->run_lotus(extra => $extra);
    $t->check_counts;
    $t->variant_counts(1);
    $t->metadata;
    my $commands = $t->cmds;
    contains($commands, '-seedSubclusters 1');
    lacks($commands, '-derepIdentity ');
};

handoff test_coarse_flag_overrides_optional_exports_and_coarse_parent_output => sub {
    my $t = shift;
    my $extra = $t->inputs(0);
    $t->options("derepStoreQuals\t1\nderepStoreDiffs\t1\nderepSubclusterFasta\t1\nderepReassign\t1\nderepCoarseClusters\t1\n");
    $t->run_lotus(extra => [@$extra, '-coarseDerep', '0.97']);
    $t->check_counts;
    $t->metadata;
    $t->variant_counts(0);
    my $base = "$t->{out}/tmpFiles";
    is(count_of(read_text("$base/derep.fas"), '>'), 2, 'two exact parents');
    is_deeply(glob_names($base, '*.diff'), [], 'no diff files');
    is_deeply(glob_names($base, '*.subclusters*.fna'), [], 'no subcluster FASTAs');
};

handoff test_below_cutoff_parent_is_recovered_once => sub {
    my $t = shift;
    $t->run_lotus(extra => [@{ $t->inputs(0) }, '-coarseDerep', '0.97', '-derepMin', '2']);
    $t->check_counts;
    $t->variant_counts(0);
    is(count_of(read_text("$t->{out}/tmpFiles/derep.fas"), '>'), 1, 'one passing parent');
    is(count_of(read_text("$t->{out}/tmpFiles/derep.fas.rest"), '>'), 1, 'one rest parent');
};

handoff test_exact_output_parity_at_97_and_100_with_multiple_workers => sub {
    my $t = shift;
    my @snapshots;
    for my $workers (1, 4, 12) {
        for my $identity ('0.97', '1.0') {
            subtest "workers=$workers identity=$identity" => sub {
                $t->{out} = "$t->{root}/parity_${workers}_$identity";
                $t->run_lotus(extra => [@{ $t->inputs(1, counts => [600, 400]) }, '-coarseDerep', $identity, '-sdmThreads', $workers]);
                $t->check_counts({ s1 => 600, s2 => 400 });
                push @snapshots, $t->variant_counts(1);
                my @parents = text_lines(read_text("$t->{out}/tmpFiles/derep.fas"));
                is_deeply([@parents[1 .. $#parents]], [$t->{seq}], 'one full-length parent');
                contains($parents[0], ';size=1000;');
            };
        }
    }
    is_deeply($_, $snapshots[0], 'variants match the first run') for @snapshots[1 .. $#snapshots];
};

handoff test_mixed_quality_admission_and_deferred_recovery => sub {
    my $t = shift;
    my $seq = $t->{seq};
    for my $identity ('0.97', '1.0') {
        for my $workers (1, 4) {
            subtest "identity=$identity workers=$workers" => sub {
                $t->{out} = "$t->{root}/mixed_${identity}_$workers";
                my $extra = $t->inputs(1);
                # One parent has only R1-passing observations. Another has
                # failed observations before its passing anchor in sample 2.
                my $variant = (substr($seq, 0, 1) ne 'T' ? 'T' : 'A') . substr($seq, 1);
                my %observations = (s1 => [([$variant, 'I', '+']) x 4, ([$seq, '+', '+']) x 2], s2 => [([$seq, 'I', 'I']) x 3]);
                for my $sample (sort keys %observations) {
                    my $rows = $observations{$sample};
                    for my $mate (1, 2) {
                        write_text("$t->{reads}/$sample.$mate.fq", join '', map {
                            my ($s, @quality) = @{ $rows->[$_] };
                            "\@${sample}_$_/$mate\n$s\n+\n" . ($quality[$mate - 1] x length $s) . "\n" } 0 .. $#$rows);
                    }
                }
                $t->options("minAvgQuality\t27\nderepSrchLen\t-1\nderepPrefix\t1\nfastqVersion\t1\n");
                my $result = $t->run_lotus(extra => [@$extra, '-coarseDerep', $identity, '-sdmThreads', $workers]);
                $t->check_counts({ s1 => 6, s2 => 3 });
                $t->variant_counts(1);
                $t->metadata;
                is(audit("$t->{out}/tmpFiles/derep.fas")->{map_total}, 9, 'nine reads in the map');
                contains(join(' ', split ' ', $result->{output}), '9 total in map; 9 passing (main + merged); 0 in rest');
            };
        }
    }
};

handoff test_later_better_variant_becomes_full_length_seed => sub {
    my $t = shift;
    my $extra = $t->inputs(0);
    my $seq = $t->{seq};
    my $better = substr($seq, 0, 300) . (substr($seq, 300, 1) ne 'A' ? 'A' : 'T') . substr($seq, 301);
    for my $c (['s1', 4, $seq, '5'], ['s2', 3, $better, 'I']) {
        my ($sample, $count, $sequence, $quality) = @$c;
        write_text("$t->{reads}/$sample.1.fq", join '', map { "\@${sample}_$_\n$sequence\n+\n" . ($quality x length $sequence) . "\n" } 0 .. $count - 1);
    }
    $t->options("TruncateSequenceLength\t200\nfastqVersion\t1\n");
    $t->{env}{ONT_TEST_CONSENSUS} = substr($seq, 0, 200);
    $t->run_lotus(extra => [@$extra, '-coarseDerep', '0.97']);
    my $state = $t->check_counts;
    my $records = $t->variant_counts(0);
    is($records->{ vkey($seq, '') }[1], '5' x length $seq, 'first variant keeps its qualities');
    is($records->{ vkey($better, '') }[1], 'I' x length $better, 'better variant keeps its qualities');
    is(length((text_lines(read_text("$t->{out}/tmpFiles/derep.fas")))[1] // ''), 200, 'truncated parent');
    my @hq = text_lines(read_text("$t->{out}/tmpFiles/derep.1.hq.fq"));
    is($hq[5], $better, 'second exported variant wins');
    my @seed = text_lines(read_text($state->{seed}));
    is(join('', @seed[1 .. $#seed]), $better, 'better variant becomes the full-length seed');
};

handoff test_metadata_records_effective_search_and_cut_settings => sub {
    my $t = shift;
    my $extra = $t->inputs(0);
    $t->options("derepSrchLen\t150\nderepSrchLen\t200\nderepPrefix\tauto\n"
        . "TruncateSequenceLength\t250\nkeepBarcodeSeq\t0\nTrimStartNTs\t5\n");
    $t->{env}{ONT_TEST_CONSENSUS} = substr($t->{seq}, 5, 200);
    $t->run_lotus(extra => [@$extra, '-coarseDerep', '0.97']);
    my $m = $t->metadata;
    is_py($m->{search_source}, 'R1', 'search_source');
    is_py($m->{search_options}{derepSrchLen}, '200', 'derepSrchLen');
    is_py($m->{search_options}{TruncateSequenceLength}, '250', 'TruncateSequenceLength');
    is_py($m->{sdm_options}{keepBarcodeSeq}, '0', 'keepBarcodeSeq');
    is_py($m->{sdm_options}{TrimStartNTs}, '5', 'TrimStartNTs');
    contains($m->{hq_sequence_policy}, 'logical tails retained');
};

handoff test_paired_hq_keeps_tails_and_removes_technical_sequences_once => sub {
    my $t = shift;
    my @primers = qw(AGGTCAGTACCGTAAC TCCAGATGCTACGTCA);
    my @barcodes = ([qw(GACTGA CTGTAC)], [qw(TACCTG AGTCGA)]);
    for my $retained (0, 1) {
        subtest "retained=$retained" => sub {
            $t->{out} = "$t->{root}/full_pair_$retained";
            my $extra = $t->inputs(1, retain => $retained);
            write_text($t->{map}, "#SampleID\tfastqFile\tBarcodeSequence\tBarcode2ndPair\tLinkerPrimerSequence\tReversePrimer\n"
                . join('', map { my $s = (qw(s1 s2))[$_]; "$s\t$s.1.fq,$s.2.fq\t" . join("\t", @{ $barcodes[$_] }, @primers) . "\n" } 0, 1));
            for my $c ([0, 's1', 4], [1, 's2', 3]) {
                my ($i, $sample, $count) = @$c;
                for my $mate (0, 1) {
                    my $sequence = $barcodes[$i][$mate] . $primers[$mate] . $t->{seq};
                    my $quality = ('I' x (length($sequence) - 20)) . ('-' x 20);
                    write_text("$t->{reads}/$sample." . ($mate + 1) . '.fq',
                        join '', map { "\@${sample}_$_/" . ($mate + 1) . "\n$sequence\n+\n$quality\n" } 0 .. $count - 1);
                }
            }
            $t->options("TruncateSequenceLength\t200\nTrimWindowThreshhold\t25\nTrimWindowWidth\t10\n"
                . "keepBarcodeSeq\t0\nkeepPrimerSeq\t0\nfastqVersion\t1\n");
            my $expected = $t->{seq};
            $t->{env}{ONT_TEST_CONSENSUS} = substr($expected, 0, 200);
            $t->run_lotus(extra => [@$extra, $retained ? ('-coarseDerep', '0.97') : ()]);
            $t->check_counts;
            for my $mate (1, 2) {
                my @hq = text_lines(read_text("$t->{out}/tmpFiles/derep.$mate.hq.fq"));
                is_deeply(every4(\@hq, 1), [$expected], "mate $mate: full sequence");
                is_deeply(every4(\@hq, 3), [('I' x (length($expected) - 20)) . ('-' x 20)], "mate $mate: tail qualities kept");
            }
            my @seeds = text_lines(read_text("$t->{out}/tmpFiles/otu_seeds.merg.fq"));
            is_deeply(every4(\@seeds, 1), [$expected], 'merged seed');
            contains(read_text("$t->{out}/LotuSLogS/SeedExtensionStats.log"), '1 were merged paired reads');
            $t->metadata(retained => $retained ? 1 : 0);
            $t->variant_counts(1) if $retained;
        };
    }
};

handoff test_dada2_sequencing_runs_keep_cumulative_variant_hq => sub {
    my $t = shift;
    my $extra = $t->inputs(1);
    my @lines = text_lines(read_text($t->{map}));
    write_text($t->{map}, "$lines[0]\tSequencingRun\n$lines[1]\trunA\n$lines[2]\trunB\n");
    chmod 0755, write_text("$t->{tools}/Rscript", $RSCRIPT) or die "chmod Rscript: $!\n";
    $t->{env}{PATH} = "$t->{tools}:" . ($t->{env}{PATH} // '');
    append_text($t->{cfg}, "dada2R $ROOT/bin/R/dada2_pip_v2.R\n");
    $t->run_lotus(extra => [@$extra, '-CL', 'dada2', '-coarseDerep', '0.97']);
    $t->check_counts;
    $t->variant_counts(1);
    is_py($t->metadata->{derep_per_sequencing_run}, 1, 'derep_per_sequencing_run');
    my $commands = $t->cmds;
    contains($commands, '-derep_format fq -derepPerSR 1');
    lacks($commands, '-derepPerSR 0');
    is_deeply(glob_names("$t->{out}/tmpFiles", 'derep.*.1.hq.fq'), [], 'no per-run HQ FASTQs');
};

handoff test_fresh_run_regenerates_older_coarse_outputs => sub {
    my $t = shift;
    $t->old_output;
    mkdir "$t->{out}/tmpFiles" or die "$t->{out}/tmpFiles: $!\n";
    write_text("$t->{out}/tmpFiles/derep.1.hq.fq", "\@old_parent;size=99;\nACGT\n+\nIIII\n");
    write_text("$t->{out}/tmpFiles/derep.subclusters.fq", "obsolete layout\n");
    mkdir "$t->{out}/primary" or die "$t->{out}/primary: $!\n";
    write_text("$t->{out}/primary/sdm_dereplication.json", qq({"contract":"old_coarse_parents"}\n));
    $t->run_lotus(extra => [@{ $t->inputs(0) }, '-coarseDerep', '0.97']);
    $t->check_counts;
    $t->variant_counts(0);
    $t->metadata;
    lacks(read_text("$t->{out}/tmpFiles/derep.1.hq.fq"), 'old_parent');
};

handoff test_native_seed_extension_preserves_sample_counts => sub {
    my $t = shift;
    for my $paired (0, 1) {
        subtest "paired=$paired" => sub {
            $t->{out} = "$t->{root}/native_$paired";
            $t->run_lotus(extra => [@{ $t->inputs($paired) }, '-coarseDerep', '0.97']);
            $t->check_counts({ s1 => 4, s2 => 3 });
            $t->variant_counts($paired);
            $t->metadata;
        };
    }
};

handoff test_storage_only_matches_ordinary_outputs_qualities_and_seed_command => sub {
    my $t = shift;
    # Capture native preprocessing outputs before LotuS can append merged
    # clustering inputs. Both preprocessing and seed selection use real SDM.
    my $wrapper = "$t->{tools}/sdm";
    chmod 0755, write_text($wrapper, $SNAPSHOT_WRAPPER) or die "chmod $wrapper: $!\n";
    write_text($t->{cfg}, replace_all(read_text($t->{cfg}), $SDM, $wrapper));
    my $clusterer = "$t->{tools}/vsearch";
    write_text($clusterer, replace_once(read_text($clusterer), q{unless basename($query) eq 'derep.fas';},
        q{unless grep { basename($query) eq $_ } qw(derep.fas derep.merg.fas);}));
    my $snapshot = "$t->{root}/preprocessing.json";
    my $seed_call = "$t->{root}/seed_call.json";
    @{ $t->{env} }{qw(COARSE_TEST_REAL_SDM COARSE_TEST_SNAPSHOT COARSE_TEST_SEED_CALL)} = ($SDM, $snapshot, $seed_call);
    # The existing checkpoint precedes FASTA assembly for paired seeds.
    my $checkpoint = 'atomic_write_text("$outdir/ont_test_state.json",';
    write_text($t->{script}, replace_once(read_text($t->{script}), $checkpoint, "mergeRds() if \$numInput == 2;\n$checkpoint"));
    my $fragment = substr($t->{seq}, 0, 420);
    for my $mode ([0, 0], [1, 0], [1, 1]) {
        my ($paired, $merge) = @$mode;
        my $baseline;
        for my $run ([undef, 1], ['0.97', 1], ['0.97', 4], ['1.0', 1]) {
            my ($identity, $workers) = @$run;
            my $label = "${paired}_${merge}_" . ($identity // 'none') . "_$workers";
            subtest "paired=$paired merge=$merge identity=" . ($identity // 'none') . " workers=$workers" => sub {
                $t->{out} = "$t->{root}/ordinary_parity_$label";
                my $extra = $t->inputs($paired, retain => 0);
                $t->{env}{ONT_TEST_CONSENSUS} = $fragment;
                for my $c (['s1', 4], ['s2', 3]) {
                    my ($sample, $count) = @$c;
                    for my $mate (1 .. ($paired ? 2 : 1)) {
                        my $fq = '';
                        for my $i (0 .. $count - 1) {
                            my $seq = $paired ? substr($fragment, 0, 260) : $fragment;
                            $seq = rc(substr($fragment, 160)) if $mate == 2;
                            # Different R2s must not expand HQ output into
                            # every full-pair variant. Unique best qualities
                            # also verify that selected mate qualities survive.
                            if ($i == 0) {
                                my $pos = $mate == 2 ? 20 : 40;
                                substr($seq, $pos, 1) = (grep { $_ ne substr($seq, $pos, 1) } qw(A C G T))[0];
                            }
                            my $quality = chr(33 + 20 + $i * 3 + ($sample eq 's2' ? 1 : 0));
                            $fq .= "\@${sample}_$i/$mate\n$seq\n+\n" . ($quality x length $seq) . "\n";
                        }
                        write_text("$t->{reads}/$sample.$mate.fq", $fq);
                    }
                }
                # Two exact keys, one below cutoff; native count recovery
                # must still conserve all seven observations.
                push @$extra, '-derepMin', ($merge ? '1' : '3'), '-mergePreClusterReads', $merge, '-sdmThreads', $workers;
                push @$extra, '-coarseDerep', $identity if defined $identity;
                $t->run_lotus(extra => $extra);
                my $state = $t->check_counts;
                my $files = json_decode(read_text($snapshot));
                is_deeply([sort keys %$files], [sort 'derep.fas', 'derep.fas.rest', 'derep.map', 'derep.1.hq.fq',
                    ($paired ? 'derep.2.hq.fq' : ()), ($merge ? 'derep.merg.fas' : ())], 'preprocessing outputs');
                lacks($files->{'derep.1.hq.fq'}, '.sub1');
                my @hq = text_lines($files->{'derep.1.hq.fq'} // '');
                is(int(@hq / 4), 2, 'two HQ representatives');
                my ($seqs, $quals) = (every4(\@hq, 1), every4(\@hq, 3));
                my $pairs = @$seqs < @$quals ? @$seqs : @$quals;
                ok(!grep({ length($seqs->[$_]) != length($quals->[$_]) || !length($quals->[$_]) } 0 .. $pairs - 1),
                    'HQ qualities are non-empty and match their sequences');
                my $metadata = json_py("$t->{out}/primary/sdm_dereplication.json");
                is_py($metadata->{quality_retention}, 0, 'quality_retention');
                is_py($metadata->{hq_record_layout}, 'representatives', 'hq_record_layout');
                my $args = [map { replace_all($_, $t->{out}, '<OUTPUT>') } @{ json_decode(read_text($seed_call)) }];
                ok(!has_arg($args, '-seedSubclusters'), 'no -seedSubclusters');
                my $outputs = { primary => $files, seed_args => $args, seeds => read_text($state->{seed}),
                    matrix => read_text("$t->{out}/OTU.txt"), seed_stats => read_text("$t->{out}/LotuSLogS/SeedExtensionStats.log") };
                ok(length $outputs->{seeds}, 'seeds written');
                if (!defined $baseline) { $baseline = $outputs }
                else { is_deeply($outputs, $baseline, 'outputs match the ordinary run') }
            };
        }
    }
};

handoff test_storage_only_preserves_main_fastq_quality_averaging => sub {
    my $t = shift;
    my $fragment = substr($t->{seq}, 0, 420);
    for my $mode ([0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 0, 27], [1, 0, 27], [1, 1, 27]) {
        my ($paired, $merge, $min_quality) = @$mode;
        $t->inputs($paired, retain => 0);
        $t->options("minAvgQuality\t$min_quality\n");
        for my $c (['s1', 4], ['s2', 3]) {
            my ($sample, $count) = @$c;
            for my $mate (1 .. ($paired ? 2 : 1)) {
                my $seq = $mate == 1 ? $fragment : rc($fragment);
                write_text("$t->{reads}/$sample.$mate.fq", join '', map {
                    "\@${sample}_$_/$mate\n$seq\n+\n" . (chr(33 + 20 + $_ * 3 + ($sample eq 's2' ? 1 : 0)) x length $seq) . "\n" } 0 .. $count - 1);
            }
        }
        my $baseline;
        for my $run ([undef, 1], ['97', 1], ['97', 4], ['100', 1]) {
            my ($identity, $workers) = @$run;
            my $label = "${paired}_${merge}_" . ($identity // 'none') . "_${workers}_$min_quality";
            subtest "paired=$paired merge=$merge identity=" . ($identity // 'none') . " workers=$workers min_quality=$min_quality" => sub {
                my $base = "$t->{root}/fq_parity_$label";
                mkdir $base or die "$base: $!\n";
                my @args = ($SDM, '-i_path', $t->{reads}, '-map', $t->{map}, '-options', "$t->{root}/options.txt",
                    '-paired', ($paired ? '2' : '1'), '-i_qual_offset', '33', '-o_qual_offset', '33', '-threads', $workers,
                    '-o_fastq', "$base/filtered.fq", '-o_dereplicate', "$base/derep.fas", '-derep_format', 'fq',
                    '-min_derep_copies', '1', '-suppressOutput', '3', '-merge_pairs_derep', $merge);
                push @args, '-derepIdentity', $identity if defined $identity;
                my ($output, $status) = run_command(\%ENV, 30, @args);
                is($status, 0, 'SDM succeeds') or diag($output);
                my %files = map { ($_ => read_text("$base/$_")) } @{ glob_names($base, 'derep.*') };
                ok(length($files{'derep.1.hq.fq'} // ''), 'HQ FASTQ written');
                my @main = text_lines($files{ $merge ? 'derep.merg.fas' : 'derep.fas' } // '');
                is(scalar @main, 4, 'one main FASTQ record');
                is(length($main[1] // ''), length($main[3] // ''), 'main quality length');
                # Selection can depend on paired/merged ranking. Preserve
                # the ordinary winner's full quality vector, not a presumed
                # maximum-quality winner; main FASTQ also collects evidence.
                my @hq = text_lines($files{'derep.1.hq.fq'} // '');
                is(length($hq[3] // ''), length $fragment, 'full HQ quality vector');
                ok(!grep({ my $q = ord($_) - 33; $q < 20 || $q > 29 } split //, $hq[3] // ''), 'HQ qualities within the input range');
                if (!defined $baseline) { $baseline = \%files }
                else { is_deeply(\%files, $baseline, 'outputs match the ordinary run') }
            };
        }
    }
};

done_testing();

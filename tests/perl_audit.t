#!/usr/bin/env perl
# Regression tests for the Perl audit; no database downloads or scientific claims.
#
# Run: prove -v tests/perl_audit.t
# Reuse the ONT fixture, with real SDM. Small probes invoke existing Perl helpers
# inside a temporary script copy; complete-flow tests retain the whole pipeline.
use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/lib";
use Fcntl qw(:flock);
use List::Util qw(sum0);
use Test::More;
use LotusTest qw($ROOT $TOOL_MAPPER_BRANCH case contains lacks read_text write_text append_text text_lines
    json_decode json_encode count_of has_arg replace_first);

# Compare via canonical JSON so the number 1 and the string "1" stay distinct.
sub same_json { local $Test::Builder::Level = $Test::Builder::Level + 1; is(json_encode($_[0]), json_encode($_[1]), $_[2]) }
sub result { my (undef, $json) = split /RESULT:/, $_[0]{output}, 2; return json_decode($json // '') }

# Stand-in USEARCH: prints a configurable version banner and, like usearch12
# (label.cpp GetSizeFromLabel, uchime3denovo.cpp), rejects de novo chimera input
# that lacks ;size= labels or is not sorted by decreasing size.
my $USEARCH = <<'PERL';
#!/usr/bin/env perl
use strict; use warnings; use JSON::PP ();
my @args = @ARGV;
my $banner = $ENV{TEST_USEARCH_BANNER} // 'usearch v11.0.667_i86linux32';
sub quit { print STDERR "$_[0]\n"; exit 1 }
if (@args == 1 && $args[0] eq '--version') {
    exit 1 if ($ENV{TEST_USEARCH_MODE} // '') eq 'banner';
    print "$banner\n"; exit 0;
}
if (!@args) { print "$banner, 4.0Gb RAM, 8 cores\n(C) Copyright 2013-18 Robert C. Edgar.\n"; exit 0 }
open my $log, '>>', $ENV{TEST_USEARCH_CALLS} or die "$ENV{TEST_USEARCH_CALLS}: $!\n";
print {$log} JSON::PP->new->encode(\@args), "\n"; close $log;
sub arg { my ($flag) = @_; for my $i (0 .. $#args - 1) { return $args[$i+1] if $args[$i] eq $flag } die "missing $flag\n" }
sub records { open my $fh, '<', $_[0] or die "$_[0]: $!\n"; local $/; my (undef, @p) = split />/, <$fh> // '', -1;
    return map { my @l = split /\n/; [$l[0], join '', @l[1 .. $#l]] } @p }
sub size { return $_[0] =~ /;size=(\d+)/ ? $1 : quit("Missing size= in >$_[0]") }
sub spew { my ($path, @recs) = @_; open my $out, '>', $path or die "$path: $!\n"; print {$out} map { ">$_->[0]\n$_->[1]\n" } @recs; close $out }
if ($args[0] eq '-sortbysize') {
    my @recs = records($args[1]); my @s = map { size($_->[0]) } @recs;
    spew(arg('-fastaout'), @recs[sort { $s[$b] <=> $s[$a] || $a <=> $b } 0 .. $#recs]);
}
elsif (grep { $args[0] eq $_ } qw(-uchime3_denovo -uchime2_denovo -uchime_denovo)) {
    my @recs = records($args[1]); my @s = map { size($_->[0]) } @recs;
    quit('Not sorted by size') if "@s" ne join ' ', sort { $b <=> $a } @s;
    spew(arg('-nonchimeras'), @recs); spew(arg('-chimeras')); spew(arg('-log'));
}
else { quit("Unexpected usearch call: @args") }
PERL


sub taxonomy {
    my ($t, $rows, %o) = @_;
    $t->{env}{AUDIT_TAX} = write_text("$t->{root}/hierarchy.tsv", "header\n$rows");
    return sprintf 'my @tax = readTaxIn($ENV{AUDIT_TAX}, %d, %d, %d); print "RESULT:", JSON::PP->new->encode(\@tax);',
        map { $o{$_} // 1 } qw(lca biom hit);
}

sub mapping_tool {
    my ($t, $name, $contents) = @_;
    $t->{env}{AUDIT_ALIGNMENTS} = write_text("$t->{root}/alignments", $contents);
    write_text("$t->{tools}/$name", <<'PERL');
#!/usr/bin/env perl
use File::Basename qw(basename); use File::Copy qw(copy);
if (grep { $_ eq '--version' } @ARGV) { print +(basename($0) eq 'minimap2' ? "2.28\n" : "vsearch v2.29.0\n"); exit 0 }
my $flag = (grep { $_ eq '-o' } @ARGV) ? '-o' : '-uc';
my ($i) = grep { $ARGV[$_] eq $flag } 0 .. $#ARGV;
die "missing $flag\n" unless defined $i;
copy($ENV{AUDIT_ALIGNMENTS}, $ARGV[$i + 1]) or die "$ENV{AUDIT_ALIGNMENTS}: $!\n";
PERL
}

sub contamination {
    return sprintf <<'PERL', $_[0] // 1;
ensure_dir($logDir); $mini2Bin = $ENV{AUDIT_MINIMAP}; $VSBin = $ENV{AUDIT_VSEARCH};
$useMini4map = %d; $doPhiX = 1; $uthreads = 1;
my $hits = contamination_rem($input, $refDBwanted, "phiX", 0);
print "RESULT:", JSON::PP->new->encode($hits);
PERL
}

sub setup_contamination {
    my ($t) = @_;
    @{ $t->{env} }{qw(AUDIT_MINIMAP AUDIT_VSEARCH)} = ("$t->{tools}/minimap2", "$t->{tools}/vsearch");
    return ['-i', $t->{ref}];
}

sub clean_table {
    my ($t, $rows) = @_;
    $t->{env}{AUDIT_FASTA} = write_text("$t->{root}/otus.fna", ">good\nACGT\n>zero\nACTG\n");
    my $table = $t->{env}{AUDIT_TABLE} = write_text("$t->{root}/table.tsv", "OTU\tsample\n$rows");
    return ('ensure_dir($logDir); $extendedLogs = 0; clean_otu_mat($ENV{AUDIT_FASTA}, $ENV{AUDIT_TABLE}, {});', $table);
}

sub fake_lambda {
    my ($t) = @_;
    my $binary = write_text("$t->{tools}/lambda3", <<'PERL');
#!/usr/bin/env perl
if (grep { $_ eq '--version' } @ARGV) { print "lambda3 version: 3.0.0\n" } else { exit 17 }
PERL
    chmod 0755, $binary or die "chmod $binary: $!\n";
    $t->{env}{AUDIT_LAMBDA} = $binary;
}

my $LAMBDA_BODY = <<'PERL';
ensure_dir($logDir); $lotus_tempDir = "$outdir/tmp"; ensure_dir($lotus_tempDir);
$doBlasting = 2; $lambda3Bin = $ENV{AUDIT_LAMBDA}; $BlastCores = 1;
doDBblasting($refDBwanted, $refDBwanted, "$lotus_tempDir/tax.out");
PERL

sub fake_usearch {
    my ($t) = @_;
    my $usearch = write_text("$t->{tools}/usearch", $USEARCH);
    chmod 0755, $usearch or die "chmod $usearch: $!\n";
    $t->{usearch_calls} = "$t->{root}/usearch_calls.jsonl";
    @{ $t->{env} }{qw(TEST_USEARCH TEST_USEARCH_CALLS)} = ($usearch, $t->{usearch_calls});
    return $usearch;
}

sub complete_pipeline {
    my ($t, %o) = @_;
    my $barbell = $o{barbell} ? 1 : 0;
    write_text($t->{script}, read_text("$ROOT/lotus3"));
    # Keep real LCA as well as real SDM; the aligner writes its documented
    # eleven-column input, so LCA/BIOM/taxonomy aggregation run unmodified.
    unlink "$t->{tools}/LCA" or die "unlink LCA: $!\n";
    symlink("$ROOT/bin/LCA", "$t->{tools}/LCA") or die "symlink LCA: $!\n";
    my $vsearch = read_text("$t->{tools}/vsearch");
    my $at = index($vsearch, $TOOL_MAPPER_BRANCH);
    die "stand-in mapper branch missing\n" if $at < 0;
    substr($vsearch, $at, 0) = <<'PERL';
elsif ($name eq 'vsearch' && has('--makeudb_usearch')) {
    spew(arg('-output'), 'synthetic test index');
}
elsif ($name eq 'vsearch' && has('-userout')) {
    spew(arg('-userout'), join '', map { my ($rid, $n) = ($_->[0], length $_->[1]);
        "$rid\tref\t99.9\t$n\t1\t0\t1\t$n\t1\t$n\t$n\n" } fasta(slurp(arg('--usearch_global'))));
}
PERL
    write_text("$t->{tools}/vsearch", $vsearch);
    $t->run_lotus(extra => $barbell ? ['-ontMinReads', '2'] : [], barbell => $barbell);
    my $citations = read_text("$t->{out}/LotuSLogS/citations.txt");
    is(count_of($citations, '10.64898/2026.05.26.727271'), 1, 'Savont cited once');
    is(count_of($citations, '10.1093/bioinformatics/btag349'), $barbell, 'Barbell citation');
    is_deeply([text_lines(read_text("$t->{out}/OTU.txt"))], ["OTU\ts1\ts2", "ASV1\t4\t3"], 'OTU table');
    my $biom = json_decode(read_text("$t->{out}/OTU.biom"));
    same_json($biom->{shape}, [1, 2], 'BIOM shape');
    same_json($biom->{data}, [[4, 3]], 'BIOM data');
    is($biom->{rows}[0]{id}, 'ASV1', 'BIOM row id');
    is($biom->{rows}[0]{metadata}{taxonomy}[0], 'k__Bacteria', 'BIOM taxonomy');
    contains(read_text("$t->{out}/OTU.fna") =~ s/\n//gr, $t->{consensus}, 'consensus in OTU.fna');
    contains(read_text("$t->{out}/higherLvl/Phylum.txt"), 'Bacteria');
    ok(-e "$t->{out}/LotuSLogS/run_manifest.txt", 'run manifest written');
    ok(!-e "$t->{out}/tmpFiles", 'temporary files removed');
    # SDM counting receives the documented -count_chimeras default.
    contains(read_text("$t->{out}/LotuSLogS/LotuS_cmds.log"), '-count_chimeras F');
}

# Make the helpers above callable as fixture methods.
{ no strict 'refs'; *{"LotusTest::$_"} = \&{"main::$_"} for qw(taxonomy mapping_tool setup_contamination clean_table fake_lambda fake_usearch complete_pipeline); }

case test_lock_failure_preserves_active_logs => sub {
    my $t = shift;
    my $log = $t->old_output;
    open my $lock, '>>', "$t->{root}/.lotus3.output.lock" or die "lock: $!\n";
    flock($lock, LOCK_EX | LOCK_NB) or die "flock: $!\n";
    my $result = $t->run_lotus(ok => 0);
    close $lock;
    contains($result->{output}, 'in use by another');
    is(read_text($log), "active run output\n", 'active log preserved');
};

case test_lock_file_keeps_same_inode_between_runs => sub {
    my $t = shift;
    my $body = 'acquire_output_lock($outdir); release_output_lock();';
    $t->probe($body);
    my $lock = "$t->{root}/.lotus3.output.lock";
    ok(-e $lock, 'lock file exists');
    my $inode = (stat $lock)[1];
    $t->probe($body);
    is((stat $lock)[1], $inode, 'lock file keeps its inode');
};

case test_output_reset_preserves_nested_input_and_config => sub {
    my $t = shift;
    for my $c (['-m', read_text($t->{map})], ['-c', read_text($t->{cfg})], ['-i', read_text($t->{raw})], ['-tax4refDB', read_text($t->{tax})]) {
        my ($flag, $data) = @$c;
        subtest "flag=$flag" => sub {
            my $log = $t->old_output;
            my $asset = write_text("$t->{out}/protected.txt", $data);
            for my $dry ([], ['--dry-run']) {
                my $result = $t->run_lotus(extra => [$flag, $asset, @$dry], ok => 0);
                contains($result->{output}, 'contains input or configuration');
                is(read_text($asset), $data, 'protected asset unchanged');
                is(read_text($log), "active run output\n", 'active log preserved');
            }
        };
    }
};

case test_output_reset_preserves_configured_reference => sub {
    my $t = shift;
    $t->old_output;
    my $asset = write_text("$t->{out}/reference.fna", ">ref\nACGT\n");
    append_text($t->{cfg}, "TAX_REFDB_KSGP $asset\n");
    my $result = $t->run_lotus(ok => 0);
    contains($result->{output}, 'contains input or configuration');
    ok(-e $asset, 'configured reference kept');
};

case test_failed_preflight_preserves_previous_run => sub {
    my $t = shift;
    my $log = $t->old_output;
    write_text("$t->{out}/OTU.txt", "previous results\n");
    write_text($t->{map}, "this is not a valid mapping file\n");
    my $result = $t->run_lotus(ok => 0);
    contains($result->{output}, 'previous output preserved');
    is(read_text($log), "active run output\n", 'active log preserved');
    is(read_text("$t->{out}/OTU.txt"), "previous results\n", 'previous results kept');
};

case test_successful_preflight_allows_output_reset => sub {
    my $t = shift;
    $t->old_output;
    write_text("$t->{out}/previous.txt", "old output\n");
    $t->run_lotus;
    ok(!-e "$t->{out}/previous.txt", 'old output removed');
    contains(read_text("$t->{out}/OTU.txt"), 'OTU_0');
};

case test_temp_cleanup_preserves_input => sub {
    my $t = shift;
    my $scratch = "$t->{root}/scratch";
    mkdir $scratch or die "$scratch: $!\n";
    my $source = write_text("$scratch/s1.fq", read_text($t->{raw}));
    write_text("$scratch/s2.fq", read_text($t->{raw}));
    write_text("$scratch/.lotus3_tmp_owned", "Output: $t->{out}\n");
    my $result = $t->run_lotus(extra => ['-i', $scratch, '-T', $scratch], ok => 0);
    contains($result->{output}, 'contains input or configuration');
    is(read_text($source), read_text($t->{raw}), 'input unchanged');
};

case test_unknown_clusterer_and_unsupported_modes_fail_before_reset => sub {
    my $t = shift;
    my $log = $t->old_output;
    for my $c ([['-CL', 'vsarch'], 'Unknown -CL'], [['-p', 'miSeq', '-CL', 'vsearch', '-highmem', '0'], 'SDM dereplication'],
               [['-exe', '2'], '-exe must'], [['-useMini4map', '2'], '-useMini4map must'], [['-saveDemultiplex', '3'], '-saveDemultiplex must']) {
        my ($args, $diagnostic) = @$c;
        subtest "args=@$args" => sub {
            contains($t->run_lotus(extra => $args, ok => 0)->{output}, $diagnostic);
            is(read_text($log), "active run output\n", 'active log preserved');
        };
    }
};

case test_shell_active_path_characters_rejected => sub {
    my $t = shift;
    for my $name ("quote'path", 'quote"path', 'glob[1]', 'glob*', 'glob?', 'back\\slash', 'paren(path)', 'semi;colon') {
        subtest "name=$name" => sub {
            my $result = $t->run_lotus(extra => ['-o', "$t->{root}/$name"], ok => 0);
            contains($result->{output}, 'Unsafe output path');
            ok(!-e "$t->{root}/$name", 'unsafe output not created');
        };
    }
};

case test_executable_crash_does_not_pass_version_check => sub {
    my $t = shift;
    write_text("$t->{tools}/vsearch", <<'PERL');
#!/usr/bin/env perl
$| = 1;
print "vsearch v2.29.0\n";
kill 'TERM', $$;
PERL
    my $result = $t->run_lotus(extra => ['--dry-run'], ok => 0);
    contains($result->{output}, 'Executable check failed (exit 143)');
    is_deeply($t->tool_calls, [], 'no tools called');
};

case test_old_minimap_version_rejected => sub {
    my $t = shift;
    write_text("$t->{tools}/minimap2", "#!/usr/bin/env perl\nprint \"2.9-r123\\n\";\n");
    my $result = $t->run_lotus(extra => ['--dry-run'], ok => 0);
    contains($result->{output}, 'too low, expected at least 2.17');
};

case test_semantic_version_comparison => sub {
    my $t = shift;
    $t->probe('die "wrong version order" if version_at_least("0.9", "0.25") || version_at_least("3.9", "3.43") || !version_at_least("2.28.1", "2.28");');
};

case test_duplicate_ids_rejected_by_default => sub {
    my $t = shift;
    my ($fasta, $taxonomy) = ("$t->{root}/validation.fasta", "$t->{root}/validation.tax");
    @{ $t->{env} }{qw(AUDIT_DUP_FASTA AUDIT_DUP_TAX)} = ($fasta, $taxonomy);
    my $body = <<'PERL';
my ($fn,$fp,$fw,$fi) = fasta_validation_scan($ENV{AUDIT_DUP_FASTA}, 10);
my ($tn,$tp,$tw) = taxonomy_validation_scan($ENV{AUDIT_DUP_TAX}, 10);
print "RESULT:", JSON::PP->new->encode({fasta => $fp, taxonomy => $tp});
PERL
    for my $duplicate (0, 1) {
        subtest 'duplicate=' . ($duplicate ? 'True' : 'False') => sub {
            my @ids = ('asv1', $duplicate ? 'asv1' : 'asv2');
            write_text($fasta, join '', map { ">$_\nACGT\n" } @ids);
            write_text($taxonomy, join '', map { "$_\tk__Bacteria\n" } @ids);
            my $errors = result($t->probe($body));
            if ($duplicate) {
                is(scalar @{ $errors->{fasta} }, 1, 'one FASTA error');
                contains($errors->{fasta}[0], 'Duplicate FASTA ID');
                is(scalar @{ $errors->{taxonomy} }, 1, 'one taxonomy error');
                contains($errors->{taxonomy}[0], 'Duplicate taxonomy ID');
            }
            else {
                is_deeply($errors, { fasta => [], taxonomy => [] }, 'no validation errors');
            }
        };
    }
};

case test_taxonomy_prefixes_missing_ranks_and_reference_ids => sub {
    my $t = shift;
    my $body = $t->taxonomy("ASV1\tk__Bacteria\tp__P\tc__C\to__O\tf__F\tg__G\ts__S\tref1\nASV2\tBacteria\tP\tC\tO\tF\tG\t\n");
    my $r = result($t->probe($body));
    is(scalar @$r, 4, 'readTaxIn returns four values');
    my ($tax, $levels, $hits, $hit_tax) = @$r;
    like($tax->{ASV1}, qr/\A\Qk__Bacteria", "p__P\E/, 'ASV1 keeps rank prefixes');
    lacks($tax->{ASV1}, 'k__k__');
    like($tax->{ASV2}, qr/\Qs__?\E\z/, 'ASV2 gets missing species rank');
    is_deeply($hits, { ref1 => ['ASV1'], ASV2 => ['ASV2'] }, 'reference hits');
    is($hit_tax->{ref1}, $tax->{ASV1}, 'reference taxonomy');
};

case test_rdp_taxonomy_keeps_otu_id => sub {
    my $t = shift;
    my $r = result($t->probe($t->taxonomy("Bacteria\tP\tC\tO\tF\tG\t\tOTU1\n", lca => 0)));
    is(scalar @$r, 4, 'readTaxIn returns four values');
    my ($tax, $levels, $hits, $hit_tax) = @$r;
    is_deeply([keys %$tax], ['OTU1'], 'OTU id kept');
    like($tax->{OTU1}, qr/\Qs__?\E\z/, 'missing species rank');
};

case test_malformed_and_duplicate_taxonomy_rejected => sub {
    my $t = shift;
    for my $c (["ASV1\tBacteria\tP\tC\tO\tF\n", 'Malformed taxonomy'], ["ASV1\tBacteria\tP\tC\tO\tF\tG\tS\n" x 2, 'Duplicate taxonomy']) {
        my ($rows, $diagnostic) = @$c;
        subtest "diagnostic=$diagnostic" => sub {
            contains($t->probe($t->taxonomy($rows), ok => 0)->{output}, $diagnostic);
        };
    }
};

case test_contamination_requires_query_coverage_and_identity => sub {
    my $t = shift;
    $t->mapping_tool('minimap2', join '', map { my ($rid, $span, $matches) = @$_; "$rid\t1000\t0\t$span\t+\tgenome\t5000000\t0\t$span\t$matches\t$span\t60\n" }
        ['short', 100, 100], ['lowid', 900, 700], ['valid', 900, 850], ['boundary', 500, 450], ['valid', 900, 850]);
    my $hits = result($t->probe(contamination(), extra => $t->setup_contamination));
    same_json($hits, { 'phiX.0' => { valid => 1, boundary => 1 } }, 'contaminant hits');
};

case test_malformed_paf_fails_with_context => sub {
    my $t = shift;
    for my $c (["bad\trow\n", 'Malformed PAF'], ["bad\t100\t0\t0\t+\tr\t100\t0\t0\t0\t0\t60\n", 'Invalid PAF']) {
        my ($line, $diagnostic) = @$c;
        $t->mapping_tool('minimap2', $line);
        contains($t->probe(contamination(), extra => $t->setup_contamination, ok => 0)->{output}, $diagnostic);
    }
};

case test_contamination_respects_vsearch_mapper => sub {
    my $t = shift;
    $t->mapping_tool('vsearch', "H\t0\t1000\t99\t+\t0\t0\t1000M\tvalid\tref\n");
    my $extra = $t->setup_contamination;
    $t->{env}{AUDIT_MINIMAP} = '/unavailable/minimap2';
    same_json(result($t->probe(contamination(0), extra => $extra)), { 'phiX.0' => { valid => 1 } }, 'contaminant hits');
};

case test_phix_hits_reach_matrix_filtering => sub {
    my $t = shift;
    my $marker = '# ////////////////////////// TAXONOMY';
    write_text($t->{script}, replace_first(read_text("$ROOT/lotus3"), $marker, "release_output_lock(); exit(0);\n$marker"));
    # The mapper uses each query header for both backmapping and PhiX search.
    append_text($t->{cfg}, "REFDB_PHIX $t->{ref}\n");
    $t->run_lotus(extra => ['-removePhiX', '1', '-keepOfftargets', '1']);
    my @rows = text_lines(read_text("$t->{out}/OTU.txt"));
    my ($id, @counts) = split /\t/, $rows[1];
    contains($id, '.phiX.');
    is(sum0(@counts), 7, 'PhiX row keeps all reads');
    contains(read_text("$t->{out}/LotuSLogS/OTU.contaminants.fa"), 'phiX');
};

case test_one_nonzero_and_one_zero_otu_is_not_empty => sub {
    my $t = shift;
    my ($body, $table) = $t->clean_table("good\t7\nzero\t0\n");
    $t->probe($body);
    is_deeply([text_lines(read_text($table))], ["OTU\tsample", "OTU1\t7"], 'zero OTU dropped');
};

case test_all_zero_otu_matrix_rejected => sub {
    my $t = shift;
    my ($body) = $t->clean_table("zero\t0\n");
    contains($t->probe($body, ok => 0)->{output}, 'Empty OTU matrix');
};

case test_lambda_index_failure_preserves_input_name_and_content => sub {
    my $t = shift;
    $t->fake_lambda;
    my $original = read_text($t->{ref});
    $t->probe($LAMBDA_BODY, ok => 0);
    ok(-e $t->{ref}, 'reference kept');
    is(read_text($t->{ref}), $original, 'reference content kept');
    ok(!-e ($t->{ref} =~ s/\.[^.\/]*\z/.fa/r), 'reference not renamed to .fa');
};

case test_lambda_index_refresh_preserves_other_database_files => sub {
    my $t = shift;
    $t->fake_lambda;
    my @sidecars = map { "$t->{ref}$_" } '.tax', '.notes', '.dna5.fm.sa.val';
    write_text($_, "preserve me\n") for @sidecars;
    my $current_index = write_text("$t->{ref}.lba.gz", 'x' x 101);
    # Encountering a legacy index must not remove arbitrary DB.* sidecars.
    $t->probe($LAMBDA_BODY, ok => 0);
    ok(-e $current_index, 'current index kept');
    is(read_text($_), "preserve me\n", "sidecar $_ kept") for @sidecars;
    # Explicit rebuild clears only the current Lambda3 index.
    $t->probe($LAMBDA_BODY, extra => ['-recalcTaxDB', '1'], ok => 0);
    ok(!-e $current_index, 'current index cleared');
    is(read_text($_), "preserve me\n", "sidecar $_ kept after rebuild") for @sidecars;
};

case test_tree_building_without_extended_logs => sub {
    my $t = shift;
    my $mafft = write_text("$t->{tools}/mafft", <<'PERL');
#!/usr/bin/env perl
open my $fh, '<', $ARGV[-1] or die "$ARGV[-1]: $!\n"; local $/; my $text = <$fh> // ''; print $text, "\n";
PERL
    my $tree = write_text("$t->{tools}/FastTree", <<'PERL');
#!/usr/bin/env perl
my ($i) = grep { $ARGV[$_] eq '-out' } 0 .. $#ARGV;
die "missing -out\n" unless defined $i;
open my $fh, '>', $ARGV[$i + 1] or die "$ARGV[$i + 1]: $!\n"; print {$fh} "(ASV1,ASV2);\n"; close $fh or die "$!\n";
PERL
    chmod 0755, $mafft, $tree or die "chmod: $!\n";
    @{ $t->{env} }{qw(AUDIT_MAFFT AUDIT_TREE)} = ($mafft, $tree);
    $t->probe(<<'PERL', extra => ['-extendedLogs', '0']);
ensure_dir($logDir); $buildPhylo = 1; $extendedLogs = 0; $uthreads = 1;
$mafftBin = $ENV{AUDIT_MAFFT}; $fasttreeBin = $ENV{AUDIT_TREE};
$lotus_tempDir = "$outdir/tmp"; ensure_dir($lotus_tempDir);
my $tree = buildTree($refDBwanted, $outdir); die "Missing tree" unless -s $tree;
PERL
    ok(-e "$t->{out}/ExtraFiles/OTU.MSA.fna", 'MSA written');
};

case test_count_chimeras_defaults_to_false_and_accepts_documented_values => sub {
    my $t = shift;
    my $body = 'print "CHIMCNT=$chimCnt\n";';
    for my $c ([[], 'F'], [['-count_chimeras', 'F'], 'F'], [['-count_chimeras', 'T'], 'T'],
               [['-count_chimeras'], 'T'], [['-count_chimeras', '-t', '1'], 'T']) {
        my ($extra, $expected) = @$c;
        subtest "extra=@$extra" => sub {
            contains($t->probe($body, extra => $extra)->{output}, "CHIMCNT=$expected\n");
        };
    }
    my $result = $t->probe($body, extra => ['-count_chimeras', 'maybe'], ok => 0);
    contains($result->{output}, '-count_chimeras must be T or F');
};

case test_usearch_versions_parse_from_both_banner_styles => sub {
    my $t = shift;
    $t->fake_usearch;
    my $body = '$usBin = $ENV{TEST_USEARCH}; detect_usearch_version();'
        . ' print "USV=$usearchVer|$usearchsubV|$usearchVerFull|"'
        . ' . (usearch_version_at_least("8.1") ? 1 : 0) . (usearch_version_at_least("10.0.241") ? 1 : 0) . "\n";';
    for my $c (['usearch v11.0.667_i86linux32', '', '11|0.667|11.0.667|11'],
               ['usearch v10.0.240_i86linux64, 16.3Gb RAM', '', '10|0.240|10.0.240|10'],
               ['usearch v8.1.1861_i86linux32', '', '8|1.1861|8.1.1861|10'],
               ['usearch v8.0.1623_i86linux32', '', '8|0.1623|8.0.1623|00'],
               ['usearch v12.0 [bd9d6e]', '', '12|0|12.0|11'],
               ['usearch v9.2.64_i86linux32', 'banner', '9|2.64|9.2.64|10']) {  # no --version: parse the banner
        my ($banner, $mode, $expected) = @$c;
        subtest "banner=$banner" => sub {
            @{ $t->{env} }{qw(TEST_USEARCH_BANNER TEST_USEARCH_MODE)} = ($banner, $mode);
            contains($t->probe($body)->{output}, "USV=$expected\n");
        };
    }
};

case test_usearch_denovo_chimera_check_restores_plain_otu_ids => sub {
    my $t = shift;
    $t->fake_usearch;
    my %seqs = (OTU_1 => 'ACGTACGTAA', OTU_2 => 'CCGGTTAACC', OTU_3 => 'TTTTGGGGCC');
    my $body = <<'PERL';
$lotus_tempDir = "$outdir/tmp"; ensure_dir($lotus_tempDir); ensure_dir($logDir); ensure_dir($extendedLogD);
$usBin = $ENV{TEST_USEARCH}; $VSBin = $usBin; $VSused = 0; $ClusterPipe = 8; $noChimChk = 0;
$useVsearch = $ENV{TEST_USE_VSEARCH}; $usearchVer = $ENV{TEST_USEARCH_MAJOR};
chimera_denovo("$outdir/otus.fna", "$outdir/otu_matrix.txt");
PERL
    # (useVsearch, usearch major) -> expected USEARCH chimera command and extra options
    for my $c (['0', '11', '-uchime3_denovo', []], ['-1', '12', '-uchime3_denovo', []],
               ['-1', '9', '-uchime2_denovo', ['-abskew', '16']], ['-1', '8', '-uchime_denovo', ['-abskew', '2']]) {
        my ($use_vsearch, $major, $command, $opts) = @$c;
        subtest "useVsearch=$use_vsearch usearch=$major" => sub {
            -d $t->{out} or mkdir $t->{out} or die "$t->{out}: $!\n";
            write_text("$t->{out}/otus.fna", join '', map { ">$_\n$seqs{$_}\n" } sort keys %seqs);
            write_text("$t->{out}/otu_matrix.txt", "OTU\ts1\ts2\nOTU_1\t1\t2\nOTU_2\t10\t5\nOTU_3\t7\t0\n");
            unlink $t->{usearch_calls};
            @{ $t->{env} }{qw(TEST_USE_VSEARCH TEST_USEARCH_MAJOR)} = ($use_vsearch, $major);
            $t->probe($body);
            my $text = read_text("$t->{out}/otus.fna");
            lacks($text, ';size=');
            my (undef, @entries) = split />/, $text, -1;
            my %got = map { my @l = text_lines($_); ($l[0] => join '', @l[1 .. $#l]) } @entries;
            is_deeply(\%got, \%seqs, 'plain OTU ids and sequences');
            my @calls = map { json_decode($_) } text_lines(read_text($t->{usearch_calls}));
            is_deeply([map { $_->[0] } @calls], ['-sortbysize', $command], 'USEARCH commands');
            is($calls[1][1], "$t->{out}/otus.fna.srt", 'chimera check reads the sorted FASTA');
            ok(has_arg($calls[1], $_), "chimera option $_") for @$opts;
        };
    }
};

case test_complete_ont_pipeline_outputs => sub {
    my $t = shift;
    $t->complete_pipeline;
};

case test_complete_barbell_pipeline_outputs_and_cleanup => sub {
    my $t = shift;
    $t->write_map(barbell => 1);
    my $original = read_text($t->{raw});
    my $original_map = read_text($t->{map});
    $t->complete_pipeline(barbell => 1);
    is(read_text($t->{raw}), $original, 'raw input unchanged');
    is(read_text($t->{map}), $original_map, 'user map unchanged');
    contains(read_text("$t->{out}/LotuSLogS/run_manifest.txt"), "Original input: $t->{raw}");
};

done_testing();

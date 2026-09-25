package LotusTest;
# Shared fixture for the LotuS3 regression tests: a temporary installation with
# the real bundled SDM, a lotus3 copy that exits once SDM has built the abundance
# matrix, and controlled Barbell/Savont/mapper stand-ins. The stand-ins make no
# scientific claims; they only exercise LotuS3's wiring.
use strict;
use warnings;
use Cwd qw(abs_path);
use Exporter 'import';
use File::Basename qw(basename dirname);
use File::Temp qw(tempdir);
use JSON::PP ();
use POSIX ();
use Test::More;

our @EXPORT_OK = qw($ROOT $SDM $TOOL $TOOL_MAPPER_BRANCH case contains lacks read_text write_text
    append_text text_lines json_decode json_encode rc run_command count_of replace_all replace_first
    after has_arg every4);

our $ROOT = abs_path(dirname(__FILE__) . '/../..');
our $SDM  = abs_path($ENV{LOTUS_TEST_SDM} // "$ROOT/bin/sdm") // ($ENV{LOTUS_TEST_SDM} // "$ROOT/bin/sdm");

# Fixed 1200-bp amplicon body (the historical fixture, random.Random(42) over ACGT).
my $SEQ = 'AAGCCCAATAAACCACTCTGACTGGCCGAATAGGGATATAGGCAACGACATGTGCGGCGACCCTTGCGACAGTGACGCTTTCGCCGTTGCCTAAACCTATTTGAAGGAGTCTAGCAGCCGCAGTAAGGCACAATACCTCGTCCGTGTTACCAGACCAAACAAGACGTCCTCTTCAATGTTTAAATGACCCTCTCGTCATAAAACCTTTCTACTATGTGTTCCGCAAGAATCAACAACTACAATGGCGCGTCGTGAATAACGCGACGGCTGAGACGAACGGCGCGTGAATGAAGCGCTTAAACAGCTCAGGAGCCAGTCCCCTACGTCGCATATCCTGGCCACTGGAGGTGAAGCGAATGGTATCGATACGTAGGAGGTGTGCCTTCGTAGGCTGTTTCTCAGGACGCCCAACTATTCTTTCCAATCCTACATCTGTTTCTTGCGTCGTAGCGGGACCCTCCATTGTTACTTATTAGGTTCTCGTTATGTCTCATAATCTCAGTGCTGGTGTGATAAGCAAACCACCCTACTGGCACGAAGTTCACAGAAGTGAGATTATGTCTCGTTTGGCAGTCTTGATGCTCGGGGGACACTTCTTTAAGCTCGGTGTGGTGGGCACGACCCTGGACGCGCGACGAAGCTAAGTTTGCAGTAATTAACCGACATCTTTGTGAACCGACCCACATTTGACGGTACGCTACCGCAACGGTATGTGTTAATGGAACAGACTTGCTTATGTGGACGTTGTATAGGGATATTACGTTACGCGTTAACCGATACATACTGGTTTCTCTCCAGTGGAGGTCTTGGTTGCCTCTAGTTTCTACGATATACTCATGGTAGTGTAACGCATAATCGAAGAGGGTCCTCCCATCTCCTGTGATGCATGGTGTGCTTACTGGGATGAATGCGCCGCAAGTAGCAGGTCCCGGCGTGGATACCTGATAGATGGTGACTAGCATGTACAAGTAACCTTGTCTATTGAGCTTCGAGGATGCATACAAGCCCACCCGCAGCCGCAACAGCGACGACTAATTGATCAGTAATTTATTAAGCACGGTGTTAACTTCTGTTTAGTGGGCTAAAATAGCAGATGTAGGGACCTCAGGAGCTAGACGGGGACCTACAACTTTGCGGGAACCAAGTTTTTGCAGTAGTGACTAACGCCGGGAATTCCTCGATATATAGTTTGATAGCTGA';

# Tests may add stand-in behaviour by inserting "elsif (...) {...}" branches
# in front of this line of $TOOL.
our $TOOL_MAPPER_BRANCH = "elsif (\$name eq 'minimap2' || \$name eq 'vsearch') {\n";

our $TOOL = <<'PERL' . $TOOL_MAPPER_BRANCH . <<'PERL';
#!/usr/bin/env perl
use strict; use warnings;
use File::Basename qw(basename); use File::Copy qw(copy); use JSON::PP ();
my $name = basename($0);
my @a = @ARGV;
sub arg { my ($flag) = @_; for my $i (0 .. $#a - 1) { return $a[$i+1] if $a[$i] eq $flag } die "$name: missing $flag\n" }
sub has { my ($flag) = @_; return scalar grep { $_ eq $flag } @a }
sub slurp { my ($f) = @_; open my $fh, '<', $f or die "$f: $!\n"; local $/; my $t = <$fh>; return $t // '' }
sub spew { my ($f, $t) = @_; open my $fh, '>', $f or die "$f: $!\n"; print {$fh} $t; close $fh or die "$f: $!\n" }
sub fasta { my @p = split />/, $_[0], -1; shift @p;
    return map { my @l = split /\n/; my ($id) = split ' ', $l[0]; [$id, join('', @l[1 .. $#l])] } @p }
if (has('--version') || has('-v')) {
    my %v = (minimap2 => '2.28', vsearch => 'vsearch v2.29.0', LCA => '0.28');
    print(($v{$name} // '1.0'), "\n"); exit 0;
}
open my $log, '>>', $ENV{ONT_TEST_CALLS} or die "$ENV{ONT_TEST_CALLS}: $!\n";
print {$log} JSON::PP->new->canonical->encode([$name, @a]), "\n";
close $log;
if ($name eq 'barbell') {
    my $out = arg('-o'); mkdir $out or die "$out: $!\n";
    copy($_, "$out/" . basename($_)) or die "$_: $!\n" for glob("$ENV{ONT_TEST_BARCODES}/*.fastq");
}
elsif ($name eq 'savont') {
    exit 17 if $ENV{ONT_TEST_FAIL};
    my $out = arg('-o'); mkdir $out or die "$out: $!\n";
    copy($a[1], "$out/input_snapshot.fq") or die "$a[1]: $!\n";
    unless ($ENV{ONT_TEST_EMPTY}) {
        my $seqs = $ENV{ONT_TEST_CONSENSUSES} ? JSON::PP::decode_json($ENV{ONT_TEST_CONSENSUSES}) : undef;
        $seqs = [$ENV{ONT_TEST_CONSENSUS}] unless $seqs && @$seqs;
        spew("$out/final_asvs.fasta", join '', map { ">final_consensus_$_ debug_id:$_ chimera_score:0\n$seqs->[$_]\n" } 0 .. $#$seqs);
    }
}
elsif ($name eq 'vsearch' && has('--cluster_size')) {
    my $query = arg('--cluster_size');
    die "unexpected clustering input $query\n" unless basename($query) eq 'derep.fas';
    spew(arg('--consout'), ">OTU0\n$ENV{ONT_TEST_CONSENSUS}\n");
    spew(arg('--uc'), join '', map { my $n = length $_->[1]; "H\t0\t$n\t99.9\t+\t0\t0\t${n}M\t$_->[0]\tOTU0\n" } fasta(slurp($query)));
}
PERL
    my ($query, $out, $db) = $name eq 'minimap2' ? ($a[-1], arg('-o'), $a[-2]) : (arg('--usearch_global'), arg('-uc'), arg('-db'));
    my $text = slurp($query);
    my @entries;
    if ($text =~ /^\@/) {
        my @l = split /\n/, $text;
        for (my $i = 0; $i < @l; $i += 4) { my ($id) = split ' ', substr($l[$i], 1); push @entries, [$id, $l[$i+1]] }
    } else {
        die "unexpected mapping input format\n" unless $text =~ /^>/;
        @entries = fasta($text);
    }
    my @targets = fasta(slurp($db));
    open my $fh, '>', $out or die "$out: $!\n";
    my $hits = 0;
    for my $e (@entries) {
        my ($id, $seq) = @$e; my $n = length $seq;
        next if $ENV{ONT_TEST_SKIP_PREFIX} && index($id, $ENV{ONT_TEST_SKIP_PREFIX}) == 0;
        $hits++;
        # closest target: mismatches over the shared prefix plus the length difference
        my ($target, $best);
        for my $t (@targets) {
            my $d = ($seq ^ $t->[1]) =~ tr/\0//c;
            ($target, $best) = ($t->[0], $d) if !defined $best || $d < $best;
        }
        print {$fh} $name eq 'minimap2'
            ? "$id\t$n\t0\t$n\t+\t$target\t$n\t0\t$n\t" . ($n - 1) . "\t$n\t60\tcg:Z:${n}M\n"
            : "H\t0\t$n\t99.9\t+\t0\t0\t${n}M\t$id\t$target\n";
    }
    close $fh or die "$out: $!\n";
    # The run summaries the real tools print to stderr.
    my $queries = @entries;
    print STDERR $name eq 'minimap2'
        ? "[M::worker_pipeline::0.001*1.00] mapped $queries sequences\n"
        : sprintf("Matching unique query sequences: %d of %d (%.2f%%)\n", $hits, $queries, $queries ? 100 * $hits / $queries : 0);
}
else {
    print STDERR "Unexpected tool call: $name\n"; exit 1;
}
PERL

sub read_text {
    my ($file) = @_;
    open my $fh, '<:raw', $file or die "Cannot read $file: $!\n";
    local $/; my $text = <$fh>;
    return $text // '';
}

sub write_text {
    my ($file, $text) = @_;
    open my $fh, '>:raw', $file or die "Cannot write $file: $!\n";
    print {$fh} $text;
    close $fh or die "Cannot close $file: $!\n";
    return $file;
}

sub append_text {
    my ($file, $text) = @_;
    open my $fh, '>>:raw', $file or die "Cannot append to $file: $!\n";
    print {$fh} $text;
    close $fh or die "Cannot close $file: $!\n";
}

sub text_lines { return split /\r?\n/, $_[0] // '' }
sub json_decode { return JSON::PP->new->decode($_[0]) }
sub json_encode { return JSON::PP->new->canonical->encode($_[0]) }
sub rc { my $s = reverse $_[0]; $s =~ tr/ACGT/TGCA/; return $s }
sub count_of { my ($text, $needle) = @_; my $n = () = $text =~ /\Q$needle\E/g; return $n }
sub replace_all { my ($text, $from, $to) = @_; $text =~ s/\Q$from\E/$to/g; return $text }
sub replace_first { my ($text, $from, $to) = @_; $text =~ s/\Q$from\E/$to/; return $text }
# Value following $flag in an argument list (undef when absent).
sub after { my ($args, $flag) = @_; for my $i (0 .. $#$args - 1) { return $args->[$i+1] if $args->[$i] eq $flag } return undef }
sub has_arg { my ($args, $flag) = @_; return scalar grep { $_ eq $flag } @$args }
# Every fourth line from $offset, e.g. FASTQ sequences with offset 1.
sub every4 { my ($lines, $offset) = @_; return [map { $lines->[$_] } grep { $_ % 4 == $offset } 0 .. $#$lines] }

sub contains {
    my ($haystack, $needle, $name) = @_;
    local $Test::Builder::Level = $Test::Builder::Level + 1;
    my $ok = ok(index($haystack // '', $needle) >= 0, $name // "contains '$needle'");
    diag('in: ' . substr($haystack // '', -3000)) unless $ok;
    return $ok;
}

sub lacks {
    my ($haystack, $needle, $name) = @_;
    local $Test::Builder::Level = $Test::Builder::Level + 1;
    my $ok = ok(index($haystack // '', $needle) < 0, $name // "lacks '$needle'");
    diag('in: ' . substr($haystack // '', -3000)) unless $ok;
    return $ok;
}

# Run a command with the given environment; stdout and stderr are merged.
sub run_command {
    my ($env, $timeout, @cmd) = @_;
    my $pid = open(my $fh, '-|') // die "Cannot fork: $!\n";
    if (!$pid) {
        %ENV = %$env;
        open(STDERR, '>&', \*STDOUT) or POSIX::_exit(127);
        exec { $cmd[0] } @cmd or do { print "Cannot execute $cmd[0]: $!\n"; POSIX::_exit(127) };
    }
    my $output = '';
    my $finished = eval {
        local $SIG{ALRM} = sub { die "timeout\n" };
        alarm $timeout;
        local $/; $output = <$fh> // '';
        alarm 0;
        1;
    };
    kill 'KILL', $pid unless $finished;
    close $fh;
    my $status = $?;
    $output .= "\n[timed out after ${timeout}s]\n" unless $finished;
    return ($output, $status);
}

# One subtest per test case, each with a fresh fixture (the former setUp).
sub case {
    my ($name, $code) = @_;
    subtest $name => sub {
        # A case passes its own run checks (e.g. run_lotus/probe) even without explicit assertions.
        eval { $code->(LotusTest->new); 1 } ? pass("$name completed") : fail("$name died: $@");
    };
}

sub new {
    my ($class) = @_;
    my $t = bless {}, $class;
    my $root = $t->{root} = tempdir('lotus3-ont-test-XXXXXXXX', TMPDIR => 1, CLEANUP => 1);
    my $install = $t->{install} = "$root/install";
    mkdir $install or die "$install: $!\n";
    for my $name (qw(configs helpers Example)) {
        symlink("$ROOT/$name", "$install/$name") or die "symlink $name: $!\n";
    }
    mkdir "$install/bin" or die "$install/bin: $!\n";
    symlink($SDM, "$install/bin/sdm") or die "symlink sdm: $!\n";
    my $src = read_text("$ROOT/lotus3");
    my $checkpoint = 'undef $tmpOTU;';
    my $at = index($src, $checkpoint);
    die "lotus3 checkpoint missing or ambiguous\n" if $at < 0 || index($src, $checkpoint, $at + 1) >= 0;
    substr($src, $at, 0) = <<'PERL';
atomic_write_text("$outdir/ont_test_state.json", JSON::PP->new->encode({
    seed => $OTUSEED, map => $mapHref, combined => $combHref,
    sdm_options => $sdmOpt, cluster => $ClusterPipe, preset => $mini2RdPreset, dereplication => $sdmDerepDo}));
write_repro_manifest(); release_output_lock(); exit(0);
PERL
    $t->{script} = write_text("$install/lotus3", $src);
    my $tools = $t->{tools} = "$root/tools";
    mkdir $tools or die "$tools: $!\n";
    my @names = qw(barbell savont minimap2 vsearch LCA);
    for my $name (@names) { chmod 0755, write_text("$tools/$name", $TOOL) or die "chmod $name: $!\n" }
    $t->{cfg} = write_text("$root/lotus.cfg", "sdm $SDM\n" . join('', map { "$_ $tools/$_\n" } @names) . "CheckForUpdates 0\n");
    $t->{seq} = $SEQ;
    $t->{consensus} = (substr($SEQ, 0, 1) ne 'T' ? 'T' : 'A') . substr($SEQ, 1);
    $t->{fwd} = 'AGAGTTTGATCCTGGCTCAG';
    $t->{rev} = 'TACGGYTACCTTGTTACGACTT';
    ($t->{revcomp} = rc($t->{rev} =~ tr/Y/T/r));
    $t->{reads} = "$root/reads";       mkdir $t->{reads} or die "$!\n";
    $t->{barcodes} = "$root/barcodes"; mkdir $t->{barcodes} or die "$!\n";
    for my $spec (['s1', 'BC01', 4], ['s2', 'BC02', 3], ['low', 'BC03', 1], ['noise', 'BC99', 2]) {
        my ($sample, $barcode, $count) = @$spec;
        my $read = $t->{fwd} . $SEQ . $t->{revcomp};
        my $fq = join '', map { "\@${sample}_$_\n$read\n+\n" . ('I' x length $read) . "\n" } 0 .. $count - 1;
        write_text("$t->{reads}/$sample.fq", $fq);
        write_text("$t->{barcodes}/$barcode.trimmed.fastq", $fq);
    }
    $t->{raw} = write_text("$root/raw.fastq", read_text("$t->{reads}/s1.fq"));
    $t->{map} = "$root/map.tsv";
    $t->write_map;
    $t->{ref} = write_text("$root/ref.fna", ">ref\n$SEQ\n");
    $t->{tax} = write_text("$root/ref.tax", "ref\tBacteria;P;C;O;F;G;S\n");
    $t->{out} = "$root/output";
    $t->{calls} = "$root/calls.jsonl";
    $t->{env} = { %ENV, ONT_TEST_CALLS => $t->{calls}, ONT_TEST_BARCODES => $t->{barcodes}, ONT_TEST_CONSENSUS => $t->{consensus} };
    return $t;
}

sub write_map {
    my ($t, %o) = @_;
    my $header = $o{header} // ($o{barbell} ? 'ONTBarcode' : 'fastqFile');
    my $rows = $o{rows} // ($o{barbell}
        ? [['s1', 'BC01'], ['s2', 'BC02'], ['low', 'BC03'], ['missing', 'BC04']]
        : [['s1', 's1.fq'], ['s2', 's2.fq']]);
    write_text($t->{map}, "#SampleID\t$header\tForwardPrimer\tReversePrimer\n"
        . join('', map { "$_->[0]\t$_->[1]\t$t->{fwd}\t$t->{rev}\n" } @$rows));
}

# Returns {output, status}. With ok (default) a failed run dies with its log tail;
# with ok => 0 a successful run is recorded as a test failure.
sub run_lotus {
    my ($t, %o) = @_;
    my $ok = $o{ok} // 1;
    my @args = ('perl', $t->{script}, '-i', ($o{barbell} ? $t->{raw} : $t->{reads}), '-m', $t->{map},
        '-o', $t->{out}, '-c', $t->{cfg}, '-p', 'ONT', '-t', '1', '-lulu', '0', '-removePhiX', '0',
        '-buildPhylo', '0', '-deactivateChimeraCheck', '1', '-refDB', $t->{ref}, '-tax4refDB', $t->{tax},
        '-taxAligner', 'vsearch');
    push @args, '-ontDemux', 'barbell', '-ontKit', 'SQK-RBK114-96' if $o{barbell};
    my ($output, $status) = run_command($t->{env}, 30, @args, @{ $o{extra} // [] });
    write_text("$t->{root}/run.log", $output);
    if ($ok && $status != 0) {
        my $prog = "$t->{out}/LotuSLogS/LotuS_progout.log";
        die "lotus3 failed (status $status):\n" . substr($output, -6000)
            . (-e $prog ? substr(read_text($prog), -6000) : '') . "\n";
    }
    isnt($status, 0, 'lotus3 fails as expected') or diag(substr($output, -3000)) if !$ok;
    return { output => $output, status => $status };
}

# Replace the lotus3 copy with the original, running $body instead of prepLtsOptions().
sub probe {
    my ($t, $body, %o) = @_;
    my $src = read_text("$ROOT/lotus3");
    my $at = index($src, 'prepLtsOptions();');
    die "prepLtsOptions() call not found\n" if $at < 0;
    substr($src, $at, length 'prepLtsOptions();') = "$body\nexit(0);";
    write_text($t->{script}, $src);
    return $t->run_lotus(%o);
}

# Simulate a previous LotuS3 run in the output directory; returns its program log.
sub old_output {
    my ($t) = @_;
    -d $_ or mkdir $_ or die "$_: $!\n" for $t->{out}, "$t->{out}/LotuSLogS";
    write_text("$t->{out}/.lotus3_created_by_this_run", "previous run\n");
    return write_text("$t->{out}/LotuSLogS/LotuS_progout.log", "active run output\n");
}

sub tool_calls {
    my ($t, $name) = @_;
    return [] unless -e $t->{calls};
    return [grep { !defined $name || $_->[0] eq $name } map { json_decode($_) } text_lines(read_text($t->{calls}))];
}

sub sdm_commands {
    my ($t) = @_;
    return [grep { index($_, $SDM) >= 0 && index($_, '-sample_sep') >= 0 }
        text_lines(read_text("$t->{out}/LotuSLogS/LotuS_cmds.log"))];
}

sub check_counts {
    my ($t, $expected) = @_;
    local $Test::Builder::Level = $Test::Builder::Level + 1;
    my @rows = text_lines(read_text("$t->{out}/OTU.txt"));
    is(scalar @rows, 2, 'one OTU row') or diag(join "\n", @rows);
    my @head = split /\t/, $rows[0];
    my @row = split /\t/, $rows[1] // '';
    is($row[0], 'OTU_0', 'OTU id');
    my %got; @got{ @head[1 .. $#head] } = map { 0 + $_ } @row[1 .. $#row];
    is_deeply(\%got, $expected // { s1 => 4, s2 => 3 }, 'sample counts');
    my $state = json_decode(read_text("$t->{out}/ont_test_state.json"));
    return $state;
}

1;

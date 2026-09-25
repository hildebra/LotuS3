#!/usr/bin/env perl
# Regressions for the 2026-09-25 pipeline audit fixes: command splitting and failure
# reporting, tree-building failures, mapping-file checks, SINTAX ranks, reference
# database checks, -xtalk, path depth checks and -create_map.
use strict;
use warnings;
use FindBin;
use Test::More;
use lib "$FindBin::Bin/lib";
use LotusTest qw($ROOT case contains lacks read_text write_text text_lines run_command);

# probe() bodies run in place of prepLtsOptions(), so the log directory is not set up yet.
my $LOGS = 'ensure_dir($logDir); open(LOG, ">>", $mainLogFile) or die; open(cmdLOG, ">>", $cmdLogFile) or die;';

case test_newlines_separate_checked_commands => sub {
    my $t = shift;
    $t->probe(<<'PERL');
my @c = split_shell_commands("a 1\nb 2; c 3\n'd\n4'");
die "got @{[scalar @c]} commands: @c" unless @c == 4 && $c[0] eq "a 1" && $c[3] eq "'d\n4'";
die "status text" unless describe_system_status(3 << 8) eq "exit code 3"
    && describe_system_status(11) eq "killed by signal 11" && describe_system_status(-1) eq "could not be started";
PERL
};

case test_failing_iqtree_aborts_tree_building => sub {
    my $t = shift;
    my $mafft = write_text("$t->{root}/mafft", "#!/bin/sh\nfor a; do last=\$a; done\ncat \"\$last\"\n");
    my $iqtree = write_text("$t->{root}/iqtree2", "#!/bin/sh\nexit 3\n");
    chmod 0755, $mafft, $iqtree;
    my $r = $t->probe(<<"PERL" . <<'PERL', ok => 0);
$LOGS
(\$mafftBin, \$iqTreeBin) = ("$mafft", "$iqtree");
PERL
$buildPhylo = 2; $uthreads = 1;
$lotus_tempDir = "$outdir/tmpFiles"; ensure_dir($lotus_tempDir);
my $fa = "$outdir/otus.fna"; open(my $o, ">", $fa) or die; print $o ">OTU_1\nACGT\n"; close $o;
buildTree($fa, $outdir);
print "TREE BUILT\n";
PERL
    is($r->{status} >> 8, 9, 'exit status 9');
    contains($r->{output}, 'Tree building: CMD failed (exit code 3)');
    lacks($r->{output}, 'TREE BUILT');
    contains(read_text("$t->{out}/LotuSLogS/LotuS_run.log"), 'CMD failed (exit code 3)');
};

case test_empty_sequencing_run_aborts => sub {
    my $t = shift;
    write_text($t->{map}, "#SampleID\tfastqFile\tForwardPrimer\tReversePrimer\tSequencingRun\n"
        . "s1\ts1.fq\t$t->{fwd}\t$t->{rev}\trun1\ns2\ts2.fq\t$t->{fwd}\t$t->{rev}\n");
    my $r = $t->run_lotus(ok => 0);
    is($r->{status} >> 8, 66, 'exit status 66');
    contains($r->{output}, "SampleID 's2' has an empty or whitespace-padded \"SequencingRun\" value");
};

case test_copied_map_starts_with_header => sub {
    my $t = shift;
    $t->probe(<<'PERL');
my $f = "$outdir/copy.map"; ensure_dir($outdir);
writeMap({ '!early' => { s1 => ['a.fq', '!early'] }, '#SampleID' => { '#SampleID' => ['fastqFile', 'SequencingRun'] },
    run2 => { s2 => ['b.fq', 'run2'] } }, $f);
open(my $h, "<", $f) or die; my @l = <$h>;
die "header not first: @l" unless $l[0] eq "#SampleID\tfastqFile\tSequencingRun\n" && @l == 3;
PERL
};

case test_whitespace_and_dashes_in_sample_names_abort => sub {
    my $t = shift;
    for my $spec (["s 1\ts1.fq", "SampleID 's 1' contains whitespace"], ["s1\ts1.fq\tgrp-1", "CombineSamples value 'grp-1'"]) {
        my ($row, $message) = @$spec;
        my $combine = $row =~ /grp/ ? "\tCombineSamples" : "";
        write_text($t->{map}, "#SampleID\tfastqFile$combine\n$row\n");
        my $r = $t->run_lotus(ok => 0, extra => ['-forwardPrimer', $t->{fwd}, '-reversePrimer', $t->{rev}]);
        is($r->{status} >> 8, 5, "exit status 5 for '$row'");
        contains($r->{output}, $message);
    }
};

case test_sintax_ranks_follow_rank_letters => sub {
    my $t = shift;
    $t->probe(<<'PERL');
ensure_dir($outdir); $utaxConf = 0.8;
my $ut = "$outdir/sintax.txt"; open(my $o, ">", $ut) or die;
print $o "OTU_1\td:Bacteria(1.00),p:Firmicutes(0.99),g:Lactobacillus(0.95),s:L_casei(0.50)\t+\n";
print $o "OTU_2\tk:Fungi(1.00),d:Other(1.00),p:Ascomycota(0.90)\t+\n";
close $o;
writeUTAXhiera($ut, ['OTU_1', 'OTU_2'], {});
open(my $h, "<", $SIM_hierFile) or die; my @l = map { chomp; [split /\t/, $_, -1] } <$h>;
die "OTU_1: @{$l[1]}" unless join(",", @{$l[1]}) eq "OTU_1,Bacteria,Firmicutes,?,?,?,Lactobacillus,?";
die "OTU_2: @{$l[2]}" unless join(",", @{$l[2]}) eq "OTU_2,Fungi,Ascomycota,?,?,?,?,?";
PERL
};

case test_reference_taxonomy_without_matching_ids_aborts => sub {
    my $t = shift;
    write_text($t->{tax}, "unrelated\tBacteria;P;C;O;F;G;S\n");
    my $r = $t->run_lotus(ok => 0);
    is($r->{status} >> 8, 55, 'exit status 55');
    contains($r->{output}, 'None of the first reference IDs');
};

case test_xtalk_needs_real_usearch => sub {
    my $t = shift;
    my $r = $t->run_lotus(ok => 0, extra => ['-xtalk', '1']);
    is($r->{status} >> 8, 83, 'exit status 83');
    contains($r->{output}, '-xtalk needs USEARCH 11');
};

case test_option_values_are_validated => sub {
    my $t = shift;
    for my $spec ([['-redoTaxOnly', '2'], '-redoTaxOnly must be 0 or 1'], [['-VXtr', '2'], '-VXtr must be 0 or 1'],
                  [['-endRem', 'ACG(T'], '-endRem must be']) {
        my ($extra, $message) = @$spec;
        my $r = $t->run_lotus(ok => 0, extra => $extra);
        contains($r->{output}, $message);
    }
};

case test_path_depth_not_length_guards_removal => sub {
    my $t = shift;
    $t->probe(<<'PERL');
die "root" unless _path_too_shallow("/") && _path_too_shallow("/tmp") && _path_too_shallow("/tmp/");
die "short but deep" if _path_too_shallow("/tmp/o1") || _path_too_shallow("/a/b");
PERL
};

case test_create_map_pairs_by_file_stem => sub {
    my $t = shift;
    my $dir = "$t->{root}/fastq";
    mkdir $dir or die "$dir: $!\n";
    write_text("$dir/$_", '') for qw(A_S1_L001_R1_001.fastq.gz A_S1_L001_R2_001.fastq.gz B-2_R1.fastq.gz B-2_R2.fastq.gz C_1.fq C_2.fq);
    my $map = "$t->{root}/auto.map";
    my ($output, $status) = run_command($t->{env}, 30, 'perl', "$ROOT/lotus3", '-create_map', $map, '-i', $dir);
    is($status, 0, 'map created') or diag($output);
    is(read_text($map), "#SampleID\tfastqFile\nA\tfastq/A_S1_L001_R1_001.fastq.gz,fastq/A_S1_L001_R2_001.fastq.gz\n"
        . "B_2\tfastq/B-2_R1.fastq.gz,fastq/B-2_R2.fastq.gz\nC\tfastq/C_1.fq,fastq/C_2.fq\n", 'pairs and names');
    write_text("$dir/D_R1.fastq.gz", '');
    ($output, $status) = run_command($t->{env}, 30, 'perl', "$ROOT/lotus3", '-create_map', $map, '-i', $dir);
    isnt($status, 0, 'a file without its mate fails');
    contains($output, 'these files have no mate: D_R1.fastq.gz');
};

done_testing();

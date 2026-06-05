#!/usr/bin/env perl
# savont2uc.pl  --  convert savont ASV output into the DADA2-style .uc read->ASV
#                   mapping + uniqueSeqs.fna that LotuS3's sdm step consumes.
#
# WHY THIS EXISTS
#   In the LotuS3 ONT path, savont is the ASV caller (analogous to DADA2). Like the
#   DADA2 R script (bin/R/dada2_pip_v2.R, which writes dada2.uc + uniqueSeqs.fna),
#   savont must hand LotuS3 two things:
#     (1) a .uc file mapping every assigned READ to its ASV  -> used by sdm to build
#         the per-sample ASV abundance matrix (sdm already knows read->sample from demux)
#     (2) a fasta of the ASV sequences (uniqueSeqs.fna)       -> used for taxonomy & tree
#   savont does NOT emit per-sample counts itself (it concatenates all inputs); the
#   per-sample tally is sdm's job, exactly as for DADA2. This script only reshapes
#   savont's native read->ASV assignments into the .uc format sdm expects.
#
# INPUTS
#   --asvs   final_asvs.fasta                      (savont output dir, top level)
#   --map    temp/read_to_asv_mappings.tsv         (savont output dir, temp/)
# OUTPUTS
#   --ucout  dada2-style .uc (read->ASV)
#   --fnaout uniqueSeqs.fna (ASV sequences, one representative-read header per ASV)
#
# .uc FORMAT (mirrors bin/R/dada2_pip_v2.R lines 547/551 exactly):
#   seed read   :  <readID>\t otu<i> \t *
#   member read :  <readID>\t match  \t dqt=1;top=<seedReadID>(99%);
#   ASV internal tag is "otu"+<i> (i=1..N, in fasta order) to match DADA2 (ASVname="otu").
#
# read_to_asv_mappings.tsv FORMAT (savont, tab-separated, NO header):
#   <readID>\t asv:<debug_id>\t <mismatches>\t <alignment_score>
#   - readID may itself contain a space (e.g. "SRR..7110 7110/1"); we keep the FULL
#     first tab-delimited field as the read ID, then normalise to its first token to
#     match how sdm/demux name reads (configurable via --idmode).
#   - a read may appear on MULTIPLE rows (aligns to several ASVs); we keep ONE best hit.
#
# BEST-HIT SELECTION (deterministic):
#   per read, choose: min mismatches, then max alignment_score, then first-seen.
#
# SURVIVOR FILTERING:
#   final_asvs.fasta defines which debug_ids survived quality/chimera filtering. Reads
#   mapping only to non-surviving debug_ids are dropped (correctly unassigned).

use strict;
use warnings;
use Getopt::Long;

my ($asvFile, $mapFile, $ucOut, $fnaOut);
my $idmode  = "firsttoken";   # firsttoken | full  -- how to render the read ID
my $verbose = 0;

GetOptions(
    "asvs=s"   => \$asvFile,
    "map=s"    => \$mapFile,
    "ucout=s"  => \$ucOut,
    "fnaout=s" => \$fnaOut,
    "idmode=s" => \$idmode,
    "verbose"  => \$verbose,
) or die "Bad options\n";

die "Usage: $0 --asvs final_asvs.fasta --map read_to_asv_mappings.tsv --ucout out.uc --fnaout out.fna [--idmode firsttoken|full]\n"
    unless ($asvFile && $mapFile && $ucOut && $fnaOut);

die "ASV fasta not found: $asvFile\n" unless -e $asvFile;
die "mapping tsv not found: $mapFile\n" unless -e $mapFile;

# ---------------------------------------------------------------------------
# 1. Parse final_asvs.fasta: learn surviving debug_ids, their order (-> otu index),
#    and their consensus sequences.
# ---------------------------------------------------------------------------
my %dbg2otu;     # debug_id -> otu index (1-based, in fasta order)
my %dbg2seq;     # debug_id -> consensus sequence
my @dbgOrder;    # debug_ids in fasta order
my $otuIdx = 0;
{
    open(my $fh, "<", $asvFile) or die "open $asvFile: $!\n";
    my $curDbg;
    while (my $line = <$fh>) {
        chomp $line;
        if ($line =~ /^>/) {
            # header e.g. ">final_consensus_0_depth_1433 debug_id:7 chimera_score:0 ..."
            if ($line =~ /debug_id:(\d+)/) {
                $curDbg = $1;
                $otuIdx++;
                $dbg2otu{$curDbg} = $otuIdx;
                push @dbgOrder, $curDbg;
                $dbg2seq{$curDbg} = "";
            } else {
                die "ASV header without debug_id (cannot link to mapping):\n$line\n";
            }
        } else {
            $dbg2seq{$curDbg} .= $line if defined $curDbg;
        }
    }
    close($fh);
}
my $nASV = scalar @dbgOrder;
die "No ASVs parsed from $asvFile (0 surviving consensi)\n" if $nASV == 0;
print STDERR "[savont2uc] parsed $nASV surviving ASVs (debug_ids: ".join(",",@dbgOrder).")\n" if $verbose;

# ---------------------------------------------------------------------------
# 2. Walk read_to_asv_mappings.tsv, keep best hit per read (survivors only).
# ---------------------------------------------------------------------------
my %bestDbg;    # readID -> chosen debug_id
my %bestMM;     # readID -> mismatches of chosen
my %bestSc;     # readID -> score of chosen
my $rowsTotal = 0; my $rowsDropped = 0;
{
    open(my $fh, "<", $mapFile) or die "open $mapFile: $!\n";
    while (my $line = <$fh>) {
        chomp $line;
        next if $line eq "";
        my @f = split(/\t/, $line);
        # expected: readID, asv:<dbg>, mismatches, score
        next unless (@f >= 4);
        my ($readRaw, $asvTag, $mm, $sc) = @f[0,1,2,3];
        next unless ($asvTag =~ /^asv:(\d+)$/);
        my $dbg = $1;
        $rowsTotal++;
        # survivor filter
        unless (exists $dbg2otu{$dbg}) { $rowsDropped++; next; }
        $mm += 0; $sc += 0;   # numeric

        my $readID = $readRaw;
        if ($idmode eq "firsttoken") {
            $readID =~ s/\s.*$//;   # keep up to first whitespace ("SRR..7110 7110/1" -> "SRR..7110")
        }

        if (!exists $bestDbg{$readID}
            || $mm < $bestMM{$readID}
            || ($mm == $bestMM{$readID} && $sc > $bestSc{$readID})) {
            $bestDbg{$readID} = $dbg;
            $bestMM{$readID}  = $mm;
            $bestSc{$readID}  = $sc;
        }
    }
    close($fh);
}
my $nReads = scalar keys %bestDbg;
print STDERR "[savont2uc] read $rowsTotal mapping rows; dropped $rowsDropped (non-survivor); assigned $nReads unique reads\n" if $verbose;
die "No reads assigned to surviving ASVs - aborting (check inputs)\n" if $nReads == 0;

# ---------------------------------------------------------------------------
# 3. Group reads by ASV, emit .uc (seed first, then members) + fna.
# ---------------------------------------------------------------------------
my %reads_by_dbg;   # debug_id -> [readIDs]
for my $r (sort keys %bestDbg) {    # sort -> deterministic seed selection & output
    push @{ $reads_by_dbg{ $bestDbg{$r} } }, $r;
}

open(my $uc,  ">", $ucOut)  or die "open $ucOut: $!\n";
open(my $fna, ">", $fnaOut) or die "open $fnaOut: $!\n";

my $written_asv = 0; my $written_reads = 0;
for my $dbg (@dbgOrder) {                  # fasta order -> otu1, otu2, ...
    my $otu = "otu" . $dbg2otu{$dbg};
    my $reads = $reads_by_dbg{$dbg};
    if (!defined $reads || !@$reads) {
        # ASV survived but had no best-hit reads assigned (possible if all its reads
        # were won by another ASV on tie-break). Still emit the ASV sequence so it is
        # not lost, using its own debug tag as a synthetic seed label.
        print STDERR "[savont2uc] WARNING: ASV $otu (debug_id:$dbg) has 0 assigned reads after best-hit selection\n" if $verbose;
        my $seedSynthetic = "savontASV_${dbg}";
        print $uc  "$seedSynthetic\t$otu\t*\n";
        print $fna ">$seedSynthetic\n$dbg2seq{$dbg}\n";
        $written_asv++;
        next;
    }
    my $seed = $reads->[0];                 # first (sorted) read is the seed label
    # seed line + ASV consensus sequence in fna (header = seed read ID, seq = savont consensus)
    print $uc  "$seed\t$otu\t*\n";
    print $fna ">$seed\n$dbg2seq{$dbg}\n";
    $written_reads++;
    # member lines
    for my $i (1 .. $#$reads) {
        my $r = $reads->[$i];
        print $uc "$r\tmatch\tdqt=1;top=$seed(99%);\n";
        $written_reads++;
    }
    $written_asv++;
}
close($uc); close($fna);

print STDERR "[savont2uc] wrote $written_asv ASVs and $written_reads read records to:\n  $ucOut\n  $fnaOut\n";

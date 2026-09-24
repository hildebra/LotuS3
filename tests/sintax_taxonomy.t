#!/usr/bin/env perl
# SINTAX (-taxAligner sintax|utax) wiring: stand-in usearch, real SDM.
#
# Run: prove -v tests/sintax_taxonomy.t
# The stand-in usearch writes sintax -tabbedout rows chosen by each case; no
# classifier accuracy is implied. Full runs reuse the ONT fixture to reach the
# taxonomy tables; the LCA stand-in fails if SINTAX mode ever calls LCA.
use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/lib";
use File::Spec;
use Test::More;
use LotusTest qw($ROOT case contains read_text write_text text_lines json_decode json_encode after run_command);

sub canon { return File::Spec->canonpath($_[0]) }

my $USEARCH = <<'PERL';
#!/usr/bin/env perl
use strict; use warnings; use JSON::PP ();
my @a = @ARGV;
if (!@a || grep { $_ eq '--version' } @a) { print +($ENV{SINTAX_TEST_BANNER} // 'usearch v11.0.667_i86linux32'), "\n"; exit 0 }
open my $log, '>>', $ENV{ONT_TEST_CALLS} or die "$ENV{ONT_TEST_CALLS}: $!\n";
print {$log} JSON::PP->new->encode(['usearch', @a]), "\n"; close $log;
sub arg { my ($flag) = @_; for my $i (0 .. $#a - 1) { return $a[$i+1] if $a[$i] eq $flag } die "usearch: missing $flag\n" }
if ($a[0] eq '-makeudb_usearch') {
    die "usearch: missing reference $a[1]\n" unless -f $a[1];
    open my $out, '>', arg('-output') or die "$!\n"; print {$out} "stand-in UDB\n"; close $out;
}
elsif ($a[0] eq '-sintax') {
    die "usearch: missing UDB\n" unless -f arg('-db');
    my $predictions = JSON::PP->new->decode($ENV{SINTAX_TEST_PREDICTIONS});
    open my $in, '<', $a[1] or die "$a[1]: $!\n";
    open my $out, '>', arg('-tabbedout') or die "$!\n";
    while (my $line = <$in>) {
        next unless $line =~ /^>(.*?)\r?$/;
        my $label = $1; my ($id) = split ' ', $label;  # usearch reports the full header
        my $prediction = $predictions->{$id};
        next unless defined $prediction;               # absent: the query gets no row
        print {$out} "$label\t$prediction\t" . ($prediction eq '' ? '' : '+') . "\t\n";
    }
    close $out;
}
else { print STDERR "Unexpected usearch call: @a\n"; exit 1 }
PERL

# Checks what bin/R/l2phyloseq.R reads, then saves what it would save.
my $RSCRIPT = <<'PERL';
#!/usr/bin/env perl
use strict; use warnings; use File::Basename qw(basename dirname); use JSON::PP ();
open my $log, '>>', $ENV{ONT_TEST_CALLS} or die "$ENV{ONT_TEST_CALLS}: $!\n";
print {$log} JSON::PP->new->encode(['Rscript', @ARGV]), "\n"; close $log;
my ($script, $table, $tax) = grep { $_ ne '--vanilla' } @ARGV;
die "unexpected R script $script\n" unless basename($script) eq 'l2phyloseq.R';
sub lines { open my $fh, '<', $_[0] or die "$_[0]: $!\n"; chomp(my @l = <$fh>); return @l }
my @rows = map { [split /\t/, $_, -1] } lines($tax);
for (@rows) { die "hiera row needs an OTU ID and seven ranks: @$_\n" unless @$_ == 8 }
my %ids = map { $_->[0] => 1 } @rows[1 .. $#rows];
my (undef, @otus) = lines($table);
for (map { (split /\t/)[0] } @otus) { die "$_ missing from $tax\n" unless $ids{$_} }
open my $out, '>', dirname($table) . '/phyloseq.Rdata' or die "$!\n"; print {$out} "stand-in phyloseq object\n"; close $out;
PERL

my $FIRMICUTES = 'd:Bacteria(1.0000),p:"Firmicutes"(0.9900),c:Bacilli(0.9700),o:Lactobacillales(0.9500),'
    . 'f:Streptococcaceae(0.9000),g:Streptococcus(0.7000),s:Streptococcus_mitis(0.4000)';
my $HEADER = join "\t", qw(Domain Phylum Class Order Family Genus Species);
my $UNASSIGNED = "\t?" x 7;

sub setup_sintax {
    my ($t) = @_;
    write_text($t->{script}, read_text("$ROOT/lotus3"));  # whole pipeline, no checkpoint
    $t->{usearch} = write_text("$t->{tools}/usearch", $USEARCH);
    chmod 0755, $t->{usearch} or die "chmod usearch: $!\n";
    my $rbin = "$t->{root}/rbin";
    mkdir $rbin or die "$rbin: $!\n";
    chmod 0755, write_text("$rbin/Rscript", $RSCRIPT) or die "chmod Rscript: $!\n";
    $t->{env}{PATH} = "$rbin:$t->{env}{PATH}";
    for my $marker (qw(SSU ITS)) {
        my $db = $t->{db}{$marker} = "$t->{root}/sintax_$marker";
        mkdir $_ or die "$_: $!\n" for $db, "$db/fasta";
        write_text("$db/fasta/refdb.fa", ">ref;tax=d:Bacteria,p:Firmicutes;\n$t->{seq}\n");
    }
    $t->{base_cfg} = read_text($t->{cfg});
    $t->configure;
    $t->predict({ ASV1 => $FIRMICUTES });
    return $t;
}

sub configure {
    my ($t, %o) = @_;
    # R helper paths as in a standard install (-lulu defaults to 1, even for -taxOnly)
    my $cfg = $t->{base_cfg} . "phyloLnk $ROOT/bin/R/l2phyloseq.R\nLULUR $ROOT/bin/R/LULU.R\n";
    $cfg .= "usearch $t->{usearch}\n" unless $o{no_usearch};
    $cfg .= 'TAX_REFDB_SSU_UTAX ' . ($o{ssu_db} // $t->{db}{SSU}) . "\nTAX_REFDB_ITS_UTAX $t->{db}{ITS}\n";
    write_text($t->{cfg}, $cfg);
}

sub predict { my ($t, $predictions) = @_; $t->{env}{SINTAX_TEST_PREDICTIONS} = json_encode($predictions) }

sub lotus {
    my ($t, $args, %o) = @_;
    my $ok = $o{ok} // 1;
    my ($output, $status) = run_command($t->{env}, 60, 'perl', $t->{script}, '-c', $t->{cfg}, '-t', '1', @$args);
    write_text("$t->{root}/run.log", $output);
    die "lotus3 failed (status $status):\n" . substr($output, -6000) . "\n" if $ok && $status != 0;
    isnt($status, 0, 'lotus3 fails as expected') or diag(substr($output, -3000)) if !$ok;
    return { output => $output, status => $status, code => $status >> 8 };
}

sub run_sintax {
    my ($t, %o) = @_;
    return $t->lotus(['-i', $t->{reads}, '-m', $t->{map}, '-o', $o{out} // $t->{out}, '-p', 'ONT', '-lulu', '0',
        '-removePhiX', '0', '-buildPhylo', '0', '-deactivateChimeraCheck', '1', '-taxAligner', 'sintax',
        @{ $o{extra} // [] }], ok => $o{ok});
}

sub tax_only {
    my ($t, $query, $out, $extra, %o) = @_;
    return $t->lotus(['-taxOnly', $query, '-o', $out, @$extra], %o);
}

sub hiera { return [text_lines(read_text("$_[0]/hiera_BLAST.txt"))] }

# Make the helpers above callable as fixture methods.
{ no strict 'refs'; *{"LotusTest::$_"} = \&{"main::$_"} for qw(setup_sintax configure predict lotus run_sintax tax_only); }

case test_sintax_run_writes_lca_layout_tables_biom_and_phyloseq_input => sub {
    my $t = setup_sintax(shift);
    # An unparsed banner would fall back to the vsearch-implied version 11.
    $t->{env}{SINTAX_TEST_BANNER} = 'usearch v10.0.240_i86linux32';
    $t->run_sintax;
    is_deeply(hiera($t->{out}), ["ASV\t$HEADER", "ASV1\tBacteria\tFirmicutes\tBacilli\tLactobacillales\tStreptococcaceae\t?\t?"],
        'OTU-first hierarchy');
    ok(!-e "$t->{out}/hiera_RDP.txt", 'no RDP hierarchy');
    my @phylum = text_lines(read_text("$t->{out}/higherLvl/Phylum.txt"));
    is_deeply([@phylum[0, 1]], ["Phylum\ts1\ts2", "Bacteria;Firmicutes\t4\t3"], 'phylum table');
    contains(read_text("$t->{out}/higherLvl/Genus.txt"), "Bacteria;Firmicutes;Bacilli;Lactobacillales;Streptococcaceae;?\t4\t3");
    my $biom = json_decode(read_text("$t->{out}/OTU.biom"));
    is(json_encode($biom->{data}), '[[4,3]]', 'BIOM counts');
    is_deeply($biom->{rows}[0]{metadata}{taxonomy},
        [qw(k__Bacteria p__Firmicutes c__Bacilli o__Lactobacillales f__Streptococcaceae g__? s__?)], 'BIOM taxonomy');
    my $rscript = $t->tool_calls('Rscript');
    is(scalar @$rscript, 1, 'phyloseq helper called once');
    is(canon($rscript->[0][4]), canon("$t->{out}/hiera_BLAST.txt"), 'phyloseq reads hiera_BLAST.txt');
    ok(-e "$t->{out}/phyloseq.Rdata", 'phyloseq object written');
    my $usearch = $t->tool_calls('usearch');
    is_deeply([map { $_->[1] } @$usearch], ['-makeudb_usearch', '-sintax'], 'SINTAX, not -utax');
    is(canon(after($usearch->[1], '-db')), canon("$t->{db}{SSU}/fasta/refdb.120.udb"), 'SSU SINTAX database');
    is_deeply($t->tool_calls('LCA'), [], 'LCA not called');
    contains(read_text("$t->{out}/LotuSLogS/citations.txt"), 'SINTAX (USEARCH v10)');
};

case test_missing_usearch_or_database_fails_before_clustering => sub {
    my $t = setup_sintax(shift);
    for my $c ([{ no_usearch => 1 }, 93, 'no usearch binary found'],
               [{ ssu_db => "$t->{root}/missing" }, 55, 'Could not find SINTAX/UTAX database directory']) {
        my ($cfg, $code, $diagnostic) = @$c;
        subtest $diagnostic => sub {
            $t->configure(%$cfg);
            my $result = $t->run_sintax(ok => 0, out => "$t->{root}/out_$code");
            is($result->{code}, $code, "exit $code");
            contains($result->{output}, $diagnostic);
            is_deeply($t->tool_calls, [], 'no tool ran');
        };
    }
};

case test_taxonomy_only_writes_otu_first_hierarchy => sub {
    my $t = setup_sintax(shift);
    my $query = write_text("$t->{root}/query.fna", join '', map { ">$_->[0]\n$_->[1]\n" }
        ['seqA', 'ACGT'], ['seqB described sequence', 'ACGG'], ['seqC;size=3;', 'ACCT'], ['seqD', 'AGGT'], ['seqE', 'TTTT'], ['seqF', 'GGGG']);
    my $original = read_text($query);
    $t->predict({ seqA => $FIRMICUTES, seqB => 'd:Bacteria(0.9000),p:Proteobacteria(0.8500),c:Gammaproteobacteria(0.5000)',
        'seqC;size=3;' => 'd:Archaea(1.0000)', seqD => 'd:Bacteria(0.5000),p:Firmicutes(0.4000)', seqE => '' });
    my @expected = ("OTU\t$HEADER", "seqA\tBacteria\tFirmicutes\tBacilli\tLactobacillales\tStreptococcaceae\t?\t?",
        "seqB\tBacteria\tProteobacteria" . ("\t?" x 5), "seqC\tArchaea" . ("\t?" x 6),
        "seqD$UNASSIGNED", "seqE$UNASSIGNED", "seqF$UNASSIGNED");
    # -refDB SINTAX alone selects the SINTAX classifier.
    for my $c (['aligner', ['-taxAligner', 'sintax']], ['refdb', ['-refDB', 'SINTAX']]) {
        my ($name, $extra) = @$c;
        subtest $name => sub {
            my $out = "$t->{root}/taxonly_$name";
            $t->tax_only($query, $out, $extra);
            is_deeply(hiera($out), \@expected, 'hierarchy');
            contains(read_text("$out/LotuSLogS/LotuS_run.TO.log"), "1 OTU's had no entry");
            is(canon($t->tool_calls('usearch')->[-1][2]), canon($query), 'classifies the input FASTA');
        };
    }
    is(read_text($query), $original, 'input unchanged');
    is_deeply($t->tool_calls('LCA'), [], 'LCA not called');
};

case test_its_sintax_uses_only_the_sintax_database_in_strict_mode => sub {
    my $t = setup_sintax(shift);
    my $query = write_text("$t->{root}/its.fna", ">its1\nACGT\n");
    $t->predict({ its1 => 'k:Fungi(1.0000),p:Ascomycota(0.9500)' });
    my $out = "$t->{root}/taxonly_its";
    $t->tax_only($query, $out, ['-taxAligner', 'sintax', '-amplicon_type', 'ITS', '-ITSx', '0', '--strict']);
    is(hiera($out)->[1], "its1\tFungi\tAscomycota" . ("\t?" x 5), 'ITS hierarchy');
    my $manifest = read_text("$out/LotuSLogS/run_manifest.txt");
    contains($manifest, "DB[0]: $t->{db}{ITS}");
    ok(index($manifest, 'DB[1]') < 0, 'no UNITE FASTA added');
    is(canon(after($t->tool_calls('usearch')->[-1], '-db')), canon("$t->{db}{ITS}/fasta/refdb.250.udb"), 'ITS SINTAX database');
};

case test_sintax_refdb_rejects_similarity_aligners => sub {
    my $t = setup_sintax(shift);
    my $query = write_text("$t->{root}/query.fna", ">seqA\nACGT\n");
    my $result = $t->tax_only($query, "$t->{root}/taxonly_bad", ['-refDB', 'SINTAX', '-taxAligner', 'vsearch'], ok => 0);
    contains($result->{output}, 'requires -taxAligner sintax');
    is_deeply($t->tool_calls, [], 'no tool ran');
};

done_testing();

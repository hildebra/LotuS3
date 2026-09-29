#!/usr/bin/env perl
# lotus3 macOS mode (-macMode): programs that cannot run on the Mac count as not installed and
# clustering, taxonomy, mapping and phylogeny switch to available alternatives with a warning.
# Runs on any host: -macMode 1 forces the mode, and MAC_BAD (a regex) marks configured paths
# that the stubbed binary check reports as Linux executables.
use strict;
use warnings;
use FindBin;
use Test::More;
use lib "$FindBin::Bin/lib";
use LotusTest qw(read_text write_text contains lacks);

my $STUB = <<'PERL';
{ no warnings 'redefine';
  *mac_binary_problem = sub { my $p = shift // ''; return ($ENV{MAC_BAD} && $p =~ /$ENV{MAC_BAD}/) ? "is a Linux executable" : "" }; }
PERL

# The fixture's lotus3 with the binary check stubbed; MAC_BAD set to $bad.
sub mac_fixture {
    my ($t, $bad) = @_;
    my $src = read_text($t->{script});
    my $at = index($src, 'prepLtsOptions();');
    die "prepLtsOptions() call not found\n" if $at < 0;
    substr($src, $at, 0) = $STUB;
    write_text($t->{script}, $src);
    $t->{env}{MAC_BAD} = $bad // '';
}

LotusTest::case('test_linux_minimap2_switches_mapping_to_vsearch', sub {
    my ($t) = @_;
    mac_fixture($t, 'minimap2$');
    my $out = $t->run_lotus(extra => ['-macMode', '1', '--dry-run'])->{output};
    contains($out, 'minimap2 at');
    contains($out, 'is a Linux executable; treated as not installed');
    contains($out, 'mapping reads with VSEARCH instead (-useMini4map 0)');
    contains($out, 'macOS mode              active, 2 adjustment(s)');
});

LotusTest::case('test_missing_swarm_switches_clustering', sub {
    my ($t) = @_;
    mac_fixture($t);
    my $out = $t->run_lotus(extra => ['-p', 'miSeq', '-CL', 'swarm', '-macMode', '1', '--dry-run'])->{output};
    contains($out, 'SWARM clustering is not available on this Mac (brew install swarm); clustering with VSEARCH instead');
    contains($out, 'Clustering algorithm    Vsearch');
});

LotusTest::case('test_missing_lambda_switches_taxonomy', sub {
    my ($t) = @_;
    mac_fixture($t);
    my $out = $t->run_lotus(extra => ['-taxAligner', 'lambda', '-macMode', '1', '--dry-run'])->{output};
    contains($out, 'Lambda is not available on this Mac');
    contains($out, 'using VSEARCH for taxonomy instead (-taxAligner 4)');
    contains($out, 'Tax assignment          VSEARCH');
});

LotusTest::case('test_missing_mafft_disables_phylogeny', sub {
    my ($t) = @_;
    mac_fixture($t);
    my $out = $t->run_lotus(extra => ['-buildPhylo', '1', '-macMode', '1', '--dry-run'])->{output};
    contains($out, 'the phylogeny needs MAFFT (brew install mafft); no phylogenetic tree is built (-buildPhylo 0)');
    contains($out, 'phylogeny           No');
});

LotusTest::case('test_linux_sdm_stops_with_installer_hint', sub {
    my ($t) = @_;
    mac_fixture($t, '/sdm$');
    my $out = $t->run_lotus(extra => ['-macMode', '1', '--dry-run'], ok => 0)->{output};
    contains($out, 'is a Linux executable. sdm is required');
    contains($out, 'helpers/autoInstall.pl');
});

LotusTest::case('test_mode_off_and_auto_on_linux', sub {
    my ($t) = @_;
    mac_fixture($t, 'minimap2$');
    for my $mode ('0', ($^O =~ /darwin/i ? () : ('auto'))) {
        my $out = $t->run_lotus(extra => ['-macMode', $mode, '--dry-run'])->{output};
        lacks($out, 'macOS mode', "-macMode $mode leaves the tools alone");
    }
    my $out = $t->run_lotus(extra => ['-macMode', '2', '--dry-run'], ok => 0)->{output};
    contains($out, '-macMode must be auto, 0 or 1');
});

LotusTest::case('test_binary_problem_detection', sub {
    my ($t) = @_;
    my $dir = $t->{root};
    write_text("$dir/elf", "\x7fELF\x02\x01\x01\x00rest");
    write_text("$dir/intel", pack('V V', 0xfeedfacf, 0x01000007) . 'rest');
    write_text("$dir/arm", pack('V V', 0xfeedfacf, 0x0100000c) . 'rest');
    write_text("$dir/script", "#!/bin/sh\necho hi\n");
    my $out = $t->probe(<<"PERL")->{output};
\$macNoRosetta = 1;
print "\$_=[" . mac_binary_problem("$dir/\$_") . "]\\n" for qw(elf intel arm script missing);
\$macNoRosetta = 0;
print "rosetta=[" . mac_binary_problem("$dir/intel") . "]\\n";
PERL
    contains($out, 'elf=[is a Linux executable]');
    contains($out, 'intel=[is an Intel-only build and Rosetta 2 is not available');
    contains($out, 'arm=[]');
    contains($out, 'script=[]');
    contains($out, 'missing=[]');
    contains($out, 'rosetta=[]');
});

done_testing();

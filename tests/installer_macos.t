#!/usr/bin/env perl
# macOS paths in helpers/autoInstall.pl, exercised on any host by switching $isMac on inside
# a probe: programs without a macOS build come from PATH or are skipped with a warning,
# Barbell is skipped on Intel Macs without Rust, and sdm/LCA sources are fetched at a pinned commit.
use strict;
use warnings;
use Cwd qw(abs_path);
use FindBin;
use Test::More;
use lib "$FindBin::Bin/lib";
use InstallerTest qw(perl_script read_file write_file capture contains_ok lacks_ok);

my $ROOT = abs_path("$FindBin::Bin/..");
my $SOURCE = read_file("$ROOT/helpers/autoInstall.pl");

# A copy of the installer that runs $snippet (as macOS) right after its setup, then exits.
sub mac_probe {
    my ($snippet, %opt) = @_;
    my $t = InstallerTest->new($ROOT);
    my $anchor = "#usearch binary linking is handled by GetOptions above\n";
    my $source = $SOURCE;
    ok($source =~ s/\Q$anchor\E/$anchor\$isMac = 1; \$macRosetta = 0;\n$snippet\nexit(0);\n/, 'probe anchor found') or return;
    write_file($t->{installer}, $source);
    $opt{setup}->($t) if $opt{setup};
    my ($status, $output) = capture([$^X, $t->{installer}], env => $t->{env}, merge => 1, timeout => 30);
    return ($t, $status, $output);
}

subtest 'test_program_from_path_is_registered' => sub {
    my ($t, $status, $output) = mac_probe(<<'PERL', setup => sub { $_[0]->tool(mafft => perl_script("exit 0;\n")) });
my $p = mac_program_from_path("mafft", ["mafft"], "brew install mafft", "No alignments.");
print "RETURNED:$p\n";
print "WARNINGS:$finalWarning\n";
PERL
    is($status, 0, 'probe ran') or diag($output);
    contains_ok($output, "RETURNED:$t->{tools}/mafft");
    like(read_file($t->{cfg}), qr{^mafft \Q$t->{tools}\E/mafft$}m, 'PATH program written to lOTUs.cfg');
    lacks_ok($output, 'was not installed');
};

subtest 'test_missing_program_is_skipped_with_warning' => sub {
    my ($t, $status, $output) = mac_probe(<<'PERL');
my $p = mac_program_from_path("iqtree", ["iqtree2", "iqtree"], "brew install iqtree2", "No IQ-TREE trees.");
print "RETURNED:", defined($p) ? $p : "undef", "\n";
print "WARNINGS:$finalWarning\n";
PERL
    is($status, 0, 'probe ran') or diag($output);
    contains_ok($output, 'RETURNED:undef');
    contains_ok($output, 'WARNINGS:iqtree2 was not installed');
    contains_ok($output, 'brew install iqtree2');
    ok(!-e $t->{cfg} || read_file($t->{cfg}) !~ /^iqtree /m, 'no iqtree entry written');
};

subtest 'test_barbell_skipped_on_intel_mac_without_rust' => sub {
    my ($t, $status, $output) = mac_probe(<<'PERL');
$installONT = 1;
no warnings 'redefine';
*ont_architecture = sub { 'x86_64' };
my @first = ont_programs_to_install();
my @second = ont_programs_to_install();
print "PROGRAMS:@first|@second\n";
my $n = () = $finalWarning =~ /barbell was not installed/g;
print "WARNED:$n\n";
PERL
    is($status, 0, 'probe ran') or diag($output);
    contains_ok($output, 'PROGRAMS:minimap2 savont|minimap2 savont');
    contains_ok($output, 'WARNED:1');
};

subtest 'test_barbell_kept_on_apple_silicon' => sub {
    my ($t, $status, $output) = mac_probe(<<'PERL');
$installONT = 1;
no warnings 'redefine';
*ont_architecture = sub { 'aarch64' };
print "PROGRAMS:@{[ont_programs_to_install()]}\n";
PERL
    is($status, 0, 'probe ran') or diag($output);
    contains_ok($output, 'PROGRAMS:minimap2 savont barbell');
};

# Stand-in git: init creates the directory, checkout adds a Makefile, rev-parse reports $GIT_HEAD.
my $GIT = perl_script(<<'PERL');
use File::Path qw(make_path);
my @a = @ARGV;
my $dir = '';
if ($a[0] eq '-C') { $dir = $a[1]; splice @a, 0, 2 }
splice @a, 0, 2 if $a[0] eq '-c';
if ($a[0] eq 'init') { make_path($a[-1]) }
elsif ($a[0] eq 'checkout') { open my $m, '>', "$dir/Makefile" or die; print {$m} "all:\n"; close $m }
elsif ($a[0] eq 'rev-parse') { print "$ENV{GIT_HEAD}\n" }
PERL

subtest 'test_pinned_source_is_fetched' => sub {
    my ($commit) = $SOURCE =~ /sdm => \['[^']+', '([0-9a-f]{40})'\]/;
    ok(defined $commit, 'sdm commit is pinned') or return;
    my ($t, $status, $output) = mac_probe(<<'PERL', setup => sub { $_[0]->tool(git => $GIT); $_[0]{env}{GIT_HEAD} = $commit });
fetch_pinned_source("sdm", "$ldir/sdm_src");
print "MAKEFILE\n" if -f "$ldir/sdm_src/Makefile";
PERL
    is($status, 0, 'probe ran') or diag($output);
    contains_ok($output, 'MAKEFILE');
    contains_ok($output, "commit $commit");
};

subtest 'test_wrong_source_commit_is_refused' => sub {
    my ($t, $status, $output) = mac_probe(<<'PERL', setup => sub { $_[0]->tool(git => $GIT); $_[0]{env}{GIT_HEAD} = '0' x 40 });
eval { fetch_pinned_source("LCA", "$ldir/LCA_src") };
print "ERROR:$@";
print "EXISTS\n" if -e "$ldir/LCA_src";
PERL
    is($status, 0, 'probe ran') or diag($output);
    contains_ok($output, 'ERROR:Fetched LCA source is at commit ' . ('0' x 40));
    lacks_ok($output, 'EXISTS');
};

subtest 'test_existing_source_dir_is_not_replaced' => sub {
    my ($t, $status, $output) = mac_probe(<<'PERL', setup => sub { $_[0]->tool(git => $GIT); mkdir "$_[0]{install}/sdm_src"; write_file("$_[0]{install}/sdm_src/notes", "mine\n") });
eval { fetch_pinned_source("sdm", "$ldir/sdm_src") };
print "ERROR:$@";
print "KEPT\n" if -f "$ldir/sdm_src/notes";
PERL
    is($status, 0, 'probe ran') or diag($output);
    contains_ok($output, 'exists but has no Makefile');
    contains_ok($output, 'KEPT');
};

done_testing();

#!/usr/bin/env perl
# Download verification in helpers/autoInstall.pl: every download is pinned to a SHA-256,
# mismatches are never installed, bundled archives are reused, and the retired online
# updater and the --help text behave as documented. No network access is used.
use strict;
use warnings;
use Cwd qw(abs_path);
use Digest::SHA qw(sha256_hex);
use FindBin;
use Test::More;
use lib "$FindBin::Bin/lib";
use InstallerTest qw(perl_script read_file write_file capture contains_ok lacks_ok);

my $ROOT = abs_path("$FindBin::Bin/..");
my $SOURCE = read_file("$ROOT/helpers/autoInstall.pl");

subtest 'test_every_download_url_is_pinned_and_https' => sub {
    my ($table) = $SOURCE =~ /\nmy %PINNED_SHA256 = \((.*?)\n\);\n/s;
    ok(defined $table, 'pinned checksum table found') or return;
    my %pinned;
    while ($table =~ /^\s*'([^']+)'\s*=>\s*'([^']*)',/mg) {
        ok(!exists $pinned{$1}, "listed once: $1");
        $pinned{$1} = $2;
    }
    cmp_ok(scalar keys %pinned, '>', 40, 'table has the installer downloads');
    for my $url (sort keys %pinned) {
        like($url, qr{^https://}, "https: $url");
        like($pinned{$url}, qr/^[0-9a-f]{64}$/, "SHA-256 for $url");
    }
    # literal URLs passed to getS2 in live (uncommented) code
    my @literal;
    for my $line (split /\n/, $SOURCE) {
        next if $line =~ /^\s*#/;
        push @literal, $1 while $line =~ /getS2\(\s*"([^"\$]+)"\s*,/g; #the whole URL is one literal
    }
    cmp_ok(scalar @literal, '>', 20, 'found the literal download URLs');
    ok(exists $pinned{$_}, "pinned: $_") for @literal;
};

# A copy of the installer that runs $snippet right after its setup, then exits.
sub installer_probe {
    my ($snippet, %opt) = @_;
    my $t = InstallerTest->new($ROOT);
    my $anchor = "#usearch binary linking is handled by GetOptions above\n";
    my $source = $SOURCE;
    ok($source =~ s/\Q$anchor\E/$anchor$snippet\nexit(0);\n/, 'probe anchor found') or return;
    write_file($t->{installer}, $source);
    $t->tool(wget => perl_script(<<'PERL'));
my ($i) = grep { $ARGV[$_] eq '-O' } 0 .. $#ARGV;
open my $o, '>', $ARGV[$i + 1] or die "wget: $!\n"; print {$o} "downloaded payload\n"; close $o;
open my $log, '>>', "$ENV{WGET_LOG}" or die; print {$log} "$ARGV[-1]\n"; close $log;
PERL
    $t->{env}{WGET_LOG} = "$t->{root}/wget.log";
    $opt{setup}->($t) if $opt{setup};
    my ($status, $output) = capture([$^X, $t->{installer}], env => $t->{env}, merge => 1, timeout => 30);
    return ($t, $status, $output);
}

subtest 'test_unpinned_download_is_refused' => sub {
    my ($t, $status, $output) = installer_probe(<<'PERL');
eval { getS2("https://example.org/unknown.tar.gz", "$bdir/unknown.tar.gz") };
print "ERROR:$@";
print "EXISTS\n" if -e "$bdir/unknown.tar.gz";
PERL
    is($status, 0, 'probe ran') or diag($output);
    contains_ok($output, 'ERROR:No pinned SHA-256 checksum for https://example.org/unknown.tar.gz');
    lacks_ok($output, 'EXISTS');
    ok(!-e "$t->{root}/wget.log", 'nothing was downloaded');
};

subtest 'test_checksum_mismatch_is_not_installed' => sub {
    my ($t, $status, $output) = installer_probe(<<'PERL');
eval { getS2("https://example.org/pkg.tar.gz", "$bdir/pkg.tar.gz", "0" x 64) };
print "ERROR:$@";
print "EXISTS\n" if -e "$bdir/pkg.tar.gz";
print "LEFTOVER\n" if grep { /pkg\.tar\.gz\.tmp/ } glob("$bdir/*");
PERL
    is($status, 0, 'probe ran') or diag($output);
    contains_ok($output, 'failed verification; nothing was installed');
    contains_ok($output, 'Checksum mismatch');
    lacks_ok($output, 'EXISTS');
    lacks_ok($output, 'LEFTOVER');
};

subtest 'test_verified_download_and_bundled_copy' => sub {
    my $payload = "downloaded payload\n";
    my $bundle = "bundled archive\n";
    my ($t, $status, $output) = installer_probe(sprintf(<<'PERL', sha256_hex($payload), sha256_hex($bundle)),
getS2("https://example.org/fresh.tar.gz", "$bdir/fresh.tar.gz", "%s");
getS2("https://example.org/pkg.tgz", "$bdir/pkg.tgz", "%s");
print "FRESH:", -s "$bdir/fresh.tar.gz", " PKG:", -s "$bdir/pkg.tgz", "\n";
PERL
        setup => sub {
            my $t = shift;
            mkdir "$t->{install}/bin"; mkdir "$t->{install}/bin/installs";
            write_file("$t->{install}/bin/installs/pkg.tgz", $bundle);
        });
    is($status, 0, 'probe ran') or diag($output);
    contains_ok($output, 'FRESH:' . length($payload) . ' PKG:' . length($bundle));
    contains_ok($output, 'Using bundled copy');
    is(read_file("$t->{root}/wget.log"), "https://example.org/fresh.tar.gz\n", 'only the missing file was downloaded');
};

subtest 'test_help_and_retired_updater' => sub {
    my $t = InstallerTest->new($ROOT);
    my ($status, $output) = capture([$^X, $t->{installer}, '--help'], env => $t->{env}, merge => 1, timeout => 30);
    is($status, 0, '--help exits 0');
    contains_ok($output, $_) for '--ont-only', '--no-telemetry', '-lambdaIndex', 'pinned SHA-256';
    ($status, $output) = capture([$^X, $t->{installer}, '-forceUpdate'], env => $t->{env}, merge => 1, timeout => 30);
    isnt($status, 0, '-forceUpdate fails');
    contains_ok($output, '-forceUpdate is no longer supported');
    ok(!-e $t->{cfg}, 'no config written');
};

done_testing();

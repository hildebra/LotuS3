#!/usr/bin/env perl
# Extract one executable from a checksum-verified .conda package, without Conda
# or Python. Uses Perl's ZIP/tar readers and the zstd command-line decompressor.
# The installer verifies the package checksum before invoking this helper.
use strict;
use warnings;
use Archive::Tar;
use Cwd qw(abs_path);
use Fcntl qw(O_WRONLY O_CREAT O_EXCL);
use File::Spec;
use File::Temp qw(tempdir);
use IO::Uncompress::Unzip qw($UnzipError);

use constant MAX_ARCHIVE_BYTES => 256 * 1024 * 1024;

sub find_zstd {
    for my $dir (File::Spec->path()) {
        my $path = File::Spec->catfile($dir, 'zstd');
        return abs_path($path) if -f $path && -x $path;
    }
    die "Zstandard reader unavailable: install the zstd command-line tool.\n";
}

sub unpack_payload {
    my ($package, $destination) = @_;
    die "Missing or oversized .conda package\n"
        unless -f $package && -s $package && -s $package <= MAX_ARCHIVE_BYTES;
    my $zip = IO::Uncompress::Unzip->new($package, Strict => 1, Transparent => 0, MultiStream => 0)
        or die "Cannot read .conda ZIP archive: $UnzipError\n";
    my ($payloads, $total, $status) = (0, 0, 1);
    while ($status > 0) {
        my $name = $zip->getHeaderInfo()->{Name};
        my $wanted = defined($name) && $name =~ /\Apkg-[^\/\\]+\.tar\.zst\z/;
        my $out;
        if ($wanted) {
            die "Expected exactly one pkg-*.tar.zst payload\n" if ++$payloads > 1;
            open $out, '>:raw', $destination or die "Cannot write payload: $!\n";
        }
        while (($status = $zip->read(my $buffer, 65536)) > 0) {
            $total += $status;
            die "Oversized ZIP contents\n" if $total > MAX_ARCHIVE_BYTES;
            print {$out} $buffer or die "Cannot write payload: $!\n" if $wanted;
        }
        die "Invalid ZIP payload: $UnzipError\n" if $status < 0;
        close $out or die "Cannot close payload: $!\n" if $wanted;
        $status = $zip->nextStream();
    }
    die "Invalid ZIP archive: $UnzipError\n" if $status < 0;
    $zip->close() or die "Cannot close ZIP archive: $UnzipError\n";
    die "Expected one nonempty pkg-*.tar.zst payload\n" unless $payloads == 1 && -s $destination;
}

sub decompress_payload {
    my ($zstd, $payload, $destination) = @_;
    open my $pipe, '-|', $zstd, '-d', '-q', '-c', '-M256MB', '--', $payload
        or die "Cannot run zstd: $!\n";
    binmode $pipe;
    open my $out, '>:raw', $destination or die "Cannot write tar archive: $!\n";
    my $total = 0;
    while (1) {
        my $n = read($pipe, my $buffer, 65536);
        die "Cannot read zstd output: $!\n" unless defined($n);
        last unless $n;
        $total += $n;
        die "Oversized decompressed package\n" if $total > MAX_ARCHIVE_BYTES;
        print {$out} $buffer or die "Cannot write tar archive: $!\n";
    }
    close $out or die "Cannot close tar archive: $!\n";
    close $pipe or die "Zstandard decompression failed (process status $?)\n";
    die "Empty decompressed package\n" unless $total;
}

sub extract_executable {
    my ($package, $member, $destination) = @_;
    my $zstd = find_zstd();
    my $stage = tempdir('lotus-conda-XXXXXXXX', TMPDIR => 1, CLEANUP => 1);
    unpack_payload($package, "$stage/payload.tar.zst");
    decompress_payload($zstd, "$stage/payload.tar.zst", "$stage/payload.tar");
    my $tar = Archive::Tar->new();
    $tar->read("$stage/payload.tar", 0, { filter_cb => sub {
        die "Oversized tar member\n" if $_[0]->size > MAX_ARCHIVE_BYTES;
        return 1;
    } }) or die "Cannot read tar archive: " . $tar->error . "\n";
    die "Invalid tar archive: " . $tar->error . "\n" if $tar->error;
    my @matches = grep { $_->full_path eq $member } $tar->get_files();
    die "Expected one nonempty regular executable: $member\n"
        unless @matches == 1 && $matches[0]->is_file && $matches[0]->type =~ /\A0?\z/
            && $matches[0]->size > 0;
    my $content = $matches[0]->get_content_by_ref();
    die "Incomplete executable: $member\n" unless length($$content) == $matches[0]->size;

    # Write only the requested regular file. Never extract package paths/links.
    sysopen my $out, $destination, O_WRONLY | O_CREAT | O_EXCL, 0600
        or die "Cannot create $destination (must not already exist): $!\n";
    my $ok = eval {
        binmode $out;
        print {$out} $$content or die "Cannot write $destination: $!\n";
        close $out or die "Cannot close $destination: $!\n";
        chmod 0755, $destination or die "Cannot make $destination executable: $!\n";
        1;
    };
    if (!$ok) {
        my $error = $@;
        close $out if defined(fileno($out));
        unlink $destination;
        die $error;
    }
}

my $ok = eval {
    if (@ARGV == 1 && $ARGV[0] eq '--check') {
        my $zstd = find_zstd();
        system($zstd, '--version') == 0 or die "Cannot execute zstd (process status $?)\n";
    } elsif (@ARGV == 3) {
        extract_executable(@ARGV);
    } else {
        die "Usage: extract_conda_executable.pl PACKAGE MEMBER DESTINATION | --check\n";
    }
    1;
};
if (!$ok) {
    print STDERR "Cannot extract Bioconda executable: $@";
    exit(1);
}

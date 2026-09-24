#!/usr/bin/env perl
# Exercise the installer CLI in temporary installations, without network/builds.
# Download tests substitute fixture digests in a temporary script copy, leaving
# checksum verification, unpacking, executable probing, and config writes intact.
use strict;
use warnings;
use Archive::Tar;
use Archive::Tar::Constant qw(HARDLINK SYMLINK);
use Cwd qw(abs_path);
use Digest::SHA qw(sha256_hex);
use File::Basename qw(dirname);
use File::Temp qw(tempdir);
use FindBin;
use IO::Compress::Bzip2 qw(bzip2 $Bzip2Error);
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Compress::Zip qw(:zip_method $ZipError);
use JSON::PP qw(decode_json);
use Test::More;
use lib "$FindBin::Bin/lib";
use InstallerTest qw(program perl_script write_executable which read_file write_file replaced
    capture contains_ok lacks_ok);

my $ROOT = abs_path("$FindBin::Bin/..");
use constant ZSTD_SKIP => 'zstd is needed for extraction tests; set LOTUS_TEST_ZSTD or add it to PATH';

sub zstd_tool {
    my $path = $ENV{LOTUS_TEST_ZSTD} || which('zstd');
    skip(ZSTD_SKIP, 1) unless defined $path && -f $path && -x $path;
    return abs_path($path);
}

# Tar archive of [name, data, \%properties] members, padded to 10240-byte records
# like tar(1); the decompressed size-limit test relies on that padding.
sub tar_archive {
    my $tar = Archive::Tar->new;
    $tar->add_data(@$_) or die "Cannot add $_->[0]: " . $tar->error . "\n" for @_;
    my $bytes = $tar->write;
    return $bytes . "\0" x (-length($bytes) % 10240);
}

# Uncompressed ZIP of (name => data) pairs, in order.
sub zip_members {
    my (@members) = @_;
    my ($zip, $bytes);
    while (my ($name, $data) = splice @members, 0, 2) {
        my %entry = (Name => $name, Method => ZIP_CM_STORE, Stream => 0);
        if ($zip) { $zip->newStream(%entry) or die "Cannot add $name: $ZipError\n" }
        else { $zip = IO::Compress::Zip->new(\$bytes, %entry) or die "Cannot create ZIP: $ZipError\n" }
        $zip->print($data);
    }
    $zip->close or die "Cannot finish ZIP: $ZipError\n";
    return $bytes;
}

# Members of a small real ZIP/Zstandard/tar .conda package; no Conda or Rust needed.
# Options: contents, member (default bin/savont), link (1 or 'hard'), duplicate.
# The payload is compressed with the zstd tool, so fixture creation skips like extraction.
sub conda_members {
    my (%opt) = @_;
    my $member = $opt{member} // 'bin/savont';
    my $entry = $opt{link}
        ? [$member, '', { type => $opt{link} eq 'hard' ? HARDLINK : SYMLINK, linkname => '/bin/sh' }]
        : [$member, $opt{contents} // program('savont')];
    my $dir = tempdir(CLEANUP => 1);
    write_file("$dir/payload.tar", tar_archive(($entry) x ($opt{duplicate} ? 2 : 1)));
    my ($status, $compressed) = capture([zstd_tool(), '-q', '-3', '-c', '--', "$dir/payload.tar"]);
    die "Unable to compress test fixture\n" if $status;
    return ('metadata.json' => '{"conda_pkg_format_version": 2}', 'pkg-savont-fixture.tar.zst' => $compressed);
}

sub conda_package { return zip_members(conda_members(@_)) }

sub setup_downloads {
    my ($t, %opt) = @_;
    my $fixtures = "$t->{root}/fixtures";
    mkdir $fixtures or die "Cannot create $fixtures: $!\n";
    my $savont = tar_archive(['savont-0.7.0/Cargo.toml', qq{[package]\nname="savont"\n}, { mode => 0755 }],
                             ['savont-0.7.0/Cargo.lock', "version = 3\n", { mode => 0755 }]);
    my $minimap2 = tar_archive(['minimap2-2.28_x64-linux/minimap2', program('minimap2'), { mode => 0755 }]);
    my %data = (barbell => program('barbell'));
    gzip(\$savont => \$data{savont}) or die "Cannot gzip fixture: $GzipError\n";
    bzip2(\$minimap2 => \$data{minimap2}) or die "Cannot bzip2 fixture: $Bzip2Error\n";
    write_file("$fixtures/$_", $data{$_}) for keys %data;
    my $source = read_file($t->{installer});
    for (['savont', 'a60141fd4d4e83cdcbc3601220dc41487f4cdf7bd81558c7125e876117ed756c'],
         ['barbell', 'f23a599eceb8b27211178facf5504935a98880c031eb9f2b4f136382193f2081'],
         ['minimap2', '51f2cf0e486d0f9f88ace1aa58fdc56571382a676ea0889ae607301c60693377']) {
        $source = replaced($source, $_->[1], sha256_hex($data{$_->[0]}));
    }
    write_file($t->{installer}, $source);
    write_file("$fixtures/barbell", 'corrupt download') if $opt{corrupt};
    $t->{env}{ONT_FIXTURES} = $fixtures;
    $t->{env}{ONT_BUILD_LOG} = "$t->{root}/build.json";
    $t->tool(wget => perl_script(<<'PERL'));
use File::Copy qw(copy);
my $url = $ARGV[-1];
my $name = index($url, '/savont/') >= 0 ? 'savont' : index($url, '/barbell/') >= 0 ? 'barbell' : 'minimap2';
my ($i) = grep { $ARGV[$_] eq '-O' } 0 .. $#ARGV;
defined $i or die "wget: missing -O\n";
copy("$ENV{ONT_FIXTURES}/$name", $ARGV[$i + 1]) or die "wget: cannot copy $name fixture: $!\n";
PERL
    $t->tool(rustc => perl_script(<<'PERL'));
print "rustc 1.88.0\n";
PERL
    $t->tool($_ => "#!/bin/sh\nexit 0\n") for qw(cc c++ cmake);
    # Logs its arguments as a JSON list and "builds" savont from the __DATA__ section.
    $t->tool(cargo => perl_script(<<'PERL') . program('savont'));
use File::Path qw(make_path);
use JSON::PP;
my @args = @ARGV;
my $json = JSON::PP->new->allow_nonref->ascii;
open my $log, '>', $ENV{ONT_BUILD_LOG} or die "cargo: cannot write build log: $!\n";
print {$log} '[' . join(', ', map { $json->encode($_) } @args) . ']';
close $log or die "cargo: cannot write build log: $!\n";
my ($i) = grep { $args[$_] eq '--target-dir' } 0 .. $#args;
defined $i or die "cargo: missing --target-dir\n";
my $dir = "$args[$i + 1]/release";
make_path($dir);
open my $out, '>', "$dir/savont" or die "cargo: cannot write $dir/savont: $!\n";
print {$out} do { local $/; <DATA> };
close $out or die "cargo: cannot write $dir/savont: $!\n";
chmod 0755, "$dir/savont" or die "cargo: cannot chmod $dir/savont: $!\n";
__DATA__
PERL
}

sub setup_bioconda {
    my ($t, $payload) = @_;
    setup_downloads($t);
    unlink "$t->{tools}/$_" or die "Cannot remove $_: $!\n" for qw(cargo rustc cc c++ cmake);
    symlink(zstd_tool(), "$t->{tools}/zstd") or die "Cannot link zstd: $!\n";
    $payload //= conda_package();
    write_file("$t->{root}/fixtures/savont", $payload);
    write_file($t->{installer}, replaced(read_file($t->{installer}),
        'e7ea28b084d176379d9fa273a3cb349c9e58e54436e2efd02c295acf91edbdd7', sha256_hex($payload)));
    write_file($t->{cfg}, read_file($t->{default}) . "savont previous-install\n");
}

sub assert_bioconda_failure {
    my ($t, $output) = @_;
    local $Test::Builder::Level = $Test::Builder::Level + 1;
    contains_ok($output, 'Savont Bioconda binary fallback failed');
    contains_ok($output, 'Install Rust >= 1.88 (including Cargo)');
    is($t->config->{savont}, 'previous-install', 'savont entry preserved');
    ok(!-e "$t->{install}/bin/savont", 'no savont installed');
    ok(!-e "$t->{install}/bin/barbell", 'no barbell installed');
    ok(!-e "$t->{root}/build.json", 'nothing built');
}

sub fixture { return InstallerTest->new($ROOT) }

# --- ONTInstaller ---------------------------------------------------------------

subtest 'test_new_install_registers_absolute_paths' => sub {
    my $t = fixture();
    $t->installed_tools; $t->run_installer;
    my $config = $t->config;
    is($config->{$_}, "$t->{tools}/$_", "$_ registered") for qw(minimap2 savont barbell);
    is($config->{TAX_REFDB_KSGP}, '/my/reference.fasta', 'other entries kept');
    ok(!-e "$t->{install}/DB", 'no DB directory');
};

subtest 'test_default_install_downloads_barbell_without_rust' => sub {
    my $t = fixture();
    setup_downloads($t);
    $t->tool('minimap2'); $t->tool('savont');
    unlink "$t->{tools}/$_" or die "Cannot remove $_: $!\n" for qw(cargo rustc cc c++ cmake);
    $t->run_installer;
    is($t->config->{barbell}, "$t->{install}/bin/barbell", 'barbell registered');
    ok(-x "$t->{install}/bin/barbell", 'barbell executable');
    ok(!-e "$t->{root}/build.json", 'nothing built');
};

subtest 'test_default_install_reuses_configured_barbell' => sub {
    my $t = fixture();
    $t->installed_tools; $t->tool(barbell => program('barbell', crash => 1));
    my $configured = write_executable("$t->{install}/custom-barbell", program('barbell'));
    write_file($t->{cfg}, read_file($t->{default}) . "barbell custom-barbell\n");
    my $output = $t->run_installer;
    is($t->config->{barbell}, $configured, 'configured barbell kept');
    lacks_ok($output, "$t->{tools}/barbell");
};

subtest 'test_old_config_duplicates_backup_and_rerun' => sub {
    my $t = fixture();
    $t->installed_tools;
    my $original = read_file($t->{default}) . "savont missing\nsavont duplicate\nbarbell missing\n# no final newline";
    write_file($t->{cfg}, $original);
    $t->run_installer;
    is(read_file("$t->{cfg}.bak"), $original, 'backup holds the original config');
    my $saved = read_file($t->{cfg});
    $t->run_installer;
    is(read_file($t->{cfg}), $saved, 'rerun leaves config unchanged');
    for my $name (qw(savont barbell minimap2)) {
        is(scalar(grep { $_->[0] eq $name } $t->entries), 1, "one $name entry");
    }
    contains_ok($saved, "# no final newline\n", 'unterminated last line is kept and terminated');
};

subtest 'test_configured_relative_paths_are_reused' => sub {
    my $t = fixture();
    $t->installed_tools; mkdir "$t->{install}/custom" or die "Cannot create custom: $!\n";
    my $savont = write_executable("$t->{install}/custom/my-savont", program('savont'));
    write_file($t->{cfg}, read_file($t->{default}) . "savont custom/my-savont\n");
    $t->run_installer;
    is($t->config->{savont}, $savont, 'relative savont entry resolved');
};

subtest 'test_missing_build_tools_preserves_config' => sub {
    my $t = fixture();
    $t->tool('minimap2'); $t->tool('barbell');
    $t->tool(rustc => qq{#!/bin/sh\necho "rustc 1.88.0"\n});
    $t->tool('cargo');
    write_file($t->{cfg}, read_file($t->{default})); my $before = read_file($t->{cfg});
    my $output = $t->run_installer(ok => 0);
    contains_ok($output, 'ONT source builds require');
    is(read_file($t->{cfg}), $before, 'config unchanged');
    ok(!-e "$t->{cfg}.bak", 'no backup written');
};

subtest 'test_crashing_program_is_not_registered' => sub {
    my $t = fixture();
    $t->installed_tools; $t->tool(savont => program('savont', crash => 1));
    my $output = $t->run_installer(ok => 0);
    contains_ok($output, 'install the zstd command-line tool');
    contains_ok($output, 'Install Rust >= 1.88');
    ok(!-e $t->{cfg}, 'no config written');
};

subtest 'test_rejects_conflicting_modes' => sub {
    my $t = fixture();
    $t->installed_tools;
    contains_ok($t->run_installer(extra => ['--condaDBinstall'], ok => 0), 'cannot be combined');
    ok(!-e $t->{cfg}, 'no config written');
};

subtest 'test_download_build_verify_and_register' => sub {
    my $t = fixture();
    setup_downloads($t); $t->run_installer;
    my $config = $t->config;
    for my $name (qw(minimap2 savont barbell)) {
        is($config->{$name}, "$t->{install}/bin/$name", "$name registered");
        ok(defined $config->{$name} && -x $config->{$name}, "$name executable");
    }
    my $args = decode_json(read_file("$t->{root}/build.json"));
    for my $flag (qw(--locked --release)) {
        ok(scalar(grep { $_ eq $flag } @$args), "cargo called with $flag");
    }
    opendir my $bin, "$t->{install}/bin" or die "Cannot list bin: $!\n";
    is_deeply([sort grep { !/\A\.\.?\z/ } readdir $bin], [qw(barbell minimap2 savont)], 'bin contents');
};

subtest 'test_checksum_failure_keeps_old_tool_entry' => sub {
    my $t = fixture();
    setup_downloads($t, corrupt => 1);
    write_file($t->{cfg}, read_file($t->{default}) . "barbell previous-install\n");
    my $output = $t->run_installer(ok => 0);
    contains_ok($output, 'Checksum mismatch');
    is($t->config->{barbell}, 'previous-install', 'barbell entry preserved');
    ok(!-e "$t->{install}/bin/barbell", 'no barbell installed');
};

subtest 'test_bioconda_without_rust_installs_and_reruns' => sub { SKIP: {
    my $t = fixture();
    setup_bioconda($t);
    is(which($_, $t->{env}{PATH}), undef, "no $_ on PATH") for qw(python python3 conda cargo rustc);
    $t->run_installer;
    is($t->config->{$_}, "$t->{install}/bin/$_", "$_ registered") for qw(minimap2 savont barbell);
    ok(!-e "$t->{root}/build.json", 'nothing built');
    my $saved = read_file($t->{cfg});
    unlink "$t->{tools}/$_" or die "Cannot remove $_: $!\n" for qw(wget zstd);
    $t->run_installer;
    is(read_file($t->{cfg}), $saved, 'rerun leaves config unchanged');
} };

subtest 'test_bioconda_with_old_rust' => sub { SKIP: {
    my $t = fixture();
    setup_bioconda($t);
    $t->tool(cargo => "#!/bin/sh\nexit 1\n");
    $t->tool(rustc => qq{#!/bin/sh\necho "rustc 1.87.0"\n});
    $t->run_installer;
    is($t->config->{savont}, "$t->{install}/bin/savont", 'savont registered');
} };

subtest 'test_bioconda_without_cargo' => sub { SKIP: {
    my $t = fixture();
    setup_bioconda($t);
    $t->tool(rustc => qq{#!/bin/sh\necho "rustc 1.88.0"\n});
    $t->run_installer;
    is($t->config->{savont}, "$t->{install}/bin/savont", 'savont registered');
} };

subtest 'test_bioconda_corrupt_download_aborts' => sub { SKIP: {
    my $t = fixture();
    setup_bioconda($t);
    write_file("$t->{root}/fixtures/savont", 'corrupt download');
    my $output = $t->run_installer(ok => 0);
    assert_bioconda_failure($t, $output);
    contains_ok($output, 'Checksum mismatch');
} };

subtest 'test_bioconda_failed_download_aborts' => sub { SKIP: {
    my $t = fixture();
    setup_bioconda($t);
    $t->tool('minimap2');
    $t->tool(wget => "#!/bin/sh\nexit 1\n");
    assert_bioconda_failure($t, $t->run_installer(ok => 0));
} };

subtest 'test_bioconda_bad_archive_aborts' => sub { SKIP: {
    my $t = fixture();
    setup_bioconda($t, 'not a ZIP archive');
    my $output = $t->run_installer(ok => 0);
    assert_bioconda_failure($t, $output);
    contains_ok($output, 'Cannot extract Bioconda executable');
} };

subtest 'test_bioconda_crashing_executable_aborts' => sub { SKIP: {
    my $t = fixture();
    setup_bioconda($t, conda_package(contents => program('savont', crash => 1)));
    assert_bioconda_failure($t, $t->run_installer(ok => 0));
} };

subtest 'test_bioconda_missing_cli_flags_aborts' => sub { SKIP: {
    my $t = fixture();
    setup_bioconda($t, conda_package(contents => qq{#!/bin/sh\necho "savont 0.7.0"\n}));
    assert_bioconda_failure($t, $t->run_installer(ok => 0));
} };

subtest 'test_bioconda_missing_zstd_aborts_before_config_update' => sub { SKIP: {
    my $t = fixture();
    setup_bioconda($t); unlink "$t->{tools}/zstd" or die "Cannot remove zstd: $!\n";
    my $before = read_file($t->{cfg});
    assert_bioconda_failure($t, $t->run_installer(ok => 0));
    is(read_file($t->{cfg}), $before, 'config unchanged');
} };

subtest 'test_bioconda_unsupported_platform_aborts' => sub { SKIP: {
    my $t = fixture();
    setup_bioconda($t);
    write_file($t->{installer}, replaced(read_file($t->{installer}),
        'my $arch = lc($host[4]);', q{my $arch = 'unsupported';}));
    my $output = $t->run_installer(ok => 0);
    assert_bioconda_failure($t, $output);
    contains_ok($output, 'No pinned Savont Bioconda binary');
} };

# --- CondaExtraction ------------------------------------------------------------

sub extract {
    my ($package, $destination, $script) = @_;
    my %env = (%ENV, PATH => dirname(zstd_tool()));
    my ($status, $stdout, $stderr) = capture([$^X, $script // "$ROOT/helpers/extract_conda_executable.pl",
        $package, 'bin/savont', $destination], env => \%env);
    return { status => $status, stdout => $stdout, stderr => $stderr };
}

sub extraction_fails {
    my ($result) = @_;
    local $Test::Builder::Level = $Test::Builder::Level + 1;
    isnt($result->{status}, 0, 'extraction fails') or diag($result->{stdout});
}

subtest 'test_rejects_missing_or_linked_executable' => sub {
    for my $kwargs ([member => 'bin/other'], [link => 1], [link => 'hard'], [contents => ''], [duplicate => 1]) {
        subtest "kwargs: $kwargs->[0] => '$kwargs->[1]'" => sub { SKIP: {
            my $temp = tempdir(CLEANUP => 1);
            my ($package, $destination) = ("$temp/savont.conda", "$temp/savont");
            write_file($package, conda_package(@$kwargs));
            extraction_fails(extract($package, $destination));
            ok(!-e $destination, 'no executable written');
        } };
    }
};

subtest 'test_refuses_to_overwrite_existing_file' => sub {
    for my $linked (0, 1) {
        subtest "linked=$linked" => sub { SKIP: {
            my $temp = tempdir(CLEANUP => 1);
            my ($package, $destination) = ("$temp/savont.conda", "$temp/savont");
            write_file($package, conda_package());
            my $target = "$temp/original"; write_file($target, 'preserve');
            if ($linked) { symlink($target, $destination) or die "Cannot link destination: $!\n" }
            else { write_file($destination, 'preserve') }
            extraction_fails(extract($package, $destination));
            is(read_file($destination), 'preserve', 'destination preserved');
            is(read_file($target), 'preserve', 'link target preserved');
        } };
    }
};

subtest 'test_rejects_duplicate_zip_payloads' => sub { SKIP: {
    my $temp = tempdir(CLEANUP => 1);
    my ($package, $destination) = ("$temp/savont.conda", "$temp/savont");
    # The regular package plus a second copy of its payload member.
    my @members = conda_members();
    my %member = @members;
    write_file($package, zip_members(@members, 'pkg-second.tar.zst' => $member{'pkg-savont-fixture.tar.zst'}));
    my $result = extract($package, $destination);
    extraction_fails($result);
    contains_ok($result->{stderr}, 'exactly one pkg-');
    ok(!-e $destination, 'no executable written');
} };

subtest 'test_zstd_failure_leaves_no_executable' => sub { SKIP: {
    my $temp = tempdir(CLEANUP => 1);
    my ($package, $destination) = ("$temp/savont.conda", "$temp/savont");
    write_file($package, zip_members('pkg-broken.tar.zst' => 'not a zstd frame'));
    my $result = extract($package, $destination);
    extraction_fails($result);
    contains_ok($result->{stderr}, 'Zstandard decompression failed');
    ok(!-e $destination, 'no executable written');
} };

subtest 'test_decompressed_size_limit' => sub { SKIP: {
    my $temp = tempdir(CLEANUP => 1);
    my ($package, $destination) = ("$temp/savont.conda", "$temp/savont");
    write_file($package, conda_package());
    my $script = "$temp/extract.pl";
    write_file($script, replaced(read_file("$ROOT/helpers/extract_conda_executable.pl"), '256 * 1024 * 1024', '4096'));
    my $result = extract($package, $destination, $script);
    extraction_fails($result);
    contains_ok($result->{stderr}, 'Oversized decompressed package');
    ok(!-e $destination, 'no executable written');
} };

done_testing();

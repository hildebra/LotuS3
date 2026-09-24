package InstallerTest;
# Temporary LotuS3 installation for the installer tests: a copy of the installer
# scripts, a default config and a tools directory that is the only PATH entry, so
# every program the installer sees is a stand-in written by the test.
use strict;
use warnings;
use Exporter qw(import);
use File::Copy qw(copy);
use File::Temp qw(tempdir);
use POSIX ();
use Test::More;

our @EXPORT_OK = qw(program perl_script write_executable which read_file write_file replaced
    capture contains_ok lacks_ok);

sub read_file {
    my ($path) = @_;
    open my $fh, '<:raw', $path or die "Cannot read $path: $!\n";
    local $/;
    my $data = <$fh>;
    return $data // '';
}

sub write_file {
    my ($path, $data) = @_;
    open my $fh, '>:raw', $path or die "Cannot write $path: $!\n";
    print {$fh} $data or die "Cannot write $path: $!\n";
    close $fh or die "Cannot close $path: $!\n";
}

sub write_executable {
    my ($path, $contents) = @_;
    write_file($path, $contents);
    chmod 0755, $path or die "Cannot chmod $path: $!\n";
    return $path;
}

# $text with every literal occurrence of $old replaced by $new.
sub replaced {
    my ($text, $old, $new) = @_;
    $text =~ s/\Q$old\E/$new/g;
    return $text;
}

sub which {
    my ($name, $path) = @_;
    for my $dir (grep { length } split /:/, $path // $ENV{PATH} // '') {
        return "$dir/$name" if -f "$dir/$name" && -x _;
    }
    return;
}

# Stand-ins run with the test's perl by absolute path: PATH holds only the tools dir.
sub perl_script { return "#!$^X\n" . shift }

# Stand-in ONT tool: reports a version for --version and otherwise prints every CLI
# flag the installer probes for. With crash => 1 it kills itself with SIGTERM.
sub program {
    my ($name, %opt) = @_;
    return perl_script("kill 'TERM', \$\$;\n") if $opt{crash};
    my $version = $name eq 'minimap2' ? '2.28-r1209' : "$name 0.7.0";
    (my $body = <<'PERL') =~ s/VERSION/$version/;
if (grep { $_ eq '--version' } @ARGV) {
    print "VERSION\n";
} else {
    print "--kit --input --output --maximize --threads --quality-value-cutoff --minimum-base-quality --chimera-allowable-errors --single-strand\n";
}
PERL
    return perl_script($body);
}

# Run @$cmd with optional env => \%env, input => $text (STDIN is inherited when
# undef), merge => 1 (STDERR into STDOUT) and timeout => $seconds.
# Returns ($wait_status, $stdout, $stderr).
sub capture {
    my ($cmd, %opt) = @_;
    my $dir = tempdir(CLEANUP => 1);
    write_file("$dir/stdin", $opt{input}) if defined $opt{input};
    my $pid = fork() // die "Cannot fork: $!\n";
    if (!$pid) {    # child: never return into the test code
        eval {
            %ENV = %{$opt{env}} if $opt{env};
            if (defined $opt{input}) { open STDIN, '<', "$dir/stdin" or die "Cannot redirect STDIN: $!\n" }
            open STDOUT, '>', "$dir/stdout" or die "Cannot redirect STDOUT: $!\n";
            if ($opt{merge}) { open STDERR, '>&', \*STDOUT or die "Cannot redirect STDERR: $!\n" }
            else { open STDERR, '>', "$dir/stderr" or die "Cannot redirect STDERR: $!\n" }
            exec { $cmd->[0] } @$cmd or die "Cannot run $cmd->[0]: $!\n";
        };
        print STDERR $@;
        POSIX::_exit(127);
    }
    my $status;
    my $finished = eval {
        local $SIG{ALRM} = sub { die "timeout\n" };
        alarm($opt{timeout} // 0);
        waitpid($pid, 0);
        $status = $?;
        alarm(0);
        1;
    };
    if (!$finished) {
        kill 'KILL', $pid;
        waitpid($pid, 0);
        die "Timed out after $opt{timeout} seconds: @$cmd\n" . read_file("$dir/stdout");
    }
    return ($status, read_file("$dir/stdout"), $opt{merge} ? '' : read_file("$dir/stderr"));
}

sub contains_ok {
    my ($text, $expected, $name) = @_;
    local $Test::Builder::Level = $Test::Builder::Level + 1;
    return like($text, qr/\Q$expected\E/, $name // "output contains '$expected'");
}

sub lacks_ok {
    my ($text, $unexpected, $name) = @_;
    local $Test::Builder::Level = $Test::Builder::Level + 1;
    return unlike($text, qr/\Q$unexpected\E/, $name // "output lacks '$unexpected'");
}

sub new {
    my ($class, $repo) = @_;
    my $root = tempdir('lotus-ont-installer-test-XXXXXXXX', TMPDIR => 1, CLEANUP => 1);
    my $install = "$root/lotus";
    my $self = bless {
        root => $root, install => $install, tools => "$root/tools",
        installer => "$install/helpers/autoInstall.pl",
        default => "$install/configs/LotuS.cfg.def", cfg => "$install/lOTUs.cfg",
    }, $class;
    for my $dir ($install, "$install/helpers", "$install/configs", $self->{tools}) {
        mkdir $dir or die "Cannot create $dir: $!\n";
    }
    write_file("$install/lotus3", "# test installation root\n");
    for my $helper (qw(autoInstall.pl extract_conda_executable.pl)) {
        copy("$repo/helpers/$helper", "$install/helpers/$helper") or die "Cannot copy $helper: $!\n";
    }
    write_file($self->{default}, "UID ??\nusearch unavailable\n# preserve this\nTAX_REFDB_KSGP /my/reference.fasta\n");
    for my $name (qw(tar chmod gzip bzip2)) {
        my $path = which($name) // die "$name is needed on PATH\n";
        symlink($path, "$self->{tools}/$name") or die "Cannot link $name: $!\n";
    }
    $self->{env} = { %ENV, PATH => $self->{tools} };
    return $self;
}

sub tool {
    my ($self, $name, $contents) = @_;
    return write_executable("$self->{tools}/$name", $contents || program($name));
}

sub installed_tools {
    my ($self) = @_;
    $self->tool($_) for qw(minimap2 savont barbell);
}

# Options: extra => [...], ok => 0 (expect failure), ont_only => 0, answers => $stdin.
sub run_installer {
    my ($self, %opt) = @_;
    my @mode = ($opt{ont_only} // 1) ? ('--ont-only') : ();
    my ($status, $output) = capture([$^X, $self->{installer}, @mode, @{$opt{extra} || []}],
        env => $self->{env}, input => $opt{answers}, merge => 1, timeout => 30);
    local $Test::Builder::Level = $Test::Builder::Level + 1;
    if ($opt{ok} // 1) { is($status, 0, 'installer exits with status 0') or diag($output) }
    else { isnt($status, 0, 'installer exits with an error') or diag($output) }
    return $output;
}

# Config lines as [key, value] pairs, skipping blank lines and comments.
sub entries {
    my ($self) = @_;
    return map { [split ' ', $_, 2] } grep { $_ ne '' && !/^#/ } split /\r?\n/, read_file($self->{cfg});
}

sub config {
    my ($self) = @_;
    return +{ map { @$_ == 2 ? @$_ : die "Malformed config entry: '@$_'\n" } $self->entries };
}

1;

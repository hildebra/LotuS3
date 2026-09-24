package DerepAudit;
# Reconcile finalized SDM map, main, merged and rest counts (derepPerSR=0).
#
# Perl port of hildebra/sdm tests/audit_derep_counts.py, 2026-09-15 integration
# contract. Run on fresh SDM output before a consumer appends/rewrites records.
# audit() dies with a message on any inconsistency. Returned hashes keep
# insertion order, which fixes the key order of the CLI's JSON output.
use strict;
use warnings;
use Exporter 'import';

our @EXPORT_OK = qw(audit records pure_path py_int);

my $WS = qr/[\s\x1c-\x1f]/;    # whitespace, including the \x1c-\x1f separators
sub rstrip { my ($s) = @_; $s =~ s/$WS+\z//; return $s }
sub strip  { my ($s) = @_; $s =~ s/\A$WS+//; $s =~ s/$WS+\z//; return $s }

# Integer: optional sign, digits with single underscores, surrounding whitespace.
sub py_int {
    my ($s) = @_;
    die "invalid integer: '$s'\n"
        unless defined $s && $s =~ /\A$WS*([+-]?[0-9]+(?:_[0-9]+)*)$WS*\z/;
    (my $n = $1) =~ tr/_//d;
    return 0 + $n;
}

# Normalised POSIX path: repeated slashes and "." parts dropped, "//" root kept.
sub pure_path {
    my ($p) = @_;
    my $root = $p =~ m{\A//(?!/)} ? '//' : $p =~ m{\A/} ? '/' : '';
    my $s = $root . join '/', grep { $_ ne '' && $_ ne '.' } split m{/}, $p;
    return $s eq '' ? '.' : $s;
}

# (path without its last suffix, suffix); leading-dot names have no suffix.
sub split_suffix {
    my ($p) = @_;
    my ($dir, $name) = $p =~ m{\A(.*/)?([^/]*)\z}s;
    my $i = rindex($name, '.');
    return ($p, '') unless $i > 0 && $i < length($name) - 1;
    return (($dir // '') . substr($name, 0, $i), substr($name, $i));
}

sub odict { tie my %h, 'DerepAudit::Ordered'; %h = @_; return \%h }

# Stream FASTA (including wrapped sequences) or SDM four-line FASTQ.
# Returns an iterator yielding (name, length) until empty.
sub records {
    my ($path) = @_;
    die "Is a directory: $path\n" if -d $path;
    open my $fh, '<:encoding(UTF-8)', $path or die "Cannot open $path: $!\n";
    my $first = <$fh>;
    return sub { return } unless defined $first;
    if ($first =~ /\A>/) {
        my ($name, $length, $done) = (rstrip(substr($first, 1)), 0, 0);
        return sub {
            return if $done;
            while (defined(my $line = <$fh>)) {
                if ($line =~ /\A>/) {
                    my @record = ($name, $length);
                    ($name, $length) = (rstrip(substr($line, 1)), 0);
                    return @record;
                }
                $length += length strip($line);
            }
            $done = 1;
            return ($name, $length);
        };
    }
    if ($first =~ /\A\@/) {
        my $name = $first;
        return sub {
            return unless length $name;
            my ($seq, $plus, $qual) = map { scalar(<$fh>) // '' } 1 .. 3;
            die "Malformed FASTQ: $path\n" unless $name =~ /\A\@/ && $plus =~ /\A\+/ && length $qual;
            die "Quality length: $path\n" unless length(rstrip($seq)) == length(rstrip($qual));
            my @record = (rstrip(substr($name, 1)), length rstrip($seq));
            $name = <$fh> // '';
            return @record;
        };
    }
    die "Expected FASTA or FASTQ: $path\n";
}

sub add_counts { my ($to, $from) = @_; $to->{$_} += $from->{$_} for keys %$from }
sub total { my $n = 0; $n += $_ for values %{ $_[0] }; return $n }

sub audit {
    my $main = pure_path($_[0]);
    my ($stem, $suffix) = split_suffix($main);
    my (%parents, %sample_names);
    my $total = odict();
    open my $map, '<:encoding(UTF-8)', "$stem.map" or die "Cannot open $stem.map: $!\n";
    while (my $line = <$map>) {
        my @fields = split /\t/, rstrip($line), -1;
        @fields = ('') unless @fields;
        if ($fields[0] eq '#SMPLS') {
            for my $entry (@fields[1 .. $#fields]) {
                die "Malformed #SMPLS entry: $entry\n" unless $entry =~ /:/;
                my ($key, $value) = split /:/, $entry, 2;
                $sample_names{ py_int($key) } = $value;
            }
        }
        elsif ($line !~ /\A#/ && strip($line) ne '') {
            my ($name, $values) = ($fields[0], odict());
            die "Duplicate map ID: $name\n" if exists $parents{$name};
            for my $entry (@fields[1 .. $#fields]) {
                next unless $entry =~ /:/;
                my ($key, $count) = split /:/, $entry, 2;
                die "$entry\n" unless py_int($count) >= 0 && exists $sample_names{ py_int($key) };
                $values->{ $sample_names{ py_int($key) } } += py_int($count);
            }
            $parents{$name} = $values;
            add_counts($total, $values);
        }
    }
    close $map;
    my @files = (main => $main, merged => "$stem.merg$suffix", rest => "$main.rest");
    my $partitions = odict();
    my %seen;
    while (my ($label, $path) = splice @files, 0, 2) {
        my ($counts, $lengths, $number) = (odict(), odict(), 0);
        if (-e $path) {
            my $next = records($path);
            while (my ($name, $length) = $next->()) {
                die "Unknown/duplicate output parent: $name\n" unless exists $parents{$name} && !$seen{$name};
                $seen{$name} = 1;
                if (index($name, ';size=') >= 0) {
                    my $size = (split /;/, (split /;size=/, $name, -1)[1], -1)[0];
                    die "$name\n" unless py_int($size // '') == total($parents{$name});
                }
                add_counts($counts, $parents{$name});
                $lengths->{$length} += 1;
                $number += 1;
            }
        }
        $partitions->{$label} = odict(records => $number, counts => total($counts), samples => $counts, lengths => $lengths);
    }
    my $absent = grep { !$seen{$_} } keys %parents;
    die "$absent map parents absent from passing/rest outputs\n" if $absent;
    my %reconciled;
    add_counts(\%reconciled, $_->{samples}) for values %$partitions;
    for my $sample (keys %reconciled, keys %$total) {
        die "Per-sample map != passing + rest\n" unless ($reconciled{$sample} // 0) == ($total->{$sample} // 0);
    }
    return odict(main => $main, map_total => total($total), samples => $total, parents => scalar(keys %parents),
        passing => $partitions->{main}{counts} + $partitions->{merged}{counts},
        rest => $partitions->{rest}{counts}, outputs => $partitions);
}

# Insertion-ordered hash.
package DerepAudit::Ordered;
sub TIEHASH  { return bless { keys => [], values => {} }, shift }
sub STORE    { my ($s, $k, $v) = @_; push @{ $s->{keys} }, $k unless exists $s->{values}{$k}; $s->{values}{$k} = $v }
sub FETCH    { return $_[0]{values}{ $_[1] } }
sub EXISTS   { return exists $_[0]{values}{ $_[1] } }
sub DELETE   { my ($s, $k) = @_; $s->{keys} = [grep { $_ ne $k } @{ $s->{keys} }]; return delete $s->{values}{$k} }
sub CLEAR    { my ($s) = @_; $s->{keys} = []; $s->{values} = {} }
sub FIRSTKEY { $_[0]{at} = 0; return $_[0]->NEXTKEY }
sub NEXTKEY  { my ($s) = @_; return $s->{at} < @{ $s->{keys} } ? $s->{keys}[ $s->{at}++ ] : undef }
sub SCALAR   { return scalar @{ $_[0]{keys} } }

1;

#!/usr/bin/env perl
# Reconcile finalized SDM map, main, merged and rest counts (derepPerSR=0).
#
# Run: perl tests/audit_derep_counts.pl <-o_dereplicate path>
# Thin CLI over tests/lib/DerepAudit.pm with the argument handling, help text,
# JSON output (indent 2, ASCII) and exit codes of the former
# audit_derep_counts.py: 0 success, 1 audit failure, 2 usage error.
use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/lib";
use File::Basename qw(basename);
use JSON::PP ();
use DerepAudit qw(audit);

my $PROG = basename($0);
my $USAGE = "usage: $PROG [-h] output\n";
my $DESCRIPTION = 'Reconcile finalized SDM map, main, merged and rest counts (derepPerSR=0). '
    . 'Copied from hildebra/sdm tests/audit_derep_counts.py, 2026-09-15 integration '
    . 'contract. Run on fresh SDM output before a consumer appends/rewrites records.';

# Help width: COLUMNS, else the terminal width (Linux), else 80.
sub columns {
    my $c = $ENV{COLUMNS};
    return 0 + $c if defined $c && $c =~ /\A\s*\+?[0-9]+\s*\z/ && $c > 0;
    my $size = "\0" x 8;
    if ($^O eq 'linux' && -t STDOUT && ioctl(STDOUT, 0x5413, $size)) {    # TIOCGWINSZ
        my (undef, $cols) = unpack 'S2', $size;
        return $cols if $cols;
    }
    return 80;
}

sub wrap {
    my ($text, $width, $indent) = @_;
    my (@lines, $line);
    for my $word (split ' ', $text) {
        if (defined $line && length($line) + 1 + length($word) <= $width) { $line .= " $word" }
        else { push @lines, $line if defined $line; $line = $word }
    }
    push @lines, $line if defined $line;
    return join "\n$indent", @lines;
}

# Help layout of the original CLI.
sub help {
    my $width = columns() - 2;
    my $max_position = $width - 20 > 4 ? $width - 20 : 4;
    $max_position = 24 if $max_position > 24;
    my $position = 14 < $max_position ? 14 : $max_position;
    my $help_width = $width - $position > 11 ? $width - $position : 11;
    my $row = sub {
        my ($name, $text) = @_;
        my $body = wrap($text, $help_width, ' ' x $position) . "\n";
        return length($name) <= $position - 4 ? sprintf('%-*s', $position, "  $name") . $body
            : "  $name\n" . (' ' x $position) . $body;
    };
    print $USAGE, "\n", wrap($DESCRIPTION, $width > 11 ? $width : 11, ''), "\n\n",
        "positional arguments:\n", $row->('output', 'Literal -o_dereplicate path; derepPerSR=0'), "\n",
        "options:\n", $row->('-h, --help', 'show this help message and exit');
    exit 0;
}

sub usage_error { print STDERR $USAGE, "$PROG: error: $_[0]\n"; exit 2 }

my (@positional, @extras, $literal);
for my $arg (@ARGV) {
    if (!$literal && $arg eq '--') { $literal = 1; next }
    if (!$literal && $arg =~ /\A-./s && $arg !~ / / && $arg !~ /\A-[0-9]+\z|\A-[0-9]*\.[0-9]+\z/) {
        help() if $arg eq '-h' || (length($arg) > 2 && index('--help', $arg) == 0);
        push @extras, $arg;
        next;
    }
    @positional ? push @extras, $arg : push @positional, $arg;
}
usage_error('the following arguments are required: output') unless @positional;
usage_error("unrecognized arguments: @extras") if @extras;

my $result = eval { audit($positional[0]) };
if (!$result) {
    binmode STDERR, ':encoding(UTF-8)';
    print STDERR $@ =~ /\n\z/ ? $@ : "$@\n";
    exit 1;
}
utf8::decode(my $main = $result->{main});
$result->{main} = $main;
(my $json = JSON::PP->new->ascii->indent->indent_length(2)->space_after->encode($result)) =~ s/\x7f/\\u007f/g;
print $json;

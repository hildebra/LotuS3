#!/usr/bin/env perl
# LotuS3 uses Linux (LF) line endings only: no tracked text file may contain CR, and
# .gitattributes has to keep enforcing eol=lf (also on Windows, whatever core.autocrlf says).
use strict;
use warnings;
use Cwd qw(abs_path);
use FindBin;
use Test::More;

my $ROOT = abs_path("$FindBin::Bin/..");

open(my $ga, '<:raw', "$ROOT/.gitattributes") or BAIL_OUT("no .gitattributes: $!");
my $attrs = do { local $/; <$ga> };
close $ga;
unlike($attrs, qr/\r/, '.gitattributes itself uses LF');
like($attrs, qr/^\*[ \t]+text=auto[ \t]+eol=lf[ \t]*$/m, '.gitattributes enforces eol=lf for all files');

# Only meaningful inside a git checkout (skipped for a plain copy of the tree).
my @eol = `git -C "$ROOT" ls-files --eol 2>/dev/null`;
SKIP: {
    skip 'not a git checkout', 1 unless @eol;
    my @bad = grep { m{^i/(crlf|mixed)\b|\bw/(crlf|mixed)\b} && !m{^i/-text} } @eol;
    is(scalar(@bad), 0, 'no tracked text file has CRLF/mixed line endings') or diag(@bad);
}
done_testing();

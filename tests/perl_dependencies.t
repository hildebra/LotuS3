#!/usr/bin/env perl
# LotuS3 needs no Perl modules beyond those shipped with Perl itself: every module that
# lotus3, the helpers, bin/ scripts and the tests load must be core in Perl 5.14
# (the version lotus3 declares) and still core in the running Perl.
use strict;
use warnings;
use FindBin;
use Module::CoreList;
use Test::More;

my $ROOT = "$FindBin::Bin/..";
my $MIN_PERL = 5.014;
my %LOCAL = map { $_ => 1 } qw(LotusTest InstallerTest DerepAudit); # tests/lib
my @files = ("$ROOT/lotus3", glob("$ROOT/helpers/*.pl"), glob("$ROOT/bin/*.pl"),
    glob("$ROOT/tests/*.t"), glob("$ROOT/tests/*.pl"), glob("$ROOT/tests/lib/*.pm"));

subtest test_only_core_modules => sub {
    my %users;
    for my $file (@files) {
        open my $fh, '<', $file or die "$file: $!\n";
        while (my $line = <$fh>) {
            next if $line =~ /^\s*#/;
            # statements only (line start, or after ; or {), not prose inside quoted strings
            (my $code = $line) =~ s/"(?:[^"\\]|\\.)*"|'(?:[^'\\]|\\.)*'/""/g;
            while ($code =~ /(?:^|[;{]\s*)(?:use|require)\s+([A-Z][\w:]*)(?=[\s;(])/g) {
                push @{ $users{$1} }, $file =~ s{^\Q$ROOT\E/}{}r;
            }
        }
    }
    ok(scalar keys %users > 20, 'found the module list') or diag(join ' ', sort keys %users);
    for my $module (sort keys %users) {
        next if $LOCAL{$module};
        my $first = Module::CoreList->first_release($module);
        my $where = join(', ', sort { $a cmp $b } keys %{{ map { $_ => 1 } @{ $users{$module} } }});
        ok(defined $first && $first <= $MIN_PERL, "$module is core in Perl $MIN_PERL")
            or diag("$module (used by $where) first core release: " . ($first // 'never'));
        ok(!Module::CoreList->removed_from($module), "$module is still core") or diag("used by $where");
    }
};

done_testing();

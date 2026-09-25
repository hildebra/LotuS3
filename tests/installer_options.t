#!/usr/bin/env perl
# Exercise interactive choices and real ONT preflight/registration in temp installs.
#
# Stop at the database-download boundary, recording selected database/program modes.
# No reference databases or non-ONT packages are downloaded or built by these tests.
use strict;
use warnings;
use Cwd qw(abs_path);
use FindBin;
use JSON::PP qw(decode_json);
use Test::More;
use lib "$FindBin::Bin/lib";
use InstallerTest qw(program perl_script read_file write_file write_executable replaced contains_ok lacks_ok);

my $ROOT = abs_path("$FindBin::Bin/..");

# Replaces the get_DBs() call: report the selected options, register ONT programs, stop.
my $CHECKPOINT = "\n" . <<'PERL';
require JSON::PP;
print "CHOICES:" . JSON::PP->new->encode({
    search => $installBlast, databases => \@refDBinstall, its => $ITSready,
    utax => $getUTAX, ont => $installONT, programs => [ont_programs_to_install()],
    db_only => ($onlyDbinstall || $condaDBinstall) ? 1 : 0,
    r_packages => $install_dada,
}) . "\n";
install_ont_programs() unless $onlyDbinstall || $condaDBinstall;
finishAI("none");
exit(0);
PERL

sub fixture { return InstallerTest->new($ROOT) }

sub run_choices {
    my ($t, $answers, %opt) = @_;
    my $source = read_file($t->{installer});
    ok($source =~ s/\nget_DBs\(\);\n/$CHECKPOINT/, 'installer calls get_DBs()');
    write_file($t->{installer}, $source);
    my $output = $t->run_installer(%opt, ont_only => 0, answers => $answers);
    my ($line) = grep { /^CHOICES:/ } split /\r?\n/, $output;
    return ($output, defined $line ? decode_json(substr($line, length 'CHOICES:')) : undef);
}

# Compare JSON-encoded, so numbers and strings stay distinct.
sub is_json {
    my ($got, $expected, $name) = @_;
    local $Test::Builder::Level = $Test::Builder::Level + 1;
    my $json = JSON::PP->new->canonical->allow_nonref;
    is($json->encode($got), $json->encode($expected), $name);
}

sub prepare_programs {
    my ($t) = @_;
    $t->installed_tools;
    $t->tool(Rscript => qq{#!/bin/sh\necho "R scripting front-end version 4.4.0"\n});
    # checked by the full-install preflight before any download
    $t->tool($_ => "#!/bin/sh\nexit 0\n") for qw(unzip make cc xz);
    -d "$t->{install}/bin" or mkdir "$t->{install}/bin" or die "$t->{install}/bin: $!\n";
    write_executable("$t->{install}/bin/sdm", perl_script(qq{print "sdm 3.53 beta\\n";\n}));
    write_executable("$t->{install}/bin/LCA", perl_script(qq{print "LCA 0.29\\n";\n}));
}

sub assert_all {
    my ($t, $answer, $prefix) = @_;
    local $Test::Builder::Level = $Test::Builder::Level + 1;
    prepare_programs($t);
    my ($output, $state) = run_choices($t, ($prefix // '') . "$answer\ny\n");
    is($state->{search}, 3, 'search');
    is_json($state->{databases}, [(0) x 8, 1, 0], 'databases');
    is($state->{$_}, 1, $_) for qw(its utax ont r_packages);
    is_json($state->{programs}, [qw(minimap2 savont barbell)], 'programs');
    is($t->config->{$_}, "$t->{tools}/$_", "$_ registered") for @{$state->{programs}};
    lacks_ok($output, $_) for 'For similarity based', 'Do you want to install a reference database',
        '-- ITS --', '-- UTAX --', '-- ONT --';
    contains_ok($output, 'Do you accept (y/n)?');
    contains_ok($output, 'Install LotuS3 with all possible dependencies');
    lacks_ok($output, 'Continue (y/n)?');
}

subtest 'test_enter_selects_every_dependency' => sub {
    assert_all(fixture(), '');
};

subtest 'test_one_selects_every_dependency' => sub {
    assert_all(fixture(), '1');
};

subtest 'test_refresh_programs_offers_all_dependencies' => sub {
    my $t = fixture();
    write_file($t->{cfg}, replaced(read_file($t->{default}), 'UID ??', 'UID 123'));
    assert_all($t, '1', "1\n");
};

subtest 'test_detailed_without_ont_preserves_entries_and_skips_rust_checks' => sub {
    my $t = fixture();
    prepare_programs($t);
    $t->tool(savont => program('savont', crash => 1)); $t->tool(barbell => program('barbell', crash => 1));
    write_file($t->{cfg}, read_file($t->{default}) . "savont previous-savont\nbarbell previous-barbell\n");
    my ($output, $state) = run_choices($t, "0\n2\n0\n0\n0\n0\n");
    is_json($state->{programs}, ['minimap2'], 'programs');
    is($state->{search}, 2, 'search');
    is_json($state->{databases}, [1, (0) x 9], 'databases');
    is($state->{$_}, 0, $_) for qw(its utax ont);
    for my $name (qw(savont barbell)) {
        is($t->config->{$name}, "previous-$name", "$name entry preserved");
        lacks_ok($output, "$t->{tools}/$name");
        ok(!-e "$t->{install}/bin/$name", "no $name installed");
    }
    is($t->config->{minimap2}, "$t->{tools}/minimap2", 'minimap2 registered');
    contains_ok($output, '-- ONT -- Install ONT tools');
    lacks_ok($output, 'Bioconda');
};

sub assert_detailed_ont {
    my ($t, $answer) = @_;
    local $Test::Builder::Level = $Test::Builder::Level + 1;
    prepare_programs($t);
    my ($output, $state) = run_choices($t, "0\n1\n1\n1\n0\n$answer\n");
    is($state->{search}, 1, 'search');
    is_json($state->{databases}, [0, 1, (0) x 8], 'databases');
    is($state->{its}, 1, 'its');
    is($state->{utax}, 0, 'utax');
    is($state->{ont}, 1, 'ont');
    is_json($state->{programs}, [qw(minimap2 savont barbell)], 'programs');
    is($t->config->{$_}, "$t->{tools}/$_", "$_ registered") for @{$state->{programs}};
    contains_ok($output, '-- ONT -- Install ONT tools');
    lacks_ok($output, 'Do you accept (y/n)?');
}

subtest 'test_detailed_with_ont_registers_all_three_tools' => sub {
    assert_detailed_ont(fixture(), '1');
};

subtest 'test_detailed_ont_defaults_to_yes' => sub {
    assert_detailed_ont(fixture(), '');
};

subtest 'test_invalid_choices_are_reprompted' => sub {
    my $t = fixture();
    prepare_programs($t);
    my ($output, $state) = run_choices($t, "2\nbad\n0\n2\n1\n1\n0\n2\n\n");
    is(scalar(() = $output =~ /Invalid answer;/g), 3, 'three reprompts');
    is($state->{search}, 2, 'search');
    is($state->{ont}, 1, 'ont');
};

subtest 'test_all_dependencies_still_requires_silva_acceptance' => sub {
    my $t = fixture();
    my ($output, $state) = run_choices($t, "\nn\n", ok => 0);
    is($state, undef, 'stopped before the checkpoint');
    contains_ok($output, 'You need to accept the SILVA license');
    ok(!-e $t->{cfg}, 'no config written');
    ok(!-e "$t->{install}/DB", 'no DB directory');
};

subtest 'test_end_of_input_is_not_default_acceptance' => sub {
    my $t = fixture();
    my $original = read_file($t->{installer});
    for (["", 'all-dependencies choice'], ["\n", 'SILVA license response'], ["0\n2\n0\n0\n0\n", 'ONT tools choice']) {
        my ($answers, $context) = @$_;
        subtest "context=$context" => sub {
            write_file($t->{installer}, $original);
            my ($output, $state) = run_choices($t, $answers, ok => 0);
            is($state, undef, 'stopped before the checkpoint');
            contains_ok($output, "End of input while waiting for the $context");
            ok(!-e $t->{cfg}, 'no config written');
        };
    }
};

subtest 'test_database_only_refresh_has_no_program_questions' => sub {
    my $t = fixture();
    write_file($t->{cfg}, replaced(read_file($t->{default}), 'UID ??', 'UID 123'));
    my ($output, $state) = run_choices($t, "2\n8\ny\n1\n1\n");
    is($state->{db_only}, 1, 'db_only');
    is_json($state->{databases}, [(0) x 8, 1, 0], 'databases');
    lacks_ok($output, 'all possible dependencies');
    lacks_ok($output, '-- ONT --');
    ok(!exists $t->config->{barbell}, 'no barbell entry');
    ok(!-e "$t->{install}/bin/barbell", 'no barbell installed');
};

subtest 'test_rscript_version_on_stderr_is_recognised' => sub {
    my $t = fixture();
    prepare_programs($t);
    # R before 4.2 prints "Rscript --version" to stderr only
    $t->tool(Rscript => qq{#!/bin/sh\necho "R scripting front-end version 4.1.2 (2021-11-01)" >&2\n});
    my ($output, $state) = run_choices($t, "\ny\n");
    lacks_ok($output, 'older than 4');
    is($state->{r_packages}, 1, 'R packages still installed');
};

subtest 'test_missing_build_tools_stop_before_downloads' => sub {
    my $t = fixture();
    prepare_programs($t);
    unlink "$t->{tools}/unzip" or die "Cannot remove unzip: $!\n";
    write_executable("$t->{install}/bin/sdm", perl_script("kill 'TERM', \$\$;\n"));
    my ($output, $state) = run_choices($t, "\ny\n", ok => 0);
    is($state, undef, 'stopped before the database downloads');
    contains_ok($output, 'The full installation needs:');
    contains_ok($output, ' - unzip');
    contains_ok($output, 'a working sdm');
    contains_ok($output, 'Nothing has been downloaded yet.');
};

subtest 'test_conda_database_mode_stays_noninteractive' => sub {
    my $t = fixture();
    my ($output, $state) = run_choices($t, '', extra => ['--condaDBinstall']);
    is($state->{db_only}, 1, 'db_only');
    lacks_ok($output, 'Answer:');
    lacks_ok($output, 'all possible dependencies');
    ok(!exists $t->config->{barbell}, 'no barbell entry');
};

done_testing();

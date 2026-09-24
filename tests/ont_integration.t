#!/usr/bin/env perl
# ONT wiring regressions: real bundled SDM, controlled Barbell/Savont/mappers.
#
# Run: prove -v tests/ont_integration.t
# No scientific validation of the stand-in tools is implied. A temporary copy of
# lotus3 exits after SDM builds the abundance matrix, before unrelated taxonomy.
use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/lib";
use File::Basename qw(basename);
use File::Spec;
use File::Temp qw(tempdir);
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use JSON::PP ();
use Test::More;
use LotusTest qw($ROOT $SDM case contains lacks read_text write_text append_text text_lines
    json_decode json_encode rc run_command after count_of every4 has_arg replace_all replace_first);

sub canon { return File::Spec->canonpath($_[0]) }
sub entries { my ($dir) = @_; opendir my $dh, $dir or die "$dir: $!\n"; return [sort grep { !/^\.\.?$/ } readdir $dh] }
sub savont_opts { return read_text("$ROOT/configs/sdm_ONT_SAVONT_opt.txt") }
sub gunzip_text { my ($file) = @_; gunzip($file => \my $text) or die "gunzip $file: $GunzipError\n"; return $text }

sub append_read {
    my ($t, $name, %o) = @_;
    my $sequence = $o{sequence} // $t->{fwd} . $t->{seq} . $t->{revcomp};
    my $quality = $o{quality} // ('I' x length $sequence);
    is(length $sequence, length $quality, "$name: sequence and quality lengths agree");
    append_text(($o{folder} // $t->{reads}) . '/' . ($o{sample} // 's1') . '.fq', "\@$name\n$sequence\n+\n$quality\n");
}

sub ont_read {
    my ($t, %o) = @_;
    my $body = $o{body} || $t->{seq};
    my $left = ('G' x ($o{offset} // 21)) . ($o{front} // '') . 'CTAGCATGATC';
    my $seq = $left . $t->{fwd} . $body . $t->{revcomp} . 'CATTGAC' . rc($o{rear} // '') . ('T' x 19);
    my $qual = join '', map { chr(33 + 25 + $_ % 15) } 0 .. length($seq) - 1;
    my $body_qual = substr($qual, length($left) + length($t->{fwd}), length $body);
    return $o{reverse} ? [rc($seq), scalar reverse($qual), $body_qual] : [$seq, $qual, $body_qual];
}

sub pooled_input {
    my ($t, $reads) = @_;
    write_text($t->{map}, "#SampleID\tBarcodeSequence\tForwardPrimer\tReversePrimer\n"
        . join('', map { "$_->[0]\t$_->[1]\t$t->{fwd}\t$t->{rev}\n" } ['s1', 'ACGTCAGTGCTAGACG'], ['s2', 'TGCATCGACAGTTCGA']));
    my $pooled = "$t->{root}/pooled.fastq.gz";
    my $fq = join '', map { "\@pooled_$_\n$reads->[$_][0]\n+\n$reads->[$_][1]\n" } 0 .. $#$reads;
    gzip(\$fq => $pooled) or die "gzip failed: $GzipError\n";
    return ['-i', $pooled];
}

sub sdm_version_wrapper {
    my ($t, $banner) = @_;
    my $wrapper = "$t->{tools}/sdm";
    my $literal = JSON::PP->new->allow_nonref->encode($banner);
    write_text($wrapper, "#!/usr/bin/env perl\nif (\@ARGV && (\$ARGV[0] eq \"-version\" || \$ARGV[0] eq \"-v\")) "
        . "{ print $literal; exit 0; }\nexec \"$SDM\", \@ARGV;\n");
    chmod 0755, $wrapper;
    write_text($t->{cfg}, replace_all(read_text($t->{cfg}), "sdm $SDM\n", "sdm $wrapper\n"));
}

sub strip_fixture_primers {
    my ($t, $folder) = @_;
    for my $name (@{ entries($folder) }) {
        my @lines = text_lines(read_text("$folder/$name"));
        for (my $i = 0; $i < @lines; $i += 4) {
            for my $j ($i + 1, $i + 3) {
                $lines[$j] = substr($lines[$j], length $t->{fwd}, length($lines[$j]) - length($t->{fwd}) - length($t->{revcomp}));
            }
        }
        write_text("$folder/$name", join("\n", @lines) . "\n");
    }
}

sub check_removed_primers {
    my ($t, %o) = @_;
    my $original_map = read_text($t->{map});
    $t->strip_fixture_primers($o{barbell} ? $t->{barcodes} : $t->{reads});
    $t->run_lotus(extra => ['-ontPrimerState', 'removed', '-ontMinReads', '2', @{ $o{extra} // [] }], barbell => $o{barbell});
    $t->check_counts;
    is(read_text($t->{map}), $original_map, 'user map unchanged');
    contains(read_text("$t->{out}/primary/in.map"), 'ForwardPrimer');
    my $effective_map = read_text("$t->{out}/primary/sdm_input.map");
    lacks($effective_map, $_) for qw(ForwardPrimer ReversePrimer LinkerPrimerSequence);
    contains($effective_map, 'fastqFile');
    my @lines = text_lines(read_text("$t->{out}/tmpFiles/savont_out/input_snapshot.fq"));
    is_deeply(every4(\@lines, 1), [($t->{seq}) x 7], 'full sequences');
    is_deeply(every4(\@lines, 3), [('I' x length $t->{seq}) x 7], 'full qualities');
    contains(read_text("$t->{out}/LotuSLogS/run_manifest.txt"), 'ONT primers entering SDM: removed');
}

# Make the helpers above callable as fixture methods.
{ no strict 'refs'; *{"LotusTest::$_"} = \&{"main::$_"} for qw(append_read ont_read pooled_input sdm_version_wrapper strip_fixture_primers check_removed_primers); }

case test_savont_counts_and_filtered_reads_without_dereplication => sub {
    my $t = shift;
    $t->run_lotus;
    my $state = $t->check_counts;
    contains($state->{preset}, 'map-ont');
    is($state->{cluster}, 9, 'Savont clusterer');
    my $snapshot = read_text("$t->{out}/tmpFiles/savont_out/input_snapshot.fq");
    is(count_of($snapshot, "\n+\n"), 7, 'seven reads reach Savont');
    lacks($snapshot, $t->{fwd});
    is($state->{dereplication}, 0, 'no SDM dereplication');
    is_deeply([glob("$t->{out}/tmpFiles/derep*")], [], 'no derep files');
    ok(!-e "$t->{out}/tmpFiles/savont_in.fq", 'no separate Savont input copy');
    my $prepared = "$t->{out}/tmpFiles/savont_reads.fq";
    is($snapshot, read_text($prepared), 'Savont reads are the prepared FASTQ');
    is(canon($t->tool_calls('savont')->[0][2]), canon($prepared), 'Savont input');
    is(canon($t->tool_calls('minimap2')->[0][-1]), canon($prepared), 'minimap2 input');
    like($snapshot, qr/\A\@s1___/, 'sample-prefixed read IDs');
    is(basename($state->{sdm_options}), 'sdm_ONT_SAVONT_effective.txt', 'effective options file');
    contains(read_text("$t->{out}/primary/sdm_original_options.txt"), "maxAmbiguousNT\t5%");
    contains(read_text($state->{sdm_options}), "maxAmbiguousNT\t-1");
    contains(read_text("$t->{out}/tmpFiles/tmp_otu.fa"), $t->{consensus});
    my @seed = text_lines(read_text($state->{seed}));
    is(join('', @seed[1 .. $#seed]), $t->{consensus}, 'seed keeps Savont consensus');
    ok(!-e "$t->{out}/tmpFiles/demultiplexed", 'no demultiplexed copy');
    my $args = $t->tool_calls('savont')->[0];
    is(after($args, '--quality-value-cutoff'), '80');
    is(after($args, '--min-read-length'), '1000');
    is(after($args, '--max-read-length'), '2000');
    ok(has_arg($args, '--single-strand'), 'single-strand mode');
};

case test_savont_both_strand_override => sub {
    my $t = shift;
    $t->run_lotus(extra => ['-savontSingleStrand', '0']);
    $t->check_counts;
    ok(!has_arg($t->tool_calls('savont')->[0], '--single-strand'), 'both strands');
};

case test_savont_main_quality_filters_exclude_reads_from_asvs_and_counts => sub {
    my $t = shift;
    my $full = $t->{fwd} . $t->{seq} . $t->{revcomp};
    $t->append_read('low_average', quality => '+' x length $full);  # Q10
    $t->append_read('low_window', quality => ('I' x 219) . ('&' x 150) . ('I' x (length($full) - 369)));
    $t->append_read('ambiguous', sequence => substr($full, 0, 200) . ('N' x 61) . substr($full, 261));
    $t->append_read('missing_reverse', sequence => $t->{fwd} . $t->{seq});
    $t->run_lotus;
    $t->check_counts;
    my $snapshot = read_text("$t->{out}/tmpFiles/savont_out/input_snapshot.fq");
    is(count_of($snapshot, "\n+\n"), 7, 'seven reads reach Savont');
    lacks($snapshot, $_) for qw(low_average low_window ambiguous missing_reverse);
    my $commands = $t->sdm_commands;
    lacks($commands->[0], '-o_dereplicate');
    lacks($commands->[1], '-derep_map');
    lacks($commands->[1], '-options');  # Prepared reads are not trimmed again.
};

case test_savont_percentage_and_length_boundaries => sub {
    my $t = shift;
    my %passing;
    my @cases = (['short_boundary', 1000, 50, 'N', 1], ['long_boundary', 2000, 100, 'N', 1],
        ['iupac_boundary', 1000, 50, 'R', 1], ['short_over', 1000, 51, 'N', 0],
        ['long_over', 2000, 101, 'N', 0], ['too_short', 999, 0, 'N', 0], ['too_long', 2001, 0, 'N', 0]);
    for my $c (@cases) {
        my ($name, $length, $count, $base, $accepted) = @$c;
        my $sequence = substr($t->{seq} x 2, 0, $length);
        substr($sequence, 100 + 2 * $_, 1) = $base for 0 .. $count - 1;
        $t->append_read($name, sequence => $t->{fwd} . $sequence . $t->{revcomp});
        $passing{$name} = $sequence if $accepted;
    }
    $t->run_lotus(extra => ['-saveDemultiplex', '2']);
    $t->check_counts({ s1 => 7, s2 => 3 });
    my @lines = text_lines(read_text("$t->{out}/tmpFiles/savont_reads.fq"));
    my %actual;
    for (my $i = 0; $i < @lines; $i += 4) {
        my ($header) = split ' ', $lines[$i];
        my (undef, $read) = split /___/, $header, 2;
        $actual{$read} = $lines[$i + 1];
    }
    for my $c (@cases) {
        my ($name, undef, undef, undef, $accepted) = @$c;
        if ($accepted) { is($actual{$name}, $passing{$name}, "$name kept intact") }
        else { ok(!exists $actual{$name}, "$name rejected") }
    }
    # Saved per-sample copies are filtered first, then compressed.
    is_deeply([glob("$t->{out}/demultiplexed/*.fq")], [], 'no uncompressed per-sample FASTQ');
    my $saved = join '', map { gunzip_text($_) } glob("$t->{out}/demultiplexed/*.fq.gz");
    contains($saved, 'short_boundary');
    lacks($saved, 'short_over');
    lacks($saved, 'long_over');
    contains(read_text("$t->{out}/LotuSLogS/savont_ambiguity_filter.log"), 'kept 10; rejected 2');
};

case test_savont_average_quality_fifteen => sub {
    my $t = shift;
    my $custom = write_text("$t->{root}/quality-only.txt",
        replace_all(savont_opts(), "QualWindowThreshhold\t14", "QualWindowThreshhold\t-1"));
    my $full = $t->{fwd} . $t->{seq} . $t->{revcomp};
    $t->append_read('q15', quality => '0' x length $full);
    $t->append_read('q14', quality => '/' x length $full);
    $t->run_lotus(extra => ['-s', $custom]);
    $t->check_counts({ s1 => 5, s2 => 3 });
    my $prepared = read_text("$t->{out}/tmpFiles/savont_reads.fq");
    contains($prepared, '___q15 ');
    lacks($prepared, '___q14 ');
};

case test_savont_integer_ambiguity_override => sub {
    my $t = shift;
    my $custom = write_text("$t->{root}/no-ambiguity.txt", replace_all(savont_opts(), "maxAmbiguousNT\t5%", "maxAmbiguousNT\t0"));
    my $seq = $t->{fwd} . $t->{seq} . $t->{revcomp};
    $t->append_read('one_n', sequence => substr($seq, 0, 200) . 'N' . substr($seq, 201));
    $t->run_lotus(extra => ['-s', $custom]);
    $t->check_counts;
    lacks(read_text("$t->{out}/tmpFiles/savont_reads.fq"), 'one_n');
    ok(!-e "$t->{out}/LotuSLogS/savont_ambiguity_filter.log", 'no percentage filter log');
};

case test_savont_invalid_ambiguity_percentage => sub {
    my $t = shift;
    my $custom = "$t->{root}/invalid-ambiguity.txt";
    for my $value ('-5%', '101%', '5.5', 'wrong') {
        subtest "value=$value" => sub {
            write_text($custom, replace_all(savont_opts(), "maxAmbiguousNT\t5%", "maxAmbiguousNT\t$value"));
            my $result = $t->run_lotus(extra => ['-s', $custom, '--dry-run'], ok => 0);
            contains($result->{output}, 'integer count or a percentage from 0% to 100%');
            is_deeply($t->tool_calls, [], 'no tools called');
        };
    }
};

case test_savont_preserves_high_quality_homopolymer_tail => sub {
    my $t = shift;
    for my $sample (qw(s1 s2)) {
        my $path = "$t->{reads}/$sample.fq";
        my @lines = text_lines(read_text($path));
        for (my $i = 0; $i < @lines; $i += 4) {
            $lines[$i + 1] = $t->{seq} . ('G' x 15);
            $lines[$i + 3] = 'I' x length $lines[$i + 1];
        }
        write_text($path, join("\n", @lines) . "\n");
    }
    $t->run_lotus(extra => ['-ontPrimerState', 'removed']);
    $t->check_counts;
    my @prepared = text_lines(read_text("$t->{out}/tmpFiles/savont_reads.fq"));
    is_deeply(every4(\@prepared, 1), [($t->{seq} . ('G' x 15)) x 7], 'tail kept');
    is_deeply(every4(\@prepared, 3), [('I' x (length($t->{seq}) + 15)) x 7], 'tail qualities kept');
};

case test_savont_does_not_recover_secondary_quality_reads => sub {
    my $t = shift;
    my $opts = replace_all(savont_opts(), "minAvgQuality\t15", "minAvgQuality\t25");
    my $custom = write_text("$t->{root}/custom.txt", replace_all($opts, "*minAvgQuality\t25", "*minAvgQuality\t0"));
    my $full = $t->{fwd} . $t->{seq} . $t->{revcomp};
    $t->append_read('secondary_only', quality => '6' x length $full);  # Q21
    $t->run_lotus(extra => ['-s', $custom, '-saveDemultiplex', '2']);
    $t->check_counts;
    lacks(read_text("$t->{out}/tmpFiles/savont_out/input_snapshot.fq"), 'secondary_only');
    is_deeply([glob("$t->{out}/tmpFiles/derep*")], [], 'no derep files');
};

case test_savont_all_reads_filtered_stops_before_clustering => sub {
    my $t = shift;
    for my $sample (qw(s1 s2)) {
        my $file = "$t->{reads}/$sample.fq";
        my @lines = text_lines(read_text($file));
        for (my $i = 3; $i < @lines; $i += 4) { $lines[$i] = '+' x length $lines[$i] }
        write_text($file, join("\n", @lines) . "\n");
    }
    my $result = $t->run_lotus(ok => 0);
    contains($result->{output}, 'No reads passed SDM filtering for Savont');
    is_deeply($t->tool_calls, [], 'no tools called');
    ok(!-e "$t->{out}/OTU.txt", 'no OTU table');
};

case test_savont_sample_prefixes_preserve_shared_read_names => sub {
    my $t = shift;
    for my $sample (qw(s1 s2)) {
        my $file = "$t->{reads}/$sample.fq";
        write_text($file, replace_all(read_text($file), "\@${sample}_", '@shared_'));
    }
    $t->run_lotus;
    $t->check_counts;
};

case test_savont_zero_hit_sample_is_preserved => sub {
    my $t = shift;
    $t->{env}{ONT_TEST_SKIP_PREFIX} = 's2___';
    $t->run_lotus;
    $t->check_counts({ s1 => 4, s2 => 0 });
};

case test_savont_combine_samples => sub {
    my $t = shift;
    my @rows = text_lines(read_text($t->{map}));
    write_text($t->{map}, "$rows[0]\tCombineSamples\n" . join("\n", map { "$_\tpooled" } @rows[1 .. $#rows]) . "\n");
    for my $sample (qw(s1 s2)) {
        my $file = "$t->{reads}/$sample.fq";
        write_text($file, replace_all(read_text($file), "\@${sample}_", '@shared_'));
    }
    $t->run_lotus;
    $t->check_counts({ pooled => 7 });
    lacks(read_text("$t->{out}/primary/sdm_input.map"), 'CombineSamples');
    contains(read_text("$t->{out}/primary/in.map"), 'CombineSamples');
    my $snapshot = read_text("$t->{out}/tmpFiles/savont_out/input_snapshot.fq");
    contains($snapshot, '@s1___shared_0');
    contains($snapshot, '@s2___shared_0');
};

case test_savont_mixed_groups_and_ungrouped_samples => sub {
    my $t = shift;
    $t->write_map(rows => [['s1', 's1.fq'], ['s2', 's2.fq'], ['low', 'low.fq']]);
    my @rows = text_lines(read_text($t->{map}));
    write_text($t->{map}, "$rows[0]\tCombineSamples\n"
        . join("\n", map { $_ . (/^s[12]\t/ ? "\tpooled" : "\t") } @rows[1 .. $#rows]) . "\n");
    $t->run_lotus;
    $t->check_counts({ pooled => 7, low => 1 });
};

case test_other_clusterers_keep_sdm_dereplication => sub {
    my $t = shift;
    for my $platform (qw(ONT miSeq)) {
        subtest "platform=$platform" => sub {
            $t->{out} = "$t->{root}/output_$platform";
            $t->run_lotus(extra => ['-p', $platform, '-CL', 'vsearch', '-derepMin', '0', '-s', "$ROOT/configs/sdm_ONT_LSSU.txt"]);
            my $state = $t->check_counts;
            isnt($state->{dereplication}, 0, 'SDM dereplication kept');
            my $derep = read_text("$t->{out}/tmpFiles/derep.fas");
            is(count_of($derep, '>'), 1, 'one dereplicate');
            contains($derep, ';size=7;');
            ok(-f "$t->{out}/tmpFiles/derep.map", 'derep map');
            contains($t->sdm_commands->[0], '-o_dereplicate');
            contains($t->sdm_commands->[1], '-derep_map');
            is_deeply($t->tool_calls('savont'), [], 'Savont not called');
            ok(!-e "$t->{out}/tmpFiles/savont_reads.fq", 'no Savont reads');
        };
    }
};

case test_pooled_offset_barcodes_indels_orientation_and_full_qualities => sub {
    my $t = shift;
    my $bc = 'ACGTCAGTGCTAGACG';
    my @variants = ($bc, substr($bc, 0, 6) . 'A' . substr($bc, 7), substr($bc, 0, 7) . 'T' . substr($bc, 7), substr($bc, 0, 7) . substr($bc, 8));
    my @reads = map { my $x = $_; map { $t->ont_read(front => $x, reverse => $_) } 0, 1 } @variants;
    push @reads, map { $t->ont_read(rear => 'TGCATCGACAGTTCGA', reverse => $_) } 0, 1, 0;
    $t->run_lotus(extra => $t->pooled_input(\@reads));
    $t->check_counts({ s1 => 8, s2 => 3 });
    my @fq = text_lines(read_text("$t->{out}/tmpFiles/savont_out/input_snapshot.fq"));
    is_deeply([sort @{ every4(\@fq, 1) }], [($t->{seq}) x @reads], 'trimmed, oriented sequences');
    is_deeply([sort @{ every4(\@fq, 3) }], [sort map { $_->[2] } @reads], 'full qualities');
    my $commands = $t->sdm_commands;
    is(scalar @$commands, 2, 'two SDM commands') or diag(join "\n", @$commands);
    contains($commands->[0], '-ontMode 1 -barcodeSearchWindow 200 -ontBarcodeEnds either');
    contains($commands->[0], '-o_fastq ');
    lacks($commands->[0], '-o_dereplicate');
    lacks($commands->[1], $_) for qw(-derep_map -ontMode -barcodeSearchWindow -ontBarcodeEnds);
    my $manifest = read_text("$t->{out}/LotuSLogS/run_manifest.txt");
    contains($manifest, 'SDM barcode/primer search window: 200');
    contains($manifest, 'SDM barcode ends: either');
};

case test_barcode_end_requirement_changes_counts_and_rejects_conflicts => sub {
    my $t = shift;
    my ($a, $b) = ('ACGTCAGTGCTAGACG', 'TGCATCGACAGTTCGA');
    my @reads = ($t->ont_read(front => $a), $t->ont_read(front => $a, rear => $a), $t->ont_read(rear => $b),
        $t->ont_read(front => $b, rear => $b, reverse => 1), $t->ont_read(front => $a, rear => $b));
    my $extra = $t->pooled_input(\@reads);
    for my $c (['either', 2], ['both', 1]) {
        my ($ends, $count) = @$c;
        subtest "ends=$ends" => sub {
            $t->{out} = "$t->{root}/output_$ends";
            $t->run_lotus(extra => [@$extra, '-ontBarcodeEnds', $ends]);
            $t->check_counts({ s1 => $count, s2 => $count });
            contains($t->sdm_commands->[0], "-ontBarcodeEnds $ends");
        };
    }
};

case test_barcode_search_window_changes_retained_reads => sub {
    my $t = shift;
    my @reads = map { my $bc = $_; map { $t->ont_read(front => $bc, offset => $_) } 0, 230 } 'ACGTCAGTGCTAGACG', 'TGCATCGACAGTTCGA';
    my $extra = $t->pooled_input(\@reads);
    for my $c ([200, 1], [300, 2]) {
        my ($window, $count) = @$c;
        subtest "window=$window" => sub {
            $t->{out} = "$t->{root}/output_$window";
            $t->run_lotus(extra => [@$extra, '-barcodeSearchWindow', $window]);
            $t->check_counts({ s1 => $count, s2 => $count });
        };
    }
};

case test_both_ends_ignored_for_filename_assignment_after_barbell => sub {
    my $t = shift;
    $t->write_map(barbell => 1);
    $t->run_lotus(extra => ['-ontMinReads', '2', '-ontBarcodeEnds', 'both'], barbell => 1);
    $t->check_counts;
    contains($t->sdm_commands->[0], '-ontMode 1');
};

case test_non_ont_primary_command_omits_ont_options => sub {
    my $t = shift;
    $t->run_lotus(extra => ['-p', 'miSeq', '-CL', 'vsearch', '-saveDemultiplex', '1']);
    my $commands = $t->sdm_commands;
    is(scalar @$commands, 1, 'one SDM command') or diag(join "\n", @$commands);
    lacks($commands->[0], $_) for qw(-ontMode -barcodeSearchWindow -ontBarcodeEnds);
};

case test_old_sdm_rejected_before_replacing_output => sub {
    my $t = shift;
    $t->sdm_version_wrapper("sdm 3.50 beta\n");
    mkdir $t->{out} or die "$!\n";
    my $marker = write_text("$t->{out}/LotuS_output_schema_version.txt", 'keep');
    my $result = $t->run_lotus(ok => 0);
    contains($result->{output}, 'ONT preprocessing requires SDM >= 3.51');
    is(read_text($marker), 'keep', 'previous output kept');
    is_deeply($t->tool_calls, [], 'no tools called');
};

case test_capability_marker_allows_older_numbered_build => sub {
    my $t = shift;
    $t->sdm_version_wrapper("sdm 3.50 beta\nONT amplicon end matching: enabled (-ontMode 1)\n");
    $t->run_lotus(extra => ['--dry-run']);
    ok(!-e $t->{out}, 'dry run creates no output');
};

case test_non_ont_can_use_older_sdm => sub {
    my $t = shift;
    $t->sdm_version_wrapper("sdm 3.43\n");
    $t->run_lotus(extra => ['-p', 'miSeq', '-CL', 'vsearch', '--dry-run']);
    ok(!-e $t->{out}, 'dry run creates no output');
};

case test_ont_rejects_incompatible_map_columns_even_empty => sub {
    my $t = shift;
    my $original = read_text($t->{map});
    for my $field (qw(Barcode2ndPair MIDfqFile SampleIDinHead alignmentFile fnaFile qualFile)) {
        subtest "field=$field" => sub {
            my @rows = text_lines($original);
            write_text($t->{map}, "$rows[0]\t$field\n" . join("\n", map { "$_\t" } @rows[1 .. $#rows]) . "\n");
            my $result = $t->run_lotus(extra => ['--dry-run'], ok => 0);
            contains($result->{output}, $field);
            ok(!-e $t->{out}, 'no output created');
        };
    }
};

case test_ont_barcode_budget_and_alphabet_rejected_early => sub {
    my $t = shift;
    $t->pooled_input([$t->ont_read(front => 'ACGTCAGTGCTAGACG')]);
    my $original = read_text($t->{map});
    my $custom = "$t->{root}/sdm_custom.txt";
    for my $c (['4', 'ACGTCAGTGCTAGACG'], ['-1', 'ACGTCAGTGCTAGACG'], ['1.5', 'ACGTCAGTGCTAGACG'], ['', 'ACGTCAGTGCTAGACG'],
               ['1 #comment', 'ACGTCAGTGCTAGACG'], ['1', 'ACGTNCGT'], ['1', 'A']) {
        my ($errors, $barcode) = @$c;
        subtest "errors='$errors' barcode=$barcode" => sub {
            write_text($custom, replace_all(read_text("$ROOT/configs/sdm_ONT_savont.txt"), "maxBarcodeErrs\t1", "maxBarcodeErrs\t$errors"));
            write_text($t->{map}, replace_all($original, 'ACGTCAGTGCTAGACG', $barcode));
            my $result = $t->run_lotus(extra => ['-i', "$t->{root}/pooled.fastq.gz", '-s', $custom, '--dry-run'], ok => 0);
            contains($result->{output}, 'maxBarcodeErrs');
            ok(!-e $t->{out}, 'no output created');
        };
    }
};

# Reads for two Savont ASVs with distinct counts; returns the counts keyed by ASV sequence.
sub two_asv_reads {
    my ($t) = @_;
    (my $second = $t->{seq}) =~ tr/ACGT/TGCA/;
    my $consensus2 = (substr($second, 0, 1) ne 'T' ? 'T' : 'A') . substr($second, 1);
    $t->{env}{ONT_TEST_CONSENSUSES} = json_encode([$t->{consensus}, $consensus2]);
    for my $c (['s1', [4, 2]], ['s2', [3, 5]]) {
        my ($sample, $counts) = @$c;
        my $fq = '';
        for my $j (0, 1) {
            my $seq = $t->{fwd} . ($t->{seq}, $second)[$j] . $t->{revcomp};
            $fq .= "\@${sample}_${j}_$_\n$seq\n+\n" . ('I' x length $seq) . "\n" for 0 .. $counts->[$j] - 1;
        }
        write_text("$t->{reads}/$sample.fq", $fq);
    }
    return { $t->{consensus} => { s1 => 4, s2 => 3 }, $consensus2 => { s1 => 2, s2 => 5 } };
}

sub check_counts_by_sequence {
    my ($t, $expected) = @_;
    local $Test::Builder::Level = $Test::Builder::Level + 1;
    my $state = json_decode(read_text("$t->{out}/ont_test_state.json"));
    my %seeds;
    for my $entry (split />/, read_text($state->{seed})) {
        next if $entry eq '';
        my @lines = text_lines($entry);
        $seeds{ $lines[0] } = join '', @lines[1 .. $#lines];
    }
    my @table = text_lines(read_text("$t->{out}/OTU.txt"));
    is_deeply([sort values %seeds], [sort keys %$expected], 'both consensuses kept as seeds');
    my @head = split /\t/, $table[0];
    for my $row (@table[1 .. $#table]) {
        my ($key, @counts) = split /\t/, $row;
        my %got; @got{ @head[1 .. $#head] } = map { 0 + $_ } @counts;
        is_deeply(\%got, $expected->{ $seeds{$key} // '' }, "counts for $key");
    }
}

# Real SDM, then edit the seed FASTA written by the Savont counting call (-otu_matrix).
my $SEED_EDIT_SDM = <<'PERL';
#!/usr/bin/env perl
use strict; use warnings;
system($ENV{ONT_TEST_REAL_SDM}, @ARGV);
my $status = $? == -1 ? 127 : $? >> 8;
if ($status == 0 && grep { $_ eq '-otu_matrix' } @ARGV) {
    my ($i) = grep { $ARGV[$_] eq '-o_fna' } 0 .. $#ARGV;
    my $file = $ARGV[$i + 1];
    open my $in, '<', $file or die "$file: $!\n";
    my @records = grep { length } split />/, do { local $/; <$in> };
    close $in;
    if ($ENV{ONT_TEST_SEED_EDIT} eq 'reverse') { @records = reverse @records }
    else { $records[0] =~ s/\n(.)/"\n" . ($1 eq 'A' ? 'C' : 'A')/e }  # the first seed no longer equals any ASV
    open my $out, '>', $file or die "$file: $!\n";
    print {$out} map { ">$_" } @records;
    close $out;
}
exit $status;
PERL

sub seed_editing_sdm {
    my ($t, $mode) = @_;
    my $wrapper = write_text("$t->{tools}/sdm", $SEED_EDIT_SDM);
    chmod 0755, $wrapper or die "chmod $wrapper: $!\n";
    @{ $t->{env} }{qw(ONT_TEST_REAL_SDM ONT_TEST_SEED_EDIT)} = ($SDM, $mode);
    write_text($t->{cfg}, replace_all(read_text($t->{cfg}), "sdm $SDM\n", "sdm $wrapper\n"));
}

{ no strict 'refs'; *{"LotusTest::$_"} = \&{"main::$_"} for qw(two_asv_reads check_counts_by_sequence seed_editing_sdm); }

case test_multiple_savont_consensuses_keep_correct_counts => sub {
    my $t = shift;
    my $expected = $t->two_asv_reads;
    $t->run_lotus;
    $t->check_counts_by_sequence($expected);
};

case test_savont_pairing_does_not_depend_on_sdm_seed_order => sub {
    my $t = shift;
    my $expected = $t->two_asv_reads;
    $t->seed_editing_sdm('reverse');
    $t->run_lotus;
    $t->check_counts_by_sequence($expected);
};

case test_savont_seed_that_is_not_an_asv_stops_the_run => sub {
    my $t = shift;
    $t->two_asv_reads;
    $t->seed_editing_sdm('mutate');
    contains($t->run_lotus(ok => 0)->{output}, 'is not one of the Savont ASV sequences');
    ok(!-e "$t->{out}/ont_test_state.json", 'stopped before the checkpoint');
};

case test_ont_backmapping_identity => sub {
    my $t = shift;
    my $lssu = ['-CL', 'vsearch', '-derepMin', '0', '-s', "$ROOT/configs/sdm_ONT_LSSU.txt"];
    my $n = 0;
    for my $c ([[], '95', 0], [$lssu, '95', 0], [['-backmap_id', '0.97'], '97', 0], [['-backmap_id', '0.9'], '90', 1]) {
        my ($extra, $identity, $warned) = @$c;
        subtest "extra=@$extra" => sub {
            $t->{out} = "$t->{root}/output_" . $n++;
            my $result = $t->run_lotus(extra => $extra);
            $t->check_counts;
            contains($t->sdm_commands->[1], " -id $identity -minQueryCov 0.8 ");
            is(index($result->{output}, 'lower than OTU clustering threshhold') >= 0 ? 1 : 0, $warned, 'backmap-below-id warning');
        };
    }
};

case test_savont_explicit_options_and_vsearch_mapping => sub {
    my $t = shift;
    $t->run_lotus(extra => ['-CL', '9', '-useMini4map', '0', '-saveDemultiplex', '2', '-savontSingleStrand', '1',
        '-savontQualCutoff', '95.5', '-savontMinBaseQual', '20', '-savontChimeraErrors', '2']);
    $t->check_counts;
    my $args = $t->tool_calls('savont')->[0];
    is(after($args, $_->[0]), $_->[1], $_->[0])
        for ['--quality-value-cutoff', '95.5'], ['--minimum-base-quality', '20'], ['--chimera-allowable-errors', '2'];
    ok(has_arg($args, '--single-strand'), 'single-strand mode');
    ok(scalar(() = glob("$t->{out}/demultiplexed/*.fq.gz")), 'demultiplexed FASTQ saved compressed');
    is_deeply([glob("$t->{out}/demultiplexed/*.fq")], [], 'no uncompressed copies');
    my $mapper = $t->tool_calls('vsearch')->[0];
    my $mapping_reads = after($mapper, '--usearch_global');
    is(canon($mapping_reads), canon("$t->{out}/tmpFiles/savont_mapping.fna"), 'VSEARCH maps converted reads');
    my @fasta = text_lines(read_text($mapping_reads));
    my @fastq = text_lines(read_text("$t->{out}/tmpFiles/savont_reads.fq"));
    is_deeply([map { $fasta[$_] } grep { $_ % 2 == 0 } 0 .. $#fasta], [map { '>' . substr($_, 1) } @{ every4(\@fastq, 0) }], 'headers');
    is_deeply([map { $fasta[$_] } grep { $_ % 2 == 1 } 0 .. $#fasta], every4(\@fastq, 1), 'sequences');
    is(int(@fasta / 2), 7, 'every read mapped');
};

case test_barbell_drops_samples_and_rewrites_map => sub {
    my $t = shift;
    $t->write_map(barbell => 1);
    my $original = read_text($t->{map});
    $t->run_lotus(extra => ['-ontMinReads', '2'], barbell => 1);
    ok(!has_arg($t->tool_calls('barbell')->[0], '--maximize'), 'conservative Barbell');
    my $state = $t->check_counts;
    is_deeply([sort keys %{ $state->{map} }], ['#SampleID', 's1', 's2'], 'low/missing samples dropped');
    is_deeply($state->{combined}, { s1 => 's1', s2 => 's2' }, 'combined samples');
    my $working_map = read_text("$t->{out}/primary/in.map");
    contains($working_map, 'fastqFile');
    lacks($working_map, 'missing');
    is(read_text($t->{map}), $original, 'user map unchanged');
    is_deeply(entries("$t->{out}/tmpFiles/ont_demux"), ['s1.fq', 's2.fq'], 'demultiplexed files');
};

case test_barbell_maximize_is_opt_in => sub {
    my $t = shift;
    $t->write_map(barbell => 1);
    $t->run_lotus(extra => ['-ontMinReads', '2', '-ontBarbellMaximize', '1'], barbell => 1);
    $t->check_counts;
    ok(has_arg($t->tool_calls('barbell')->[0], '--maximize'), 'Barbell --maximize');
};

case test_primary_ont_path_does_not_need_barbell => sub {
    my $t = shift;
    unlink "$t->{tools}/barbell" or die "$!\n";
    $t->run_lotus;
    $t->check_counts;
    my $citations = read_text("$t->{out}/LotuSLogS/citations.txt");
    contains($citations, '10.64898/2026.05.26.727271');
    lacks($citations, '10.1093/bioinformatics/btag349');
};

case test_removed_primers_with_barbell => sub {
    my $t = shift;
    $t->write_map(barbell => 1);
    $t->check_removed_primers(barbell => 1);
};

case test_removed_primers_preserves_custom_preset => sub {
    my $t = shift;
    my $custom = write_text("$t->{root}/sdm_custom.txt", replace_all(savont_opts(), "minAvgQuality\t15", "minAvgQuality\t25"));
    my $before = read_text($custom);
    $t->check_removed_primers(extra => ['-s', $custom]);
    is(read_text($custom), $before, 'custom preset unchanged');
    my $effective = read_text("$t->{out}/primary/sdm_ONT_removed.txt");
    contains($effective, $_) for "minAvgQuality\t25", "RejectSeqWithoutFwdPrim\tF", "RejectSeqWithoutRevPrim\tF", "ExtensivePrimerChecks\tF";
};

case test_removed_primers_without_map_primer_columns => sub {
    my $t = shift;
    $t->strip_fixture_primers($t->{reads});
    write_text($t->{map}, "#SampleID\tfastqFile\ns1\ts1.fq\ns2\ts2.fq\n");
    my $result = $t->run_lotus(extra => ['-ontPrimerState', 'removed']);
    $t->check_counts;
    lacks($result->{output}, 'No forward PCR primer');
};

case test_barbell_demultiplex_only_cites_only_barbell => sub {
    my $t = shift;
    $t->write_map(barbell => 1);
    $t->run_lotus(extra => ['-saveDemultiplex', '1'], barbell => 1);
    my $citations = read_text("$t->{out}/LotuSLogS/citations.txt");
    is(count_of($citations, '10.1093/bioinformatics/btag349'), 1, 'Barbell cited once');
    lacks($citations, '10.64898/2026.05.26.727271');
    is_deeply($t->tool_calls('savont'), [], 'Savont not called');
};

case test_barbell_without_new_fastq_column => sub {
    my $t = shift;
    $t->write_map(barbell => 1);
    my $map = replace_first(read_text($t->{map}), "ONTBarcode\t", "ONTBarcode\tfastqFile\t");
    $map = replace_all($map, "BC0$_\t", "BC0$_\told$_.fq\t") for 1 .. 4;
    write_text($t->{map}, $map);
    $t->run_lotus(extra => ['-ontMinReads', '2', '-ontWriteFastqCol', '0'], barbell => 1);
    $t->check_counts;
    my $working_map = read_text("$t->{out}/primary/in.map");
    contains($working_map, 'fastqFile');
    lacks($working_map, 'old1.fq');
    lacks($working_map, 'missing');
};

case test_barbell_missing_fastq_column_in_mode_zero_rejected => sub {
    my $t = shift;
    $t->write_map(barbell => 1);
    my $result = $t->run_lotus(extra => ['-ontWriteFastqCol', '0', '--dry-run'], barbell => 1, ok => 0);
    contains($result->{output}, 'requires an existing fastqFile');
    is_deeply($t->tool_calls, [], 'no tools called');
};

case test_ont_custom_sdm_preset_is_respected => sub {
    my $t = shift;
    $t->run_lotus(extra => ['-s', "$ROOT/configs/sdm_ONT_LSSU.txt", '--dry-run']);
    is_deeply($t->tool_calls, [], 'no tools called');
};

case test_missing_savont_is_reported => sub {
    my $t = shift;
    unlink "$t->{tools}/savont" or die "$!\n";
    my $result = $t->run_lotus(extra => ['--dry-run'], ok => 0);
    contains($result->{output}, 'No valid savont binary');
    ok(!-e $t->{out}, 'no output created');
};

case test_missing_barbell_is_reported => sub {
    my $t = shift;
    $t->write_map(barbell => 1);
    unlink "$t->{tools}/barbell" or die "$!\n";
    my $result = $t->run_lotus(extra => ['--dry-run'], barbell => 1, ok => 0);
    contains($result->{output}, 'barbell executable unavailable');
    is_deeply($t->tool_calls, [], 'no tools called');
};

case test_paired_reads_are_rejected => sub {
    my $t = shift;
    my $map = replace_all(read_text($t->{map}), 's1.fq', 's1.fq,s2.fq');
    write_text($t->{map}, replace_all($map, "s2\ts2.fq", "s2\ts2.fq,s1.fq"));
    my $result = $t->run_lotus(extra => ['--dry-run'], ok => 0);
    contains($result->{output}, 'single-end FASTQ');
};

case test_savont_accepts_highmem_zero => sub {
    my $t = shift;
    $t->run_lotus(extra => ['-highmem', '0']);
    $t->check_counts;
    is_deeply([glob("$t->{out}/tmpFiles/derep*")], [], 'no derep files');
};

case test_dry_run_is_non_destructive => sub {
    my $t = shift;
    $t->write_map(barbell => 1);
    mkdir $t->{out} or die "$!\n";
    write_text("$t->{out}/LotuS_output_schema_version.txt", 'keep');
    my $report = write_text("$t->{root}/report.txt", 'keep report');
    my $manifest = write_text("$t->{root}/manifest.txt", 'keep manifest');
    my $result = $t->run_lotus(extra => ['--dry-run', '-ontPrimerState', 'removed', '-ontBarbellMaximize', '1',
        '--dependencyReport', $report, '--manifest', $manifest], barbell => 1);
    contains($result->{output}, 'Dry-run validation complete');
    is_deeply(entries($t->{out}), ['LotuS_output_schema_version.txt'], 'output untouched');
    is(read_text($report), 'keep report', 'report untouched');
    is(read_text($manifest), 'keep manifest', 'manifest untouched');
    is_deeply($t->tool_calls, [], 'no tools called');
};

case test_barbell_duplicate_barcode_rejected_in_dry_run => sub {
    my $t = shift;
    $t->write_map(barbell => 1, rows => [['s1', 'BC01'], ['s2', 'BC01']]);
    my $result = $t->run_lotus(extra => ['--dry-run'], barbell => 1, ok => 0);
    contains($result->{output}, 'assigned to more than one sample');
    is_deeply($t->tool_calls, [], 'no tools called');
    ok(!-e $t->{out}, 'no output created');
};

case test_barbell_unsafe_sample_rejected => sub {
    my $t = shift;
    $t->write_map(barbell => 1, rows => [['../outside', 'BC01']]);
    my $result = $t->run_lotus(extra => ['--dry-run'], barbell => 1, ok => 0);
    contains($result->{output}, 'Unsafe ONT SampleID');
    is_deeply($t->tool_calls, [], 'no tools called');
};

case test_savont_failure_stops_before_backmapping => sub {
    my $t = shift;
    $t->{env}{ONT_TEST_FAIL} = '1';
    $t->run_lotus(ok => 0);
    is_deeply($t->tool_calls('minimap2'), [], 'no backmapping');
    ok(!-e "$t->{out}/OTU.txt", 'no OTU table');
    ok(-e "$t->{out}/tmpFiles/savont_reads.fq", 'prepared reads kept');
};

case test_empty_savont_output_stops_before_backmapping => sub {
    my $t = shift;
    $t->{env}{ONT_TEST_EMPTY} = '1';
    my $result = $t->run_lotus(ok => 0);
    contains($result->{output}, 'savont produced no ASVs');
    is_deeply($t->tool_calls('minimap2'), [], 'no backmapping');
};

case test_invalid_flags_rejected => sub {
    my $t = shift;
    for my $args (['-barcodeSearchWindow', '0'], ['-barcodeSearchWindow', '10001'], ['-barcodeSearchWindow', '2.5'],
        ['-ontBarcodeEnds', 'bad'], ['-p', 'miSeq', '-barcodeSearchWindow', '200'], ['-p', 'miSeq', '-ontBarcodeEnds', 'either'],
        ['-ontDemux', 'invalid'], ['-ontBarbellMaximize', '2'], ['-ontBarbellMaximize', '1'], ['-ontPrimerState', 'invalid'],
        ['-p', 'miSeq', '-ontPrimerState', 'removed'], ['-p', 'miSeq', '-ontDemux', 'barbell'], ['-p', 'miSeq', '-CL', 'savont'],
        ['-ontMinReads', '-1'], ['-ontWriteFastqCol', '2'], ['-savontSingleStrand', '2'], ['-savontQualCutoff', '101'],
        ['-savontMinBaseQual', '-1'], ['-savontChimeraErrors', '-1']) {
        subtest "args=@$args" => sub {
            $t->run_lotus(extra => $args, ok => 0);
            is_deeply($t->tool_calls, [], 'no tools called');
        };
    }
};

case test_miseq_dry_run_remains_supported => sub {
    my $t = shift;
    $t->run_lotus(extra => ['-p', 'miSeq', '-CL', 'vsearch', '--dry-run']);
    is_deeply($t->tool_calls, [], 'no tools called');
};

case test_removed_dnaclust_stays_removed => sub {
    my $t = shift;
    my $result = $t->run_lotus(extra => ['-CL', '4'], ok => 0);
    contains($result->{output}, 'has been removed');
    ok(!-e $t->{out}, 'no output created');
};

# savont2uc.pl converter
sub converter {
    my ($root, @extra) = @_;
    return run_command(\%ENV, 30, 'perl', "$ROOT/bin/savont2uc.pl", '--asvs', "$root/asvs.fa", '--map', "$root/map.tsv",
        '--ucout', "$root/out.uc", '--fnaout', "$root/out.fa", @extra);
}

subtest test_best_hit_and_survivor_filter => sub {
    my $root = tempdir('savont-converter-test-XXXXXXXX', TMPDIR => 1, CLEANUP => 1);
    write_text("$root/asvs.fa", ">final_0 debug_id:7\nACGT\n>final_1 debug_id:8\nTGCA\n");
    write_text("$root/map.tsv", "r1 description\tasv:7\t1\t90\nr1 description\tasv:8\t0\t80\nr2\tasv:7\t0\t90\nr3\tasv:7\t0\t90\nr4\tasv:99\t0\t100\n");
    my ($output, $status) = converter($root);
    is($status, 0, 'converter succeeds') or diag($output);
    is(read_text("$root/out.fa"), ">r2\nACGT\n>r1\nTGCA\n", 'representative FASTA');
    is(read_text("$root/out.uc"), "r2\totu1\t*\nr3\tmatch\tdqt=1;top=r2(99%);\nr1\totu2\t*\n", 'UC mapping');
};

subtest test_invalid_id_mode_and_duplicate_asv_rejected => sub {
    my $root = tempdir('savont-converter-test-XXXXXXXX', TMPDIR => 1, CLEANUP => 1);
    write_text("$root/asvs.fa", ">a debug_id:7\nACGT\n>b debug_id:7\nTGCA\n");
    write_text("$root/map.tsv", "r1\tasv:7\t0\t90\n");
    for my $c ([['--idmode', 'invalid'], '--idmode must be'], [[], 'Duplicate ASV debug_id']) {
        my ($extra, $message) = @$c;
        subtest "extra=@$extra" => sub {
            my ($output, $status) = converter($root, @$extra);
            isnt($status, 0, 'converter fails');
            contains($output, $message);
            ok(!-e "$root/out.uc", 'no UC output');
        };
    }
};

done_testing();

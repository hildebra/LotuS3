#!/usr/bin/env python3
"""Coarse dereplication preflight and FASTQ handoff regressions.

Reuse the ONT fixture and its controlled clusterer/mappers. Preprocessing and
native seed/count checks use the installed SDM. A separate capture wrapper
exercises failed-seed handling. Standard retained-quality HQ support is required.

LOTUS_TEST_SDM=/path/to/sdm python3 -m unittest discover -s tests -p coarse_derep.py -v
"""
from collections import Counter
import json
import hashlib
import os
from pathlib import Path
import re
import subprocess
import unittest
import ont_integration as ont
import perl_audit as audit
from audit_derep_counts import audit as audit_counts

VERSION = subprocess.run([str(ont.SDM), '-v'], text=True, capture_output=True)
version = re.search(r'sdm\s+(\d+)\.(\d+)', VERSION.stdout + VERSION.stderr)
HAS_COARSE = bool(version and tuple(map(int, version.groups())) >= (3, 52))
FLAGS = subprocess.run([str(ont.SDM), '-help_flags'], text=True, capture_output=True)
HAS_SEEDS = '-seedSubclusters' in FLAGS.stdout + FLAGS.stderr
HAS_STANDARD_HQ = HAS_SEEDS and 'Standard HQ output: .1.hq.fq; paired additionally .2.hq.fq' in FLAGS.stdout + FLAGS.stderr

WRAPPER = r'''#!/usr/bin/env python3
import json, os, pathlib, subprocess, sys
args = sys.argv[1:]
if args and args[0] in ('-v', '-version'):
    print(os.environ.get('COARSE_TEST_VERSION', 'sdm 3.52 beta'))
elif args == ['-help_flags']:
    print('-seedSubclusters <0|1>' if not os.environ.get('COARSE_TEST_NO_CAPABILITY') else '-derepStoreQuals <0|1>')
    if not os.environ.get('COARSE_TEST_LEGACY_HQ'):
        print('Standard HQ output: .1.hq.fq; paired additionally .2.hq.fq\n-derepCoarseClusters <0|1>')
elif '-optimalRead2Cluster' in args:
    pathlib.Path(os.environ['COARSE_TEST_SEED_CALL']).write_text(json.dumps(args))
    sys.exit(17)  # Intentional failure: never claim to emulate native seed selection.
else:
    result = subprocess.run([os.environ['COARSE_TEST_REAL_SDM']] + args)
    if result.returncode == 0 and os.environ.get('COARSE_TEST_REMOVE_MATE') and '-o_dereplicate' in args:
        base = pathlib.Path(args[args.index('-o_dereplicate')+1]).with_suffix('')
        pathlib.Path(str(base)+'.2.hq.fq').unlink()
    sys.exit(result.returncode)
'''


class CoarsePreflight(unittest.TestCase):
    setUp = ont.ONTIntegration.setUp
    write_map = ont.ONTIntegration.write_map
    run_lotus = ont.ONTIntegration.run_lotus
    probe = audit.PerlAudit.probe

    def capture_wrapper(self):
        wrapper = self.tools/'sdm'
        wrapper.write_text(WRAPPER); wrapper.chmod(0o755)
        self.cfg.write_text(self.cfg.read_text().replace(str(ont.SDM), str(wrapper)))
        self.seed_call = self.root/'seed_call.json'
        self.env.update(COARSE_TEST_REAL_SDM=str(ont.SDM), COARSE_TEST_SEED_CALL=str(self.seed_call))

    def test_invalid_ani_rejected_before_output_creation_even_without_strict(self):
        for value in ('0', '0.9499', '1.0001', '95', '-1', 'NaN', 'Inf', 'abc', '0.97;echo', '1e999', ''):
            with self.subTest(value=value):
                result = self.run_lotus(['-coarseDerep='+value, '--no-strict'], ok=False)
                self.assertIn('-coarseDerep must be a numeric ANI' if value else 'coarseDerep requires an argument', result.stdout)
                self.assertFalse(self.out.exists())

    def test_inclusive_ani_bounds_and_fraction_formats(self):
        for value in ('0.95', '.975', '9.9e-1', '1.0'):
            with self.subTest(value=value):
                result = self.run_lotus(['-coarseDerep', value, '-v'])
                self.assertIn('LotuS 3.', result.stdout)

    def test_old_sdm_rejected_before_processing(self):
        self.capture_wrapper()
        self.env['COARSE_TEST_VERSION'] = 'sdm 3.51 beta'
        result = self.run_lotus(['-CL', 'vsearch', '-coarseDerep', '0.97'], ok=False)
        self.assertIn('requires SDM >= 3.52', result.stdout)
        self.assertFalse(self.seed_call.exists())
        self.assertFalse((self.out/'tmpFiles').exists())

    def test_unpatched_352_rejected_before_processing(self):
        self.capture_wrapper()
        self.env['COARSE_TEST_NO_CAPABILITY'] = '1'
        result = self.run_lotus(['-CL', 'vsearch', '-coarseDerep', '0.97'], ok=False)
        self.assertIn('does not advertise -seedSubclusters', result.stdout)
        self.assertFalse(self.seed_call.exists())
        self.assertFalse((self.out/'tmpFiles').exists())

    def test_old_retained_quality_layout_rejected_despite_seed_flag(self):
        self.capture_wrapper()
        self.env['COARSE_TEST_LEGACY_HQ'] = '1'
        result = self.run_lotus(['-CL', 'vsearch', '-coarseDerep', '0.97'], ok=False)
        self.assertIn('does not advertise the standard retained-quality HQ layout', result.stdout)
        self.assertFalse(self.seed_call.exists())
        self.assertFalse((self.out/'tmpFiles').exists())

    def test_incompatible_modes_rejected_even_without_strict(self):
        cases = [
            ('$sdmDerepDo = 0;', 'full run with SDM dereplication'),
            ('$mergePreCluster = 1;', '-mergePreClusterReads 0'),
            ('$saveDemulti = 1;', 'demultiplex-only'),
            ('$TaxOnly = "1";', 'taxonomy-only'),
            ('$onlyTaxRedo = 1;', 'taxonomy-only'),
        ]
        for settings, diagnostic in cases:
            with self.subTest(settings=settings):
                result = self.probe(settings+' validate_option_combinations();', ['-coarseDerep', '0.97', '--no-strict'], ok=False)
                self.assertIn(diagnostic, result.stdout)

    def test_finalized_abundance_and_pair_bound_are_reported(self):
        result = self.probe(r'''print parse_sdm_short_report(
            "Maximum pairs with both reads accepted: 900 (90.0%) n/a\n",
            "Dereplication: 10 unique sequences (avg size 78; 780 counts, 10 merged)\n"
            . "Dereplication abundance: 1200 total in map; 780 passing; 420 in rest\n"
            . "Retention exclusions: 40 pairs with empty R2 skipped before admission\n");''')
        self.assertIn('Pairs accepted on both ends (upper bound): 900 (90.0%)', result.stdout)
        self.assertIn('10 passing unique sequences from 780 reads', result.stdout)
        self.assertIn('1,200 total in map; 780 passing (main + merged); 420 in rest', result.stdout)
        self.assertIn('Retention exclusions: 40 pairs with physically empty R2', result.stdout)

    def test_count_auditor_includes_merged_and_rest_and_rejects_missing_parents(self):
        main = self.root/'derep.fas'
        (self.root/'derep.map').write_text('#SMPLS\t0:s1\t1:s2\n'
            'p;size=4;\t0:3\t1:1\nm;size=2;\t0:1\t1:1\nr;size=3;\t1:3\n')
        (self.root/'derep.fas.rest').write_text('>r;size=3;\nACGT\n')
        for fastq in (False, True):
            with self.subTest(fastq=fastq):
                main.write_text('@p;size=4;\nACGT\n+\nIIII\n' if fastq else '>p;size=4;\nAC\nGT\n')
                merged = self.root/'derep.merg.fas'
                merged.write_text('@m;size=2;\nAGGT\n+\nIIII\n' if fastq else '>m;size=2;\nAGGT\n')
                counts = audit_counts(main)
                self.assertEqual((counts['map_total'], counts['passing'], counts['rest']), (9, 6, 3))
                self.assertEqual(counts['samples'], {'s1': 4, 's2': 5})
                merged.unlink()
                with self.assertRaisesRegex(AssertionError, 'map parents absent'):
                    audit_counts(main)

    def test_ordinary_preprocessing_merge_policy_is_preserved(self):
        result = self.probe(r'''$mergePreCluster = 1; $ClusterPipe = 7;
            print JSON::PP->new->canonical->encode(sdm_preprocessing_merge_options());''')
        self.assertIn('{"merge_pairs_demulti":1,"merge_pairs_derep":1}', result.stdout)



@unittest.skipUnless(HAS_COARSE and HAS_STANDARD_HQ, 'requires SDM with standard retained-quality HQ support')
class CoarseHandoff(unittest.TestCase):
    setUp = ont.ONTIntegration.setUp
    write_map = ont.ONTIntegration.write_map
    run_lotus = ont.ONTIntegration.run_lotus
    capture_wrapper = CoarsePreflight.capture_wrapper
    check_counts = ont.ONTIntegration.check_counts

    def inputs(self, paired, counts=(4, 3)):
        self.map.write_text('#SampleID\tfastqFile\n' + ''.join(
            f'{s}\t{s}.1.fq' + (f',{s}.2.fq' if paired else '') + '\n' for s in ('s1', 's2')))
        for sample, count in zip(('s1', 's2'), counts):
            for mate in range(1, 3 if paired else 2):
                records = []
                for i in range(count):
                    seq = self.seq
                    if sample == 's1' and i == 0 and mate == (2 if paired else 1):
                        seq = ('A' if seq[0] != 'A' else 'T') + seq[1:]
                    records.append(f'@{sample}_{i}/{mate}\n{seq}\n+\n'+'I'*len(seq)+'\n')
                (self.reads/f'{sample}.{mate}.fq').write_text(''.join(records))
        options = self.root/'options.txt'
        options.write_text('minSeqLength\t100\nmaxSeqLength\t2000\nminAvgQuality\t0\n'
                           'RejectSeqWithoutFwdPrim\tF\nRejectSeqWithoutRevPrim\tF\n'
                           'TrimWindowThreshhold\t0\nmaxHomonucleotide\t100\n')
        return ['-p', 'miSeq', '-CL', 'vsearch', '-s', str(options), '-derepMin', '1', '-sdmThreads', '1']

    def check_handoff(self, paired, identity):
        self.capture_wrapper()
        result = self.run_lotus(self.inputs(paired)+['-coarseDerep', identity], ok=False)
        self.assertIn('SDM retained-variant seed extension failed', result.stdout)
        self.assertNotIn('Fallback to', result.stdout)
        args = json.loads(self.seed_call.read_text())
        self.assertEqual(args[args.index('-seedSubclusters')+1], '1')
        if paired:
            self.assertEqual(args[args.index('-merge_pairs_seed')+1], '1')
        else:
            self.assertNotIn('-merge_pairs_seed', args)
        self.assertEqual(args[args.index('-i_qual_offset')+1], '33')
        files = [Path(x) for x in args[args.index('-i_fastq')+1].split(',')]
        self.assertEqual([f.name for f in files], ['derep.1.hq.fq', 'derep.2.hq.fq'] if paired else ['derep.1.hq.fq'])
        records = [f.read_text().splitlines() for f in files]
        self.assertEqual(len(records[0])//4, 2)
        self.assertEqual(sum(int(re.search(r';size=(\d+);', h).group(1)) for h in records[0][::4]), 7)
        if paired:
            self.assertEqual(records[0][::4], records[1][::4])
            self.assertEqual(len(set(records[0][1::4])), 1)
            self.assertEqual(len(set(records[1][1::4])), 2)
        commands = (self.out/'LotuSLogS/LotuS_cmds.log').read_text()
        self.assertIn('-derepIdentity '+format(float(identity)*100, '.12g'), commands)
        for flag in ('-derepStoreQuals 1', '-derepStoreDiffs 0', '-derepSubclusterFasta 0', '-derepReassign 0', '-derepCoarseClusters 0', '-merge_pairs_derep 0', '-merge_pairs_filter 0', '-merge_pairs_demulti 0'):
            self.assertIn(flag, commands)
        self.assertEqual(Path(args[args.index('-derep_map')+1]).name, 'derep.map')
        self.assertEqual(list((self.out/'tmpFiles').glob('*.subclusters*.fq')), [])
        self.assertEqual(list((self.out/'tmpFiles').glob('*.diff')), [])

    def test_single_end_reconstructed_fastq_and_percent_conversion(self):
        self.check_handoff(False, '0.975')

    def test_paired_reconstructed_fastqs_preserve_r2_only_variants(self):
        self.check_handoff(True, '0.95')

    def test_explicit_identity_one_still_exports_subclusters(self):
        self.check_handoff(True, '1.0')

    def test_missing_reconstructed_mate_fails_before_seed_command(self):
        self.capture_wrapper()
        self.env['COARSE_TEST_REMOVE_MATE'] = '1'
        result = self.run_lotus(self.inputs(True)+['-coarseDerep', '0.97'], ok=False)
        self.assertIn('Missing or empty retained-variant seed FASTQ', result.stdout)
        self.assertFalse(self.seed_call.exists())
        self.assertFalse((self.out/'primary/sdm_dereplication.json').exists())

    def test_default_hq_path_still_counts_without_coarse_flags(self):
        self.run_lotus(self.inputs(False))
        self.check_counts({'s1': 4, 's2': 3})
        commands = (self.out/'LotuSLogS/LotuS_cmds.log').read_text()
        self.assertIn('derep.1.hq.fq', commands)
        self.assertNotIn('-seedSubclusters', commands)
        self.assertNotIn('-derepIdentity', commands)
        self.assertNotIn('-derepStoreQuals', commands)
        self.metadata(retained=0)

    def metadata(self, retained=1, granularity='exact search keys with configured prefix consolidation'):
        metadata = json.loads((self.out/'primary/sdm_dereplication.json').read_text())
        self.assertEqual(metadata['status'], 'complete')
        self.assertEqual(metadata['contract'], 'standard_hq_v2')
        self.assertEqual(metadata['quality_retention'], retained)
        self.assertEqual(metadata['output_granularity'], granularity)
        self.assertEqual(metadata['hq_record_layout'], 'exact variants' if retained else 'representatives')
        self.assertEqual(metadata['sdm_binary_sha256'], hashlib.sha256(ont.SDM.read_bytes()).hexdigest())
        self.assertEqual(metadata['sdm_options_sha256'], hashlib.sha256((self.root/'options.txt').read_bytes()).hexdigest())
        self.assertEqual(metadata['noSearchWithMerge'], {'present': False, 'value': None})
        self.assertEqual(metadata['preprocessing_merge_flags'], {f'merge_pairs_{stage}': 0 for stage in ('derep', 'filter', 'demulti')})
        if not metadata['derep_per_sequencing_run']:
            counts = audit_counts(self.out/'tmpFiles/derep.fas')
            self.assertEqual(counts['map_total'], counts['passing'] + counts['rest'])
        manifest = (self.out/'LotuSLogS/run_manifest.txt').read_text()
        self.assertIn('SDM quality retention: '+str(retained), manifest)
        self.assertIn('SDM dereplication output: '+granularity, manifest)
        return metadata

    def variant_counts(self, paired):
        base = self.out/'tmpFiles'
        parents = {}
        for line in (base/'derep.map').read_text().splitlines():
            if line.startswith('#'): continue
            header, *samples = line.split('\t')
            count = int(re.search(r';size=(\d+);', header).group(1))
            self.assertEqual(count, sum(int(x.split(':')[1]) for x in samples))
            parents[header.split(';size=')[0]] = count
        r1 = (base/'derep.1.hq.fq').read_text().splitlines()
        r2 = (base/'derep.2.hq.fq').read_text().splitlines() if paired else ['']*len(r1)
        if paired: self.assertEqual(r1[::4], r2[::4])
        sums = Counter()
        records = {}
        for i in range(0, len(r1), 4):
            header = r1[i][1:]
            parent, number = header.split(';size=')[0].rsplit('.sub', 1)
            self.assertTrue(number.isdigit())
            count = int(re.search(r';size=(\d+);', header).group(1))
            sums[parent] += count
            key = (r1[i+1], r2[i+1])
            old = records.get(key, (0, r1[i+3], r2[i+3]))
            records[key] = (old[0]+count, r1[i+3], r2[i+3])
        self.assertEqual(dict(sums), parents)
        self.assertEqual(list(base.glob('*.subclusters*.fq')), [])
        return records

    def test_options_file_retention_selects_variant_reader_at_100_percent(self):
        extra = self.inputs(True)
        with (self.root/'options.txt').open('a') as out:
            out.write('derepStoreQuals\t0\nderepStoreQuals\t1\nderepIdentity\t100\n')
        self.run_lotus(extra)
        self.check_counts()
        self.variant_counts(True)
        self.metadata()
        commands = (self.out/'LotuSLogS/LotuS_cmds.log').read_text()
        self.assertIn('-seedSubclusters 1', commands)
        self.assertNotIn('-derepIdentity ', commands)

    def test_coarse_flag_overrides_optional_exports_and_coarse_parent_output(self):
        extra = self.inputs(False)
        with (self.root/'options.txt').open('a') as out:
            out.write('derepStoreQuals\t0\nderepStoreDiffs\t1\nderepSubclusterFasta\t1\nderepReassign\t1\nderepCoarseClusters\t1\n')
        self.run_lotus(extra+['-coarseDerep', '0.97'])
        self.check_counts()
        self.metadata()
        self.variant_counts(False)
        base = self.out/'tmpFiles'
        self.assertEqual((base/'derep.fas').read_text().count('>'), 2)
        self.assertEqual(list(base.glob('*.diff')), [])
        self.assertEqual(list(base.glob('*.subclusters*.fna')), [])

    def test_below_cutoff_parent_is_recovered_once(self):
        self.run_lotus(self.inputs(False)+['-coarseDerep', '0.97', '-derepMin', '2'])
        self.check_counts()
        self.variant_counts(False)
        self.assertEqual((self.out/'tmpFiles/derep.fas').read_text().count('>'), 1)
        self.assertEqual((self.out/'tmpFiles/derep.fas.rest').read_text().count('>'), 1)

    def test_exact_output_parity_at_97_and_100_with_multiple_workers(self):
        snapshots = []
        for workers in (1, 4, 12):
            for identity in ('0.97', '1.0'):
                with self.subTest(workers=workers, identity=identity):
                    self.out = self.root/f'parity_{workers}_{identity}'
                    self.run_lotus(self.inputs(True, (600, 400))+['-coarseDerep', identity, '-sdmThreads', str(workers)])
                    self.check_counts({'s1': 600, 's2': 400})
                    snapshots.append(self.variant_counts(True))
                    parents = (self.out/'tmpFiles/derep.fas').read_text().splitlines()
                    self.assertEqual(parents[1:], [self.seq])
                    self.assertIn(';size=1000;', parents[0])
        for snapshot in snapshots[1:]: self.assertEqual(snapshot, snapshots[0])

    def test_mixed_quality_admission_and_deferred_recovery(self):
        for identity in ('0.97', '1.0'):
            for workers in (1, 4):
                with self.subTest(identity=identity, workers=workers):
                    self.out = self.root/f'mixed_{identity}_{workers}'
                    extra = self.inputs(True)
                    # One parent has only R1-passing observations. Another has
                    # failed observations before its passing anchor in sample 2.
                    variant = ('T' if self.seq[0] != 'T' else 'A') + self.seq[1:]
                    observations = {'s1': [(variant, 'I', '+')]*4 + [(self.seq, '+', '+')]*2,
                                    's2': [(self.seq, 'I', 'I')]*3}
                    for sample, rows in observations.items():
                        for mate in (1, 2):
                            (self.reads/f'{sample}.{mate}.fq').write_text(''.join(
                                f'@{sample}_{i}/{mate}\n{seq}\n+\n'+quality[mate-1]*len(seq)+'\n'
                                for i, (seq, *quality) in enumerate(rows)))
                    with (self.root/'options.txt').open('a') as out:
                        out.write('minAvgQuality\t27\nderepSrchLen\t-1\nderepPrefix\t1\nfastqVersion\t1\n')
                    result = self.run_lotus(extra+['-coarseDerep', identity, '-sdmThreads', str(workers)])
                    self.check_counts({'s1': 6, 's2': 3})
                    self.variant_counts(True)
                    self.metadata()
                    counts = audit_counts(self.out/'tmpFiles/derep.fas')
                    self.assertEqual(counts['map_total'], 9)
                    self.assertIn('9 total in map; 9 passing (main + merged); 0 in rest', ' '.join(result.stdout.split()))

    def test_later_better_variant_becomes_full_length_seed(self):
        extra = self.inputs(False)
        better = self.seq[:300]+('A' if self.seq[300] != 'A' else 'T')+self.seq[301:]
        for sample, count, sequence, quality in [('s1', 4, self.seq, '5'), ('s2', 3, better, 'I')]:
            (self.reads/f'{sample}.1.fq').write_text(''.join(f'@{sample}_{i}\n{sequence}\n+\n'+quality*len(sequence)+'\n' for i in range(count)))
        with (self.root/'options.txt').open('a') as out:
            out.write('TruncateSequenceLength\t200\nfastqVersion\t1\n')
        self.env['ONT_TEST_CONSENSUS'] = self.seq[:200]
        self.run_lotus(extra+['-coarseDerep', '0.97'])
        state = self.check_counts()
        records = self.variant_counts(False)
        self.assertEqual(records[(self.seq, '')][1], '5'*len(self.seq))
        self.assertEqual(records[(better, '')][1], 'I'*len(better))
        self.assertEqual(len((self.out/'tmpFiles/derep.fas').read_text().splitlines()[1]), 200)
        hq = (self.out/'tmpFiles/derep.1.hq.fq').read_text().splitlines()
        self.assertEqual(hq[5], better)  # The second exported variant wins.
        self.assertEqual(''.join(Path(state['seed']).read_text().splitlines()[1:]), better)

    def test_metadata_records_effective_search_and_cut_settings(self):
        extra = self.inputs(False)
        with (self.root/'options.txt').open('a') as out:
            out.write('derepSrchLen\t150\nderepSrchLen\t200\nderepPrefix\tauto\n'
                      'TruncateSequenceLength\t250\nkeepBarcodeSeq\t0\nTrimStartNTs\t5\n')
        self.env['ONT_TEST_CONSENSUS'] = self.seq[5:205]
        self.run_lotus(extra+['-coarseDerep', '0.97'])
        metadata = self.metadata()
        self.assertEqual(metadata['search_source'], 'R1')
        self.assertEqual(metadata['search_options']['derepSrchLen'], '200')
        self.assertEqual(metadata['search_options']['TruncateSequenceLength'], '250')
        self.assertEqual(metadata['sdm_options']['keepBarcodeSeq'], '0')
        self.assertEqual(metadata['sdm_options']['TrimStartNTs'], '5')
        self.assertIn('logical tails retained', metadata['hq_sequence_policy'])

    def test_paired_hq_keeps_tails_and_removes_technical_sequences_once(self):
        primers = ('AGGTCAGTACCGTAAC', 'TCCAGATGCTACGTCA')
        barcodes = (('GACTGA', 'CTGTAC'), ('TACCTG', 'AGTCGA'))
        for retained in (False, True):
            with self.subTest(retained=retained):
                self.out = self.root/f'full_pair_{retained}'
                extra = self.inputs(True)
                self.map.write_text('#SampleID\tfastqFile\tBarcodeSequence\tBarcode2ndPair\tLinkerPrimerSequence\tReversePrimer\n'
                    + ''.join(f'{sample}\t{sample}.1.fq,{sample}.2.fq\t'
                        + '\t'.join((*barcodes[i], *primers))+'\n' for i,sample in enumerate(('s1', 's2'))))
                for i, (sample, count) in enumerate((('s1', 4), ('s2', 3))):
                    for mate in range(2):
                        prefix = barcodes[i][mate] + primers[mate]
                        sequence = prefix + self.seq
                        quality = 'I'*(len(sequence)-20)+'-'*20
                        (self.reads/f'{sample}.{mate+1}.fq').write_text(''.join(
                            f'@{sample}_{j}/{mate+1}\n{sequence}\n+\n{quality}\n' for j in range(count)))
                with (self.root/'options.txt').open('a') as out:
                    out.write('TruncateSequenceLength\t200\nTrimWindowThreshhold\t25\nTrimWindowWidth\t10\n'
                              'keepBarcodeSeq\t0\nkeepPrimerSeq\t0\nfastqVersion\t1\n')
                expected = self.seq
                self.env['ONT_TEST_CONSENSUS'] = expected[:200]
                self.run_lotus(extra+(['-coarseDerep', '0.97'] if retained else []))
                self.check_counts()
                for mate in (1, 2):
                    hq = (self.out/f'tmpFiles/derep.{mate}.hq.fq').read_text().splitlines()
                    self.assertEqual(hq[1::4], [expected])
                    self.assertEqual(hq[3::4], ['I'*(len(expected)-20)+'-'*20])
                seeds = (self.out/'tmpFiles/otu_seeds.merg.fq').read_text().splitlines()
                self.assertEqual(seeds[1::4], [expected])
                self.assertIn('1 were merged paired reads', (self.out/'LotuSLogS/SeedExtensionStats.log').read_text())
                self.metadata(retained=int(retained))
                if retained: self.variant_counts(True)

    def test_dada2_sequencing_runs_keep_cumulative_variant_hq(self):
        extra = self.inputs(True)
        lines = self.map.read_text().splitlines()
        self.map.write_text(lines[0]+'\tSequencingRun\n'+lines[1]+'\trunA\n'+lines[2]+'\trunB\n')
        rscript = self.tools/'Rscript'
        rscript.write_text("""#!/usr/bin/env python3
import os, pathlib, sys
args = sys.argv[1:]
assert args[0] == '--vanilla'
out = pathlib.Path(args[3])
files = list(out.glob('derep.*.fas'))
assert len(files) == 2, files
assert all(f.read_text().startswith('@') for f in files)
(out/'dada2.uc').write_text('')
(out/'dada2_p1_errF.pdf').write_text('test')
(out/'uniqueSeqs.fna').write_text('>OTU0\\n'+os.environ['ONT_TEST_CONSENSUS']+'\\n')
""")
        rscript.chmod(0o755)
        self.env['PATH'] = str(self.tools)+os.pathsep+self.env['PATH']
        self.cfg.write_text(self.cfg.read_text()+f'dada2R {ont.ROOT}/bin/R/dada2_pip_v2.R\n')
        self.run_lotus(extra+['-CL', 'dada2', '-coarseDerep', '0.97'])
        self.check_counts()
        self.variant_counts(True)
        self.assertEqual(self.metadata()['derep_per_sequencing_run'], 1)
        commands = (self.out/'LotuSLogS/LotuS_cmds.log').read_text()
        self.assertIn('-derep_format fq -derepPerSR 1', commands)
        self.assertNotIn('-derepPerSR 0', commands)
        self.assertEqual(list((self.out/'tmpFiles').glob('derep.*.1.hq.fq')), [])

    def test_fresh_run_regenerates_older_coarse_outputs(self):
        audit.PerlAudit.old_output(self)
        (self.out/'tmpFiles').mkdir()
        (self.out/'tmpFiles/derep.1.hq.fq').write_text('@old_parent;size=99;\nACGT\n+\nIIII\n')
        (self.out/'tmpFiles/derep.subclusters.fq').write_text('obsolete layout\n')
        (self.out/'primary').mkdir()
        (self.out/'primary/sdm_dereplication.json').write_text('{"contract":"old_coarse_parents"}\n')
        self.run_lotus(self.inputs(False)+['-coarseDerep', '0.97'])
        self.check_counts()
        self.variant_counts(False)
        self.metadata()
        self.assertNotIn('old_parent', (self.out/'tmpFiles/derep.1.hq.fq').read_text())

    def test_native_seed_extension_preserves_sample_counts(self):
        for paired in (False, True):
            with self.subTest(paired=paired):
                self.out = self.root/f'native_{paired}'
                self.run_lotus(self.inputs(paired)+['-coarseDerep', '0.97'])
                self.check_counts({'s1': 4, 's2': 3})
                self.variant_counts(paired)
                self.metadata()


if __name__ == '__main__':
    unittest.main()

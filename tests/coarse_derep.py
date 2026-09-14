#!/usr/bin/env python3
"""Coarse dereplication preflight and FASTQ handoff regressions.

Reuse the ONT fixture and its controlled clusterer/mappers. With SDM >= 3.52,
preprocessing is real; the capture wrapper deliberately fails seed extension,
whose required upstream support is tracked in the worker brief. Once that fix
is installed, the native seed/count regression runs as well.

LOTUS_TEST_SDM=/path/to/sdm python3 -m unittest discover -s tests -p coarse_derep.py -v
"""
import json
import os
from pathlib import Path
import re
import subprocess
import unittest
import ont_integration as ont
import perl_audit as audit

VERSION = subprocess.run([str(ont.SDM), '-v'], text=True, capture_output=True)
version = re.search(r'sdm\s+(\d+)\.(\d+)', VERSION.stdout + VERSION.stderr)
HAS_COARSE = bool(version and tuple(map(int, version.groups())) >= (3, 52))
FLAGS = subprocess.run([str(ont.SDM), '-help_flags'], text=True, capture_output=True)
HAS_SEEDS = '-seedSubclusters' in FLAGS.stdout + FLAGS.stderr

WRAPPER = r'''#!/usr/bin/env python3
import json, os, pathlib, subprocess, sys
args = sys.argv[1:]
if args and args[0] in ('-v', '-version'):
    print(os.environ.get('COARSE_TEST_VERSION', 'sdm 3.52 beta'))
elif args == ['-help_flags']:
    print('-seedSubclusters <0|1>' if not os.environ.get('COARSE_TEST_NO_CAPABILITY') else '-derepStoreQuals <0|1>')
elif '-optimalRead2Cluster' in args:
    pathlib.Path(os.environ['COARSE_TEST_SEED_CALL']).write_text(json.dumps(args))
    sys.exit(17)  # Intentional failure: never claim to emulate native seed selection.
else:
    result = subprocess.run([os.environ['COARSE_TEST_REAL_SDM']] + args)
    if result.returncode == 0 and os.environ.get('COARSE_TEST_REMOVE_MATE') and '-o_dereplicate' in args:
        base = pathlib.Path(args[args.index('-o_dereplicate')+1]).with_suffix('')
        pathlib.Path(str(base)+'.subclusters.2.fq').unlink()
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

    def test_incompatible_modes_rejected_even_without_strict(self):
        cases = [
            ('$sdmDerepDo = 0;', 'full run with SDM dereplication'),
            ('$ClusterPipe = 7;', 'DADA2'),
            ('$mergePreCluster = 1;', '-mergePreClusterReads 0'),
            ('$saveDemulti = 1;', 'demultiplex-only'),
            ('$TaxOnly = "1";', 'taxonomy-only'),
            ('$onlyTaxRedo = 1;', 'taxonomy-only'),
        ]
        for settings, diagnostic in cases:
            with self.subTest(settings=settings):
                result = self.probe(settings+' validate_option_combinations();', ['-coarseDerep', '0.97', '--no-strict'], ok=False)
                self.assertIn(diagnostic, result.stdout)


@unittest.skipUnless(HAS_COARSE, 'requires SDM >= 3.52 via LOTUS_TEST_SDM')
class CoarseHandoff(unittest.TestCase):
    setUp = ont.ONTIntegration.setUp
    write_map = ont.ONTIntegration.write_map
    run_lotus = ont.ONTIntegration.run_lotus
    capture_wrapper = CoarsePreflight.capture_wrapper
    check_counts = ont.ONTIntegration.check_counts

    def inputs(self, paired):
        self.map.write_text('#SampleID\tfastqFile\n' + ''.join(
            f'{s}\t{s}.1.fq' + (f',{s}.2.fq' if paired else '') + '\n' for s in ('s1', 's2')))
        for sample, count in [('s1', 4), ('s2', 3)]:
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
        self.assertIn('SDM coarse-subcluster seed extension failed', result.stdout)
        self.assertNotIn('Fallback to', result.stdout)
        args = json.loads(self.seed_call.read_text())
        self.assertEqual(args[args.index('-seedSubclusters')+1], '1')
        self.assertEqual(args[args.index('-merge_pairs_seed')+1], '1')
        files = [Path(x) for x in args[args.index('-i_fastq')+1].split(',')]
        self.assertEqual([f.name for f in files], ['derep.subclusters.1.fq', 'derep.subclusters.2.fq'] if paired else ['derep.subclusters.fq'])
        records = [f.read_text().splitlines() for f in files]
        self.assertEqual(len(records[0])//4, 2)
        self.assertEqual(sum(int(re.search(r';size=(\d+);', h).group(1)) for h in records[0][::4]), 7)
        if paired:
            self.assertEqual(records[0][::4], records[1][::4])
            self.assertEqual(len(set(records[0][1::4])), 1)
            self.assertEqual(len(set(records[1][1::4])), 2)
        commands = (self.out/'LotuSLogS/LotuS_cmds.log').read_text()
        self.assertIn('-derepIdentity '+format(float(identity)*100, '.12g'), commands)
        for flag in ('-derepStoreQuals 1', '-derepPerSR 0', '-merge_pairs_derep 0', '-merge_pairs_filter 0', '-merge_pairs_demulti 0'):
            self.assertIn(flag, commands)
        self.assertEqual(Path(args[args.index('-derep_map')+1]).name, 'derep.map')

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
        self.assertIn('Missing or empty coarse-dereplication seed FASTQ', result.stdout)
        self.assertFalse(self.seed_call.exists())

    def test_default_hq_path_still_counts_without_coarse_flags(self):
        self.run_lotus(self.inputs(False))
        self.check_counts({'s1': 4, 's2': 3})
        commands = (self.out/'LotuSLogS/LotuS_cmds.log').read_text()
        self.assertIn('derep.1.hq.fq', commands)
        self.assertNotIn('-seedSubclusters', commands)
        self.assertNotIn('-derepIdentity', commands)
        self.assertNotIn('-derepStoreQuals', commands)

    @unittest.skipUnless(HAS_SEEDS, 'requires upstream -seedSubclusters implementation')
    def test_native_seed_extension_preserves_sample_counts(self):
        for paired in (False, True):
            with self.subTest(paired=paired):
                self.out = self.root/f'native_{paired}'
                self.run_lotus(self.inputs(paired)+['-coarseDerep', '0.97'])
                self.check_counts({'s1': 4, 's2': 3})


if __name__ == '__main__':
    unittest.main()

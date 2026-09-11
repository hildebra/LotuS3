#!/usr/bin/env python3
"""ONT wiring regressions: real bundled SDM, controlled Barbell/Savont/mappers.

Run: python3 -m unittest discover -s tests -p 'ont_integration.py' -v
No scientific validation of the stand-in tools is implied. A temporary copy of
lotus3 exits after SDM builds the abundance matrix, before unrelated taxonomy.
"""
import json
import os
from pathlib import Path
import random
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
TOOL = r'''#!/usr/bin/env python3
import json, os, pathlib, shutil, sys
name = pathlib.Path(sys.argv[0]).name
a = sys.argv[1:]
def arg(flag): return a[a.index(flag)+1]
if '--version' in a or '-v' in a:
    print({'minimap2':'2.28', 'vsearch':'vsearch v2.29.0', 'LCA':'0.28'}.get(name, '1.0'))
    sys.exit(0)
with open(os.environ['ONT_TEST_CALLS'], 'a') as log:
    log.write(json.dumps([name] + a) + '\n')
if name == 'barbell':
    out = pathlib.Path(arg('-o')); out.mkdir()
    for f in pathlib.Path(os.environ['ONT_TEST_BARCODES']).glob('*.fastq'):
        shutil.copyfile(f, out/f.name)
elif name == 'savont':
    if os.environ.get('ONT_TEST_FAIL'): sys.exit(17)
    out = pathlib.Path(arg('-o')); out.mkdir()
    shutil.copyfile(a[1], out/'input_snapshot.fq')
    if not os.environ.get('ONT_TEST_EMPTY'):
        (out/'final_asvs.fasta').write_text('>final_consensus_0_depth_7 debug_id:7 chimera_score:0\n' + os.environ['ONT_TEST_CONSENSUS'] + '\n')
elif name in ('minimap2', 'vsearch'):
    if name == 'minimap2':
        query = pathlib.Path(a[-1]); out = pathlib.Path(arg('-o'))
    else:
        query = pathlib.Path(arg('--usearch_global')); out = pathlib.Path(arg('-uc'))
    # Requiring FASTA catches accidental use of SDM's FASTQ dereplication mode.
    text = query.read_text(); assert text.startswith('>'), text[:100]
    entries = [x.splitlines() for x in text.split('>')[1:]]
    with out.open('w') as f:
        for lines in entries:
            rid = lines[0].split()[0]; n = len(''.join(lines[1:]))
            if name == 'minimap2':
                f.write(f'{rid}\t{n}\t0\t{n}\t+\tASV0\t{n}\t0\t{n}\t{n-1}\t{n}\t60\tcg:Z:{n}M\n')
            else:
                f.write(f'H\t0\t{n}\t99.9\t+\t0\t0\t{n}M\t{rid}\tASV0\n')
else:
    sys.exit('Unexpected tool call: ' + name)
'''


class ONTIntegration(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='lotus3-ont-test-')
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.install = self.root/'install'; self.install.mkdir()
        for name in ('configs', 'helpers', 'Example'):
            (self.install/name).symlink_to(ROOT/name, target_is_directory=True)
        (self.install/'bin').mkdir()
        (self.install/'bin'/'sdm').symlink_to(ROOT/'bin'/'sdm')
        src = (ROOT/'lotus3').read_text()
        checkpoint = 'undef $tmpOTU;'
        assert src.count(checkpoint) == 1
        src = src.replace(checkpoint, '''
atomic_write_text("$outdir/ont_test_state.json", JSON::PP->new->encode({
    seed => $OTUSEED, seed_extension => $seedExtDone, map => $mapHref, combined => $combHref,
    sdm_options => $sdmOpt, cluster => $ClusterPipe, preset => $mini2RdPreset}));
write_repro_manifest(); release_output_lock(); exit(0);
''' + checkpoint)
        self.script = self.install/'lotus3'; self.script.write_text(src)
        self.tools = self.root/'tools'; self.tools.mkdir()
        for name in ('barbell', 'savont', 'minimap2', 'vsearch', 'LCA'):
            f = self.tools/name; f.write_text(TOOL); f.chmod(0o755)
        self.cfg = self.root/'lotus.cfg'
        self.cfg.write_text(f'sdm {ROOT}/bin/sdm\n' + ''.join(f'{n} {self.tools/n}\n' for n in ('barbell', 'savont', 'minimap2', 'vsearch', 'LCA')) + 'CheckForUpdates 0\n')
        rng = random.Random(42)
        self.seq = ''.join(rng.choice('ACGT') for _ in range(800))
        self.consensus = ('T' if self.seq[0] != 'T' else 'A') + self.seq[1:]
        self.fwd = 'AGAGTTTGATCCTGGCTCAG'
        self.rev = 'TACGGYTACCTTGTTACGACTT'
        self.revcomp = self.rev.replace('Y', 'T').translate(str.maketrans('ACGT', 'TGCA'))[::-1]
        self.reads = self.root/'reads'; self.reads.mkdir()
        self.barcodes = self.root/'barcodes'; self.barcodes.mkdir()
        for sample, barcode, count in [('s1','BC01',4), ('s2','BC02',3), ('low','BC03',1), ('noise','BC99',2)]:
            seq = self.fwd + self.seq + self.revcomp
            fq = ''.join(f'@{sample}_{i}\n{seq}\n+\n' + 'I'*len(seq) + '\n' for i in range(count))
            (self.reads/f'{sample}.fq').write_text(fq)
            (self.barcodes/f'{barcode}.trimmed.fastq').write_text(fq)
        self.raw = self.root/'raw.fastq'; self.raw.write_text((self.reads/'s1.fq').read_text())
        self.map = self.root/'map.tsv'
        self.write_map()
        self.ref = self.root/'ref.fna'; self.ref.write_text('>ref\n'+self.seq+'\n')
        self.tax = self.root/'ref.tax'; self.tax.write_text('ref\tBacteria;P;C;O;F;G;S\n')
        self.out = self.root/'output'
        self.calls = self.root/'calls.jsonl'
        self.env = dict(os.environ, ONT_TEST_CALLS=str(self.calls), ONT_TEST_BARCODES=str(self.barcodes), ONT_TEST_CONSENSUS=self.consensus)

    def write_map(self, barbell=False, rows=None, header=None):
        header = header or ('ONTBarcode' if barbell else 'fastqFile')
        rows = rows or ([('s1','BC01'),('s2','BC02'),('low','BC03'),('missing','BC04')] if barbell else [('s1','s1.fq'),('s2','s2.fq')])
        self.map.write_text(f'#SampleID\t{header}\tForwardPrimer\tReversePrimer\n' + ''.join(f'{s}\t{x}\t{self.fwd}\t{self.rev}\n' for s,x in rows))

    def run_lotus(self, extra=(), barbell=False, ok=True):
        args = ['perl', str(self.script), '-i', str(self.raw if barbell else self.reads), '-m', str(self.map), '-o', str(self.out), '-c', str(self.cfg), '-p', 'ONT', '-t', '1', '-lulu', '0', '-removePhiX', '0', '-buildPhylo', '0', '-deactivateChimeraCheck', '1', '-refDB', str(self.ref), '-tax4refDB', str(self.tax), '-taxAligner', 'vsearch']
        if barbell: args += ['-ontDemux','barbell','-ontKit','SQK-RBK114-96']
        result = subprocess.run(args+list(extra), env=self.env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=30)
        (self.root/'run.log').write_text(result.stdout)
        if ok and result.returncode != 0:
            prog = self.out/'LotuSLogS'/'LotuS_progout.log'
            self.fail(result.stdout[-6000:] + (prog.read_text()[-6000:] if prog.exists() else ''))
        if not ok: self.assertNotEqual(result.returncode, 0, result.stdout)
        return result

    def tool_calls(self, name=None):
        calls = [json.loads(l) for l in self.calls.read_text().splitlines()] if self.calls.exists() else []
        return [x for x in calls if name is None or x[0] == name]

    def check_counts(self):
        rows = (self.out/'OTU.txt').read_text().splitlines()
        self.assertEqual(len(rows), 2, rows)
        self.assertEqual(rows[1].split('\t')[0], 'OTU_0')
        self.assertEqual(dict(zip(rows[0].split('\t')[1:], map(int, rows[1].split('\t')[1:]))), {'s1':4,'s2':3})
        state = json.loads((self.out/'ont_test_state.json').read_text())
        self.assertEqual(state['seed_extension'], 1)
        return state

    def test_savont_backmap_counts_and_full_reads(self):
        self.run_lotus()
        state = self.check_counts()
        self.assertIn('map-ont', state['preset'])
        self.assertEqual(state['cluster'], 9)
        snapshot = self.out/'tmpFiles/savont_out/input_snapshot.fq'
        self.assertEqual(snapshot.read_text().count('\n+\n'), 7)
        self.assertNotIn(self.fwd, snapshot.read_text())
        self.assertEqual((self.out/'tmpFiles/derep.fas').read_text()[0], '>')
        self.assertIn(self.consensus, (self.out/'tmpFiles/tmp_otu.fa').read_text())
        self.assertEqual(self.consensus, ''.join(Path(state['seed']).read_text().splitlines()[1:]))
        self.assertFalse((self.out/'tmpFiles/demultiplexed').exists())
        args = self.tool_calls('savont')[0]
        self.assertEqual(args[args.index('--quality-value-cutoff')+1], '80')
        self.assertNotIn('--single-strand', args)

    def test_savont_explicit_options_and_vsearch_mapping(self):
        self.run_lotus(['-CL','9','-useMini4map','0','-saveDemultiplex','2','-savontSingleStrand','1','-savontQualCutoff','95.5','-savontMinBaseQual','20','-savontChimeraErrors','2'])
        self.check_counts()
        args = self.tool_calls('savont')[0]
        for flag,value in [('--quality-value-cutoff','95.5'),('--minimum-base-quality','20'),('--chimera-allowable-errors','2')]:
            self.assertEqual(args[args.index(flag)+1], value)
        self.assertIn('--single-strand', args)
        self.assertTrue(list((self.out/'demultiplexed').glob('*.fq')))

    def test_barbell_drops_samples_and_rewrites_map(self):
        self.write_map(barbell=True)
        original = self.map.read_bytes()
        self.run_lotus(['-ontMinReads','2'], barbell=True)
        state = self.check_counts()
        self.assertEqual(set(state['map']), {'#SampleID','s1','s2'})
        self.assertEqual(state['combined'], {'s1':'s1','s2':'s2'})
        working_map = (self.out/'primary/in.map').read_text()
        self.assertIn('fastqFile', working_map)
        self.assertNotIn('missing', working_map)
        self.assertEqual(self.map.read_bytes(), original)
        self.assertEqual({x.name for x in (self.out/'tmpFiles/ont_demux').iterdir()}, {'s1.fq','s2.fq'})

    def test_barbell_without_new_fastq_column(self):
        self.write_map(barbell=True)
        self.map.write_text(self.map.read_text().replace('ONTBarcode\t', 'ONTBarcode\tfastqFile\t', 1).replace('BC01\t', 'BC01\told1.fq\t').replace('BC02\t', 'BC02\told2.fq\t').replace('BC03\t', 'BC03\told3.fq\t').replace('BC04\t', 'BC04\told4.fq\t'))
        self.run_lotus(['-ontMinReads','2','-ontWriteFastqCol','0'], barbell=True)
        self.check_counts()
        working_map = (self.out/'primary/in.map').read_text()
        self.assertIn('fastqFile', working_map)
        self.assertNotIn('old1.fq', working_map)
        self.assertNotIn('missing', working_map)

    def test_barbell_missing_fastq_column_in_mode_zero_rejected(self):
        self.write_map(barbell=True)
        result = self.run_lotus(['-ontWriteFastqCol','0','--dry-run'], barbell=True, ok=False)
        self.assertIn('requires an existing fastqFile', result.stdout)
        self.assertEqual(self.tool_calls(), [])

    def test_ont_custom_sdm_preset_is_respected(self):
        custom = ROOT/'configs/sdm_ONT_LSSU.txt'
        self.run_lotus(['-s',str(custom),'--dry-run'])
        self.assertEqual(self.tool_calls(), [])

    def test_missing_savont_is_reported(self):
        (self.tools/'savont').unlink()
        result = self.run_lotus(['--dry-run'], ok=False)
        self.assertIn('No valid savont binary', result.stdout)
        self.assertFalse(self.out.exists())

    def test_missing_barbell_is_reported(self):
        self.write_map(barbell=True)
        (self.tools/'barbell').unlink()
        result = self.run_lotus(['--dry-run'], barbell=True, ok=False)
        self.assertIn('barbell executable unavailable', result.stdout)
        self.assertEqual(self.tool_calls(), [])

    def test_paired_reads_are_rejected(self):
        self.map.write_text(self.map.read_text().replace('s1.fq', 's1.fq,s2.fq').replace('s2\ts2.fq', 's2\ts2.fq,s1.fq'))
        result = self.run_lotus(['--dry-run'], ok=False)
        self.assertIn('single-end FASTQ', result.stdout)

    def test_savont_requires_dereplication(self):
        result = self.run_lotus(['-highmem','0','--dry-run'], ok=False)
        self.assertIn('requires SDM dereplication', result.stdout)

    def test_dry_run_is_non_destructive(self):
        self.write_map(barbell=True)
        self.out.mkdir(); sentinel = self.out/'LotuS_output_schema_version.txt'; sentinel.write_text('keep')
        report = self.root/'report.txt'; report.write_text('keep report')
        manifest = self.root/'manifest.txt'; manifest.write_text('keep manifest')
        result = self.run_lotus(['--dry-run','--dependencyReport',str(report),'--manifest',str(manifest)], barbell=True)
        self.assertIn('Dry-run validation complete', result.stdout)
        self.assertEqual(list(self.out.iterdir()), [sentinel])
        self.assertEqual(report.read_text(), 'keep report')
        self.assertEqual(manifest.read_text(), 'keep manifest')
        self.assertEqual(self.tool_calls(), [])

    def test_barbell_duplicate_barcode_rejected_in_dry_run(self):
        self.write_map(True, [('s1','BC01'),('s2','BC01')])
        r = self.run_lotus(['--dry-run'], barbell=True, ok=False)
        self.assertIn('assigned to more than one sample', r.stdout)
        self.assertEqual(self.tool_calls(), [])
        self.assertFalse(self.out.exists())

    def test_barbell_unsafe_sample_rejected(self):
        self.write_map(True, [('../outside','BC01')])
        r = self.run_lotus(['--dry-run'], barbell=True, ok=False)
        self.assertIn('Unsafe ONT SampleID', r.stdout)
        self.assertEqual(self.tool_calls(), [])

    def test_savont_failure_stops_before_backmapping(self):
        self.env['ONT_TEST_FAIL'] = '1'
        self.run_lotus(ok=False)
        self.assertEqual(self.tool_calls('minimap2'), [])
        self.assertFalse((self.out/'OTU.txt').exists())
        self.assertTrue((self.out/'tmpFiles/demultiplexed').exists())

    def test_empty_savont_output_stops_before_backmapping(self):
        self.env['ONT_TEST_EMPTY'] = '1'
        result = self.run_lotus(ok=False)
        self.assertIn('savont produced no ASVs', result.stdout)
        self.assertEqual(self.tool_calls('minimap2'), [])

    def test_invalid_flags_rejected(self):
        for args in [('-ontDemux','invalid'),('-p','miSeq','-ontDemux','barbell'),('-p','miSeq','-CL','savont'),('-ontMinReads','-1'),('-ontWriteFastqCol','2'),('-savontSingleStrand','2'),('-savontQualCutoff','101'),('-savontMinBaseQual','-1'),('-savontChimeraErrors','-1')]:
            with self.subTest(args=args):
                self.run_lotus(args, ok=False)
                self.assertEqual(self.tool_calls(), [])

    def test_miseq_dry_run_remains_supported(self):
        self.run_lotus(['-p','miSeq','-CL','vsearch','--dry-run'])
        self.assertEqual(self.tool_calls(), [])

    def test_removed_dnaclust_stays_removed(self):
        result = self.run_lotus(['-CL','4'], ok=False)
        self.assertIn('has been removed', result.stdout)
        self.assertFalse(self.out.exists())


class SavontConverter(unittest.TestCase):
    def test_best_hit_and_survivor_filter(self):
        with tempfile.TemporaryDirectory(prefix='savont-converter-test-') as tmp:
            root = Path(tmp)
            asvs = root/'asvs.fa'; mappings = root/'map.tsv'
            uc = root/'out.uc'; fasta = root/'out.fa'
            asvs.write_text('>final_0 debug_id:7\nACGT\n>final_1 debug_id:8\nTGCA\n')
            mappings.write_text('r1 description\tasv:7\t1\t90\nr1 description\tasv:8\t0\t80\nr2\tasv:7\t0\t90\nr3\tasv:7\t0\t90\nr4\tasv:99\t0\t100\n')
            run = subprocess.run(['perl',str(ROOT/'bin/savont2uc.pl'),'--asvs',str(asvs),'--map',str(mappings),'--ucout',str(uc),'--fnaout',str(fasta)], capture_output=True, text=True)
            self.assertEqual(run.returncode, 0, run.stderr)
            self.assertEqual(fasta.read_text(), '>r2\nACGT\n>r1\nTGCA\n')
            self.assertEqual(uc.read_text(), 'r2\totu1\t*\nr3\tmatch\tdqt=1;top=r2(99%);\nr1\totu2\t*\n')

    def test_invalid_id_mode_and_duplicate_asv_rejected(self):
        with tempfile.TemporaryDirectory(prefix='savont-converter-test-') as tmp:
            root = Path(tmp)
            asvs = root/'asvs.fa'; mappings = root/'map.tsv'
            asvs.write_text('>a debug_id:7\nACGT\n>b debug_id:7\nTGCA\n')
            mappings.write_text('r1\tasv:7\t0\t90\n')
            args = ['perl',str(ROOT/'bin/savont2uc.pl'),'--asvs',str(asvs),'--map',str(mappings),'--ucout',str(root/'out.uc'),'--fnaout',str(root/'out.fa')]
            for extra, message in [(['--idmode','invalid'],'--idmode must be'), ([], 'Duplicate ASV debug_id')]:
                with self.subTest(extra=extra):
                    run = subprocess.run(args+extra, capture_output=True, text=True)
                    self.assertNotEqual(run.returncode, 0)
                    self.assertIn(message, run.stderr)
                    self.assertFalse((root/'out.uc').exists())


if __name__ == '__main__':
    unittest.main()

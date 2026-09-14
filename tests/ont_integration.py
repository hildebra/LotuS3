#!/usr/bin/env python3
"""ONT wiring regressions: real bundled SDM, controlled Barbell/Savont/mappers.

Run: python3 -m unittest discover -s tests -p 'ont_integration.py' -v
No scientific validation of the stand-in tools is implied. A temporary copy of
lotus3 exits after SDM builds the abundance matrix, before unrelated taxonomy.
"""
from collections import Counter
import gzip
import json
import os
from pathlib import Path
import random
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SDM = Path(os.environ.get('LOTUS_TEST_SDM', str(ROOT/'bin/sdm'))).resolve()
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
        seqs = json.loads(os.environ.get('ONT_TEST_CONSENSUSES', 'null')) or [os.environ['ONT_TEST_CONSENSUS']]
        (out/'final_asvs.fasta').write_text(''.join(f'>final_consensus_{i} debug_id:{i} chimera_score:0\n{seq}\n' for i,seq in enumerate(seqs)))
elif name == 'vsearch' and '--cluster_size' in a:
    query = pathlib.Path(arg('--cluster_size'))
    assert query.name == 'derep.fas'
    entries = [x.splitlines() for x in query.read_text().split('>')[1:]]
    pathlib.Path(arg('--consout')).write_text('>OTU0\n'+os.environ['ONT_TEST_CONSENSUS']+'\n')
    with pathlib.Path(arg('--uc')).open('w') as out:
        for lines in entries:
            rid = lines[0].split()[0]; n = len(''.join(lines[1:]))
            out.write(f'H\t0\t{n}\t99.9\t+\t0\t0\t{n}M\t{rid}\tOTU0\n')
elif name in ('minimap2', 'vsearch'):
    if name == 'minimap2':
        query = pathlib.Path(a[-1]); out = pathlib.Path(arg('-o')); db = pathlib.Path(a[-2])
    else:
        query = pathlib.Path(arg('--usearch_global')); out = pathlib.Path(arg('-uc')); db = pathlib.Path(arg('-db'))
    text = query.read_text()
    if text.startswith('@'):
        lines = text.splitlines()
        entries = [[lines[i][1:], lines[i+1]] for i in range(0, len(lines), 4)]
    else:
        assert text.startswith('>'), text[:100]
        entries = [x.splitlines() for x in text.split('>')[1:]]
    targets = [(x.splitlines()[0].split()[0], ''.join(x.splitlines()[1:])) for x in db.read_text().split('>')[1:]]
    with out.open('w') as f:
        for lines in entries:
            rid = lines[0].split()[0]; seq = ''.join(lines[1:]); n = len(seq)
            if os.environ.get('ONT_TEST_SKIP_PREFIX') and rid.startswith(os.environ['ONT_TEST_SKIP_PREFIX']):
                continue
            target = min(targets, key=lambda t: sum(x != y for x,y in zip(seq, t[1])) + abs(n-len(t[1])))[0]
            if name == 'minimap2':
                f.write(f'{rid}\t{n}\t0\t{n}\t+\t{target}\t{n}\t0\t{n}\t{n-1}\t{n}\t60\tcg:Z:{n}M\n')
            else:
                f.write(f'H\t0\t{n}\t99.9\t+\t0\t0\t{n}M\t{rid}\t{target}\n')
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
        (self.install/'bin'/'sdm').symlink_to(SDM)
        src = (ROOT/'lotus3').read_text()
        checkpoint = 'undef $tmpOTU;'
        assert src.count(checkpoint) == 1
        src = src.replace(checkpoint, '''
atomic_write_text("$outdir/ont_test_state.json", JSON::PP->new->encode({
    seed => $OTUSEED, seed_extension => $seedExtDone, map => $mapHref, combined => $combHref,
    sdm_options => $sdmOpt, cluster => $ClusterPipe, preset => $mini2RdPreset, dereplication => $sdmDerepDo}));
write_repro_manifest(); release_output_lock(); exit(0);
''' + checkpoint)
        self.script = self.install/'lotus3'; self.script.write_text(src)
        self.tools = self.root/'tools'; self.tools.mkdir()
        for name in ('barbell', 'savont', 'minimap2', 'vsearch', 'LCA'):
            f = self.tools/name; f.write_text(TOOL); f.chmod(0o755)
        self.cfg = self.root/'lotus.cfg'
        self.cfg.write_text(f'sdm {SDM}\n' + ''.join(f'{n} {self.tools/n}\n' for n in ('barbell', 'savont', 'minimap2', 'vsearch', 'LCA')) + 'CheckForUpdates 0\n')
        rng = random.Random(42)
        self.seq = ''.join(rng.choice('ACGT') for _ in range(1200))
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

    def check_counts(self, expected=None):
        rows = (self.out/'OTU.txt').read_text().splitlines()
        self.assertEqual(len(rows), 2, rows)
        self.assertEqual(rows[1].split('\t')[0], 'OTU_0')
        self.assertEqual(dict(zip(rows[0].split('\t')[1:], map(int, rows[1].split('\t')[1:]))), expected or {'s1':4,'s2':3})
        state = json.loads((self.out/'ont_test_state.json').read_text())
        self.assertEqual(state['seed_extension'], 1)
        return state

    def test_savont_counts_and_filtered_reads_without_dereplication(self):
        self.run_lotus()
        state = self.check_counts()
        self.assertIn('map-ont', state['preset'])
        self.assertEqual(state['cluster'], 9)
        snapshot = self.out/'tmpFiles/savont_out/input_snapshot.fq'
        self.assertEqual(snapshot.read_text().count('\n+\n'), 7)
        self.assertNotIn(self.fwd, snapshot.read_text())
        self.assertEqual(state['dereplication'], 0)
        self.assertEqual(list((self.out/'tmpFiles').glob('derep*')), [])
        self.assertFalse((self.out/'tmpFiles/savont_in.fq').exists())
        prepared = self.out/'tmpFiles/savont_reads.fq'
        self.assertEqual(snapshot.read_bytes(), prepared.read_bytes())
        self.assertEqual(Path(self.tool_calls('savont')[0][2]), prepared)
        self.assertEqual(Path(self.tool_calls('minimap2')[0][-1]), prepared)
        self.assertTrue(snapshot.read_text().startswith('@s1___'))
        self.assertEqual(Path(state['sdm_options']).name, 'sdm_ONT_SAVONT_effective.txt')
        self.assertIn('maxAmbiguousNT\t5%', (self.out/'primary/sdm_original_options.txt').read_text())
        self.assertIn('maxAmbiguousNT\t-1', Path(state['sdm_options']).read_text())
        self.assertIn(self.consensus, (self.out/'tmpFiles/tmp_otu.fa').read_text())
        self.assertEqual(self.consensus, ''.join(Path(state['seed']).read_text().splitlines()[1:]))
        self.assertFalse((self.out/'tmpFiles/demultiplexed').exists())
        args = self.tool_calls('savont')[0]
        self.assertEqual(args[args.index('--quality-value-cutoff')+1], '80')
        self.assertEqual(args[args.index('--min-read-length')+1], '1000')
        self.assertEqual(args[args.index('--max-read-length')+1], '2000')
        self.assertIn('--single-strand', args)

    def test_savont_both_strand_override(self):
        self.run_lotus(['-savontSingleStrand','0'])
        self.check_counts()
        self.assertNotIn('--single-strand', self.tool_calls('savont')[0])

    @staticmethod
    def rc(seq):
        return seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]

    def append_read(self, name, sequence=None, quality=None, sample='s1', folder=None):
        sequence = sequence if sequence is not None else self.fwd + self.seq + self.revcomp
        quality = quality if quality is not None else 'I' * len(sequence)
        self.assertEqual(len(sequence), len(quality))
        with ((folder or self.reads) / (sample + '.fq')).open('a') as out:
            out.write(f'@{name}\n{sequence}\n+\n{quality}\n')

    def test_savont_main_quality_filters_exclude_reads_from_asvs_and_counts(self):
        full = self.fwd + self.seq + self.revcomp
        self.append_read('low_average', quality='+' * len(full))  # Q10
        self.append_read('low_window', quality='I'*219 + '&'*150 + 'I'*(len(full)-369))
        self.append_read('ambiguous', sequence=full[:200]+'N'*61+full[261:])
        self.append_read('missing_reverse', sequence=self.fwd+self.seq)
        self.run_lotus()
        self.check_counts()
        snapshot = (self.out/'tmpFiles/savont_out/input_snapshot.fq').read_text()
        self.assertEqual(snapshot.count('\n+\n'), 7)
        for name in ('low_average','low_window','ambiguous','missing_reverse'):
            self.assertNotIn(name, snapshot)
        commands = self.sdm_commands()
        self.assertNotIn('-o_dereplicate', commands[0])
        self.assertNotIn('-derep_map', commands[1])
        self.assertNotIn('-options', commands[1])  # Prepared reads are not trimmed again.

    def test_savont_percentage_and_length_boundaries(self):
        passing = {}
        cases = [('short_boundary',1000,50,'N',True), ('long_boundary',2000,100,'N',True),
                 ('iupac_boundary',1000,50,'R',True), ('short_over',1000,51,'N',False),
                 ('long_over',2000,101,'N',False), ('too_short',999,0,'N',False),
                 ('too_long',2001,0,'N',False)]
        for name,length,count,base,accepted in cases:
            seq = list((self.seq*2)[:length])
            for i in range(count): seq[100+2*i] = base
            sequence = ''.join(seq)
            self.append_read(name, sequence=self.fwd+sequence+self.revcomp)
            if accepted: passing[name] = sequence
        self.run_lotus(['-saveDemultiplex','2'])
        self.check_counts({'s1':7,'s2':3})
        lines = (self.out/'tmpFiles/savont_reads.fq').read_text().splitlines()
        actual = {lines[i].split()[0].split('___',1)[1]: lines[i+1] for i in range(0,len(lines),4)}
        for name,length,count,base,accepted in cases:
            if accepted: self.assertEqual(actual[name], passing[name])
            else: self.assertNotIn(name, actual)
        saved = ''.join(p.read_text() for p in (self.out/'demultiplexed').glob('*.fq'))
        self.assertNotIn('short_over', saved)
        self.assertNotIn('long_over', saved)
        self.assertIn('kept 10; rejected 2', (self.out/'LotuSLogS/savont_ambiguity_filter.log').read_text())

    def test_savont_average_quality_fifteen(self):
        custom = self.root/'quality-only.txt'
        custom.write_text((ROOT/'configs/sdm_ONT_SAVONT_opt.txt').read_text()
                          .replace('QualWindowThreshhold\t14','QualWindowThreshhold\t-1'))
        full = self.fwd+self.seq+self.revcomp
        self.append_read('q15', quality='0'*len(full))
        self.append_read('q14', quality='/'*len(full))
        self.run_lotus(['-s',str(custom)])
        self.check_counts({'s1':5,'s2':3})
        prepared = (self.out/'tmpFiles/savont_reads.fq').read_text()
        self.assertIn('___q15 ', prepared)
        self.assertNotIn('___q14 ', prepared)

    def test_savont_integer_ambiguity_override(self):
        custom = self.root/'no-ambiguity.txt'
        custom.write_text((ROOT/'configs/sdm_ONT_SAVONT_opt.txt').read_text()
                          .replace('maxAmbiguousNT\t5%', 'maxAmbiguousNT\t0'))
        seq = self.fwd+self.seq+self.revcomp
        self.append_read('one_n', sequence=seq[:200]+'N'+seq[201:])
        self.run_lotus(['-s',str(custom)])
        self.check_counts()
        self.assertNotIn('one_n', (self.out/'tmpFiles/savont_reads.fq').read_text())
        self.assertFalse((self.out/'LotuSLogS/savont_ambiguity_filter.log').exists())

    def test_savont_invalid_ambiguity_percentage(self):
        custom = self.root/'invalid-ambiguity.txt'
        for value in ('-5%', '101%', '5.5', 'wrong'):
            with self.subTest(value=value):
                custom.write_text((ROOT/'configs/sdm_ONT_SAVONT_opt.txt').read_text()
                                  .replace('maxAmbiguousNT\t5%', 'maxAmbiguousNT\t'+value))
                result = self.run_lotus(['-s',str(custom),'--dry-run'], ok=False)
                self.assertIn('integer count or a percentage from 0% to 100%', result.stdout)
                self.assertEqual(self.tool_calls(), [])

    def test_savont_preserves_high_quality_homopolymer_tail(self):
        for sample in ('s1','s2'):
            path = self.reads/(sample+'.fq')
            lines = path.read_text().splitlines()
            for i in range(0,len(lines),4):
                lines[i+1] = self.seq + 'G'*15
                lines[i+3] = 'I'*len(lines[i+1])
            path.write_text('\n'.join(lines)+'\n')
        self.run_lotus(['-ontPrimerState','removed'])
        self.check_counts()
        prepared = (self.out/'tmpFiles/savont_reads.fq').read_text().splitlines()
        self.assertEqual(prepared[1::4], [self.seq+'G'*15]*7)
        self.assertEqual(prepared[3::4], ['I'*(len(self.seq)+15)]*7)

    def test_savont_does_not_recover_secondary_quality_reads(self):
        custom = self.root/'custom.txt'
        opts = (ROOT/'configs/sdm_ONT_SAVONT_opt.txt').read_text()
        custom.write_text(opts.replace('minAvgQuality\t15', 'minAvgQuality\t25')
                          .replace('*minAvgQuality\t25', '*minAvgQuality\t0'))
        full = self.fwd + self.seq + self.revcomp
        self.append_read('secondary_only', quality='6'*len(full))  # Q21
        self.run_lotus(['-s',str(custom),'-saveDemultiplex','2'])
        self.check_counts()
        self.assertNotIn('secondary_only', (self.out/'tmpFiles/savont_out/input_snapshot.fq').read_text())
        self.assertEqual(list((self.out/'tmpFiles').glob('derep*')), [])

    def test_savont_all_reads_filtered_stops_before_clustering(self):
        for sample in ('s1','s2'):
            file = self.reads/(sample+'.fq'); lines = file.read_text().splitlines()
            for i in range(3,len(lines),4): lines[i] = '+'*len(lines[i])
            file.write_text('\n'.join(lines)+'\n')
        result = self.run_lotus(ok=False)
        self.assertIn('No reads passed SDM filtering for Savont', result.stdout)
        self.assertEqual(self.tool_calls(), [])
        self.assertFalse((self.out/'OTU.txt').exists())

    def test_savont_sample_prefixes_preserve_shared_read_names(self):
        for sample in ('s1','s2'):
            file = self.reads/(sample+'.fq')
            file.write_text(file.read_text().replace('@'+sample+'_', '@shared_'))
        self.run_lotus()
        self.check_counts()

    def test_savont_zero_hit_sample_is_preserved(self):
        self.env['ONT_TEST_SKIP_PREFIX'] = 's2___'
        self.run_lotus()
        self.check_counts({'s1':4,'s2':0})

    def test_savont_combine_samples(self):
        rows = self.map.read_text().splitlines()
        self.map.write_text(rows[0]+'\tCombineSamples\n' + '\n'.join(x+'\tpooled' for x in rows[1:])+'\n')
        for sample in ('s1','s2'):
            file = self.reads/(sample+'.fq')
            file.write_text(file.read_text().replace('@'+sample+'_', '@shared_'))
        self.run_lotus()
        self.check_counts({'pooled':7})
        self.assertNotIn('CombineSamples', (self.out/'primary/sdm_input.map').read_text())
        self.assertIn('CombineSamples', (self.out/'primary/in.map').read_text())
        snapshot = (self.out/'tmpFiles/savont_out/input_snapshot.fq').read_text()
        self.assertIn('@s1___shared_0', snapshot)
        self.assertIn('@s2___shared_0', snapshot)

    def test_savont_mixed_groups_and_ungrouped_samples(self):
        self.write_map(rows=[('s1','s1.fq'),('s2','s2.fq'),('low','low.fq')])
        rows = self.map.read_text().splitlines()
        self.map.write_text(rows[0]+'\tCombineSamples\n' +
                           '\n'.join(x+('\tpooled' if x.startswith(('s1\t','s2\t')) else '\t') for x in rows[1:])+'\n')
        self.run_lotus()
        self.check_counts({'pooled':7,'low':1})

    def test_other_clusterers_keep_sdm_dereplication(self):
        for platform in ('ONT','miSeq'):
            with self.subTest(platform=platform):
                self.out = self.root/('output_'+platform)
                self.run_lotus(['-p',platform,'-CL','vsearch','-derepMin','0',
                                '-s',str(ROOT/'configs/sdm_ONT_LSSU.txt')])
                state = self.check_counts()
                self.assertNotEqual(state['dereplication'], 0)
                derep = self.out/'tmpFiles/derep.fas'
                self.assertEqual(derep.read_text().count('>'), 1)
                self.assertIn(';size=7;', derep.read_text())
                self.assertTrue((self.out/'tmpFiles/derep.map').is_file())
                self.assertIn('-o_dereplicate', self.sdm_commands()[0])
                self.assertIn('-derep_map', self.sdm_commands()[1])
                self.assertEqual(self.tool_calls('savont'), [])
                self.assertFalse((self.out/'tmpFiles/savont_reads.fq').exists())

    def ont_read(self, front='', rear='', reverse=False, offset=21, body=None):
        body = body or self.seq
        left = 'G'*offset + front + 'CTAGCATGATC'
        seq = left + self.fwd + body + self.revcomp + 'CATTGAC' + self.rc(rear) + 'T'*19
        qual = ''.join(chr(33 + 25 + i % 15) for i in range(len(seq)))
        body_qual = qual[len(left)+len(self.fwd):len(left)+len(self.fwd)+len(body)]
        return (self.rc(seq), qual[::-1], body_qual) if reverse else (seq, qual, body_qual)

    def pooled_input(self, reads):
        self.map.write_text('#SampleID\tBarcodeSequence\tForwardPrimer\tReversePrimer\n'
                            + ''.join(f'{s}\t{bc}\t{self.fwd}\t{self.rev}\n' for s,bc in
                                      [('s1','ACGTCAGTGCTAGACG'),('s2','TGCATCGACAGTTCGA')]))
        pooled = self.root/'pooled.fastq.gz'
        with gzip.open(pooled, 'wt') as out:
            for i,(seq,qual,_) in enumerate(reads):
                out.write(f'@pooled_{i}\n{seq}\n+\n{qual}\n')
        return ['-i',str(pooled)]

    def sdm_commands(self):
        return [x for x in (self.out/'LotuSLogS/LotuS_cmds.log').read_text().splitlines()
                if str(SDM) in x and '-sample_sep' in x]

    def test_pooled_offset_barcodes_indels_orientation_and_full_qualities(self):
        bc = 'ACGTCAGTGCTAGACG'
        variants = [bc, bc[:6]+'A'+bc[7:], bc[:7]+'T'+bc[7:], bc[:7]+bc[8:]]
        reads = [self.ont_read(front=x, reverse=rev) for x in variants for rev in (False,True)]
        reads += [self.ont_read(rear='TGCATCGACAGTTCGA', reverse=rev) for rev in (False,True,False)]
        self.run_lotus(self.pooled_input(reads))
        self.check_counts({'s1':8,'s2':3})
        fq = (self.out/'tmpFiles/savont_out/input_snapshot.fq').read_text().splitlines()
        self.assertEqual(Counter(fq[1::4]), Counter([self.seq]*len(reads)))
        self.assertEqual(Counter(fq[3::4]), Counter(x[2] for x in reads))
        commands = self.sdm_commands()
        self.assertEqual(len(commands), 2, commands)
        self.assertIn('-ontMode 1 -barcodeSearchWindow 200 -ontBarcodeEnds either', commands[0])
        self.assertIn('-o_fastq ', commands[0])
        self.assertNotIn('-o_dereplicate', commands[0])
        self.assertNotIn('-derep_map', commands[1])
        self.assertNotIn('-ontMode', commands[1])
        self.assertNotIn('-barcodeSearchWindow', commands[1])
        self.assertNotIn('-ontBarcodeEnds', commands[1])
        manifest = (self.out/'LotuSLogS/run_manifest.txt').read_text()
        self.assertIn('SDM barcode/primer search window: 200', manifest)
        self.assertIn('SDM barcode ends: either', manifest)

    def test_barcode_end_requirement_changes_counts_and_rejects_conflicts(self):
        a,b = 'ACGTCAGTGCTAGACG','TGCATCGACAGTTCGA'
        reads = [self.ont_read(front=a), self.ont_read(front=a,rear=a),
                 self.ont_read(rear=b), self.ont_read(front=b,rear=b,reverse=True),
                 self.ont_read(front=a,rear=b)]
        extra = self.pooled_input(reads)
        for ends,count in [('either',2),('both',1)]:
            with self.subTest(ends=ends):
                self.out = self.root/('output_'+ends)
                self.run_lotus([*extra,'-ontBarcodeEnds',ends])
                self.check_counts({'s1':count,'s2':count})
                self.assertIn('-ontBarcodeEnds '+ends, self.sdm_commands()[0])

    def test_barcode_search_window_changes_retained_reads(self):
        reads = [self.ont_read(front=bc, offset=offset) for bc in
                 ('ACGTCAGTGCTAGACG','TGCATCGACAGTTCGA') for offset in (0,230)]
        extra = self.pooled_input(reads)
        for window,count in [(200,1),(300,2)]:
            with self.subTest(window=window):
                self.out = self.root/f'output_{window}'
                self.run_lotus([*extra,'-barcodeSearchWindow',str(window)])
                self.check_counts({'s1':count,'s2':count})

    def test_both_ends_ignored_for_filename_assignment_after_barbell(self):
        self.write_map(barbell=True)
        self.run_lotus(['-ontMinReads','2','-ontBarcodeEnds','both'],barbell=True)
        self.check_counts()
        self.assertIn('-ontMode 1', self.sdm_commands()[0])

    def test_non_ont_primary_command_omits_ont_options(self):
        self.run_lotus(['-p','miSeq','-CL','vsearch','-saveDemultiplex','1'])
        commands = self.sdm_commands()
        self.assertEqual(len(commands), 1, commands)
        for flag in ('-ontMode','-barcodeSearchWindow','-ontBarcodeEnds'):
            self.assertNotIn(flag, commands[0])

    def sdm_version_wrapper(self, banner):
        wrapper = self.tools/'sdm'
        wrapper.write_text('#!/usr/bin/env perl\nif (@ARGV && ($ARGV[0] eq "-version" || $ARGV[0] eq "-v")) '
                           + '{ print '+json.dumps(banner)+'; exit 0; }\nexec "'+str(SDM)+'", @ARGV;\n')
        wrapper.chmod(0o755)
        self.cfg.write_text(self.cfg.read_text().replace('sdm '+str(SDM)+'\n','sdm '+str(wrapper)+'\n'))

    def test_old_sdm_rejected_before_replacing_output(self):
        self.sdm_version_wrapper('sdm 3.50 beta\n')
        self.out.mkdir()
        marker = self.out/'LotuS_output_schema_version.txt'; marker.write_text('keep')
        result = self.run_lotus(ok=False)
        self.assertIn('ONT preprocessing requires SDM >= 3.51', result.stdout)
        self.assertEqual(marker.read_text(),'keep')
        self.assertEqual(self.tool_calls(), [])

    def test_capability_marker_allows_older_numbered_build(self):
        self.sdm_version_wrapper('sdm 3.50 beta\nONT amplicon end matching: enabled (-ontMode 1)\n')
        self.run_lotus(['--dry-run'])
        self.assertFalse(self.out.exists())

    def test_non_ont_can_use_older_sdm(self):
        self.sdm_version_wrapper('sdm 3.43\n')
        self.run_lotus(['-p','miSeq','-CL','vsearch','--dry-run'])
        self.assertFalse(self.out.exists())

    def test_ont_rejects_incompatible_map_columns_even_empty(self):
        original = self.map.read_text()
        for field in ('Barcode2ndPair','MIDfqFile','SampleIDinHead','alignmentFile','fnaFile','qualFile'):
            with self.subTest(field=field):
                rows = original.splitlines()
                self.map.write_text(rows[0]+'\t'+field+'\n'+'\n'.join(x+'\t' for x in rows[1:])+'\n')
                result = self.run_lotus(['--dry-run'],ok=False)
                self.assertIn(field, result.stdout)
                self.assertFalse(self.out.exists())

    def test_ont_barcode_budget_and_alphabet_rejected_early(self):
        self.pooled_input([self.ont_read(front='ACGTCAGTGCTAGACG')])
        original = self.map.read_text()
        custom = self.root/'sdm_custom.txt'
        for errors,barcode in [('4','ACGTCAGTGCTAGACG'),('-1','ACGTCAGTGCTAGACG'),
                               ('1.5','ACGTCAGTGCTAGACG'),('','ACGTCAGTGCTAGACG'),('1 #comment','ACGTCAGTGCTAGACG'),('1','ACGTNCGT'),('1','A')]:
            with self.subTest(errors=errors,barcode=barcode):
                custom.write_text((ROOT/'configs/sdm_ONT_savont.txt').read_text().replace('maxBarcodeErrs\t1','maxBarcodeErrs\t'+errors))
                self.map.write_text(original.replace('ACGTCAGTGCTAGACG',barcode))
                result = self.run_lotus(['-i',str(self.root/'pooled.fastq.gz'),'-s',str(custom),'--dry-run'],ok=False)
                self.assertIn('maxBarcodeErrs',result.stdout)
                self.assertFalse(self.out.exists())

    def test_multiple_savont_consensuses_keep_correct_counts(self):
        second = self.seq.translate(str.maketrans('ACGT','TGCA'))
        consensus2 = ('T' if second[0]!='T' else 'A')+second[1:]
        self.env['ONT_TEST_CONSENSUSES'] = json.dumps([self.consensus,consensus2])
        for sample,counts in [('s1',(4,2)),('s2',(3,5))]:
            with (self.reads/(sample+'.fq')).open('w') as out:
                for j,(body,count) in enumerate(zip((self.seq,second),counts)):
                    seq = self.fwd+body+self.revcomp
                    for i in range(count):
                        out.write(f'@{sample}_{j}_{i}\n{seq}\n+\n'+ 'I'*len(seq)+'\n')
        self.run_lotus()
        state = json.loads((self.out/'ont_test_state.json').read_text())
        seeds = {x.splitlines()[0]:''.join(x.splitlines()[1:]) for x in Path(state['seed']).read_text().split('>')[1:]}
        table = (self.out/'OTU.txt').read_text().splitlines()
        expected = {self.consensus:{'s1':4,'s2':3},consensus2:{'s1':2,'s2':5}}
        self.assertEqual(set(seeds.values()),set(expected))
        for row in table[1:]:
            key,*counts = row.split('\t')
            self.assertEqual(dict(zip(table[0].split('\t')[1:],map(int,counts))),expected[seeds[key]])

    def test_savont_explicit_options_and_vsearch_mapping(self):
        self.run_lotus(['-CL','9','-useMini4map','0','-saveDemultiplex','2','-savontSingleStrand','1','-savontQualCutoff','95.5','-savontMinBaseQual','20','-savontChimeraErrors','2'])
        self.check_counts()
        args = self.tool_calls('savont')[0]
        for flag,value in [('--quality-value-cutoff','95.5'),('--minimum-base-quality','20'),('--chimera-allowable-errors','2')]:
            self.assertEqual(args[args.index(flag)+1], value)
        self.assertIn('--single-strand', args)
        self.assertTrue(list((self.out/'demultiplexed').glob('*.fq')))
        mapper = self.tool_calls('vsearch')[0]
        mapping_reads = Path(mapper[mapper.index('--usearch_global')+1])
        self.assertEqual(mapping_reads, self.out/'tmpFiles/savont_mapping.fna')
        fasta = mapping_reads.read_text().splitlines()
        fastq = (self.out/'tmpFiles/savont_reads.fq').read_text().splitlines()
        self.assertEqual(fasta[::2], ['>'+h[1:] for h in fastq[::4]])
        self.assertEqual(fasta[1::2], fastq[1::4])
        self.assertEqual(len(fasta)//2, 7)

    def test_barbell_drops_samples_and_rewrites_map(self):
        self.write_map(barbell=True)
        original = self.map.read_bytes()
        self.run_lotus(['-ontMinReads','2'], barbell=True)
        self.assertNotIn('--maximize', self.tool_calls('barbell')[0])
        state = self.check_counts()
        self.assertEqual(set(state['map']), {'#SampleID','s1','s2'})
        self.assertEqual(state['combined'], {'s1':'s1','s2':'s2'})
        working_map = (self.out/'primary/in.map').read_text()
        self.assertIn('fastqFile', working_map)
        self.assertNotIn('missing', working_map)
        self.assertEqual(self.map.read_bytes(), original)
        self.assertEqual({x.name for x in (self.out/'tmpFiles/ont_demux').iterdir()}, {'s1.fq','s2.fq'})

    def test_barbell_maximize_is_opt_in(self):
        self.write_map(barbell=True)
        self.run_lotus(['-ontMinReads','2','-ontBarbellMaximize','1'], barbell=True)
        self.check_counts()
        self.assertIn('--maximize', self.tool_calls('barbell')[0])

    def test_primary_ont_path_does_not_need_barbell(self):
        (self.tools/'barbell').unlink()
        self.run_lotus()
        self.check_counts()
        citations = (self.out/'LotuSLogS/citations.txt').read_text()
        self.assertIn('10.64898/2026.05.26.727271', citations)
        self.assertNotIn('10.1093/bioinformatics/btag349', citations)

    def strip_fixture_primers(self, folder):
        for file in folder.iterdir():
            lines = file.read_text().splitlines()
            for i in range(0, len(lines), 4):
                lines[i+1] = lines[i+1][len(self.fwd):-len(self.revcomp)]
                lines[i+3] = lines[i+3][len(self.fwd):-len(self.revcomp)]
            file.write_text('\n'.join(lines)+'\n')

    def check_removed_primers(self, extra=(), barbell=False):
        original_map = self.map.read_bytes()
        if barbell: self.strip_fixture_primers(self.barcodes)
        else: self.strip_fixture_primers(self.reads)
        self.run_lotus(['-ontPrimerState','removed','-ontMinReads','2',*extra], barbell=barbell)
        self.check_counts()
        self.assertEqual(self.map.read_bytes(), original_map)
        self.assertIn('ForwardPrimer', (self.out/'primary/in.map').read_text())
        effective_map = (self.out/'primary/sdm_input.map').read_text()
        for field in ('ForwardPrimer','ReversePrimer','LinkerPrimerSequence'):
            self.assertNotIn(field, effective_map)
        self.assertIn('fastqFile', effective_map)
        lines = (self.out/'tmpFiles/savont_out/input_snapshot.fq').read_text().splitlines()
        self.assertEqual(lines[1::4], [self.seq]*7)
        self.assertEqual(lines[3::4], ['I'*len(self.seq)]*7)
        self.assertIn('ONT primers entering SDM: removed', (self.out/'LotuSLogS/run_manifest.txt').read_text())

    def test_removed_primers_with_barbell(self):
        self.write_map(barbell=True)
        self.check_removed_primers(barbell=True)

    def test_removed_primers_preserves_custom_preset(self):
        custom = self.root/'sdm_custom.txt'
        custom.write_text((ROOT/'configs/sdm_ONT_SAVONT_opt.txt').read_text().replace('minAvgQuality\t15','minAvgQuality\t25'))
        before = custom.read_bytes()
        self.check_removed_primers(['-s',str(custom)])
        self.assertEqual(custom.read_bytes(), before)
        effective = (self.out/'primary/sdm_ONT_removed.txt').read_text()
        self.assertIn('minAvgQuality\t25', effective)
        self.assertIn('RejectSeqWithoutFwdPrim\tF', effective)
        self.assertIn('RejectSeqWithoutRevPrim\tF', effective)
        self.assertIn('ExtensivePrimerChecks\tF', effective)

    def test_removed_primers_without_map_primer_columns(self):
        self.strip_fixture_primers(self.reads)
        self.map.write_text('#SampleID\tfastqFile\ns1\ts1.fq\ns2\ts2.fq\n')
        result = self.run_lotus(['-ontPrimerState','removed'])
        self.check_counts()
        self.assertNotIn('No forward PCR primer', result.stdout)

    def test_barbell_demultiplex_only_cites_only_barbell(self):
        self.write_map(barbell=True)
        self.run_lotus(['-saveDemultiplex','1'], barbell=True)
        citations = (self.out/'LotuSLogS/citations.txt').read_text()
        self.assertEqual(citations.count('10.1093/bioinformatics/btag349'), 1)
        self.assertNotIn('10.64898/2026.05.26.727271', citations)
        self.assertEqual(self.tool_calls('savont'), [])

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

    def test_savont_accepts_highmem_zero(self):
        self.run_lotus(['-highmem','0'])
        self.check_counts()
        self.assertEqual(list((self.out/'tmpFiles').glob('derep*')), [])

    def test_dry_run_is_non_destructive(self):
        self.write_map(barbell=True)
        self.out.mkdir(); sentinel = self.out/'LotuS_output_schema_version.txt'; sentinel.write_text('keep')
        report = self.root/'report.txt'; report.write_text('keep report')
        manifest = self.root/'manifest.txt'; manifest.write_text('keep manifest')
        result = self.run_lotus(['--dry-run','-ontPrimerState','removed','-ontBarbellMaximize','1','--dependencyReport',str(report),'--manifest',str(manifest)], barbell=True)
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
        self.assertTrue((self.out/'tmpFiles/savont_reads.fq').exists())

    def test_empty_savont_output_stops_before_backmapping(self):
        self.env['ONT_TEST_EMPTY'] = '1'
        result = self.run_lotus(ok=False)
        self.assertIn('savont produced no ASVs', result.stdout)
        self.assertEqual(self.tool_calls('minimap2'), [])

    def test_invalid_flags_rejected(self):
        for args in [('-barcodeSearchWindow','0'),('-barcodeSearchWindow','10001'),('-barcodeSearchWindow','2.5'),('-ontBarcodeEnds','bad'),('-p','miSeq','-barcodeSearchWindow','200'),('-p','miSeq','-ontBarcodeEnds','either'),('-ontDemux','invalid'),('-ontBarbellMaximize','2'),('-ontBarbellMaximize','1'),('-ontPrimerState','invalid'),('-p','miSeq','-ontPrimerState','removed'),('-p','miSeq','-ontDemux','barbell'),('-p','miSeq','-CL','savont'),('-ontMinReads','-1'),('-ontWriteFastqCol','2'),('-savontSingleStrand','2'),('-savontQualCutoff','101'),('-savontMinBaseQual','-1'),('-savontChimeraErrors','-1')]:
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

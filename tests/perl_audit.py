#!/usr/bin/env python3
"""Regression tests for the Perl audit; no database downloads or scientific claims.

Reuse the ONT fixture, with real SDM. Small probes invoke existing Perl helpers
inside a temporary script copy; complete-flow tests retain the whole pipeline.
"""
import fcntl
import json
import os
from pathlib import Path
import subprocess
import unittest
import ont_integration as ont

ROOT = ont.ROOT


class PerlAudit(unittest.TestCase):
    setUp = ont.ONTIntegration.setUp
    write_map = ont.ONTIntegration.write_map
    run_lotus = ont.ONTIntegration.run_lotus
    tool_calls = ont.ONTIntegration.tool_calls

    def probe(self, body, extra=(), ok=True):
        source = (ROOT/'lotus3').read_text()
        source = source.replace('prepLtsOptions();', body+'\nexit(0);', 1)
        self.script.write_text(source)
        return self.run_lotus(extra, ok=ok)

    def old_output(self):
        self.out.mkdir(exist_ok=True)
        (self.out/'.lotus3_created_by_this_run').write_text('previous run\n')
        logs = self.out/'LotuSLogS'; logs.mkdir(exist_ok=True)
        program_log = logs/'LotuS_progout.log'
        program_log.write_text('active run output\n')
        return program_log

    def test_lock_failure_preserves_active_logs(self):
        log = self.old_output()
        with (self.root/'.lotus3.output.lock').open('a') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            result = self.run_lotus(ok=False)
        self.assertIn('in use by another', result.stdout)
        self.assertEqual(log.read_text(), 'active run output\n')

    def test_lock_file_keeps_same_inode_between_runs(self):
        body = 'acquire_output_lock($outdir); release_output_lock();'
        self.probe(body)
        lock = self.root/'.lotus3.output.lock'
        self.assertTrue(lock.exists())
        inode = lock.stat().st_ino
        self.probe(body)
        self.assertEqual(inode, lock.stat().st_ino)

    def test_output_reset_preserves_nested_input_and_config(self):
        for flag, data in [('-m', self.map.read_text()), ('-c', self.cfg.read_text()), ('-i', self.raw.read_text()), ('-tax4refDB', self.tax.read_text())]:
            with self.subTest(flag=flag):
                log = self.old_output()
                asset = self.out/'protected.txt'; asset.write_text(data)
                for dry in [[], ['--dry-run']]:
                    result = self.run_lotus([flag, str(asset), *dry], ok=False)
                    self.assertIn('contains input or configuration', result.stdout)
                    self.assertEqual(asset.read_text(), data)
                    self.assertEqual(log.read_text(), 'active run output\n')

    def test_output_reset_preserves_configured_reference(self):
        self.old_output()
        asset = self.out/'reference.fna'; asset.write_text('>ref\nACGT\n')
        self.cfg.write_text(self.cfg.read_text()+f'TAX_REFDB_KSGP {asset}\n')
        result = self.run_lotus(ok=False)
        self.assertIn('contains input or configuration', result.stdout)
        self.assertTrue(asset.exists())

    def test_failed_preflight_preserves_previous_run(self):
        log = self.old_output()
        (self.out/'OTU.txt').write_text('previous results\n')
        self.map.write_text('this is not a valid mapping file\n')
        result = self.run_lotus(ok=False)
        self.assertIn('previous output preserved', result.stdout)
        self.assertEqual(log.read_text(), 'active run output\n')
        self.assertEqual((self.out/'OTU.txt').read_text(), 'previous results\n')

    def test_successful_preflight_allows_output_reset(self):
        self.old_output()
        (self.out/'previous.txt').write_text('old output\n')
        self.run_lotus()
        self.assertFalse((self.out/'previous.txt').exists())
        self.assertIn('OTU_0', (self.out/'OTU.txt').read_text())

    def test_temp_cleanup_preserves_input(self):
        scratch = self.root/'scratch'; scratch.mkdir()
        source = scratch/'s1.fq'; source.write_bytes(self.raw.read_bytes())
        (scratch/'s2.fq').write_bytes(self.raw.read_bytes())
        (scratch/'.lotus3_tmp_owned').write_text(f'Output: {self.out}\n')
        result = self.run_lotus(['-i',str(scratch),'-T',str(scratch)], ok=False)
        self.assertIn('contains input or configuration', result.stdout)
        self.assertEqual(source.read_bytes(), self.raw.read_bytes())

    def test_unknown_clusterer_and_unsupported_modes_fail_before_reset(self):
        log = self.old_output()
        for args, diagnostic in [(['-CL','vsarch'],'Unknown -CL'), (['-p','miSeq','-CL','vsearch','-highmem','0'],'SDM dereplication'), (['-exe','2'],'-exe must'), (['-useMini4map','2'],'-useMini4map must'), (['-saveDemultiplex','3'],'-saveDemultiplex must')]:
            with self.subTest(args=args):
                self.assertIn(diagnostic, self.run_lotus(args, ok=False).stdout)
                self.assertEqual(log.read_text(), 'active run output\n')

    def test_shell_active_path_characters_rejected(self):
        for name in ["quote'path", 'quote"path', 'glob[1]', 'glob*', 'glob?', 'back\\slash', 'paren(path)', 'semi;colon']:
            with self.subTest(name=name):
                result = self.run_lotus(['-o',str(self.root/name)], ok=False)
                self.assertIn('Unsafe output path', result.stdout)
                self.assertFalse((self.root/name).exists())

    def test_executable_crash_does_not_pass_version_check(self):
        (self.tools/'vsearch').write_text('#!/usr/bin/python3\nimport os, signal\nprint("vsearch v2.29.0", flush=True)\nos.kill(os.getpid(), signal.SIGTERM)\n')
        result = self.run_lotus(['--dry-run'], ok=False)
        self.assertIn('Executable check failed (exit 143)', result.stdout)
        self.assertEqual(self.tool_calls(), [])

    def test_old_minimap_version_rejected(self):
        (self.tools/'minimap2').write_text('#!/usr/bin/python3\nprint("2.9-r123")\n')
        result = self.run_lotus(['--dry-run'], ok=False)
        self.assertIn('too low, expected at least 2.17', result.stdout)

    def test_semantic_version_comparison(self):
        self.probe('die "wrong version order" if version_at_least("0.9", "0.25") || version_at_least("3.9", "3.43") || !version_at_least("2.28.1", "2.28");')

    def test_duplicate_reference_ids_are_warnings(self):
        fasta = self.root/'duplicate.fasta'
        taxonomy = self.root/'duplicate.tax'
        fasta.write_text(''.join('>dup\nACGT\n' for _ in range(5)))
        taxonomy.write_text(''.join('dup\tk__Bacteria\n' for _ in range(5)))
        self.env['AUDIT_DUP_FASTA'] = str(fasta)
        self.env['AUDIT_DUP_TAX'] = str(taxonomy)
        body = '''my ($fn,$fp,$fw,$fi) = fasta_validation_scan($ENV{AUDIT_DUP_FASTA}, 10);
my ($tn,$tp,$tw) = taxonomy_validation_scan($ENV{AUDIT_DUP_TAX}, 10);
die join("\n", @$fp, @$tp) if @$fp || @$tp;
print "RESULT:", join("\n", @$fw, @$tw);'''
        result = self.probe(body)
        self.assertIn('4 duplicate FASTA IDs', result.stdout)
        self.assertIn('4 duplicate taxonomy IDs', result.stdout)

    def taxonomy(self, rows, biom=True, hit=True, lca=True):
        f = self.root/'hierarchy.tsv'; f.write_text('header\n'+rows)
        self.env['AUDIT_TAX'] = str(f)
        return f'my @tax = readTaxIn($ENV{{AUDIT_TAX}}, {int(lca)}, {int(biom)}, {int(hit)}); print "RESULT:", JSON::PP->new->encode(\\@tax);'

    def result(self, completed):
        return json.loads(completed.stdout.split('RESULT:',1)[1])

    def test_taxonomy_prefixes_missing_ranks_and_reference_ids(self):
        body = self.taxonomy('ASV1\tk__Bacteria\tp__P\tc__C\to__O\tf__F\tg__G\ts__S\tref1\nASV2\tBacteria\tP\tC\tO\tF\tG\t\n')
        tax, levels, hits, hit_tax = self.result(self.probe(body))
        self.assertTrue(tax['ASV1'].startswith('k__Bacteria", "p__P'))
        self.assertNotIn('k__k__', tax['ASV1'])
        self.assertTrue(tax['ASV2'].endswith('s__?'))
        self.assertEqual(hits, {'ref1':['ASV1'],'ASV2':['ASV2']})
        self.assertEqual(hit_tax['ref1'], tax['ASV1'])

    def test_rdp_taxonomy_keeps_otu_id(self):
        body = self.taxonomy('Bacteria\tP\tC\tO\tF\tG\t\tOTU1\n', lca=False)
        tax, levels, hits, hit_tax = self.result(self.probe(body))
        self.assertEqual(list(tax), ['OTU1'])
        self.assertTrue(tax['OTU1'].endswith('s__?'))

    def test_malformed_and_duplicate_taxonomy_rejected(self):
        for rows, diagnostic in [('ASV1\tBacteria\tP\tC\tO\tF\n','Malformed taxonomy'), ('ASV1\tBacteria\tP\tC\tO\tF\tG\tS\n'*2,'Duplicate taxonomy')]:
            with self.subTest(diagnostic=diagnostic):
                self.assertIn(diagnostic, self.probe(self.taxonomy(rows), ok=False).stdout)

    def mapping_tool(self, name, contents):
        file = self.root/'alignments'; file.write_text(contents)
        self.env['AUDIT_ALIGNMENTS'] = str(file)
        (self.tools/name).write_text('''#!/usr/bin/python3
import os, pathlib, sys
if '--version' in sys.argv:
    print('2.28' if pathlib.Path(sys.argv[0]).name == 'minimap2' else 'vsearch v2.29.0'); sys.exit(0)
args=sys.argv[1:]
flag = '-o' if '-o' in args else '-uc'
pathlib.Path(args[args.index(flag)+1]).write_text(pathlib.Path(os.environ['AUDIT_ALIGNMENTS']).read_text())
''')

    def contamination(self, mini=1):
        return f'''ensure_dir($logDir); $mini2Bin = $ENV{{AUDIT_MINIMAP}}; $VSBin = $ENV{{AUDIT_VSEARCH}};
$useMini4map = {mini}; $doPhiX = 1; $uthreads = 1;
my $hits = contamination_rem($input, $refDBwanted, "phiX", 0);
print "RESULT:", JSON::PP->new->encode($hits);'''

    def setup_contamination(self):
        self.env['AUDIT_MINIMAP'] = str(self.tools/'minimap2')
        self.env['AUDIT_VSEARCH'] = str(self.tools/'vsearch')
        return ['-i',str(self.ref)]

    def test_contamination_requires_query_coverage_and_identity(self):
        self.mapping_tool('minimap2', ''.join(f'{rid}\t1000\t0\t{span}\t+\tgenome\t5000000\t0\t{span}\t{matches}\t{span}\t60\n' for rid,span,matches in [('short',100,100),('lowid',900,700),('valid',900,850),('boundary',500,450),('valid',900,850)]))
        hits = self.result(self.probe(self.contamination(), self.setup_contamination()))
        self.assertEqual(hits, {'phiX.0':{'valid':1,'boundary':1}})

    def test_malformed_paf_fails_with_context(self):
        for line, diagnostic in [('bad\trow\n','Malformed PAF'), ('bad\t100\t0\t0\t+\tr\t100\t0\t0\t0\t0\t60\n','Invalid PAF')]:
            self.mapping_tool('minimap2', line)
            self.assertIn(diagnostic, self.probe(self.contamination(), self.setup_contamination(), ok=False).stdout)

    def test_contamination_respects_vsearch_mapper(self):
        self.mapping_tool('vsearch', 'H\t0\t1000\t99\t+\t0\t0\t1000M\tvalid\tref\n')
        extra = self.setup_contamination(); self.env['AUDIT_MINIMAP'] = '/unavailable/minimap2'
        hits = self.result(self.probe(self.contamination(mini=0), extra))
        self.assertEqual(hits, {'phiX.0':{'valid':1}})

    def test_phix_hits_reach_matrix_filtering(self):
        source = (ROOT/'lotus3').read_text()
        marker = '# ////////////////////////// TAXONOMY'
        source = source.replace(marker, 'release_output_lock(); exit(0);\n'+marker,1)
        self.script.write_text(source)
        # The mapper uses each query header for both backmapping and PhiX search.
        self.cfg.write_text(self.cfg.read_text()+f'REFDB_PHIX {self.ref}\n')
        self.run_lotus(['-removePhiX','1','-keepOfftargets','1'])
        rows=(self.out/'OTU.txt').read_text().splitlines()
        self.assertIn('.phiX.', rows[1].split('\t')[0])
        self.assertEqual(sum(map(int, rows[1].split('\t')[1:])), 7)
        self.assertIn('phiX', (self.out/'LotuSLogS/OTU.contaminants.fa').read_text())

    def clean_table(self, rows):
        fa = self.root/'otus.fna'; fa.write_text('>good\nACGT\n>zero\nACTG\n')
        table = self.root/'table.tsv'; table.write_text('OTU\tsample\n'+rows)
        self.env['AUDIT_FASTA']=str(fa); self.env['AUDIT_TABLE']=str(table)
        return '''ensure_dir($logDir); $extendedLogs = 0; clean_otu_mat($ENV{AUDIT_FASTA}, $ENV{AUDIT_TABLE}, {});''', table

    def test_one_nonzero_and_one_zero_otu_is_not_empty(self):
        body, table = self.clean_table('good\t7\nzero\t0\n')
        self.probe(body)
        self.assertEqual(table.read_text().splitlines(), ['OTU\tsample','OTU1\t7'])

    def test_all_zero_otu_matrix_rejected(self):
        body, table = self.clean_table('zero\t0\n')
        self.assertIn('Empty OTU matrix', self.probe(body, ok=False).stdout)

    def test_lambda_index_failure_preserves_input_name_and_content(self):
        binary = self.tools/'lambda3'; binary.write_text('#!/usr/bin/python3\nimport sys\nif "--version" in sys.argv: print("lambda3 version: 3.0.0")\nelse: sys.exit(17)\n'); binary.chmod(0o755)
        self.env['AUDIT_LAMBDA']=str(binary)
        original = self.ref.read_bytes()
        body = '''ensure_dir($logDir); $lotus_tempDir = "$outdir/tmp"; ensure_dir($lotus_tempDir);
$doBlasting = 2; $lambda3Bin = $ENV{AUDIT_LAMBDA}; $BlastCores = 1;
doDBblasting($refDBwanted, $refDBwanted, "$lotus_tempDir/tax.out");'''
        self.probe(body, ok=False)
        self.assertTrue(self.ref.exists())
        self.assertEqual(self.ref.read_bytes(), original)
        self.assertFalse(self.ref.with_suffix('.fa').exists())

    def test_lambda_index_refresh_preserves_other_database_files(self):
        binary = self.tools/'lambda3'; binary.write_text('#!/usr/bin/python3\nimport sys\nif "--version" in sys.argv: print("lambda3 version: 3.0.0")\nelse: sys.exit(17)\n'); binary.chmod(0o755)
        self.env['AUDIT_LAMBDA']=str(binary)
        sidecars = [Path(str(self.ref)+suffix) for suffix in ('.tax','.notes','.dna5.fm.sa.val')]
        for path in sidecars: path.write_text('preserve me\n')
        current_index=Path(str(self.ref)+'.lba.gz'); current_index.write_bytes(b'x'*101)
        body = '''ensure_dir($logDir); $lotus_tempDir = "$outdir/tmp"; ensure_dir($lotus_tempDir);
$doBlasting = 2; $lambda3Bin = $ENV{AUDIT_LAMBDA}; $BlastCores = 1;
doDBblasting($refDBwanted, $refDBwanted, "$lotus_tempDir/tax.out");'''
        # Encountering a legacy index must not remove arbitrary DB.* sidecars.
        self.probe(body, ok=False)
        self.assertTrue(current_index.exists())
        for path in sidecars: self.assertEqual(path.read_text(), 'preserve me\n')
        # Explicit rebuild clears only the current Lambda3 index.
        self.probe(body, ['-recalcTaxDB','1'], ok=False)
        self.assertFalse(current_index.exists())
        for path in sidecars: self.assertEqual(path.read_text(), 'preserve me\n')

    def test_tree_building_without_extended_logs(self):
        mafft = self.tools/'mafft'; mafft.write_text('#!/usr/bin/python3\nimport pathlib,sys\nprint(pathlib.Path(sys.argv[-1]).read_text())\n'); mafft.chmod(0o755)
        tree = self.tools/'FastTree'; tree.write_text('#!/usr/bin/python3\nimport pathlib,sys\npathlib.Path(sys.argv[sys.argv.index("-out")+1]).write_text("(ASV1,ASV2);\\n")\n'); tree.chmod(0o755)
        self.env['AUDIT_MAFFT']=str(mafft); self.env['AUDIT_TREE']=str(tree)
        body = '''ensure_dir($logDir); $buildPhylo = 1; $extendedLogs = 0; $uthreads = 1;
$mafftBin = $ENV{AUDIT_MAFFT}; $fasttreeBin = $ENV{AUDIT_TREE};
$lotus_tempDir = "$outdir/tmp"; ensure_dir($lotus_tempDir);
my $tree = buildTree($refDBwanted, $outdir); die "Missing tree" unless -s $tree;'''
        self.probe(body, ['-extendedLogs','0'])
        self.assertTrue((self.out/'ExtraFiles/OTU.MSA.fna').exists())

    def test_complete_ont_pipeline_outputs(self):
        self.complete_pipeline()

    def test_complete_barbell_pipeline_outputs_and_cleanup(self):
        self.write_map(barbell=True)
        original = self.raw.read_bytes()
        original_map = self.map.read_bytes()
        self.complete_pipeline(barbell=True)
        self.assertEqual(self.raw.read_bytes(), original)
        self.assertEqual(self.map.read_bytes(), original_map)
        self.assertIn('Original input: '+str(self.raw), (self.out/'LotuSLogS/run_manifest.txt').read_text())

    def complete_pipeline(self, barbell=False):
        self.script.write_text((ROOT/'lotus3').read_text())
        # Keep real LCA as well as real SDM; the aligner writes its documented
        # eleven-column input, so LCA/BIOM/taxonomy aggregation run unmodified.
        (self.tools/'LCA').unlink(); (self.tools/'LCA').symlink_to(ROOT/'bin/LCA')
        vsearch = (self.tools/'vsearch').read_text()
        checkpoint = "elif name in ('minimap2', 'vsearch'):"
        vsearch = vsearch.replace(checkpoint, '''elif name == 'vsearch' and '--makeudb_usearch' in a:
    pathlib.Path(arg('-output')).write_text('synthetic test index')
elif name == 'vsearch' and '-userout' in a:
    entries = pathlib.Path(arg('--usearch_global')).read_text().split('>')[1:]
    with pathlib.Path(arg('-userout')).open('w') as out:
        for entry in entries:
            lines = entry.splitlines(); rid = lines[0].split()[0]; n=len(''.join(lines[1:]))
            out.write(f'{rid}\\tref\\t99.9\\t{n}\\t1\\t0\\t1\\t{n}\\t1\\t{n}\\t{n}\\n')
''' + checkpoint)
        (self.tools/'vsearch').write_text(vsearch)
        self.run_lotus(["-ontMinReads","2"] if barbell else [], barbell=barbell)
        table=(self.out/'OTU.txt').read_text().splitlines()
        self.assertEqual(table, ['OTU\ts1\ts2','ASV1\t4\t3'])
        biom=json.loads((self.out/'OTU.biom').read_text())
        self.assertEqual(biom['shape'], [1,2]); self.assertEqual(biom['data'], [[4,3]])
        self.assertEqual(biom['rows'][0]['id'], 'ASV1')
        self.assertEqual(biom['rows'][0]['metadata']['taxonomy'][0], 'k__Bacteria')
        self.assertIn(self.consensus, (self.out/'OTU.fna').read_text().replace('\n',''))
        self.assertIn('Bacteria', (self.out/'higherLvl/Phylum.txt').read_text())
        self.assertTrue((self.out/'LotuSLogS/run_manifest.txt').exists())
        self.assertFalse((self.out/'tmpFiles').exists())


if __name__ == '__main__':
    unittest.main()

#!/usr/bin/env python3
"""Exercise the installer CLI in temporary installations, without network/builds.
Download tests substitute fixture digests in a temporary script copy, leaving
checksum verification, unpacking, executable probing, and config writes intact.
"""
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import tarfile
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


def program(name, crash=False):
    return '#!/usr/bin/python3\n' + ('import os, signal\nos.kill(os.getpid(), signal.SIGTERM)\n' if crash else f'''import sys
if '--version' in sys.argv:
    print({repr('2.28-r1209' if name == 'minimap2' else name + ' 0.7.0')})
else:
    print('--kit --input --output --maximize --threads --quality-value-cutoff --minimum-base-quality --chimera-allowable-errors --single-strand')
''')


class ONTInstaller(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='lotus-ont-installer-test-')
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.install = self.root/'lotus'; self.install.mkdir()
        (self.install/'helpers').mkdir(); (self.install/'configs').mkdir()
        (self.install/'lotus3').write_text('# test installation root\n')
        self.installer = self.install/'helpers/autoInstall.pl'
        shutil.copyfile(ROOT/'helpers/autoInstall.pl', self.installer)
        self.default = self.install/'configs/LotuS.cfg.def'
        self.default.write_text('UID ??\nusearch unavailable\n# preserve this\nTAX_REFDB_KSGP /my/reference.fasta\n')
        self.cfg = self.install/'lOTUs.cfg'
        self.tools = self.root/'tools'; self.tools.mkdir()
        for name in ('tar','chmod','gzip','bzip2'):
            (self.tools/name).symlink_to(shutil.which(name))
        self.env = dict(os.environ, PATH=str(self.tools))

    def tool(self, name, contents=None, **kwargs):
        p = self.tools/name; p.write_text(contents or program(name, **kwargs)); p.chmod(0o755)
        return p

    def installed_tools(self):
        for name in ('minimap2','savont','barbell'): self.tool(name)

    def run_installer(self, extra=(), ok=True):
        result = subprocess.run(['/usr/bin/perl',str(self.installer),'--ont-only',*extra], env=self.env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=30)
        if ok: self.assertEqual(result.returncode, 0, result.stdout)
        else: self.assertNotEqual(result.returncode, 0, result.stdout)
        return result

    def entries(self):
        return [l.split(None,1) for l in self.cfg.read_text().splitlines() if l and not l.startswith('#')]

    def test_new_install_registers_absolute_paths(self):
        self.installed_tools(); self.run_installer()
        config = dict(self.entries())
        for name in ('minimap2','savont','barbell'):
            self.assertEqual(config[name], str(self.tools/name))
        self.assertEqual(config['TAX_REFDB_KSGP'], '/my/reference.fasta')
        self.assertFalse((self.install/'DB').exists())

    def test_old_config_duplicates_backup_and_rerun(self):
        self.installed_tools()
        original = self.default.read_text()+'savont missing\nsavont duplicate\nbarbell missing\n# no final newline'
        self.cfg.write_text(original)
        self.run_installer()
        self.assertEqual(Path(str(self.cfg)+'.bak').read_text(), original)
        saved = self.cfg.read_text()
        self.run_installer()
        self.assertEqual(self.cfg.read_text(), saved)
        for name in ('savont','barbell','minimap2'):
            self.assertEqual(sum(key == name for key,value in self.entries()), 1)
        self.assertIn('# no final newline\n', saved)

    def test_configured_relative_paths_are_reused(self):
        self.installed_tools(); (self.install/'custom').mkdir()
        savont = self.install/'custom/my-savont'; savont.write_text(program('savont')); savont.chmod(0o755)
        self.cfg.write_text(self.default.read_text()+'savont custom/my-savont\n')
        self.run_installer()
        self.assertEqual(dict(self.entries())['savont'], str(savont))

    def test_missing_build_tools_preserves_config(self):
        self.tool('minimap2'); self.tool('barbell')
        self.cfg.write_bytes(self.default.read_bytes()); before = self.cfg.read_bytes()
        result = self.run_installer(ok=False)
        self.assertIn('ONT source builds require cargo', result.stdout)
        self.assertEqual(self.cfg.read_bytes(), before)
        self.assertFalse(Path(str(self.cfg)+'.bak').exists())

    def test_crashing_program_is_not_registered(self):
        self.installed_tools(); self.tool('savont', crash=True)
        result = self.run_installer(ok=False)
        self.assertIn('require cargo', result.stdout)
        self.assertFalse(self.cfg.exists())

    def test_rejects_conflicting_modes(self):
        self.installed_tools()
        self.assertIn('cannot be combined', self.run_installer(['--condaDBinstall'], ok=False).stdout)
        self.assertFalse(self.cfg.exists())

    def setup_downloads(self, corrupt=False):
        fixtures = self.root/'fixtures'; fixtures.mkdir()
        data = {}
        for name,version in [('savont','0.7.0'),('minimap2','2.28_x64-linux')]:
            buf = io.BytesIO()
            with tarfile.open(fileobj=buf, mode='w:gz' if name=='savont' else 'w:bz2') as tar:
                for filename, text in ([(f'{name}-{version}/Cargo.toml','[package]\nname="savont"\n'),(f'{name}-{version}/Cargo.lock','version = 3\n')] if name=='savont' else [(f'{name}-{version}/minimap2',program(name))]):
                    info = tarfile.TarInfo(filename); b=text.encode(); info.size=len(b); info.mode=0o755; tar.addfile(info, io.BytesIO(b))
            data[name]=buf.getvalue()
        data['barbell']=program('barbell').encode()
        for name,b in data.items(): (fixtures/name).write_bytes(b)
        source = self.installer.read_text()
        for name,old in [('savont','a60141fd4d4e83cdcbc3601220dc41487f4cdf7bd81558c7125e876117ed756c'),('barbell','f23a599eceb8b27211178facf5504935a98880c031eb9f2b4f136382193f2081'),('minimap2','51f2cf0e486d0f9f88ace1aa58fdc56571382a676ea0889ae607301c60693377')]:
            source = source.replace(old, hashlib.sha256(data[name]).hexdigest())
        self.installer.write_text(source)
        if corrupt: (fixtures/'barbell').write_text('corrupt download')
        self.env['ONT_FIXTURES'] = str(fixtures)
        self.env['ONT_BUILD_LOG'] = str(self.root/'build.json')
        self.tool('wget', '''#!/usr/bin/python3
import os, pathlib, shutil, sys
url=sys.argv[-1]
name = 'savont' if '/savont/' in url else 'barbell' if '/barbell/' in url else 'minimap2'
shutil.copyfile(pathlib.Path(os.environ['ONT_FIXTURES'])/name, sys.argv[sys.argv.index('-O')+1])
''')
        self.tool('rustc', '#!/usr/bin/python3\nprint("rustc 1.88.0")\n')
        for name in ('cc','c++','cmake'): self.tool(name, '#!/bin/sh\nexit 0\n')
        self.tool('cargo', '''#!/usr/bin/python3
import json, os, pathlib, sys
args=sys.argv[1:]
pathlib.Path(os.environ['ONT_BUILD_LOG']).write_text(json.dumps(args))
p=pathlib.Path(args[args.index('--target-dir')+1])/'release/savont'
p.parent.mkdir(parents=True); p.write_text('''+repr(program('savont'))+'''); p.chmod(0o755)
''')

    def test_download_build_verify_and_register(self):
        self.setup_downloads(); self.run_installer()
        config = dict(self.entries())
        for name in ('minimap2','savont','barbell'):
            self.assertEqual(config[name], str(self.install/'bin'/name))
            self.assertTrue(os.access(config[name], os.X_OK))
        args=json.loads((self.root/'build.json').read_text())
        self.assertIn('--locked', args); self.assertIn('--release', args)
        self.assertEqual(sorted(p.name for p in (self.install/'bin').iterdir()), ['barbell','minimap2','savont'])

    def test_checksum_failure_keeps_old_tool_entry(self):
        self.setup_downloads(corrupt=True)
        self.cfg.write_text(self.default.read_text()+'barbell previous-install\n')
        result = self.run_installer(ok=False)
        self.assertIn('Checksum mismatch', result.stdout)
        self.assertEqual(dict(self.entries())['barbell'], 'previous-install')
        self.assertFalse((self.install/'bin/barbell').exists())


if __name__ == '__main__':
    unittest.main()

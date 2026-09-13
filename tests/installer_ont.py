#!/usr/bin/env python3
"""Exercise the installer CLI in temporary installations, without network/builds.
Download tests substitute fixture digests in a temporary script copy, leaving
checksum verification, unpacking, executable probing, and config writes intact.
"""
import hashlib
import ctypes
import ctypes.util
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import tarfile
import tempfile
import unittest
import zipfile

ROOT = Path(__file__).resolve().parents[1]


def conda_package(contents=None, member='bin/savont', link=False, duplicate=False):
    """Small real ZIP/Zstandard/tar archive; no Conda or Rust needed."""
    payload = io.BytesIO()
    with tarfile.open(fileobj=payload, mode='w') as archive:
        info = tarfile.TarInfo(member)
        data = (contents if contents is not None else program('savont')).encode()
        if link:
            info.type = tarfile.LNKTYPE if link == 'hard' else tarfile.SYMTYPE
            info.linkname = '/bin/sh'
        else:
            info.size = len(data)
        for _ in range(2 if duplicate else 1):
            archive.addfile(info, io.BytesIO(data))
    data = payload.getvalue()
    name = ctypes.util.find_library('zstd')
    if not name:
        raise unittest.SkipTest('libzstd is needed to create .conda fixtures')
    lib = ctypes.CDLL(name)
    lib.ZSTD_compressBound.argtypes = [ctypes.c_size_t]
    lib.ZSTD_compressBound.restype = ctypes.c_size_t
    lib.ZSTD_compress.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p, ctypes.c_size_t, ctypes.c_int]
    lib.ZSTD_compress.restype = ctypes.c_size_t
    lib.ZSTD_isError.argtypes = [ctypes.c_size_t]
    lib.ZSTD_isError.restype = ctypes.c_uint
    capacity = lib.ZSTD_compressBound(len(data))
    compressed = ctypes.create_string_buffer(capacity)
    size = lib.ZSTD_compress(compressed, capacity, data, len(data), 3)
    if lib.ZSTD_isError(size):
        raise RuntimeError('Unable to compress test fixture')
    package = io.BytesIO()
    with zipfile.ZipFile(package, 'w') as archive:
        archive.writestr('metadata.json', '{"conda_pkg_format_version": 2}')
        archive.writestr('pkg-savont-fixture.tar.zst', compressed.raw[:size])
    return package.getvalue()


def zstd_tool():
    path = os.environ.get('LOTUS_TEST_ZSTD') or shutil.which('zstd')
    if not path or not os.path.isfile(path) or not os.access(path, os.X_OK):
        raise unittest.SkipTest('zstd is needed for extraction tests; set LOTUS_TEST_ZSTD or add it to PATH')
    return str(Path(path).resolve())


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
        shutil.copyfile(ROOT/'helpers/extract_conda_executable.pl', self.install/'helpers/extract_conda_executable.pl')
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

    def run_installer(self, extra=(), ok=True, *, ont_only=True, answers=None):
        mode = ['--ont-only'] if ont_only else []
        result = subprocess.run(['/usr/bin/perl',str(self.installer),*mode,*extra], input=answers, env=self.env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=30)
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

    def test_default_install_downloads_barbell_without_rust(self):
        self.setup_downloads()
        self.tool('minimap2'); self.tool('savont')
        for name in ('cargo', 'rustc', 'cc', 'c++', 'cmake'):
            (self.tools/name).unlink()
        self.run_installer()
        self.assertEqual(dict(self.entries())['barbell'], str(self.install/'bin/barbell'))
        self.assertTrue(os.access(self.install/'bin/barbell', os.X_OK))
        self.assertFalse((self.root/'build.json').exists())

    def test_default_install_reuses_configured_barbell(self):
        self.installed_tools(); self.tool('barbell', crash=True)
        configured = self.install/'custom-barbell'
        configured.write_text(program('barbell')); configured.chmod(0o755)
        self.cfg.write_text(self.default.read_text()+'barbell custom-barbell\n')
        result = self.run_installer()
        self.assertEqual(dict(self.entries())['barbell'], str(configured))
        self.assertNotIn(str(self.tools/'barbell'), result.stdout)

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
        self.tool('rustc', '#!/bin/sh\necho "rustc 1.88.0"\n')
        self.tool('cargo')
        self.cfg.write_bytes(self.default.read_bytes()); before = self.cfg.read_bytes()
        result = self.run_installer(ok=False)
        self.assertIn('ONT source builds require', result.stdout)
        self.assertEqual(self.cfg.read_bytes(), before)
        self.assertFalse(Path(str(self.cfg)+'.bak').exists())

    def test_crashing_program_is_not_registered(self):
        self.installed_tools(); self.tool('savont', crash=True)
        result = self.run_installer(ok=False)
        self.assertIn('install the zstd command-line tool', result.stdout)
        self.assertIn('Install Rust >= 1.88', result.stdout)
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

    def setup_bioconda(self, payload=None):
        self.setup_downloads()
        for name in ('cargo', 'rustc', 'cc', 'c++', 'cmake'):
            (self.tools/name).unlink()
        (self.tools/'zstd').symlink_to(zstd_tool())
        payload = conda_package() if payload is None else payload
        (self.root/'fixtures/savont').write_bytes(payload)
        self.installer.write_text(self.installer.read_text().replace(
            'e7ea28b084d176379d9fa273a3cb349c9e58e54436e2efd02c295acf91edbdd7',
            hashlib.sha256(payload).hexdigest()))
        self.cfg.write_text(self.default.read_text()+'savont previous-install\n')

    def assert_bioconda_failure(self, result):
        self.assertIn('Savont Bioconda binary fallback failed', result.stdout)
        self.assertIn('Install Rust >= 1.88 (including Cargo)', result.stdout)
        self.assertEqual(dict(self.entries())['savont'], 'previous-install')
        self.assertFalse((self.install/'bin/savont').exists())
        self.assertFalse((self.install/'bin/barbell').exists())
        self.assertFalse((self.root/'build.json').exists())

    def test_bioconda_without_rust_installs_and_reruns(self):
        self.setup_bioconda()
        for name in ('python', 'python3', 'conda', 'cargo', 'rustc'):
            self.assertIsNone(shutil.which(name, path=self.env['PATH']))
        self.run_installer()
        for name in ('minimap2', 'savont', 'barbell'):
            self.assertEqual(dict(self.entries())[name], str(self.install/'bin'/name))
        self.assertFalse((self.root/'build.json').exists())
        saved = self.cfg.read_bytes()
        (self.tools/'wget').unlink(); (self.tools/'zstd').unlink()
        self.run_installer()
        self.assertEqual(self.cfg.read_bytes(), saved)

    def test_bioconda_with_old_rust(self):
        self.setup_bioconda()
        self.tool('cargo', '#!/bin/sh\nexit 1\n')
        self.tool('rustc', '#!/bin/sh\necho "rustc 1.87.0"\n')
        self.run_installer()
        self.assertEqual(dict(self.entries())['savont'], str(self.install/'bin/savont'))

    def test_bioconda_without_cargo(self):
        self.setup_bioconda()
        self.tool('rustc', '#!/bin/sh\necho "rustc 1.88.0"\n')
        self.run_installer()
        self.assertEqual(dict(self.entries())['savont'], str(self.install/'bin/savont'))

    def test_bioconda_corrupt_download_aborts(self):
        self.setup_bioconda()
        (self.root/'fixtures/savont').write_bytes(b'corrupt download')
        result = self.run_installer(ok=False)
        self.assert_bioconda_failure(result)
        self.assertIn('Checksum mismatch', result.stdout)

    def test_bioconda_failed_download_aborts(self):
        self.setup_bioconda()
        self.tool('minimap2')
        self.tool('wget', '#!/bin/sh\nexit 1\n')
        self.assert_bioconda_failure(self.run_installer(ok=False))

    def test_bioconda_bad_archive_aborts(self):
        self.setup_bioconda(b'not a ZIP archive')
        result = self.run_installer(ok=False)
        self.assert_bioconda_failure(result)
        self.assertIn('Cannot extract Bioconda executable', result.stdout)

    def test_bioconda_crashing_executable_aborts(self):
        self.setup_bioconda(conda_package(program('savont', crash=True)))
        self.assert_bioconda_failure(self.run_installer(ok=False))

    def test_bioconda_missing_cli_flags_aborts(self):
        self.setup_bioconda(conda_package('#!/bin/sh\necho "savont 0.7.0"\n'))
        self.assert_bioconda_failure(self.run_installer(ok=False))

    def test_bioconda_missing_zstd_aborts_before_config_update(self):
        self.setup_bioconda(); (self.tools/'zstd').unlink()
        before = self.cfg.read_bytes()
        self.assert_bioconda_failure(self.run_installer(ok=False))
        self.assertEqual(self.cfg.read_bytes(), before)

    def test_bioconda_unsupported_platform_aborts(self):
        self.setup_bioconda()
        self.installer.write_text(self.installer.read_text().replace('my $arch = lc($host[4]);', "my $arch = 'unsupported';"))
        result = self.run_installer(ok=False)
        self.assert_bioconda_failure(result)
        self.assertIn('No pinned Savont Bioconda binary', result.stdout)


class CondaExtraction(unittest.TestCase):
    def extract(self, package, destination, script=None):
        env = dict(os.environ, PATH=str(Path(zstd_tool()).parent))
        return subprocess.run(['/usr/bin/perl', str(script or ROOT/'helpers/extract_conda_executable.pl'), str(package), 'bin/savont', str(destination)], env=env, capture_output=True, text=True)

    def test_rejects_missing_or_linked_executable(self):
        for kwargs in ({'member': 'bin/other'}, {'link': True}, {'link': 'hard'}, {'contents': ''}, {'duplicate': True}):
            with self.subTest(kwargs=kwargs), tempfile.TemporaryDirectory() as temp:
                package = Path(temp)/'savont.conda'; destination = Path(temp)/'savont'
                package.write_bytes(conda_package(**kwargs))
                result = self.extract(package, destination)
                self.assertNotEqual(result.returncode, 0, result.stdout)
                self.assertFalse(destination.exists())

    def test_refuses_to_overwrite_existing_file(self):
        for linked in (False, True):
            with self.subTest(linked=linked), tempfile.TemporaryDirectory() as temp:
                package = Path(temp)/'savont.conda'; destination = Path(temp)/'savont'
                package.write_bytes(conda_package())
                target = Path(temp)/'original'; target.write_text('preserve')
                if linked: destination.symlink_to(target)
                else: destination.write_text('preserve')
                result = self.extract(package, destination)
                self.assertNotEqual(result.returncode, 0, result.stdout)
                self.assertEqual(destination.read_text(), 'preserve')
                self.assertEqual(target.read_text(), 'preserve')

    def test_rejects_duplicate_zip_payloads(self):
        with tempfile.TemporaryDirectory() as temp:
            package = Path(temp)/'savont.conda'; destination = Path(temp)/'savont'
            package.write_bytes(conda_package())
            with zipfile.ZipFile(package, 'a') as archive:
                archive.writestr('pkg-second.tar.zst', archive.read('pkg-savont-fixture.tar.zst'))
            result = self.extract(package, destination)
            self.assertNotEqual(result.returncode, 0, result.stdout)
            self.assertIn('exactly one pkg-', result.stderr)
            self.assertFalse(destination.exists())

    def test_zstd_failure_leaves_no_executable(self):
        with tempfile.TemporaryDirectory() as temp:
            package = Path(temp)/'savont.conda'; destination = Path(temp)/'savont'
            with zipfile.ZipFile(package, 'w') as archive:
                archive.writestr('pkg-broken.tar.zst', b'not a zstd frame')
            result = self.extract(package, destination)
            self.assertNotEqual(result.returncode, 0, result.stdout)
            self.assertIn('Zstandard decompression failed', result.stderr)
            self.assertFalse(destination.exists())

    def test_decompressed_size_limit(self):
        with tempfile.TemporaryDirectory() as temp:
            package = Path(temp)/'savont.conda'; destination = Path(temp)/'savont'
            package.write_bytes(conda_package())
            script = Path(temp)/'extract.pl'
            script.write_text((ROOT/'helpers/extract_conda_executable.pl').read_text().replace('256 * 1024 * 1024', '4096'))
            result = self.extract(package, destination, script)
            self.assertNotEqual(result.returncode, 0, result.stdout)
            self.assertIn('Oversized decompressed package', result.stderr)
            self.assertFalse(destination.exists())


if __name__ == '__main__':
    unittest.main()

#!/usr/bin/env python3
"""Exercise interactive choices and real ONT preflight/registration in temp installs.

Stop at the database-download boundary, recording selected database/program modes.
No reference databases or non-ONT packages are downloaded or built by these tests.
"""
import json
import unittest
import installer_ont as ont


class InstallerOptions(unittest.TestCase):
    setUp = ont.ONTInstaller.setUp
    tool = ont.ONTInstaller.tool
    installed_tools = ont.ONTInstaller.installed_tools
    run_installer = ont.ONTInstaller.run_installer
    entries = ont.ONTInstaller.entries

    def run_choices(self, answers, extra=(), ok=True):
        source = self.installer.read_text()
        checkpoint = '''
require JSON::PP;
print "CHOICES:" . JSON::PP->new->encode({
    search => $installBlast, databases => \\@refDBinstall, its => $ITSready,
    utax => $getUTAX, ont => $installONT, programs => [ont_programs_to_install()],
    db_only => ($onlyDbinstall || $condaDBinstall) ? 1 : 0,
    r_packages => $install_dada,
}) . "\\n";
install_ont_programs() unless $onlyDbinstall || $condaDBinstall;
finishAI("none");
exit(0);
'''
        self.assertIn('\nget_DBs();\n', source)
        self.installer.write_text(source.replace('\nget_DBs();\n', checkpoint, 1))
        result = self.run_installer(extra, ok, ont_only=False, answers=answers)
        state = next((json.loads(line.partition('CHOICES:')[2]) for line in result.stdout.splitlines() if line.startswith('CHOICES:')), None)
        return result, state

    def prepare_programs(self):
        self.installed_tools()
        self.tool('Rscript', '#!/bin/sh\necho "R scripting front-end version 4.4.0"\n')

    def assert_all(self, answer, prefix=''):
        self.prepare_programs()
        result, state = self.run_choices(prefix+answer+'\ny\n')
        self.assertEqual(int(state['search']), 3)
        self.assertEqual(state['databases'], [0]*8+[1, 0])
        for field in ('its', 'utax', 'ont', 'r_packages'):
            self.assertEqual(int(state[field]), 1)
        self.assertEqual(state['programs'], ['minimap2', 'savont', 'barbell'])
        for name in state['programs']:
            self.assertEqual(dict(self.entries())[name], str(self.tools/name))
        for question in ('For similarity based', 'Do you want to install a reference database', '-- ITS --', '-- UTAX --', '-- ONT --'):
            self.assertNotIn(question, result.stdout)
        self.assertIn('Do you accept (y/n)?', result.stdout)
        self.assertIn('Install LotuS3 with all possible dependencies', result.stdout)
        self.assertNotIn('Continue (y/n)?', result.stdout)

    def test_enter_selects_every_dependency(self):
        self.assert_all('')

    def test_one_selects_every_dependency(self):
        self.assert_all('1')

    def test_refresh_programs_offers_all_dependencies(self):
        self.cfg.write_text(self.default.read_text().replace('UID ??', 'UID 123'))
        self.assert_all('1', prefix='1\n')

    def test_detailed_without_ont_preserves_entries_and_skips_rust_checks(self):
        self.prepare_programs()
        self.tool('savont', crash=True); self.tool('barbell', crash=True)
        self.cfg.write_text(self.default.read_text()+'savont previous-savont\nbarbell previous-barbell\n')
        result, state = self.run_choices('0\n2\n0\n0\n0\n0\n')
        self.assertEqual(state['programs'], ['minimap2'])
        self.assertEqual(int(state['search']), 2)
        self.assertEqual(state['databases'], [1]+[0]*9)
        for field in ('its', 'utax', 'ont'):
            self.assertEqual(int(state[field]), 0)
        for name in ('savont', 'barbell'):
            self.assertEqual(dict(self.entries())[name], 'previous-'+name)
            self.assertNotIn(str(self.tools/name), result.stdout)
            self.assertFalse((self.install/'bin'/name).exists())
        self.assertEqual(dict(self.entries())['minimap2'], str(self.tools/'minimap2'))
        self.assertIn('-- ONT -- Install ONT tools', result.stdout)
        self.assertNotIn('Bioconda', result.stdout)

    def assert_detailed_ont(self, answer):
        self.prepare_programs()
        result, state = self.run_choices('0\n1\n1\n1\n0\n'+answer+'\n')
        self.assertEqual(int(state['search']), 1)
        self.assertEqual(state['databases'], [0, 1]+[0]*8)
        self.assertEqual(int(state['its']), 1)
        self.assertEqual(int(state['utax']), 0)
        self.assertEqual(int(state['ont']), 1)
        self.assertEqual(state['programs'], ['minimap2', 'savont', 'barbell'])
        for name in state['programs']:
            self.assertEqual(dict(self.entries())[name], str(self.tools/name))
        self.assertIn('-- ONT -- Install ONT tools', result.stdout)
        self.assertNotIn('Do you accept (y/n)?', result.stdout)

    def test_detailed_with_ont_registers_all_three_tools(self):
        self.assert_detailed_ont('1')

    def test_detailed_ont_defaults_to_yes(self):
        self.assert_detailed_ont('')

    def test_invalid_choices_are_reprompted(self):
        self.prepare_programs()
        result, state = self.run_choices('2\nbad\n0\n2\n1\n1\n0\n2\n\n')
        self.assertEqual(result.stdout.count('Invalid answer;'), 3)
        self.assertEqual(int(state['search']), 2)
        self.assertEqual(int(state['ont']), 1)

    def test_all_dependencies_still_requires_silva_acceptance(self):
        result, state = self.run_choices('\nn\n', ok=False)
        self.assertIsNone(state)
        self.assertIn('You need to accept the SILVA license', result.stdout)
        self.assertFalse(self.cfg.exists())
        self.assertFalse((self.install/'DB').exists())

    def test_end_of_input_is_not_default_acceptance(self):
        original = self.installer.read_text()
        for answers, context in [('', 'all-dependencies choice'), ('\n', 'SILVA license response'), ('0\n2\n0\n0\n0\n', 'ONT tools choice')]:
            with self.subTest(context=context):
                self.installer.write_text(original)
                result, state = self.run_choices(answers, ok=False)
                self.assertIsNone(state)
                self.assertIn('End of input while waiting for the '+context, result.stdout)
                self.assertFalse(self.cfg.exists())

    def test_database_only_refresh_has_no_program_questions(self):
        self.cfg.write_text(self.default.read_text().replace('UID ??', 'UID 123'))
        result, state = self.run_choices('2\n8\ny\n1\n1\n')
        self.assertEqual(state['db_only'], 1)
        self.assertEqual(state['databases'], [0]*8+[1, 0])
        self.assertNotIn('all possible dependencies', result.stdout)
        self.assertNotIn('-- ONT --', result.stdout)
        self.assertNotIn('barbell', dict(self.entries()))
        self.assertFalse((self.install/'bin/barbell').exists())

    def test_conda_database_mode_stays_noninteractive(self):
        result, state = self.run_choices('', extra=['--condaDBinstall'])
        self.assertEqual(state['db_only'], 1)
        self.assertNotIn('Answer:', result.stdout)
        self.assertNotIn('all possible dependencies', result.stdout)
        self.assertNotIn('barbell', dict(self.entries()))


if __name__ == '__main__':
    unittest.main()

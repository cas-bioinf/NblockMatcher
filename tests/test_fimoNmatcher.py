import unittest
import subprocess
import os
import pandas as pd
import re
from src.fimoNmatcher import Diagram


class TestNblockMatcher(unittest.TestCase):
    def setUp(self):
        self.testdata_dir = os.path.join(os.path.dirname(__file__), 'testdata')
        self.motif_file = os.path.join(self.testdata_dir, 'tm.meme')
        self.fasta_file = os.path.join(self.testdata_dir, 't.fa')
        self.gff_file = os.path.join(self.testdata_dir, 'test.gff')
        self.output_file = os.path.join(self.testdata_dir, 'output.csv')
        self.vis_prefix = os.path.join(self.testdata_dir, 'test_vis')
        self.ref_output_file = os.path.join(self.testdata_dir, 'ref_output.csv')
        self.base = '../src/fimoNmatcher.py'

    def tearDown(self):
        if os.path.exists(self.output_file):
            os.remove(self.output_file)

    def cmp2files(self, f1, f2):
        with open(f1) as fd1, open(f2) as fd2:
            self.assertListEqual(list(fd1), list(fd2))

    def test_run_shell(self):
        o = subprocess.check_output(f'python3 {self.base} -h', shell=True)
        self.assertIn("usage: fimoNmatcher.py", o.decode())

    def test_run_direct(self):
        o = subprocess.check_output(['python3', self.base, '-h'])
        self.assertIn("usage: fimoNmatcher.py", o.decode())

    def test_csv_output(self):
        subprocess.run([
            'python3',
            self.base,
            '--motif-file', self.motif_file,
            '--fasta', self.fasta_file,
            '--diagrams', "[A]7-9[T]7-9[C]7-9[G]",
            '--diag-match-output-file', self.output_file,
        ])
        self.cmp2files(self.ref_output_file, self.output_file)

    def test_no_pval_opt(self):
        subprocess.run([
            'python3',
            self.base,
            '--motif-file', self.motif_file,
            '--fasta', self.fasta_file,
            '--diagrams', "[A]7-9[T]7-9[C]7-9[G]",
            '--diag-match-output-file', self.output_file,
            '--no-pval'
        ])
        # check if p-val col is empty
        t = pd.read_csv(self.output_file, sep='\t')
        self.assertTrue(t.pval.isna().all())

        r = pd.read_csv(self.ref_output_file, sep='\t')
        rr = r.drop(columns='pval')
        tt = t.drop(columns='pval')
        self.assertTrue((rr==tt).all().all())

    def test_filter_score(self):
        subprocess.run([
            'python3',
            self.base,
            '--motif-file', self.motif_file,
            '--fasta', self.fasta_file,
            '--diagrams', "[A]7-9[T]7-9[C]7-9[G]",
            '--diag-match-output-file', self.output_file,
            '--no-pval',
            '--filter-score', '20'
        ])
        t = pd.read_csv(self.output_file, sep='\t')
        r = pd.read_csv(self.ref_output_file, sep='\t')
        self.assertTrue(len(r) > len(t))
        self.assertTrue( len(t.loc[t.score < 20]) == 0)
        self.assertFalse(len(r.loc[r.score < 20]) == 0)

    def test_fimo_opts(self):
        r = subprocess.check_output([
            'python3',
            self.base,
            '--motif-file', self.motif_file,
            '--fasta', self.fasta_file,
            '--diagrams', "[A]7-9[T]7-9[C]7-9[G]",
            '--diag-match-output-file', self.output_file,
            '--no-pval',
            '--fimo-options', "--thresh 0.000001"
        ])

        t = pd.read_csv(self.output_file, sep='\t')
        self.assertEqual(len(t), 0)
        self.assertIn('Warning: No motif hits found on s1 at all.', r.decode())

    def test_fimo_opts_cli_err(self):
        r = subprocess.run(' '.join([
            'python3',
            self.base,
            '--motif-file', self.motif_file,
            '--fasta', self.fasta_file,
            '--diagrams', "[A]7-9[T]7-9[C]7-9[G]",
            '--diag-match-output-file', self.output_file,
            '--no-pval',
            '--fimo-options', "--thresh 0.000001"
        ]), shell=True, stderr=subprocess.PIPE)

        self.assertIn('Error: argument --fimo-options: expected one argument', r.stderr.decode())
        print(r.stderr.decode())

    def test_pval_filter_exception(self):
        r = subprocess.run([
            'python3',
            self.base,
            '--motif-file', self.motif_file,
            '--fasta', self.fasta_file,
            '--diagrams', "[A]7-9[T]7-9[C]7-9[G]",
            '--diag-match-output-file', self.output_file,
            '--no-pval',
            '--filter-pval', '0.01'
        ], stderr=subprocess.PIPE)
        self.assertIn(
            "error: argument --filter-pval: not allowed with argument --no-pval",
            r.stderr.decode()
        )

    def test_pval_filter(self):
        r = subprocess.run([
            'python3',
            self.base,
            '--motif-file', self.motif_file,
            '--fasta', self.fasta_file,
            '--diagrams', "[A]7-9[T]7-9[C]7-9[G]",
            '--diag-match-output-file', self.output_file,
            '--filter-pval', '0.0000000001'
        ], stderr=subprocess.PIPE)
        t = pd.read_csv(self.output_file, sep='\t')
        r = pd.read_csv(self.ref_output_file, sep='\t')
        self.assertTrue(len(r) > len(t))
        self.assertTrue( len(t.loc[t.pval > 0.0000000001]) == 0)
        self.assertFalse(len(r.loc[r.pval > 0.0000000001]) == 0)

    def test_gff_annot(self):
        r = subprocess.run([
            'python3',
            self.base,
            '--motif-file', self.motif_file,
            '--fasta', self.fasta_file,
            '--diagrams', "[A]8[T]8[C]8[G]",
            '--diag-match-output-file', self.output_file,
            '--gff', self.gff_file,
        ], stderr=subprocess.PIPE)
        t = pd.read_csv(self.output_file, sep='\t')
        self.assertIn('gff', t.columns)
        self.assertTrue(t['gff'].isna().any())

    def test_vis(self):
        r = subprocess.run([
            'python3',
            self.base,
            '--motif-file', self.motif_file,
            '--fasta', self.fasta_file,
            '--diagrams', "[A]7-9[T]7-9[C]7-9[G]",
            '--diag-match-output-file', self.output_file,
            '--vis-track-prefix', self.vis_prefix,
        ], stderr=subprocess.PIPE)
        # we should have gff, html, interact, and copy of fasta
        files = [self.vis_prefix + i for i in ['.gff', '.html', '.interact', '_genome.fa']]
        for f in files:
            self.assertTrue(os.path.exists(f))

        r = pd.read_csv(self.ref_output_file, sep='\t')
        t = pd.read_csv(self.output_file, sep='\t')
        o = pd.read_csv(self.vis_prefix + '.gff', sep='\t', header=None)
        self.assertEqual(len(r), len(t))
        # in this case of 4 motifs per diagram, we have gff with 4 rows per diagram site
        self.assertEqual(len(t), len(o)/4)

        self.cmp2files(self.fasta_file, self.vis_prefix + '_genome.fa')

        # for each diagram site, we have 3 arcs
        ii = pd.read_csv(self.vis_prefix + '.interact', header=None)
        self.assertEqual(len(ii), len(r)*3)

        with open(self.vis_prefix + '.html') as f:
            txt = f.read()
            self.assertIn('NblockMatcher result visualisation', txt)
            self.assertIn('test_vis_genome.fa', txt)
            self.assertIn('test_vis.gff', txt)
            self.assertIn('test_vis.interact', txt)
            self.assertEqual(len(re.findall(r'<tr><td><span class="clickme"', txt)), len(r))

        for f in files:
            os.remove(f)

    def test_vis_gff(self):
        r = subprocess.run([
            'python3',
            self.base,
            '--motif-file', self.motif_file,
            '--fasta', self.fasta_file,
            '--diagrams', "[A]7-9[T]7-9[C]7-9[G]",
            '--diag-match-output-file', self.output_file,
            '--vis-track-prefix', self.vis_prefix,
            '--gff', self.gff_file,
        ], stderr=subprocess.PIPE)
        # we should have gff, html, interact, and copy of fasta
        files = [self.vis_prefix + i for i in ['.gff', '.html', '.interact', '_genome.fa', '_spec.gff']]
        for f in files:
            self.assertTrue(os.path.exists(f))

        r = pd.read_csv(self.ref_output_file, sep='\t')
        t = pd.read_csv(self.output_file, sep='\t')
        o = pd.read_csv(self.vis_prefix + '.gff', sep='\t', header=None)
        self.assertEqual(len(r), len(t))
        # in this case of 4 motifs per diagram, we have gff with 4 rows per diagram site
        self.assertEqual(len(t), len(o)/4)

        self.cmp2files(self.fasta_file, self.vis_prefix + '_genome.fa')
        self.cmp2files(self.gff_file, self.vis_prefix + '_spec.gff')

        # for each diagram site, we have 3 arcs
        ii = pd.read_csv(self.vis_prefix + '.interact', header=None)
        self.assertEqual(len(ii), len(r)*3)

        with open(self.vis_prefix + '.html') as f:
            self.assertIn('test_vis_spec.gff', f.read())

        for f in files:
            os.remove(f)

    def test_zero_dist(self):
        tf = os.path.join(self.testdata_dir, 'tmp.fa')
        with open(tf, 'w') as f:
            f.write('>t1\nGGGGGGAAAAAATTTTTTGGGGG\n')
            f.write('>t2_5\nGGGAAAAAAGGGGGTTTTTTGGG\n')
            f.write('>t3_4\nGGGAAAAAAGGGGTTTTTTGGG\n')
            f.write('>t4_6\nGGGAAAAAAGGGGGGTTTTTTGGG\n')
        try:
            r = subprocess.run([
                'python3',
                self.base,
                '--motif-file', self.motif_file,
                '--fasta', tf,
                '--diagrams', "[A]0-1[T]",
                '--diag-match-output-file', self.output_file,
            ])
            t = pd.read_csv(self.output_file, sep='\t')
            self.assertEqual(len(t.dists.unique()), 1)
            self.assertEqual(t.dists.unique().tolist()[0], 0)
        finally:
            os.remove(tf)

    def test_zero_dist2(self):
        tf = os.path.join(self.testdata_dir, 'tmp.fa')
        with open(tf, 'w') as f:
            f.write('>t1\nGGGGGGAAAAAATTTTTTGGGGG\n')
            f.write('>t2_5\nGGGAAAAAAGGGGGTTTTTTGGG\n')
            f.write('>t3_4\nGGGAAAAAAGGGGTTTTTTGGG\n')
            f.write('>t4_6\nGGGAAAAAAGGGGGGTTTTTTGGG\n')
        try:
            r = subprocess.run([
                'python3',
                self.base,
                '--motif-file', self.motif_file,
                '--fasta', tf,
                '--diagrams', "[A]0[T]",
                '--diag-match-output-file', self.output_file,
            ])
            t = pd.read_csv(self.output_file, sep='\t')
            self.assertEqual(len(t.dists.unique()), 1)
            self.assertEqual(t.dists.unique().tolist()[0], 0)
        finally:
            os.remove(tf)

    def test_zero_dist_rc(self):
        tf = os.path.join(self.testdata_dir, 'tmp.fa')
        with open(tf, 'w') as f:
            f.write('>t1\nGGGGGGAAAAAATTTTTTGGGGG\n')
            f.write('>t2_5\nGGGAAAAAAGGGGGTTTTTTGGG\n')
            f.write('>t3_4\nGGGAAAAAAGGGGTTTTTTGGG\n')
            f.write('>t4_6\nGGGAAAAAAGGGGGGTTTTTTGGG\n')
        try:
            r = subprocess.run([
                'python3',
                self.base,
                '--motif-file', self.motif_file,
                '--fasta', tf,
                '--diagrams', "[-T]0[-A]",
                '--diag-match-output-file', self.output_file,
            ])
            t = pd.read_csv(self.output_file, sep='\t')
            self.assertEqual(len(t.dists.unique()), 1)
            self.assertEqual(t.dists.unique().tolist()[0], 0)
            self.assertEqual(len(t.sequence_name.unique()), 1)
            self.assertEqual(t.sequence_name.unique().tolist()[0], 't1')
        finally:
            os.remove(tf)

    def test_invalid_fasta_same_ids(self):
        tf = os.path.join(self.testdata_dir, 'tmp.fa')
        with open(tf, 'w') as f:
            f.write('>t1\nGGGGGGAAAAAATTTTTTGGGGG\n')
            f.write('>t1\nGGGAAAAAAGGGGGTTTTTTGGG\n')
        try:
            r = subprocess.run([
                'python3',
                self.base,
                '--motif-file', self.motif_file,
                '--fasta', tf,
                '--diagrams', "[-T]0[-A]",
                '--diag-match-output-file', self.output_file,
            ], stderr=subprocess.PIPE)

            self.assertIn("Input FASTA file can't contain duplicated IDs", r.stderr.decode())

        finally:
            os.remove(tf)


class TestDiagram(unittest.TestCase):
    def test_ok_plus(self):
        d = Diagram('[b]6-8[a]7[c]')
        self.assertEqual(d.models, ['b', 'a', 'c'])
        self.assertEqual(d.dists, [[7,9],[8,9]])
        self.assertEqual(d.strands, ['+']*3)

    def test_ok_minus(self):
        d = Diagram('[-b]6-8[a]7[+c]')
        self.assertEqual(d.models, ['b', 'a', 'c'])
        self.assertEqual(d.dists, [[7,9],[8,9]])
        self.assertEqual(d.strands, ['-', '+', '+'])

    def test_various_incorrect_diagrams(self):
        with self.assertRaisesRegex(ValueError, 'Diagram should start/end with square brackets'):
            Diagram('asdas')

        with self.assertRaisesRegex(ValueError, 'Please provide at least two motifs.'):
            Diagram('[abc]')

        with self.assertRaisesRegex(ValueError, 'Empty motif name in diagram.'):
            Diagram('[a][]')

        with self.assertRaisesRegex(ValueError, 'Empty motif name in diagram.'):
            Diagram('[][b]')

        with self.assertRaisesRegex(ValueError, 'Empty motif name in diagram.'):
            Diagram('[+][]')

        with self.assertRaisesRegex(ValueError, 'Empty motif name in diagram.'):
            Diagram('[][-]')

        with self.assertRaisesRegex(ValueError, 'Empty motif name in diagram.'):
            Diagram('[+][-]')

        with self.assertRaisesRegex(ValueError, 'Unbalanced diagram: unexpected closing bracket'):
            Diagram('[a]4[b]]2[c]')

        with self.assertRaisesRegex(ValueError, 'Unbalanced diagram: unexpected opening bracket'):
            Diagram('[a]3[[b]5[c]')

        with self.assertRaisesRegex(ValueError, 'Incorrect distance specification'):
            Diagram('[a]a-b[b]')

        with self.assertRaisesRegex(ValueError, 'Invalid spacer range specification'):
            Diagram('[a]-1-5[b]')

        with self.assertRaisesRegex(ValueError, 'Incorrect distance specification'):
            Diagram('[a]4.5[b]')

        with self.assertRaisesRegex(ValueError, 'Incorrect distance specification, first index greater than second'):
            Diagram('[a]5-4[b]')

        with self.assertRaisesRegex(ValueError, 'Diagram should start/end with square brackets'):
            Diagram('5[a]02-4[b]')

        with self.assertRaisesRegex(ValueError, 'Diagram should start/end with square brackets'):
            Diagram('[a]02-4[b]5')

    def test_zero_length_diagrams(self):
        d1 = Diagram('[a]0[b]')
        d2 = Diagram('[a]0-1[b]')
        self.assertEqual(d1.dists, d2.dists)
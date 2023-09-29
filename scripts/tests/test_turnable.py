import sys
sys.path.append("../")
import os

from unittest import TestCase, main
from scripts.turnable import read_mutation_candidates

class Test_tuRNAble(TestCase):
    fp_tmpfile = None

    def setUp(self):
        self.fp_prefix = 'tests/data/'

    def tearDown(self):
        pass

    def test_read_mutation_candidates(self):
        obs = read_mutation_candidates(
            os.path.join(self.fp_prefix, "denovos_edited3.tsv"))

        # test if correct number of mutations are read and filtered
        self.assertEqual(obs.shape, (1093, 7))
        # if only_pointmutations=True, reference_length must be 1 everywhere...
        self.assertTrue(obs['reference_length'].unique() == [1])
        # ... and consist of only one character
        self.assertTrue(obs['reference'].apply(len).unique() == [1])

        # same as above for mutated sequence
        self.assertEqual(obs['mutation_length'].unique(), [1])
        self.assertEqual(obs['mutation'].apply(len).unique(), [1])

        # test that chr is prefixed
        self.assertEqual(
            obs['chromosome'].apply(lambda x: x[:3]).unique(),
            ['chr'])

        with self.assertRaisesRegex(ValueError, "positions are non-negative"):
            read_mutation_candidates(
                os.path.join(self.fp_prefix, "err_negative.tsv"))

        with self.assertRaisesRegex(ValueError, "non-DNA/RNA nucleotides"):
            read_mutation_candidates(
                os.path.join(self.fp_prefix, "err_nonXNA_ref.tsv"))
        with self.assertRaisesRegex(ValueError, "non-DNA/RNA nucleotides"):
            read_mutation_candidates(
                os.path.join(self.fp_prefix, "err_nonXNA_mut.tsv"))

        with self.assertRaisesRegex(ValueError, "have more then one"):
            read_mutation_candidates(
                os.path.join(self.fp_prefix, "err_multiMut.tsv"))

        obsID = read_mutation_candidates(
            os.path.join(self.fp_prefix, "denovos_edited3.tsv"),
            only_pointmutations=False)
        self.assertEqual(obsID.shape, (1208, 7))
        self.assertTrue(obs.shape[0] <= obsID.shape[0])


if __name__ == '__main__':
    main()

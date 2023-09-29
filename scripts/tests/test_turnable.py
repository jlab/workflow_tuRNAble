import sys
from unittest import TestCase, main
from phylofiller.converter import easel_table2pd
import pandas as pd
from skbio.util import get_data_path

sys.path.append("../")
from scripts.turnable import (read_mutation_candidates,  # noqa: E402
                              overlap_mutations_annotations,
                              extract_reference_subsequences,
                              mutate_sequence,
                              create_mutated_sequence)


class Test_tuRNAble(TestCase):
    fp_tmpfile = None

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_read_mutation_candidates(self):
        obs = read_mutation_candidates(get_data_path("denovos_edited3.tsv"))

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
                get_data_path("err_negative.tsv"))

        with self.assertRaisesRegex(ValueError, "non-DNA/RNA nucleotides"):
            read_mutation_candidates(
                get_data_path("err_nonXNA_ref.tsv"))
        with self.assertRaisesRegex(ValueError, "non-DNA/RNA nucleotides"):
            read_mutation_candidates(
                get_data_path("err_nonXNA_mut.tsv"))

        with self.assertRaisesRegex(ValueError, "have more then one"):
            read_mutation_candidates(
                get_data_path("err_multiMut.tsv"))

        with self.assertRaisesRegex(ValueError, 'not contain "PASS"'):
            read_mutation_candidates(
                get_data_path("err_noPASS.tsv"))

        obsID = read_mutation_candidates(
            get_data_path("denovos_edited3.tsv"),
            only_pointmutations=False)
        self.assertEqual(obsID.shape, (1208, 7))
        self.assertTrue(obs.shape[0] <= obsID.shape[0])

    def test_overlap_mutations_annotations(self):
        with open(get_data_path("cmsearch.tab"), "r") as f:
            annotations = easel_table2pd(f.readlines(), verbose=False)
        mutations = read_mutation_candidates(
            get_data_path("denovos_edited3.tsv"),
            only_pointmutations=False)

        obs = overlap_mutations_annotations(mutations, annotations,
                                            verbose=False)
        self.assertEqual(obs.shape, (7, 25))
        self.assertEqual(list(obs['#target name'].unique()),
                         ['chr16', 'chr19', 'chr4', 'chr7'])
        self.assertEqual(list(obs['query name'].unique()),
                         ["IRES_n-myc", "mir-762", "mir-1249", "isrG",
                          "mir-207"])

        self.assertTrue(len(obs.index.unique()) == obs.shape[0])

    def test_extract_reference_subsequences(self):
        res = pd.read_csv(get_data_path("pos_extract.tsv"),
                          sep="\t")

        exp = [
            'ACTCCTGTAAACAGCTGCTGCCAGCCAGGTGGTGGCCTGGCGGGGACAGCCTAGGCTCGCAGCCT'
            'CCAGGGGCACCCCCTCACCACCCCCCTGCCCTGACTCACCTTGGCCA',
            'GCCCGGCTCCGGGTCTCGGCCCGTACAGTCCGGCCGGCCATGCTGGCGGGGCTGGGGCCGGGGCC'
            'GAGCCCGCGG',
            'GCCCGGCTCCGGGTCTCGGCCCGTACAGTCCGGCCGGCCATGCTGGCGGGGCTGGGGCCGGGGCC'
            'GAGCCCGCGG',
            'GGCCTGGGGGGGAGGGAGTGTGCTCGATCCCACTGGTGGCCAAGCCCCCTCCCTCACCCTTCC',
            'AATTTGCTGCATCAATTTCACACTACTTCCATATCTAAAGAAACAAAAAAATTACCTGCTGCATA'
            'TAAGCATCTTGAAGTAGGTGGTGGTGGTGGTGGTGGTGGTGCTGCTGCTGCTGCTGCTGCTGCTG'
            'CTGTTGCTGTTGCTGCTGCTGCTGTTGCTGTTGCTGCTGCTGCTGCTGCTGCTGCTGGTGAGGAT'
            'GACGATGCTGTAAATGGAGTTGCTGTAATCT',
            'GAAGGAGGGGCCGGGCTGGGTCAGGGGCTGGGCGGGGCCGCGGCAGCCCCTGACGCCGCTCTTCC'
            'TCTCTCT',
            'GAAGGAGGGGCCGGGCTGGGTCAGGGGCTGGGCGGGGCCGCGGCAGCCCCTGACGCCGCTCTTCC'
            'TCTCTCT']
        obs = extract_reference_subsequences(
            get_data_path("ref.fasta.gz"), res, verbose=False)

        pd.testing.assert_series_equal(pd.Series(exp), obs)

        with self.assertRaisesRegex(ValueError, 'Length of extracted'):
            res.iloc[1, 2] = 5000
            extract_reference_subsequences(get_data_path("ref.fasta.gz"), res,
                                           verbose=False)

    def test_mutate_sequence(self):
        self.assertEqual(mutate_sequence("ACTGGTATGC", 4, "G", "C"),
                         ("ACTGGTATGC", "ACTGCTATGC"))
        self.assertEqual(mutate_sequence("ACTGGTATGC", 4, "G", "CCG"),
                         ("ACTGG--TATGC", "ACTGCCGTATGC"))
        self.assertEqual(mutate_sequence("ACTGGTATGC", 4, "GTA", ""),
                         ("ACTGGTATGC", "ACTG---TGC"))
        self.assertEqual(mutate_sequence("ACTGCCTGC", 4, "C", "CCCGTA"),
                         ("ACTGC-----CTGC", "ACTGCCCGTACTGC"))
        with self.assertRaisesRegex(ValueError,
                                    'reference is not of type str'):
            mutate_sequence(["ACTGGTATGC"], 4, "G", "C")
        with self.assertRaisesRegex(ValueError, 'position is larger'):
            mutate_sequence("ACTGGTATGC", 4000, "G", "C")
        with self.assertRaisesRegex(ValueError, 'position is < 0'):
            mutate_sequence("ACTGGTATGC", -4000, "G", "C")
        with self.assertRaisesRegex(ValueError, 'does not match your'):
            mutate_sequence("ACTGGTATGC", 4, "T", "C")

    def test_create_mutated_sequence(self):
        res = pd.read_csv(get_data_path("pos_extract.tsv"),
                          sep="\t")
        res['reference_sequence'] = extract_reference_subsequences(
            get_data_path("ref.fasta.gz"), res, verbose=False)
        exp_ref = [
            "ACTCCTGTAAACAGCTGCTGCCAGCCAGGTGGTGGCCTGGCGGGGACAGCCTAGGCTCGCAGCCT"
            "CCAGGGGCACCCCCTCACCACCCCCCTGCCCTGACTCACCTTGGCCA",
            "GCCCGGCTCCGGGTCTCGGCCCGTACAGTCCGGCCGGCCATGCTGGCGGGGCTGGGGCCGGGGCC"
            "GAGCCCGCGG",
            "GCCCGGCTCCGGGTCTCGGCCCGTACAGTCCGGCCGGCCATGCTGGCGGGGCTGGGGCCGGGGCC"
            "GAGCCCGCGG",
            "GGCCTGGGGGGGAGGGAGTGTGCTCGATCCCACTGGTGGCCAAGCCCCCTCCCTCACCCTTCC",
            "AATTTGCTGCATCAATTTCACACTACTTCCATATCTAAAGAAACAAAAAAATTACCTGCTGCATA"
            "TAAGCATCTTGAAGTAGGTGGTGGTGGTGGTGGTGGTGGTGCTGCTGCTGCTGCTGCTGCTGCTG"
            "CTGTTGCTGTTGCTGCTGCTGCTGTTGCTGTTGCTGCTGCTGCTG--------------------"
            "-CTGCTGCTGCTGGTGAGGATGACGATGCTGTAAATGGAGTTGCTGTAATCT",
            "GAAGGAGGGGCCGGGCTGGGTCAGGGGCTGGGCGGGGCCGCGGCAGCCCCTGACGCCGCTCTTCC"
            "TCTCTCT",
            "GAAGGAGGGGCCGGGCTGGGTCAGGGGCTGGGCGGGGCCGCGGCAGCCCCTGACGCCGCTCTTCC"
            "TCTCTCT",
        ]
        exp_mut = [
            "ACTCCTGTAAACAGCTGCTGCCAGCCAGGTGGTGGCCTGGCGGGGACAGCCTAGGCTCGCAGCCT"
            "CCAGGGGCACCCCCTCACCACCCCCCTGCCCTGCCTCACCTTGGCCA",
            "GCCCGGCTCCGGGTCTCGGCCCGTACAGTCCGGCCGGCCATGCTGGCGGGGCTGGGGCCCGGGCC"
            "GAGCCCGCGG",
            "GCCCGGCTCCGGGTCTCGGCCCGTACAGTCCGGCCGGCCATGCTGGCGGGGCTGGGGCCGAGGCC"
            "GAGCCCGCGG",
            "GGCCTGGGGGGGAGGGAGTGTGCTCGATCCCACTGGTGGCCAAGCCCCCTCCCGCACCCTTCC",
            "AATTTGCTGCATCAATTTCACACTACTTCCATATCTAAAGAAACAAAAAAATTACCTGCTGCATA"
            "TAAGCATCTTGAAGTAGGTGGTGGTGGTGGTGGTGGTGGTGCTGCTGCTGCTGCTGCTGCTGCTG"
            "CTGTTGCTGTTGCTGCTGCTGCTGTTGCTGTTGCTGCTGCTGCTGTTGCTGTTGCTGCTGCTGCT"
            "GCTGCTGCTGCTGGTGAGGATGACGATGCTGTAAATGGAGTTGCTGTAATCT",
            "GAAGGAGGGGCCGGGCTGGGTCAGGGGCTGGGCGGGGCCGCCGCAGCCCCTGACGCCGCTCTTCC"
            "TCTCTCT",
            "GAAGGAGGGGCCGGGCTGGGTCAGGGGCTGGGCGGGGCCCCGGCAGCCCCTGACGCCGCTCTTCC"
            "TCTCTCT",
        ]
        obs = create_mutated_sequence(res, verbose=False)
        pd.testing.assert_series_equal(
            pd.Series(exp_ref, name='aln_reference'),
            obs['aln_reference'])
        pd.testing.assert_series_equal(
            pd.Series(exp_mut, name='aln_mutation'),
            obs['aln_mutation'])


if __name__ == '__main__':
    main()

import pandas as pd
from tqdm import tqdm
import skbio


def read_mutation_candidates(fp_denovos: str, only_pointmutations=True):
    """
    Read trio-WES observed de-novo mutations from file.

    Parameters
    ----------
    fp_denovos : str
        Path to a tabular file with one (de-novo) mutation per row.
        Columns must be:
        1) 'chromosome' name
        2) 'start' position of the mutation (leftmost affected nucleotide)
        3) 'reference' sub-sequence or nucleotide starting from 2)
        4) 'mutation' mutated version of 'reference' sub-sequence or nucleotide
        5) 'qc' a quality control field. Must contain 'PASS' atm (2023-09-29)
        6) 'sample' free text to list the affected sample IDs. Is ignored.
    only_pointmutations : bool
        Listed mutations are typically replacements of individual bases.
        However, insertions or deletions or multiple subsequent nucloetides can
        also be handled. These larger mutations are excluded by default.
        Set only_pointmutations=False if you want to keep them for downstream
        analysis.

    Returns
    -------
    pd.DataFrame. Note that column 'qc' will be dropped. Chromosome names will
    be prefixed by 'chr' if missing in the file.

    Raises
    ------
    a. ValueError if column 'qc' is not PASS in each row
    b. ValueError if 'start' < 0
    c. ValueError if 'reference' or 'mutation' sequence is not of DNA/RNA
       alphabet or gap, see variable ALPHABET
    d. ValueError if multiple rows describe same genomic position
    """
    hits = pd.read_csv(
        fp_denovos, sep="\t", header=None, dtype=str,
        names=['chromosome', 'start', 'reference', 'mutation', 'qc', 'sample'])

    # compute length of reference and mutation region to discriminate between
    # point mutations, insertions, deletions
    for f in ['reference', 'mutation']:
        hits['%s_length' % f] = hits[f].apply(len)

    if hits['qc'].unique() != ['PASS']:
        raise ValueError(
            'Fifth column ("qc") of your file does not contain '
            '"PASS" in rows %s' % ', '.join(hits[hits['qc'] != 'PASS'].index))
    del hits['qc']

    # preprend "chr" to chromosome names if missing
    hits['chromosome'] = hits['chromosome'].apply(
        lambda x: x if x.startswith('chr') else 'chr%s' % x)

    # only focus on point mutations
    if only_pointmutations:
        hits = hits[(hits['reference_length'] == 1) &
                    (hits['mutation_length'] == 1)]

    # convert 'start' to integer
    hits['start'] = hits['start'].astype(int)
    if not all(hits['start'] >= 0):
        raise ValueError("Not all 'start' positions are non-negative!\n%s" %
                         hits[hits['start'] >= 0][['chromosome', 'start']])

    ALPHABET = ['A', 'C', 'G', 'T', 'U', '-', '_']
    for t in ['reference', 'mutation']:
        nonXNA = hits[~hits[t].apply(lambda s: all(map(
            lambda n: n in ALPHABET,
            s.upper())))]
        if nonXNA.shape[0] > 0:
            raise ValueError(
                "You have non-DNA/RNA nucleotides in your %s column!\n%s" % (
                    t,
                    nonXNA[['chromosome', 'start', t]]))

    # assert assumption that every reference position has only one mutated from
    # (not necessarily true as different patients might have different
    # mutations)
    mut_per_pos = hits.groupby(['chromosome', 'start']).size()
    if mut_per_pos.value_counts().shape[0] != 1:
        msg = '\n  '.join(mut_per_pos[mut_per_pos > 1].reset_index().apply(
            lambda row: '%ix %s:%i' % (
                row[0], row['chromosome'], row['start']), axis=1).values)
        raise ValueError("Some of your genomic positions have more then one"
                         " mutation listed:\n  %s" % msg)

    return hits


def overlap_mutations_annotations(mutations: pd.DataFrame,
                                  annotations: pd.DataFrame,
                                  verbose: bool = True):
    """Subsets mutations to those overlapping any annotation
       (at least with one nucleotide).

    Parameters
    ----------
    mutations: pd.DataFrame
        A table that hold information about (de-novo) mutations per row in the
        format of function 'read_mutation_candidates'.
    annotations: pd.DataFrame
        A table parsed from the tabular output of cmsearch, see pyhlofiller's
        easel_table2pd
    verbose: bool
        report progress

    Returns
    -------
    Slice of mutations input, with novel index as one annotation can be hit by
    multiple mutations.
    """
    overlaps = []
    for chrom, gd in tqdm(mutations.groupby('chromosome'),
                          desc='find intersection', disable=not verbose):
        chr_rfam = annotations[annotations['#target name'] == chrom]

        for idx_denovo, d in gd.iterrows():
            h = chr_rfam[((chr_rfam['seq from'] <= d['start']) &
                          (d['start'] <= chr_rfam['seq to'])) |
                         ((chr_rfam['seq from'] >= d['start']) &
                          (d['start'] >= chr_rfam['seq to'])) |
                         ((chr_rfam['seq from'] <=
                           d['start']+d['reference_length']) &
                          (d['start']+d['reference_length'] <=
                           chr_rfam['seq to'])) |
                         ((chr_rfam['seq from'] >=
                           d['start']+d['reference_length']) &
                          (d['start']+d['reference_length'] >=
                           chr_rfam['seq to']))].copy()
            if h.shape[0] > 0:
                for c in d.index:
                    h[c] = d[c]
                overlaps.append(h)
    overlaps = pd.concat(overlaps).reset_index()
    del overlaps['index']
    return overlaps


def extract_reference_subsequences(fp_reference: str, intervals: pd.DataFrame,
                                   verbose=True):
    """Extract sub-sequences for a table of genomic intervals.

    Parameters
    ----------
    fp_reference : str
        Path to a gzipped multiple fasta file, holding the full sequences, e.g.
        the human reference assembly hg38.
    intervals : pd.DataFrame
        Tabular information about genomic intervals to extract. Must have at
        least columns 'seq from', 'seq to', 'strand' and '#target name'.
        The latter values must correspond to fasta header in fp_reference.
    verbose: bool
        report progress

    Returns
    -------
    pd.Series of extracted sub-sequences with same index as 'intervals'.
    """
    res = pd.Series(index=intervals.index, dtype=str)

    # for some reason reading the sequences as DNA fails:
    # .DNA.read(, lowercase=True), i.e. seq lengths is 1
    for seq in tqdm(skbio.io.read(fp_reference, format='fasta'),
                    desc="extract reference sub-sequences",
                    disable=not verbose):
        # genomic positions between slices in skbio: seq[start:end] and tblout
        # from cmsearch have an offset of -1, i.e. we obtain the correct
        # sub-sequence via seq[start-1:end] ...

        # ... same with Layals de-novo positions: pos = seq[pos-1]
        for idx, hit in intervals[intervals['#target name'] ==
                                  seq.metadata['id']].iterrows():
            start, end = hit['seq from'], hit['seq to']
            if hit['strand'] == '-':
                # reverse complement
                start, end = end, start
            subseq = skbio.sequence.DNA(seq[start-1:end], lowercase=True)
            if hit['strand'] == '-':
                subseq = subseq.reverse_complement()
            res.loc[idx] = str(subseq)

    odd = intervals[res.apply(len) != intervals.apply(
        lambda x: abs(x['seq from'] - x['seq to']) + 1, axis=1)]
    if odd.shape[0] > 0:
        raise ValueError("Length of extracted sub-sequences do not match with"
                         " positions from cmsearch for: %s" % odd)

    return res

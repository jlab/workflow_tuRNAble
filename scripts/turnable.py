import pandas as pd


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

    ALPHABET = ['A','C','G','T','U','-','_']
    for t in ['reference', 'mutation']:
        nonXNA = hits[~hits[t].apply(lambda s: all(map(lambda n: n in ALPHABET, s.upper())))]
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



# def read_ncRNA_hits(fp_cmsearch, verbose=False, omit_nonexcl=False):
#     with open(fp_cmsearch, "r") as f:
#         rfam = easel_table2pd(f.readlines())
#     if verbose:
#         print('%i total hits' % rfam.shape[0])
#
#     if omit_nonexcl:
#         # omit all hits with non "inc" == ! entries, i.e. e-value to low
#         rfam = rfam[rfam['inc'] == '!']
#         if verbose:
#             print('%i hits with !' % rfam.shape[0])
#
#     return rfam
#
# def get_overlap_ncRNA_denovo(fp_ncRNA, fp_denovo, verbose=True):
#     denovos = read_and_subset_denovos(fp_denovo, only_pointmutations=False)
#     rfam = read_ncRNA_hits(fp_ncRNA, verbose=verbose)
#
#     hits = []
#     for chrom, gd in tqdm(denovos.groupby('chromosome'), desc='find intersection'):
#         #if verbose:
#         #    print(chrom, gd.shape[0])
#         chr_rfam = rfam[rfam['#target name'] == chrom]
#
#         for idx_denovo, d in gd.iterrows():
#             h = chr_rfam[((chr_rfam['seq from'] <= d['start']) & (d['start'] <= chr_rfam['seq to'])) |
#                          ((chr_rfam['seq from'] >= d['start']) & (d['start'] >= chr_rfam['seq to'])) |
#                          ((chr_rfam['seq from'] <= d['start']+d['reference_length']) & (d['start']+d['reference_length'] <= chr_rfam['seq to'])) |
#                          ((chr_rfam['seq from'] >= d['start']+d['reference_length']) & (d['start']+d['reference_length'] >= chr_rfam['seq to']))].copy()
#             if h.shape[0] > 0:
#                 for c in d.index:
#                     h[c] = d[c]
#                 hits.append(h)
#     hits = pd.concat(hits).reset_index()
#     del hits['index']
#     return hits
#
# def extract_reference_subsequences(fp_reference, res):
#     res['reference_sequence'] = ""
#     # for some reason reading the sequences as DNA fails: .DNA.read(, lowercase=True), i.e. seq lengths is 1
#     for seq in tqdm(skbio.io.read(fp_reference, format='fasta'), desc="extract reference sub-sequences"):
#         # genomic positions between slices in skbio: seq[start:end] and tblout from cmsearch have an offset
#         # of -1, i.e. we obtain the correct sub-sequence via seq[start-1:end] ...
#
#         # ... same with Layals de-novo positions: pos = seq[pos-1]
#         for idx, hit in res[res['#target name'] == seq.metadata['id']].iterrows():
#             start, end = hit['seq from'], hit['seq to']
#             if hit['strand'] == '-':
#                 # reverse complement
#                 start, end = end, start
#             subseq = skbio.sequence.DNA(seq[start-1:end], lowercase=True)
#             if hit['strand'] == '-':
#                 subseq = subseq.reverse_complement()
#             res.loc[idx, 'reference_sequence'] = str(subseq)
#
#     assert all(res['reference_sequence'].apply(len) == abs(res['seq from'] - res['seq to'])+1)
#     return res
#
# def mutate_sequence(reference : str, position : int, subseq : str, mutation : str):
#     """
#     Parameters
#     ----------
#     reference : str
#         The original sequence.
#     position : int
#         The left position at which a mutation occures (zero based)
#     subseq : str
#         The sub-string of the reference affected by the mutation.
#     mutation : str
#         The mutation that shall be integrated at <position> into <reference>
#
#     Returns
#     -------
#     Two alignment rows (i.e. can contain gap chars '-' in both rows) of
#     a) the reference sequence and b) a mutated version.
#     """
#     assert type(reference) == str
#     assert position < len(reference)
#     assert position > 0
#     if reference[position:position+len(subseq)] != subseq:
#         raise ValueError("sub-string '%s' starting at position %i does not match your provided subseq '%s'" %
#                         (reference[position:position+len(subseq)], position, subseq))
#
#     prefix = reference[:position]
#     suffix = reference[position+len(subseq):]
#
#     infix_mutation = mutation
#     aln_orig = prefix + subseq + ('-' * ((len(mutation) - len(subseq)))) + suffix
#     aln_mutation = prefix + mutation + ('-' * ((len(subseq) - len(mutation)))) + suffix
#     return (aln_orig, aln_mutation)
#
# assert mutate_sequence("ACTGGTATGC", 4, "G", "C") == ("ACTGGTATGC", "ACTGCTATGC")
# assert mutate_sequence("ACTGGTATGC", 4, "G", "CCG") == ("ACTGG--TATGC", "ACTGCCGTATGC")
# assert mutate_sequence("ACTGGTATGC", 4, "GTA", "") == ("ACTGGTATGC", "ACTG---TGC")
# assert mutate_sequence("ACTGCCTGC", 4, "C", "CCCGTA") == ("ACTGC-----CTGC", "ACTGCCCGTACTGC")
#
#
# def create_mutated_sequence(res):
#     res['aln_reference'] = ""
#     res['aln_mutation'] = ""
#     #res['mutation_sequence'] = ""
#     for idx, row in tqdm(res.iterrows(), desc="generate mutated sequences"):
#         left = row['seq from']
#         seq = skbio.sequence.DNA(row['reference_sequence'], lowercase=True)
#         if row['strand'] == '-':
#             left = row['seq to']
#             seq = seq.reverse_complement()
#     #    pos = row['start'] - start
#         aln_seq, aln_mut = mutate_sequence(str(seq), row['start'] - left, row['reference'], row['mutation'])
#         if row['strand'] == '-':
#             aln_seq, aln_mut = str(skbio.sequence.DNA(aln_seq).reverse_complement()), str(skbio.sequence.DNA(aln_mut).reverse_complement())
#     #    mutseq = skbio.sequence.DNA(str(seq[:pos]) + row['mutation'] + str(seq[pos+row['reference_length']:]))
#     #    if row['strand'] == '-':
#     #        mutseq = mutseq.reverse_complement()
#         res.loc[idx, 'aln_reference'] = aln_seq
#         res.loc[idx, 'aln_mutation'] = aln_mut
#     return res
#
# def bp_prob_shift(bp1, bp2, len_seq):
#     # the probability mass that shifts between both inputs
#     shift = pd.concat([bp1,bp2], axis=1).fillna(0).apply(lambda r: abs(r.iloc[0] - r.iloc[1]), axis=1).sum()
#     norm = len_seq * (len_seq+1) / 2
#
#     return shift / norm
#
# def getShapeProbs(sequence, level=3, grammar='overdangle', name=None):
#     sequence = sequence.replace('-','').upper().replace('T','U')
#     hsum = hashlib.md5(sequence.encode()).hexdigest()
#
#     fp_res = './seq%s_level%i_grammar%s.shapes' % (hsum, level, grammar)
#     cluster_run(
#         ['/vol/jlab/bin/RNAshapes --mode sample --grammar %s --shapeLevel %s "%s" > "%s"' % (grammar, level, sequence, fp_res)],
#         'SPS', fp_res, ppn=1, dry=False, wait=False, err=None)
#     if os.path.exists(fp_res):
#         res = pd.read_csv(fp_res, skiprows=2, header=None, sep="\s+").rename(columns={3: 'shape', 2: 'shape_probability'}).set_index('shape')['shape_probability']
#         if name is not None:
#             res.name = name
#         return res
#     else:
#         return None

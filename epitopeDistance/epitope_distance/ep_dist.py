#!/usr/bin/env python

# Author: Hongjun Bai @ VGS (Rolland Lab), MHRP, HJF (hbai@hivresearch.org)

import sys
import collections
import argparse
import functools
import itertools as it
import copy

import pandas as pd
import numpy as np

first = lambda iterable_obj: next(iter(iterable_obj))

def cleaned(instream, comment='#'):
    """ Remove comments(starts with #) or empty lines """
    for line in instream:
        cleaned_line = line.strip().split(comment)[0]
        if cleaned_line:
            yield cleaned_line.strip()

def get_aa_dist():
    _aa_dist = """
    # Converted based on BLOSUM62 matrix by 1/2*(sim(AA1, AA1)+sim(AA2, AA2)) - sim(AA1, AA2)
    # O: PNGS; -: insertion/deletion; maximum distane (-13) is assigned
    #
    	R	K	Q	E	D	N	H	T	S	A	G	P	C	V	I	L	M	F	Y	W	O	-
    R	0.0	3.0	4.0	5.0	7.5	5.5	6.5	6.0	5.5	5.5	7.5	8.0	10.0	7.5	7.5	6.5	6.0	8.5	8.0	11.0	13.0	13.0
    K	3.0	0.0	4.0	4.0	6.5	5.5	7.5	6.0	4.5	5.5	7.5	7.0	10.0	6.5	7.5	6.5	6.0	8.5	8.0	11.0	13.0	13.0
    Q	4.0	4.0	0.0	3.0	5.5	5.5	6.5	6.0	4.5	5.5	7.5	7.0	10.0	6.5	7.5	6.5	5.0	8.5	7.0	10.0	13.0	13.0
    E	5.0	4.0	3.0	0.0	3.5	5.5	6.5	6.0	4.5	5.5	7.5	7.0	11.0	6.5	7.5	7.5	7.0	8.5	8.0	11.0	13.0	13.0
    D	7.5	6.5	5.5	3.5	0.0	5.0	8.0	6.5	5.0	7.0	7.0	7.5	10.5	8.0	8.0	9.0	8.5	9.0	9.5	12.5	13.0	13.0
    N	5.5	5.5	5.5	5.5	5.0	0.0	6.0	5.5	4.0	7.0	6.0	8.5	10.5	8.0	8.0	8.0	7.5	9.0	8.5	12.5	13.0	13.0
    H	6.5	7.5	6.5	6.5	8.0	6.0	0.0	8.5	7.0	8.0	9.0	9.5	11.5	9.0	9.0	9.0	8.5	8.0	5.5	11.5	13.0	13.0
    T	6.0	6.0	6.0	6.0	6.5	5.5	8.5	0.0	3.5	4.5	7.5	7.0	8.0	4.5	5.5	5.5	6.0	7.5	8.0	10.0	13.0	13.0
    S	5.5	4.5	4.5	4.5	5.0	4.0	7.0	3.5	0.0	3.0	5.0	6.5	7.5	6.0	6.0	6.0	5.5	7.0	7.5	10.5	13.0	13.0
    A	5.5	5.5	5.5	5.5	7.0	7.0	8.0	4.5	3.0	0.0	5.0	6.5	6.5	4.0	5.0	5.0	5.5	7.0	7.5	10.5	13.0	13.0
    G	7.5	7.5	7.5	7.5	7.0	6.0	9.0	7.5	5.0	5.0	0.0	8.5	10.5	8.0	9.0	9.0	8.5	9.0	9.5	10.5	13.0	13.0
    P	8.0	7.0	7.0	7.0	7.5	8.5	9.5	7.0	6.5	6.5	8.5	0.0	11.0	7.5	8.5	8.5	8.0	10.5	10.0	13.0	13.0	13.0
    C	10.0	10.0	10.0	11.0	10.5	10.5	11.5	8.0	7.5	6.5	10.5	11.0	0.0	7.5	7.5	7.5	8.0	9.5	10.0	12.0	13.0	13.0
    V	7.5	6.5	6.5	6.5	8.0	8.0	9.0	4.5	6.0	4.0	8.0	7.5	7.5	0.0	1.0	3.0	3.5	6.0	6.5	10.5	13.0	13.0
    I	7.5	7.5	7.5	7.5	8.0	8.0	9.0	5.5	6.0	5.0	9.0	8.5	7.5	1.0	0.0	2.0	3.5	5.0	6.5	10.5	13.0	13.0
    L	6.5	6.5	6.5	7.5	9.0	8.0	9.0	5.5	6.0	5.0	9.0	8.5	7.5	3.0	2.0	0.0	2.5	5.0	6.5	9.5	13.0	13.0
    M	6.0	6.0	5.0	7.0	8.5	7.5	8.5	6.0	5.5	5.5	8.5	8.0	8.0	3.5	3.5	2.5	0.0	5.5	7.0	9.0	13.0	13.0
    F	8.5	8.5	8.5	8.5	9.0	9.0	8.0	7.5	7.0	7.0	9.0	10.5	9.5	6.0	5.0	5.0	5.5	0.0	3.5	7.5	13.0	13.0
    Y	8.0	8.0	7.0	8.0	9.5	8.5	5.5	8.0	7.5	7.5	9.5	10.0	10.0	6.5	6.5	6.5	7.0	3.5	0.0	7.0	13.0	13.0
    W	11.0	11.0	10.0	11.0	12.5	12.5	11.5	10.0	10.5	10.5	10.5	13.0	12.0	10.5	10.5	9.5	9.0	7.5	7.0	0.0	13.0	13.0
    O	13.0	13.0	13.0	13.0	13.0	13.0	13.0	13.0	13.0	13.0	13.0	13.0	13.0	13.0	13.0	13.0	13.0	13.0	13.0	13.0	0.0	13.0
    -	13.0	13.0	13.0	13.0	13.0	13.0	13.0	13.0	13.0	13.0	13.0	13.0	13.0	13.0	13.0	13.0	13.0	13.0	13.0	13.0	13.0	0.0
    """
    #
    content = []
    for line in cleaned(_aa_dist.split('\n')):
        fields = line.split('\t')
        content.append(fields)
    aa_list = ''.join(content[0])
    result = {}
    for rec  in content[1:]:
        ai = rec[0] 
        for j, aj in enumerate(aa_list):
           result[(ai, aj)] = float(rec[1+j])
    return result
DistBLOSUM62 = get_aa_dist()

#
def load_wts(wt_file, corresponding_msa):
    """ Load the precomputed interaction weight between ab and ag based on resolved complex structures """
    if type(corresponding_msa[0]) is not str: # Likely called in reticulate, take care the column-major character matrix
        msa_t = [''.join(s) for s in corresponding_msa]
        msa = transpose(msa_t)
    else:
        msa = corresponding_msa
    #
    seq_raw2aln = dict((s.replace('-', ''), s) for s in msa)

    df_wts = pd.read_table(wt_file, index_col=0, sep='\t', comment='#')
    seqs_t = list(df_wts['AA'])
    aln_raw = transpose(seqs_t)
    aln = [seq_raw2aln[s.replace('-', '')] for s in aln_raw]
    iw2imsa = idmap_msa2msa(aln_raw, aln)
    wts = dict((iw2imsa[iw], w) for iw, w in enumerate(df_wts['weight']) if w > 0.0)
    wt_seqs = aln
    return (wts, wt_seqs)

def load_msa(msa_file):
    """ For interface with reticulate only, not needed for python usage """
    with open(msa_file) as infile:
        heads, msa = read_fasta(infile)

    # Feed to the column major character matrix of R
    result = np.array([[aa for aa in seq] for seq in msa])
    n_rows, n_cols = result.shape
    result.shape = (n_cols, n_rows)
    return result

def ep_dist(seq, refseqs, wts, aa_dist=DistBLOSUM62):
    """ Get the minimum epitope distance between seq and refseqs """
    # Take care interface with R
    if type(seq) is not str:
        input_seq = ''.join(seq)
    else:
        input_seq = seq
    if type(refseqs[0]) is not str:
        refseqs_t = [''.join(s) for s in refseqs]
        refs = transpose(refseqs_t)
    else:
        refs = refseqs
    if first(wts.keys()) is not int:
        input_wts = dict((int(k), v) for k, v in wts.items())
    else:
        input_wts = wts
    # Done interface

    dist2refs = []
    for refseq in refs:
        dist2refs.append(ep_dist_pngs_ins(input_seq, refseq, input_wts, aa_dist))
    #print(dist2refs)
    return min(dist2refs)

#
def ep_dist_pngs_ins(seq, refseq, weights, distmat):
    insertions = get_insertions(seq, refseq)
    weights_ins = assign_w_ins(insertions, weights, seq)
    seq_o = mark_pngs(seq)
    refseq_o = mark_pngs(refseq)
    weighted_score = w_score(seq_o, refseq_o, weights_ins, distmat)
    #print(seq_o)
    #print(refseq_o)
    #print(''.join(refseq_o[i] for i in weights))
    #print(''.join(seq_o[i] for i in weights))
    return weighted_score / sum(weights_ins.values()) # return normalized score

def w_score(seq, refseq, weights, matrix):
    """  weighted sum of score of sites has a weight assigned """
    score = sum(matrix[(seq[i], refseq[i])]*w for i, w in weights.items())
    return score

def get_insertions(seq, ref, spacers='-. '):
    insertions = []
    gap_open = False
    for i, (aa, ref_aa) in enumerate(zip(seq, ref)):
        if aa not in spacers and ref_aa in spacers:
            gap_open = True
        elif ref_aa not in spacers:
            gap_open = False
        else:
            pass
        if gap_open:
            insertions.append(i)
    grouped = group_consecutive_int(insertions)
    return [(gp[0], gp[-1]) for gp in grouped]

def group_consecutive_int(int_list):
    if int_list:
        int_list.sort()
        grouped = [[int_list[0]]]
        for i, j in zip(int_list[:-1], int_list[1:]):
            if j == i+1:
                grouped[-1].append(j)
            elif j > i+1:
                grouped.append([j])
            else:
                assert(False)
        return grouped
    else:
        return []

def assign_w_ins(insertions, weights, seq):
    weight_with_ins = copy.deepcopy(weights)
    # get inserstions
    spacers = '-. '
    for s, e in insertions:
        w = 0.0
        if s-1 not in weights and e+1 in weights:   # seq: ===iiXXXX===
            w = weights[e+1]                        # ref: ===..XXXX===
        elif s-1 in weights and e+1 in weights:     # seq: ===XXiiXX===
            w = (weights[s-1] + weights[e+1])/2     # ref: ===XX..XX===
        elif s-1 in weights and e+1 not in weights: # seq: ===XXXXii===
            w = weights[s-1]                        # ref: ===XXXX..===
        else:
            pass
        if w == 0.0: continue
        for i in range(s, e+1):
            if seq[i] not in spacers:
                weight_with_ins[i] = w
    return weight_with_ins

#
def mark_pngs(seq, merge_overlapped_sequon=False, spacer='-', marker='O'):
    pngs = set(get_pngs(seq, merge_overlapped_sequon, spacer))
    return ''.join(marker if i in pngs else aa for i, aa in enumerate(seq))

def get_pngs(seq, merge_overlapped_sequon=True, spacer='-', restrict_P3=False):
    """ resturn [i_pngs0, i_pngs1, ...]; i_pngs0, index of PNGS0 in seq """
    idmap = idmap_raw2aln(seq, spacer)
    raw_seq = seq.upper()
    for c in spacer:
        raw_seq = raw_seq.replace(c, '')
    raw_seq += 'XXX'  # Append 'XXX' to handel sequence ends with 'N'
    pngs = set()
    overlapped = set()
    for i, aa in enumerate(raw_seq):
        if (aa.upper() == 'N'
            and raw_seq[i+1] != 'P'
            and raw_seq[i+2] in 'ST'):
            if restrict_P3 and raw_seq[i+3] == 'P':  #Ref: PMID 20510933
                continue
            pngs.add(i)
            if merge_overlapped_sequon and i-1 in pngs:
                overlapped.add(i-1)
    pngs -= overlapped
    pngs = [idmap[i] for i in pngs]
    pngs.sort()
    return pngs

def idmap_raw2aln(formated_seq, spacer='-'):
    accumulator = 0
    idmap = {}
    for i, aa in enumerate(formated_seq):
        if aa not in spacer:
            idmap[accumulator] = i
            accumulator += 1
    return idmap

def read_fasta(filestream):
    headers = []
    sequences = []
    for h, s in fasta_iter(filestream):
        headers.append(h)
        sequences.append(s)
    return (headers, sequences)

def fasta_iter(filestream):
    faiter = (x[1] for x in it.groupby(filestream, lambda line: line[0] == ">"))
    for header in faiter:
        header_str = header.__next__()[1:].strip()
        seq = "".join(s.strip() for s in faiter.__next__())
        yield header_str, seq

def idmap_msa2msa(msa_a, msa_b, spacer='-'):
    """ msa_a and msa_b must be identical after removing empty/gappy columns """
    tmsa_a = transpose(msa_a)
    tmsa_b = transpose(msa_b)
    i_a, i_b = 0, 0
    idmap = {}
    while i_a < len(tmsa_a) and i_b < len(tmsa_b):
        if tmsa_a[i_a] == tmsa_b[i_b]: # identical columns
            idmap[i_a] = i_b
            i_a += 1
            i_b += 1
        else:
            # if the two columns are different, one must be gappy (thus removed in the another) 
            n_spacers_a = tmsa_a[i_a].count(spacer)
            n_spacers_b = tmsa_b[i_b].count(spacer)
            if n_spacers_a > n_spacers_b:
                i_a += 1
            else:
                i_b += 1
    return idmap

def transpose(msa):
    try:
        msa_t = []  # Transposed msa
        for j in range(len(msa[0])):
            msa_t.append(''.join(x[j] for x in msa))
    except KeyError:
        #print('Check whether the MSA is empty or not same lengh.', file=sys.stderr)
        assert(False)
    return msa_t


def parse():
    parser = argparse.ArgumentParser(description="Get epitope distance between the references and the specified sequences")
    required = parser.add_argument_group('required arguments')
    required.add_argument('-msa',
            help="The input Multiple Sequence Alignment. Sequences should be included: (1) sequence of the antigen in the weight file, (2) reference sequence (default: the seqeunce in the weight file), (3) sequences to be analyzed",
            type=argparse.FileType('rt'),
            required=True)
    required.add_argument('-w', '--weight_file',
            help="The weight file depict the ab-ag interactions",
            type=argparse.FileType('rt'),
            required=True)
    required.add_argument('-p', '--prefix',
            help="Prefix(es) of the header of sequences to be analyzied (multiple opitions can be speicified, e.g. -p 'B C 01_AE') ",
            required=True, nargs='+')
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument('-r', '--ref', help="Prefix of the header/name of reference sequence (default: None, using the sequence from the weight file)")
    optional.add_argument('-o', '--outfile', help='Output file (default: stdout)', default=sys.stdout, type=argparse.FileType('wt'))
    return parser.parse_args()

def main():
    para = parse()

    # Load sequences
    heads, msa = read_fasta(para.msa)
    sequences = collections.OrderedDict((h, s) for h, s in zip(heads, msa))

    # Read the weight of epitope sites
    wts, wt_seqs = load_wts(para.weight_file, msa)

    # Reference sequences
    if para.ref is not None:
        refseq_names = [name for name in heads if name.startswith(para.ref)]
        print(para.ref)
        assert(len(refseq_names) > 0), f"no sequence in the MSA match the given prefix for references"
        refseqs = [sequences[n] for n in refseq_names]
        refnames = ', '.join(refseq_names)
    else:
        refseqs = wt_seqs
        refname =  'Sequence(s) from the weight file'

    # Epitope distance
    selected_seqs = [name for name in heads if (any(name.startswith(p) for p in para.prefix))]
    assert( len(selected_seqs) > 0 ), "no sequence in the MSA match the given prefix(es)"

    ep_distances = collections.OrderedDict()
    for name in selected_seqs:
        seq = sequences[name]
        ep_distances[name] = ep_dist(seq, refseqs, wts, DistBLOSUM62)
    #print(ep_distances)

    # Ouput
    fprint = functools.partial(print, file=para.outfile)
    fprint('#$ {}'.format(' '.join(sys.argv)))
    fprint(f'# Reference(s): {refnames}')
    fprint("# The final ep_distance is the minimum ep_distance to reference sequences")
    fprint('')
    fprint('sequence\tep_dist')
    for seq_name, epdist in ep_distances.items():
        fprint(f'{seq_name}\t{epdist:.3f}')

if __name__ == '__main__':
    main()

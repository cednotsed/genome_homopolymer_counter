import numpy as np
import pandas as pd
from Bio import SeqIO
# from pathlib import Path
import argparse
parser = argparse.ArgumentParser(description='Script for counting the number of homopolymers'
                                             'in a fasta file containing a single reference. '
                                             'Fields: nucleotide, start (incl.), end (incl.), homopolymer length')

parser.add_argument('-i', type=str, help='Provide filepath of input fasta file')
parser.add_argument('-o', type=str, help='Provide filepath of output tsv file')
args = parser.parse_args()

# cwd = Path.cwd()
record = SeqIO.read(args.i, 'fasta')
# record = SeqIO.read(cwd / 'krakblast/wuhan-hu-1.fasta', 'fasta')
fasta = pd.DataFrame(np.array(record.seq))
fasta['start_pos'] = fasta.index + 1    # convert to 1 based indexing
fasta['end_pos'] = -2
fasta['homopolymer_len'] = np.zeros(fasta.shape[0], dtype='int')
assert len(fasta) == 29903
fasta = fasta.rename({0 : 'nucleotide'}, axis=1)
groups = fasta.groupby('nucleotide')

for nucleotide in groups.groups.keys():
    nuc_index = groups.groups[nucleotide]
    for idx in range(len(nuc_index)):
        count = 1
        dummy_idx = idx
        while dummy_idx < len(nuc_index) - 1 and nuc_index[dummy_idx] + 1 == nuc_index[dummy_idx + 1]:
            count += 1
            dummy_idx += 1

        fasta.loc[nuc_index[idx], 'homopolymer_len'] = count

for i in range(fasta.shape[0]):
    start = fasta.iloc[i, :]['start_pos']
    length = fasta.iloc[i, :]['homopolymer_len']

    if start == 1:
        fasta.iloc[i, 2] = start

    else:
        fasta.iloc[i, 2] = start + length - 1

fasta = fasta.drop_duplicates(subset='end_pos')

fasta.to_csv(args.o, sep='\t', index=False, header=True)

#!/usr/bin/env python
import scipy as sp
import pandas as pd
import glob
import os, sys, time

summarized_seq_file = sys.argv[1] #'output/results/all_summarized_seqs'
unique_seq_file = sys.argv[2] #'output/results/all_unique_seqs'
output_file = sys.argv[3] #'output/results/alignment.txt'

# Load summarized seqs df
df = pd.read_csv(summarized_seq_file,delim_whitespace=True,skiprows=1)

# Extract seqs from unique_seq_file
seq_df = pd.read_csv(unique_seq_file,delim_whitespace=True)
df['seq'] = seq_df['seq']

# Replace lack of breakpoints with breakpoints at 88bp
indices = (df['fbp']==-1) & (df['rbp']==-1)
df.loc[indices,['fbp','rbp']] = 88

# Split sequencs into left, right, and indel
left_seqs = []
right_seqs = []
center_seqs = []
for row in df.iterrows():
    seq = row[1]['seq']
    fbp = row[1]['fbp']
    rbp = row[1]['rbp']
    left_seqs.append(seq[:fbp])
    right_seqs.append(seq[rbp:])
    cseq = row[1]['indel_seq'][1:-1]
    if len(cseq) > 5:
        cseq = '[%dbp]'%len(cseq)
    center_seqs.append(cseq)

# Create alignments as sequences
lrc_seqs = zip(left_seqs, center_seqs, right_seqs)
df['align'] = ['{:<30}{:^10}{:>30}'.format(l[60:],c,r[:-146]) for l,c,r in lrc_seqs]

# Create alignment data frame
indices = df['dL'] != 0
align_cols = [col for col in df.columns if 'ct_' in col] + ['align']
align_df = df.loc[indices,align_cols]

# Sum counts for each identical alginment
new_df = align_df.groupby('align').sum()
new_df.sort('ct_all',ascending=False,inplace=True)
new_df.reset_index(inplace=True)

# Save aligment
g = open(output_file,'w')
ct_cols = [col for col in new_df.columns if 'ct_' in col]
g.write(' '.join(['{:>6}'.format(col) for col in ct_cols]) + '  alignment\n')
ct_cols_tot = df.loc[~indices,ct_cols].sum()
g.write(' '.join(['%6d'%ct for ct in ct_cols_tot]) + '  < no length change >\n')
for n in range(20):
    g.write(' '.join(['%6d'%ct for ct in new_df.loc[n,ct_cols]]) + '  ' + new_df.loc[n,'align']+'\n')
g.close()


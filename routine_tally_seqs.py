#!/usr/bin/env python
import scipy as sp
import glob
import os, sys, time

# Get arguments
#output_dir = sys.argv[1]
#output_file = '%s/all_unique_seqs'%output_dir
output_file = sys.argv[1]
output_dir = '/'.join(output_file.split('/')[:-1])

# Get list of output directories
sample_nums = range(1,9)
sample_dirs = ['%s/sample_%d'%(output_dir,n) for n in sample_nums]

# Load unique seqs into memory
seq_to_counts_dict = {}
for n in sample_nums:
    sample_dir = '%s/sample_%d'%(output_dir,n)
    unique_seqs_file = '%s/unique_seqs'%(sample_dir)
    print 'Loading seqs from %s...'%unique_seqs_file
    f = open(unique_seqs_file)
    for line in f.readlines():
        atoms = line.split()
        seq = atoms[2]
        count = int(atoms[1])
        seq_to_counts_dict[seq] = [0]*(len(sample_nums)+1)
        #if not seq in seq_to_counts_dict.keys():
        #    seq_to_counts_dict[seq] = [0]*(len(sample_nums)+1)
        #seq_to_counts_dict[seq][n] += count
        #seq_to_counts_dict[seq][0] += count
    f.close()

# Count occurances of unique seqs
for n in sample_nums:
    sample_dir = '%s/sample_%d'%(output_dir,n)
    unique_seqs_file = '%s/unique_seqs'%(sample_dir)
    print 'Counting seqs in %s...'%unique_seqs_file
    f = open(unique_seqs_file)
    for line in f.readlines():
        atoms = line.split()
        seq = atoms[2]
        count = int(atoms[1])
        seq_to_counts_dict[seq][n] += count
        seq_to_counts_dict[seq][0] += count
    f.close()

# Get list of sequences ordered by their timepoints
total_counts = [x[0] for x in seq_to_counts_dict.values()]
seqs = seq_to_counts_dict.keys()
indices = sp.argsort(total_counts)[::-1]

# Write results
f = open(output_file,'w')
line = 'ct_all\t%s\tseq\n'%('\t'.join(['ct_%d'%n for n in sample_nums]))
f.write(line)
for i in indices:
    seq = seqs[i]
    counts = seq_to_counts_dict[seq]
    line = '\t'.join([str(x) for x in counts]) + '\t' + seq + '\n'
    f.write(line)
f.close()



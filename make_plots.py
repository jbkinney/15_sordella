#!/usr/bin/env python
import pandas as pd
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import sys

# Specify input and output directories and files
summarized_data_file = 'output/results/all_summarized_seqs.txt' #sys.argv[1]
samples_file = 'data/samples.txt'
out_dir = 'output/results' #sys.argv[2]

#summarized_data_file = sys.argv[1]
#samples_file = sys.argv[2]
#out_dir = sys.argv[3]

# Enter interactive mode; prevents plots from stealing focus
plt.ion()
plt.close('all')

# Set analysis parameters
max_l = 10
min_count = 10

# Load sample information
f = open(samples_file)
f.readline()
sample_nums = []
sample_names = []
for line in f.readlines():
    atoms = line.split()
    sample_nums.append(int(atoms[0]))
    sample_names.append(atoms[1])
f.close()
sample_nums = np.array(sample_nums)
num_samples = len(sample_nums)
sample_labels = ['sample %d'%n for n in sample_nums]

# Load wt sequence
f = open(summarized_data_file)
wt_seq = f.readline().split()[-1]
f.close()
seq_length = len(wt_seq)

# Load data into data frame
df = pd.read_csv(summarized_data_file,skiprows=1,delim_whitespace=True)

# Filter sequences based on counts
df_filtered = df[df['ct_all'] > min_count]

# Keep only sequences that have an insertion and/or deletion
df_mut = df_filtered[df_filtered['dL'] != 0]
#df_mut = df_filtered[(df_filtered['L_ins'] > 0) | (df_filtered['L_del'] > 0)]
num_mut_rows = df_mut.shape[0]

# Compute the number of sequences for each sample
filtered_seq_counts = sp.zeros(num_samples)
for s in sample_nums:
    # Get counts, L_ins, L_del, fbp, and rbp
    col_name = 'ct_%d'%(s)
    filtered_seq_counts[s-1] = 1.0*sum(df_filtered[col_name].iloc[:])

# Compute mutation statistics
L_ins_hist = np.zeros([max_l+1,num_samples])
L_del_hist = np.zeros([max_l+1,num_samples])
indel_locs = np.zeros([len(wt_seq),num_samples])
mut_seq_counts = sp.zeros(num_samples)
for n in range(num_mut_rows):
    # Tally insertion
    L_ins = min(df_mut['L_ins'].iloc[n],max_l)

    # Tally deletion
    L_del = min(df_mut['L_del'].iloc[n],max_l)

    # Tally indel location
    fbp = df_mut['fbp'].iloc[n]
    rbp = df_mut['rbp'].iloc[n]    

    for s in sample_nums:
        # Get counts, L_inss, L_del, fbp, and rbp
        col_name = 'ct_%d'%(s)
        count = 1.0*df_mut[col_name].iloc[n]
        mut_seq_counts[s-1] += count
        L_ins_hist[L_ins,s-1] += count
        L_del_hist[L_del,s-1] += count
        indel_locs[fbp:(rbp+1),s-1] += count


#
# Insertion/deletion multipanel plot
#

fig = plt.figure(figsize=(6.5,9))
gs = gridspec.GridSpec(num_samples, 4)
ymax = 1.0
yticks = [0,.5,1]
ylabels = ['%d%%'%(100*y) for y in yticks]
for s in sample_nums-1:
    ax = plt.subplot(gs[s,0:2])

    # Draw location of CRISPR cut site
    crispr_x = 88
    plt.plot([crispr_x, crispr_x],[0,1],color='m',linestyle=':')

    y = indel_locs[:,s]/filtered_seq_counts[s]
    x = np.arange(seq_length)
    plt.bar(x,y,width=1,linewidth=1,facecolor='blue',edgecolor='blue')
    plt.xlim([0,seq_length])
    plt.ylim([0,ymax])
    plt.yticks(yticks,ylabels)
    plt.ylabel(sample_labels[s],rotation=90,horizontalalignment='center')
    xticks = np.arange(0,seq_length,50)

    if s < num_samples-1:
        plt.xticks(xticks,[])
    else:
        plt.xlabel('breakpoint span (bp)')

    xticks = np.linspace(0,max_l,3,endpoint=True)

    ax = plt.subplot(gs[s,2])
    y = L_del_hist[:,s]/filtered_seq_counts[s]
    x = np.arange(max_l+1)
    plt.bar(x,y,width=0.8,linewidth=0.5,color='red')
    plt.yticks(yticks,[])
    plt.ylim([0,ymax])
    plt.xlim([0,max_l])
    if s < num_samples-1:
        plt.xticks(xticks,[])
    else:
        plt.xticks(xticks)
        plt.xlabel('deletion (bp)')

    ax = plt.subplot(gs[s,3])
    y = L_ins_hist[:,s]/filtered_seq_counts[s]
    x = np.arange(max_l+1)
    plt.bar(x,y,width=0.8,linewidth=0.5,color='green')
    plt.yticks(yticks,[])
    plt.ylim([0,ymax])
    plt.xlim([0,max_l])
    if s < num_samples-1:
        plt.xticks(xticks,[])
    else:
        plt.xticks(xticks)
        plt.xlabel('insertion (bp)')

plt.subplots_adjust(left=0.13, right=0.98, top=0.98, bottom=0.07)
plt.show()
fig.savefig(out_dir+'/mutations.pdf')

#
# Mutation efficiency bar plot
#
fig = plt.figure(figsize=(6.5,6.5))
ax = fig.add_subplot(111)
yticks = np.linspace(0,1,6,endpoint=True)
ylabels = ['%d%%'%(100*y) for y in yticks]

ps = mut_seq_counts/filtered_seq_counts
dps = sp.sqrt(ps*(1-ps)/filtered_seq_counts)
error_config = {'ecolor': '0.3'}
rects1 = ax.bar(sample_nums, ps, 0.8,
                color='blue', linewidth=1,
                yerr=dps,
                error_kw=error_config)
plt.xticks(sample_nums+0.5,sample_labels, rotation=45, horizontalalignment='right')
plt.ylim([0,1])
plt.yticks(yticks,ylabels)
plt.ylabel('Indel fraction')

plt.subplots_adjust(left=0.12, right=0.98, top=0.98, bottom=0.12)
plt.show()
fig.savefig(out_dir+'/rates.pdf')

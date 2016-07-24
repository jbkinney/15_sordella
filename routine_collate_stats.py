#!/usr/bin/env python
import scipy as sp
import pandas as pd
import glob
import os, sys, time

# Get arguments
#stats_file_glob = 'output/stats/*.stats'
stats_dir = sys.argv[1]
stats_file_glob = stats_dir+'/*.stats'

#output_file = 'output/results/stats.txt'
output_dir = sys.argv[2]
output_file = '%s/stats.txt'%output_dir

# Initiate dataframe

# Load each stat dataframe
files = glob.glob(stats_file_glob)
df_list = []
for file_name in files:
    file_base = file_name.split('/')[-1]
    sample_name = file_base.split('.')[0]
    this_df = pd.read_csv(file_name, delim_whitespace=True)
    this_df['sample_name'] = sample_name
    this_df['file'] = file_base
    del this_df['r1_file']
    del this_df['r2_file']
    df_list.append(this_df)

# Clean up dataframe 
df = pd.concat(df_list)
df.reset_index(inplace=True, drop=True)
df['good_reads'] = df['tot_reads'] * df['success_pct'] / 100.
df['good_reads'] = df['good_reads'].apply(int)
del df['success_pct']

# Group by sample name and compute pct_good
final_df = df.groupby('sample_name').sum()
final_df['pct_good'] = 100.*final_df['good_reads']/final_df['tot_reads']

# Write output dataframe
final_df.to_csv(output_file, float_format='%0.2f', sep='\t')

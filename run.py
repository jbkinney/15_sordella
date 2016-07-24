#!/usr/bin/env python
pipeline_dir = '.' 
python_to_use = '/usr/bin/env python'

# Indicate whether to use multiple nodes
use_multiple_nodes = True

# Define input files and directories
data_dir = pipeline_dir + '/data'
read1_file_glob = data_dir + '/LID296024-0*_S*_L001_R1_001.fastq'
read2_file_glob = data_dir + '/LID296024-0*_S*_L001_R2_001.fastq'
regions_file = data_dir + '/regions.txt'
samples_file = data_dir + '/samples.txt'

# Define output directories
output_dir = pipeline_dir + '/output'
reads_dir = output_dir + '/reads'
scripts_dir = output_dir + '/scripts'
results_dir = output_dir + '/results'

split_file_globs = \
    [(reads_dir + '/LID296024-0%d_S%d_L001_R1_001.fastq.*'%(n,n), \
      reads_dir + '/LID296024-0%d_S%d_L001_R2_001.fastq.*'%(n,n), \
      results_dir + '/sample_%d/observed_seqs.*'%n) \
    for n in range(1,9)]

# Run the pipeline
execfile(pipeline_dir+'/pipeline.py')


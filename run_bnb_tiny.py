#!/usr/bin/env python
pipeline_dir = '/sonas-hs/kinney/hpc/home/jkinney/projects/15_sordella'
python_to_use = '/usr/bin/env python'

# Indicate whether to use multiple nodes
use_multiple_nodes = True

# Define input files and directories
data_dir = pipeline_dir + '/data_tiny'
read1_file_glob = data_dir + '/sample_*_R1.fastq'
read2_file_glob = data_dir + '/sample_*_R2.fastq'
regions_file = data_dir + '/regions.txt'

# Define output directories
output_dir = pipeline_dir + '/output'
reads_dir = output_dir + '/reads'
scripts_dir = output_dir + '/scripts'
results_dir = output_dir + '/results'

split_file_globs = \
    [(reads_dir + '/sample_%d_R1.fastq.*'%n, \
      reads_dir + '/sample_%d_R2.fastq.*'%n, \
      results_dir + '/sample_%d/observed_seqs.*'%n) \
    for n in range(1,9)]

# Run the pipeline
execfile(pipeline_dir+'/pipeline.py')

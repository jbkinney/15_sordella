#!/usr/local/bin/python

#pipeline_dir = '/data/kinney/jkinney/github/14_crispr'
#python_to_use = '/opt/hpc/pkg/python-2.7/bin/python'
#input_dir = '/data/kinney/jkinney/github/14_crispr/example'

pipeline_dir = '/Users/jkinney/github/15_sordella'
python_to_use = '/usr/bin/env python'
input_dir = '/Users/jkinney/github/15_sordella/production'

output_dir = input_dir + '/output'
reads_dir = output_dir + '/reads'
scripts_dir = output_dir + '/scripts'
results_dir = output_dir + '/results'
summary_dir = output_dir + '/summaries'
plots_dir = output_dir + '/plots'
efficiency_file = output_dir + '/efficiencies.txt'

# Indicate whether to use multiple nodes
use_multiple_nodes = False

# Give names of fastq files as a glob. 
read1_file_glob = input_dir + '/data/LID296024-0*_S*_L001_R1_001.fastq'
read2_file_glob = input_dir + '/data/LID296024-0*_S*_L001_R2_001.fastq'
read1_split_file_glob = reads_dir + '/data/LID296024-0*_S*_L001_R1_001.fastq.*'
read2_split_file_glob = reads_dir + '/data/LID296024-0*_S*_L001_R2_001.fastq.*'

#experiments_file = input_dir + '/experiments.txt'
#barcodes_file = input_dir + '/barcodes.txt'
regions_file = input_dir + '/regions.txt'
#infection_file = input_dir + '/infection_rates.txt'

# Run the pipeline
execfile(pipeline_dir+'/pipeline.py')

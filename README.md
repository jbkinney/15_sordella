# 15_sordella

Code and processed data reported in Senturk et al.(2016), Nat. Commun. [in press]

DATA 

The raw sequence data for this study is available on the Sequence Read Archive under accession number [SRP078612](http://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP078612). The first 10K reads for each of the 8 samples are provided in the FASTQ files in the ```data_tiny/``` directory. The Cas9-targeted region of p53, as well as the primers use to amplify this locus for sequencing, are given in [```data/regions.txt```](data/regions.txt). The content of each of the 8 sequencing samples is described in [```data/samples.txt```](data/samples.txt). 

PUBLISHED RESULTS

The published results are provided in the directory ```published/```.
* [```all_unique_seqs.txt```](published/all_unique_seqs.txt) shows all of the unique reconstructed sequences for the p53 locus described in this paper. 
* [```all_summarized_seqs.txt```](published/all_summarized_seqs.txt) shows summary information for each unique sequence.
* [```stats.txt```](published/stats.txt) shows the number and percentage of sequences in each sample that could be successfully reconstructed. 
* [```alignment.txt```](published/alignment.txt) shows an alignment of the most prevalent sequences, along with their observed counts in the 8 samples. 
* [```rates.pdf```](published/rates.pdf)  is the bar chart shown in Fig. 2D
* [```mutations.pdf```](published/mutations.pdf)  is the bar chart shown in Fig. 2E

RUNNING THE PIPELINE ON SAMPLE DATA

1. Download this repository and change to the top-level directory ```15_sordella/```
2. Exectue ```$ ./run_tiny.py```, which will run the pipeline on the small sample datasets in ```data_tiny/```. The results will be stored in ```output/results/```
3. Exectue ```$ ./make_plots.py```, which will create the correspondign plots and store them in ```output/results/```
4. The results of this pipeline will be in ```output/results/```

RUNNING THE PIPELINE ON THE FULL DATASETS. 

1. Download the eight paired-end read datasets from [SRP078612](http://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP078612). Save the resulting 16 files (8 for read1, 8 for read2) in fastq format in the directory ```data/```.
2. Edit ```run.py``` so that the variables ```read1_file_glob```, ```read1_file_glob```, and ```split_file_globs``` point to these fastq files. 
3. Set the variable ```use_multiple_nodes = False``` in ```run.py``` to run the analysis in single node mode. To run analysis on multiple nodes, set ```use_multiple_nodes = True```. To get the analysis on multiple nodes working in your cluster environment, however, you might have to edit the function ```submit_and_complete_jobs()```, defined in ```pipleine.py```
4. Exectue ```$ ./run.py```, which will run the pipeline on the small sample datasets in ```data_tiny/```. The results will be stored in ```output/results/```
5. Exectue ```$ ./make_plots.py```, which will create the correspondign plots and store them in ```output/results/```
6. The results of this pipeline will be in ```output/results/```

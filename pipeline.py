import os, sys, time, commands
import glob
from ruffus import *

# Set names for various files
observed_seq_files_split = results_dir+'/*/observed_seqs.*'
observed_seq_files = results_dir+'/*/observed_seqs'
unique_seq_files = results_dir+'/*/unique_seqs'
unique_seq_file = results_dir+'/all_unique_seqs.txt'
summarized_seq_file = results_dir + '/all_summarized_seqs.txt'
alignment_file = results_dir+'/alignment.txt'
stats_dir = output_dir + '/stats'
plots_dir = output_dir + '/plots'

summarize_seqs_script = python_to_use  + ' ' + pipeline_dir + '/routine_summarize_seqs.py'
parse_seqs_script = python_to_use  + ' ' + pipeline_dir + '/routine_parse_seqs.py'
tally_seqs_script = python_to_use  + ' ' + pipeline_dir + '/routine_tally_seqs.py'
make_plots_script = python_to_use  + ' ' + pipeline_dir + '/routine_make_plots.py'
collate_stats_script = python_to_use  + ' ' + pipeline_dir + '/routine_collate_stats.py'
make_alignments_script = python_to_use  + ' ' + pipeline_dir + '/routine_make_alignments.py'

# Useful way to give feedback
def give_feedback(feedback):
    func_name =sys._getframe(1).f_code.co_name
    print '\nIn '+func_name+': '+feedback,
    sys.stdout.flush()

# Submits a list of scripts to be run as separate jobs
# Waits for all jobs to complete before continuing
def submit_and_complete_jobs(scripts, use_multiple_nodes):
    
    # Clear scripts dir
    os.system('rm -rf %s'%scripts_dir)
    os.system('mkdir %s'%scripts_dir)

    if use_multiple_nodes:

        # Write and execute scripts
        for n, script in enumerate(scripts):
            script_name = scripts_dir + '/script_num_%d.sh'%n 
            f = open(script_name, 'w')
            f.write(script)
            f.close()
            os.system('chmod +x '+script_name)
            os.system('qsub -cwd -e %s -o %s %s > .junk' % (script_name+'.e', script_name+'.o', script_name))
        
        # Monitor jobs 
        give_feedback('Waiting for jobs to complete...')
        jobs_remaining = len(scripts)
        wait_time = 60
        start_time = time.time()
        while jobs_remaining > 0:
            give_feedback(str(jobs_remaining)+' jobs remaining; waiting '+str(wait_time)+' seconds...')
            time.sleep(wait_time)
            jobs_remaining = int(commands.getoutput('qstat | grep script_num | wc -l'))
            
        finish_time = time.time()
        give_feedback('All jobs finished after %.1f seconds'%(finish_time - start_time))
            
    else:
        for n, script in enumerate(scripts):
            print 'Executing script %d of %d...'%(n+1,len(scripts))
            print script
            os.system(script)

    # Announce job completion
    give_feedback('Done.\n')

###################################################################################

# Split reads into files containing 50K reads each
#@files([read1_file_glob, read2_file_glob],[read1_split_file_glob, read2_split_file_glob])
def stage_one(ins=None, outs=None):

    # Clean output directories
    os.system('rm -rf %s'%output_dir)
    os.system('mkdir %s'%output_dir)
    os.system('rm -rf %s'%reads_dir)
    os.system('mkdir %s'%reads_dir)
    
    # Get list of r1 and r2 files
    read1_files = glob.glob(read1_file_glob)
    read2_files = glob.glob(read2_file_glob)
    
    # Make sure there are the same number of such files
    assert(len(read1_files) == len(read2_files))
    num_paired_read_files = len(read1_files)
    assert num_paired_read_files > 0
    
    scripts = []
    for pair_num in range(num_paired_read_files):
        read1_file_name = read1_files[pair_num]
        read2_file_name = read2_files[pair_num]
        
        read1_split_file_prefix = reads_dir + '/' + read1_file_name.split('/')[-1].split('.')[0]+'.fastq'
        read2_split_file_prefix = reads_dir + '/' + read2_file_name.split('/')[-1].split('.')[0]+'.fastq'
        
        #give_feedback('Splitting fastq files...\n')   
        script = '''
        source ~/.bash_profile_node  
        split -l 200000 %s %s.
        split -l 200000 %s %s.
        ''' % (read1_file_name, read1_split_file_prefix, read2_file_name, read2_split_file_prefix) 
        scripts.append(script)
    
    # Submit jobs    
    submit_and_complete_jobs(scripts, use_multiple_nodes)

# Create directories for experiments 
#@follows(stage_one)
#@files([], '%s/*/*'%results_dir)
def stage_two(ins=None, outs=None):

    # Clean results directory
    os.system('rm -rf %s'%stats_dir)
    os.system('mkdir %s'%stats_dir)
    os.system('rm -rf %s'%results_dir)
    os.system('mkdir %s'%results_dir)
    for n in range(1,9):
        os.system('mkdir %s/sample_%d'%(results_dir,n))

    give_feedback('Done!\n')

# Process reads in each fastq file into observed_seqs
@follows(stage_two)
#@files([read1_split_file_glob, read2_split_file_glob],observed_seq_files_split)
def stage_three(ins=None, outs=None):
    give_feedback('Reconstructing observed_seqs from reads...\n')
    
    scripts = []
    for read1_split_file_glob, read2_split_file_glob, out_file_glob \
        in split_file_globs:

        # Get list of extensions
        files = glob.glob(read1_split_file_glob)
        extensions = ['.'.join(x.split('.')[-1:]) for x in files]

        # For each extension, farm out a parse_seqs.py job
        for extension in extensions:
            r1_file = '%s.%s'%(read1_split_file_glob[:-2],extension)
            r2_file = '%s.%s'%(read2_split_file_glob[:-2],extension)
            out_file = '%s.%s'%(out_file_glob[:-2],extension)
            sample_name = out_file.split('/')[-2]
            stats_file = '%s/%s.%s.stats'%(stats_dir,sample_name,extension)
            command = '%s %s %s %s %s %s'%(parse_seqs_script, r1_file, r2_file, regions_file, out_file, stats_file)
            script = '''
            source ~/.bash_profile_node
            %s
            '''%command
            scripts.append(script)
    submit_and_complete_jobs(scripts, use_multiple_nodes)

# Combine observed seqs at each timepoint into a single file, then count the number
# of occurances of each one
#@follows(stage_three)
#@files(observed_seq_files_split, [observed_seq_files, unique_seq_files])
def stage_four(ins=None, outs=None):
    give_feedback('Computing results/[samples]/[observed_seqs, unique_seqs]...\n')
    sample_dirs = glob.glob('%s/*'%results_dir)
    scripts = []
    for sample_dir in sample_dirs:
        script = '''
        source ~/.bash_profile_node
        cat %s/observed_seqs.* > %s/observed_seqs
        cat %s/observed_seqs | sort | uniq -c | sort -nr | nl -b t > %s/unique_seqs
        '''%(sample_dir, sample_dir, sample_dir, sample_dir)
        scripts.append(script)
    submit_and_complete_jobs(scripts, use_multiple_nodes)

# Tally and summarize seqs
#@follows(stage_four)
#@files(unique_seq_files, summarized_seq_files)
def stage_five(ins=None, outs=None):
    give_feedback('Summarizing sequences observed in each experiment...\n')
    #sample_dirs = glob.glob('%s/*'%results_dir)
    scripts = []
    #for sample_dir in sample_dirs:
    script = '''
    source ~/.bash_profile_node
    %s %s
    %s %s %s %s
    '''%(tally_seqs_script, unique_seq_file,\
        summarize_seqs_script, unique_seq_file, regions_file, summarized_seq_file)
    scripts.append(script)    
    submit_and_complete_jobs(scripts, use_multiple_nodes)

# Collate stats
def stage_six(ins=None, outs=None):
    give_feedback('Collating read mapping statistics...\n')
    script = '''
    source ~/.bash_profile_node
    %s %s %s
    '''%(collate_stats_script, stats_dir, results_dir)
    print script
    os.system(script)

# Collate stats
def stage_seven(ins=None, outs=None):
    give_feedback('Make alignments...\n')
    script = '''
    source ~/.bash_profile_node
    %s %s %s %s
    '''%(make_alignments_script, summarized_seq_file, unique_seq_file, alignment_file)
    print script
    os.system(script)

### Run pipeline
print '##################################################'
print 'Begin read mapping pipeline'
pipeline_start_time = time.time()
stage_one()
stage_two()
stage_three()
stage_four()
stage_five()
stage_six()
stage_seven()
#pipeline_run([stage_five])
interval = time.time() - pipeline_start_time
print '\nPipeline Done! Runtime: %.1f min' % (interval/60.0)
print '##################################################'


    
    

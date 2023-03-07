#%%
import os
import glob
import pandas as pd
import git

#%%
# Date for sequencing run of the library to map
DATE = 20220423
# Description to be attached to folder names
DESCRIPTION = '_mapping/'

# Find project parental folder
repo = git.Repo('./', search_parent_directories=True)
homedir = repo.working_dir

# Path to input folder
INPUT_DIR = f'{homedir}/data/sequencing/{DATE}{DESCRIPTION}'

# Path to output folder
OUTPUT_DIR = f'{homedir}/data/processed_sequencing/{DATE}{DESCRIPTION}'

# Generate output directory if it doesn't exist
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

#%%
# Read list of demultiplexed sequences
files = glob.glob(INPUT_DIR+"*fastq.gz")

# Group sequences by index to feed them into the fastp


for file in files:
    name = file.split('/')[-1]
    name = name.split('_')[0] + "_" + name.split('_')[2]
    out_name = OUTPUT_DIR + name + ".fastq.gz"
    # Define outputs
    html_report = f'{OUTPUT_DIR}{DATE}_{name}_fastp_report.html'
    json_report = f'{OUTPUT_DIR}{DATE}_{name}_fastp_report.json'
    report_title = f'{DATE}{DESCRIPTION} fastp report'
    
    # Define string to be ran on the terminal
    fastp = f'''fastp \
        --in1 {file} \
        --out1 {out_name} \
        --trim_tail1 11 \
        --verbose \
        --disable_length_filtering \
        --html {html_report} \
        --json {json_report} \
        --report_title "{html_report}" \
        --thread 6
    '''

    # Run program
    os.system(fastp)

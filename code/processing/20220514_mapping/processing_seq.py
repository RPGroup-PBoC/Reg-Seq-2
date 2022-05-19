#%%
import os
import glob
import pandas as pd
import git

#%%
# Date for sequencing run of the library to map
DATE = 20220514
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

#%%\

# Group sequences by index to feed them into the fastp


for ind in [100, 110, 201, 204]:
    in1 = glob.glob(INPUT_DIR + str(ind) + "*_R1*.fastq.gz")[0]
    in2 = glob.glob(INPUT_DIR + str(ind) + "*_R2*.fastq.gz")[0]
    out1 = OUTPUT_DIR + str(ind) + "_R1.fastq.gz"
    out2 = OUTPUT_DIR + str(ind) + "_R2.fastq.gz"
    out_name = OUTPUT_DIR + str(ind) + ".fastq.gz"
    # Define outputs
    html_report = f'{OUTPUT_DIR}{DATE}_{ind}_fastp_report.html'
    json_report = f'{OUTPUT_DIR}{DATE}_{ind}_fastp_report.json'
    report_title = f'{DATE}{DESCRIPTION} fastp report'
    
    # Define string to be ran on the terminal
    fastp = f'''fastp \
        --in1 {in1} \
        --in2 {in2} \
        --merge \
        --merged_out {out_name} \
        --out1 {out1} \
        --out2 {out2} \
        --trim_tail1 11 \
        --trim_tail2 11 \
        --verbose \
        --disable_length_filtering \
        --html {html_report} \
        --json {json_report} \
        --report_title "{html_report}" \
        --thread 6
    '''

    # Run program
    os.system(fastp)

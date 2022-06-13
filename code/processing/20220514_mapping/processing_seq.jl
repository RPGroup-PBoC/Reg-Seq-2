using Glob
# Set path
dir = @__DIR__
home_dir = joinpath(split(dir, "/")[1:end-3])

# Date for sequencing run of the library to map
DATE = 20220514
# Description to be attached to folder names
DESCRIPTION = "_mapping/"

# Path to input folder
INPUT_DIR = "/$(home_dir)/data/sequencing/$(DATE)$(DESCRIPTION)"

# Path to output folder
OUTPUT_DIR = "/$(home_dir)/data/processed_sequencing_/$(DATE)$(DESCRIPTION)"

# Generate output directory if it doesn't exist
if ~isdir(OUTPUT_DIR)
    mkdir(OUTPUT_DIR)
end


# for loop is adapted to data structure for this run
for (ind, s) in zip([100, 110, 201, 204], ["S1", "S2", "S3", "S4"])
    in1 = INPUT_DIR * string(ind) * "_" * s * "_R1_001.fastq.gz"
    in2 = INPUT_DIR * string(ind) * "_" * s * "_R2_001.fastq.gz"
    out1 = OUTPUT_DIR * string(ind) * "_R1.fastq.gz"
    out2 = OUTPUT_DIR * string(ind) * "_R2.fastq.gz"
    out_name = OUTPUT_DIR * string(ind) * ".fastq.gz"
    # Define outputs
    html_report = "$(OUTPUT_DIR)$(DATE)_$(ind)_fastp_report.html"
    json_report = "$(OUTPUT_DIR)$(DATE)_$(ind)_fastp_report.json"
    report_title = "$(DATE)$(DESCRIPTION) fastp report"
    
    # Define string to be ran on the terminal
    fastp = `fastp 
        --in1 $(in1) 
        --in2 $(in2)
        --out1 $(out1) 
        --out2 $(out2) 
        --trim_tail1 11 
        --trim_tail2 11 
        --verbose 
        --disable_length_filtering 
        --html $html_report 
        --json $json_report 
        --report_title \"$(html_report)\" 
        --thread 6
        --qualified_quality_phred 20
    `

    # Run program
    run(`conda activate wgregseq`)
    run(fastp)
    run(`conda deactivate wgregseq`)
end

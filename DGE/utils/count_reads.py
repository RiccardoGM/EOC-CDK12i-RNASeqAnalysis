# --- import libraries --- #
import os
import subprocess
import sys


# --- Retrieve inputs from command line --- #
NAME = sys.argv[1]
idx = int(sys.argv[2])
stranded = str(sys.argv[3])

# --- Set the path to GTF annotation file --- #
gencode_file = "gencode.v44.annotation.gtf.gz" # Set gencode_file name
gtf_file = f"/path/to/your/{gencode_file}" # Set path to gencode_file


# --- Create list of BAM files found in the directory --- #
bam_dir = "/path/to/your/BAM_data_folder" # Set path to your bam files directory
bam_files = [f for f in os.listdir(bam_dir) if (f.endswith(".bam") & f.startswith(NAME))]
bam_file = bam_files[idx-1]


# --- Perform HTSeq counting --- #
counts_data_dir = "/path/to/your/counts_data_dir/"  # Set output directory
# construct output filename
output_filename = (
    os.path.splitext(bam_file)[0]
    .split('.fastq')[0]
    + f"_counts_s#{stranded}.txt"
)
# full output path
output_file = os.path.join(counts_data_dir, output_filename)
# Run HTSeq count
cmd = (
    f"htseq-count -f bam -r pos -s {stranded} "
    f"-i gene_id -t exon "
    f"{os.path.join(bam_dir, bam_file)} {gtf_file} "
    f"> {output_file}"
)
print("File:", output_file)
print("Counting started")
subprocess.call(cmd, shell=True)
print("Counting ended")


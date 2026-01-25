# Utils description

## count_reads.py

This script performs gene-level read counting from RNAseq BAM files using HTSeq-count and a GENCODE annotation file. The script selects a BAM file based on a sample name and index provided via the command line, then runs `htseq-count` to assign aligned reads to genes using exon features

### User-defined inputs

The following parameters must be customized before running the script:

- **GENCODE annotation file**
  - `gencode_file_name`
  - `gencode_file_path`

- **BAM files directory**
  - `bam_dir`

- **Output directory for count tables**
  - `counts_data_dir`

### Command-line arguments

The script expects three arguments:

- `NAME`: prefix used to identify BAM files
- `idx`: index of the BAM file to process (1-based)
- `stranded`: library strandedness (`yes`, `no`, or `reverse`)

### Output

The script generates a gene-level count table in the specified output directory with the following name format:

<sample_name>_counts_s#<stranded>.txt

## submit_slurm_job.sh

This is a SLURM batch script that submits a job to perform gene-level read counting using the `count_reads.py` Python workflow

### Description

The script:

- Allocates compute resources on a cluster using SLURM directives (`#SBATCH`)  
- Activates a dedicated Conda environment (`htseq_env`) containing HTSeq and dependencies  
- Passes command-line arguments to `count_reads.py` (`NAME`, `idx`, `stranded`)  
- Runs the Python script on the selected BAM file and generates the gene-level count table  

### User-customizable parameters

- **Path to Conda initialization script** (`. "/path/to/conda.sh"`)  
- **Conda environment name** (`htseq_env`)  
- **SLURM resource requests** (CPUs, memory, walltime, partition)  

### Inputs

Refer to the `count_reads.py` description for details on:

- Sample name (`NAME`) and BAM file index (`idx`)  
- Library strandedness (`stranded`)  
- BAM files directory, GENCODE annotation file, and output directory

### Output

Produces the same gene-level count file as described in `count_reads.py`, written to the configured output directory

## samples_loop.sh

This script automates the submission of multiple SLURM jobs for HTSeq read counting by looping over a set of BAM files

### Description

The script:

- Defines the sample prefix (`NAME`), the number of BAM files to process (`max_index`), and the library strandedness (`stranded`)  
- Loops over all BAM files from 1 to `max_index` (should match the number of files starting with `NAME`)  
- Submits a SLURM job for each file by calling `submit_slurm_job.sh` with the appropriate arguments (`NAME`, index, `stranded`)  

### User-customizable parameters

- **Sample prefix** (`NAME`)  
- **Number of BAM files** (`max_index`)  
- **Library strandedness** (`stranded`)  

### Execution

Each iteration of the loop submits a separate job to the cluster, which in turn runs `count_reads.py` for that specific BAM file  

### Output

For each BAM file, a separate gene-level count table is generated in the configured output directory, as described in `count_reads.py`

## create_counts_matrix.py

This script aggregates individual gene-level HTSeq count files into a single counts matrix suitable for downstream analysis (e.g., DESeq2)

### Description

The script:

- Loads metadata for a given sample line to determine the mapping between sample names and count files  
- Reads the individual count files from a specified directory (`counts_data_dir`)  
- Filters for gene identifiers starting with `"ENSG0"`  
- Combines all count columns into a matrix with genes as rows and samples as columns  
- Exports the resulting counts matrix as a tab-delimited file

### User-customizable parameters

- **Counts files directory** (`counts_data_dir`)  
- **Metadata directory** containing `<line_name>_metadata.csv` (`metadata_dir`)  
- **Line name** (`line_name`) to select the correct metadata and count files  
- **Output directory** for the counts matrix (`counts_mat_dir`)  
- Optional batch identifier (`batch`) to include in the output file name

### Output

The script produces a counts matrix file:

<line_name>_counts_matrix.txt

or, if `batch` is specified:

<line_name>_<batch>_counts_matrix.txt
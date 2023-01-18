# SaLTy
Staphylococcus aureus Lineage Typer (SaLTy).
SaLTy assigns a lineage to S. aureus WGS data and is suitable for processing large volumes of WGS.


usage: PROG [-h] [-t THREADS] [-f] [--report] [-i GENOME_FOLDER | -r READS_FOLDER] [-o OUTPUT_FOLDER]
            [-csv] [-s] [-l LINEAGES] [-k KMA_INDEX]

options:
  -h, --help            show this help message and exit

GENERAL:
  -t THREADS, --threads THREADS
                        Number of threads (speeds up parsing raw reads).
  -f, --force           Overwite existing output folder.
  --report              Only generate summary report from previous SALTy outputs.

INPUT:
  -i GENOME_FOLDER, --genome_folder GENOME_FOLDER
                        Input folder with assembled DNA sequence file.
  -r READS_FOLDER, --reads_folder READS_FOLDER
                        Folder with forward and reverse raw reads (fastq.gz)

OUTPUT:
  -o OUTPUT_FOLDER, --output_folder OUTPUT_FOLDER
                        Output Folder to save result.
  -csv, --format_csv    Output file in csv format.
  -s, --summary         Concatenate all output assignments into single file.

DB PATHS:
  -l LINEAGES, --lineages LINEAGES
                        Path to specific alleles for each lineage.
  -k KMA_INDEX, --kma_index KMA_INDEX
                        Path to indexed KMA database.

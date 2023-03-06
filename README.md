# <i>Staphylococcus aureus</i> Lineage Typer (SaLTy)
SaLTy assigns a lineage to <i>Staphylococcus aureus</i> WGS data and is suitable for describing large-scale <i>S. aureus</i> genomic epidemiology.

SALTy typing is highly accurate and can quickly analyse large volumes of <i>S. aureus </i> illumina reads (fastq.gz) or genome assemblies (fasta).

---
# Dependencies
1. python (v3.6 or greater)
2. kma (v1.4.9 or greater)
3. pandas  (v1.5.0 or greater)

---

# Installation
```commandline
conda install -c conda-forge -c bioconda salty
```
*See below for Mac M1 installation.

[![Anaconda-Server Badge](https://anaconda.org/bioconda/salty/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/salty/badges/downloads.svg)](https://anaconda.org/bioconda/salty)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/salty/badges/version.svg)](https://anaconda.org/bioconda/salty)

# Usage
```commandline
salty.py -i <input_folder> -o <output_folder>

GENERAL:
  -t THREADS, --threads THREADS
                        Number of threads (speeds up parsing raw reads).
  -f, --force           Overwite existing output folder.
  --report              Only generate summary report from previous SALTy
                        outputs.
  --check               Checks that required dependencies are available.
  -v, --version         Returns current SaLTy version.


INPUT:
  -i INPUT_FOLDER, --input_folder INPUT_FOLDER
                        Folder of genomes (*.fasta or *.fna) and/or pair end
                        reads (each accession must have *_1.fastq.qz and
                        *_2.fastq.

OUTPUT:
  -o OUTPUT_FOLDER, --output_folder OUTPUT_FOLDER
                        Output Folder to save result.
  -c, --csv_format    Output file in csv format.
  -s, --summary         Concatenate all output assignments into single file.

DB PATHS:
  -l LINEAGES, --lineages LINEAGES
                        Path to specific alleles for each lineage.
  -k KMA_INDEX, --kma_index KMA_INDEX
                        Path to indexed KMA database.
```


# Example Usage
Run SaLTy on a folder containing pairs of fastq files and/or assemblies.

```commandline
salty --input_folder /input/reads/folder/ --output_folder /output/folder/name/ --threads INT
```
Specify SaLTy output assigned lineages in comma-separated (.csv) format. 
```commandline
salty --input_folder /input/reads/folder/ --output_folder /output/folder/name/ --csv_format
```

# Example Output
````
Genome		Lineage		SACOL0451	SACOL1908	SACOL2725
SRR992071338	15		24		20		24
SRR992044718	15		24		20		24
SRR10591328	15		24		20		69
````

## Column Descriptions
| Column | Description |
| ----------- | ----------- |
|Sample|Input strain ID (extracted from input files).|
|Lineage|The lineage assigned by the SaLTy algorithim.|
|SACOL0451|The identified allele for the SACOL0451 locus.|
|SACOL1908|The identified allele for the SACOL1908 locus.|
|SACOL2725|The identified allele for the SACOL2725 locus.|


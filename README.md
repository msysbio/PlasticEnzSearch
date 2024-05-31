# PlasticTools

## Overview

PlasticTools is a bioinformatics tool designed to quantify and report the concentration of potential plastic-degrading enzymes in metagenomic data.

## System Requirements

PlasticTools is compatible with Linux distributions and macOS.
Unfortunately, it does not support Windows at this time.

## Installation

This project depends on several external tools. With conda installed, you can install all the dependencies with the following commands:

```bash
# hmmsearch
conda install -c bioconda hmmer
# featureCounts
conda install -c bioconda subread
# prodigal
conda install -c bioconda prodigal
# samtools
conda install -c bioconda samtools
# biopython
conda install -c anaconda biopython
```

for parallel execution of prodigal you will additionaly need pprodigal:

```bash
pip install pprodigal
```

### BLAST

To start using BLAST you will need to download the BLAST+ executables are available from the NCBI website as well as a protein database, like swissprot which can be downloaded from the UniProt site.

### Optional UI Installation

For the optional UI, you will need to install `flet`. The installation process varies depending on your operating system.

#### Linux (debian-based)

```bash
sudo apt-get install libmpv1
pip install flet
```

#### macOS

```bash
brew install mpv
pip install flet
```

## Output

In the chosen output folder, an HTML file will be present with visualization of the results, a TSV file with the data used in the visualization, fasta files containing the sequences from the HMMER tool output and XML files with the BLAST output.
To normalize for depth and gene length, RPKM (reads per kilobase per million mapped reads) is used for calculating abundance of plastic degradation enzymes. 
$$ RPKM = \frac{reads\ mapped}{\left(\frac{gene\ length}{1000}\right) \times \left(\frac{total\ reads}{1,000,000}\right)} $$
The raw counts of reads mapped and the proportion (reads mapped / total reads) are also provided for each plastic type but may not be suitable when comparing the outputs of different samples.

## Usage

You can run PlasticTools with the following commands:

### With Translation

```bash
python3 main.py --contigs <path/to/contigs> --plastic <plastic_name> --output <path/to/output/folder>
```

This will run `translate_search.py` and `quantify_hmm.py`.

### Without Translation

```bash
python3 main.py --contigs <path/to/contigs> --plastic <plastic_name> --output <path/to/output/folder> --skip_translation
```

This will only run `quantify_hmm.py`.

### With Custom HMMs

```bash
python3 main.py --contigs <path/to/contigs> --plastic <plastic_name> --output <path/to/output/folder> --custom_hmm <path/to/hmms_file>
```

## Troubleshooting

If you encounter any issues while running PlasticTools, you can refer to the following case scenario for troubleshooting:

```bash
for id in . {
	python3 PlastEnz --output ./{file} --plastic PET --bam ./{file}.bam --ORF ./{file}.fa
}
```

This will create a new folder at the specified output path (e.g., `./Blauwwe`). 
Inside this folder, it will create a subfolder for the specified plastic (e.g., `./Blauwe/PET`) and a `temps` folder to save all temporary files. 
It will also generate an `Abundance_table.tsv` file.

The `Bams/Gene_counts` folder will contain `.bam` or `.tsv` files. The `ORF` argument should point to the exact ORFs file, and the `Contigs` argument should point to the exact contigs file.

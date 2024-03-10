# PlasticEnzSearch

## Overview
PlasticEnzSearch is a bioinformatics tool designed to quantify and report the concentration of potential plastic-degrading enzymes in metagenomic data.

## System Requirements

PlasticEnzSearch is compatible with Linux distributions and macOS.
Unfortunately, it does not support Windows at this time.

## Installation

This project depends on several external tools. You can install all the dependencies with the following commands:

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

for parallel execution of prodigal you will additionaly need:

```bash
pip install pprodigal
```

## SVM Installation
```bash
pip install torch
pip install fair-esm
```

## Optional UI Installation


For the optional UI, you will need to install `flet`. The installation process varies depending on your operating system.

### Linux:
```bash
sudo apt-get install libmpv1
pip install flet
```

### macOS:
```bash
brew install mpv
pip install flet
```


## Usage
You can run PlasticEnzSearch with the following commands:

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
If you encounter any issues while running PlasticEnzSearch, you can refer to the following case scenario for troubleshooting:

```bash
for id in . {
	python3 PlastEnz --output ./{file} --plastic PET --bam ./{file}.bam --ORF ./{file}.fa
}
```
This will create a new folder at the specified output path (e.g., `./Blauwwe`). 
Inside this folder, it will create a subfolder for the specified plastic (e.g., `./Blauwe/PET`) and a `temps` folder to save all temporary files. 
It will also generate an `Abundance_table.tsv` file.

The `Bams/Gene_counts` folder will contain `.bam` or `.tsv` files. The `ORF` argument should point to the exact ORFs file, and the `Contigs` argument should point to the exact contigs file. 

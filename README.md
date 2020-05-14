    #####################################################################
    ####  __   __(_) _ __  ___   ___   _ __ | |_  ___  _ __  |___ \  ####
    ####  \ \ / /| || '__|/ __| / _ \ | '__|| __|/ _ \| '__|   __) | ####
    ####   \ V / | || |   \__ \| (_) || |   | |_|  __/| |     / __/  ####
    ####    \_/  |_||_|   |___/ \___/ |_|    \__|\___||_|    |_____| ####
    #####################################################################


# VirSorter 2 

VirSorter 2 is has made major updates to VirSorter 1:

- work with more viral groups including dsDNAphage, ssDNA, RNA, NCLDV, lavidaviridae
- apply machine learning to estimate viralness using genomic and taxonomic features and hallmark gene counts
- train with high quality virus genomes from metagenomes or other sources


# Installation

Conda is the easiest way to install dependencie. Conda can install following [this link](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).  
You can create an new env by:

```bash
conda create -n vs2 python=3 scikit-learn=0.22.1 imbalanced-learn pandas seaborn hmmer prodigal screed last ncbi-genome-download ruamel.yaml snakemake=5.16.0 click
conda activate vs2
```

Then install VirSorter 2:

```bash
git clone https://github.com/jiarong/VirSorter2.git
cd VirSorter2
pip install -e .
virsorter setup -d db -j 4
```

This will download all databases and install dependencies and takes 10+ mins, but this only need to be done once. Please be patience and wait..

# Quick run

To run the viral sequence classification, you can run `virsorter run -h` to see options. Here an example using 8 threads (-j) and including provirus step (--provirus):

```bash
cd test
virsorter run -w test-out -i 4seq.fa -j 8 --provirus
```

# More details later

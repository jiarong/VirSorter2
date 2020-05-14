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

VirSorte 2 only require python 3, snakemake, python click and ruamel package to start. If you have not them installed, conda is the easiest way (conda can install following https://docs.conda.io/projects/conda/en/latest/user-guide/install/). These can be installed in current env by:

```bash
conda install -c bioconda -c conda-forge snakemake click ruamel.yaml
```

Or else you can create an new env by:

```bash
conda create -n vs2 -c bioconda -c conda-forge snakemake click ruamel.yaml
conda activate vs2
```

Other dependencies and databases are installed by setup command:

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

    #####################################################################
    ####  __   __(_) _ __  ___   ___   _ __ | |_  ___  _ __  |___ \  ####
    ####  \ \ / /| || '__|/ __| / _ \ | '__|| __|/ _ \| '__|   __) | ####
    ####   \ V / | || |   \__ \| (_) || |   | |_|  __/| |     / __/  ####
    ####    \_/  |_||_|   |___/ \___/ |_|    \__|\___||_|    |_____| ####
    #####################################################################


# VirSorter 2 

VirSorter2: A multi-classifier, expert-guided approach to detect diverse DNA and RNA virus genomes.

VirSorter 2 is has made major updates to VirSorter 1:

- work with more viral groups including dsDNAphage, ssDNA, RNA, NCLDV, lavidaviridae
- apply machine learning to estimate viralness using genomic and taxonomic features and hallmark gene counts
- train with high quality virus genomes from metagenomes or other sources


# Installation (tested on CentOS linux; should work in all linux)

Conda is the easiest way to install VirSorter2. Conda can install by following [this link](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).

```bash
conda intall virsorter
```

To install the development version:

```bash
conda create -n vs2 python=3 scikit-learn=0.22.1 imbalanced-learn pandas seaborn hmmer prodigal screed last ncbi-genome-download ruamel.yaml snakemake=5.16.0 click
conda activate vs2
git clone https://github.com/jiarong/VirSorter2.git
cd VirSorter2
pip install -e .
```

Then download all databases and install dependencies (takes 10+ mins, but this only need to be done once). Please be patience and wait..

```bash
virsorter setup -d db -j 4
```

# Quick run

To run the viral sequence classification:

```bash
# fetch testing data
git clone https://github.com/jiarong/VirSorter2.git
cd VirSorter2/test
# run classification with 4 threads (-j)
virsorter run -w test-out -i 8seq.fa -j 4
```

Due to large HMM database, this small dataset will also test a few mins to finish. In the output directory (test-out), three files are useful:

- final-viral-combined.fa      # identified viral sequences
- final-viral-score.tsv        # score of each viral sequences
- final-viral-boundary.tsv     # start and end position in origial sequences

# More options  

VirSorter2 finds all viral groups currently included (ssDNAphage, NCLDV , RNA, ssDNA virus, and lavidavirida) by default. You can use `--include-groups` to chose specific groups:
```
virsorter run -w test-out -i 8seq.fa --include-groups "dsDNAphage,ssDNA" -j 4 --forceall
```

VirSorter2 run's main computational bottleneck is the annotation step (feature extraction). In case you want to re-run with different score cutoff (--min-score), `classify` argument can skip the annotation steps, and only re-run classify step.

```
virsorter run -w test-out -i 8seq.fa --include-groups "dsDNAphage,ssDNA" -j 4 --min-score 0.8 classify --forceall
```

You can run `virsorter run -h` to see all options.

# More details later

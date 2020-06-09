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

## option 1

Conda is the easiest way to install VirSorter2. Conda can install by following [this link](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).

```bash
conda intall virsorter
```

## option 2

To install the development version:

```bash
conda create -n vs2 python=3 scikit-learn=0.22.1 imbalanced-learn pandas seaborn hmmer prodigal screed last ncbi-genome-download ruamel.yaml snakemake=5.16.0 click
conda activate vs2
git clone https://github.com/jiarong/VirSorter2.git
cd VirSorter2
pip install -e .
```

# download database and dependencies

Then download all databases and install dependencies (takes 10+ mins, but this only need to be done once). Please be patient. The following downloads databases and dependencies to `db` directory, and its location is recorded in the tool configuration as a default for `--db-dir` of other `virsorter` subcommands.

```bash
virsorter setup -d db -j 4
```

# Quick run

To run viral sequence identification:

```bash
# fetch testing data
wget -O test.fa https://raw.githubusercontent.com/jiarong/VirSorter2/master/test/8seq.fa
# run classification with 4 threads (-j)
virsorter run -w test.out -i test.fa -j 4
```

Due to large HMM database that VirSorter2 uses, this small dataset takes a few mins to finish. In the output directory (test.out), three files are useful:

> - final-viral-combined.fa      # identified viral sequences
> - final-viral-score.tsv        # score of each viral sequences across groups
> - final-viral-boundary.tsv     # start and end position in origial sequences

# More options  

VirSorter2 finds all viral groups currently included (ssDNAphage, NCLDV , RNA, ssDNA virus, and lavidavirida) by default. You can use `--include-groups` to chose specific groups:
```
virsorter run -w test.out -i test.fa --include-groups "dsDNAphage,ssDNA" -j 4 --forceall
```

VirSorter2 takes one positional argument, `all` or `classify`. The default is `all`, which means running the whole pipeline, including 1) preprocessing 2) annotation (feature extraction), and classification. The main computational bottleneck is the annotation step, taking about 95% of CPU time. In case you just want to re-run with different score cutoff (--min-score), `classify` argument can skip the annotation steps, and only re-run classify step.

```
virsorter run -w test.out -i test.fa --include-groups "dsDNAphage,ssDNA" -j 4 --min-score 0.8 classify --forceall
```

In case you need to have some results quickly, there are two options: 1) turn off provirus step with `--provirus-off`; this reduces sensitivity on sequences that are only partially virus; 2) subsample ORFs from each sequence with `--max-orf-per-seq`; This option subsamples ORFs to a cutoff if a sequence has more ORFs than that. Note that this option is only availale when `--provirus-off` is used. 

```
virsorter run -w test.out -i test.fa --provirus-off --max-orf-per-seq 20 --forceall
```

You can run `virsorter run -h` to see all options. VirSorter2 is a wrapper around [snakemake](https://snakemake.readthedocs.io/en/stable/). It's an great pipeline management tool, designed for reproducibility, and running on computer clusters. All snakemake options still works here. You just need to append those snakemake option to virsorter options (after `all` or `classify`). The `--forceall` option in the above examples is a snakemake option. It forces to re-run the pipeline.

# Training customized classifiers

VirSorter2 currently has classifiers of five viral groups (dsDNAphage, NCLDV, RNA, ssNA virus, and *lavidaviridae*). It's designed for easy addition of more classifiers. The information of classifiers are store in the database (`-d`) specified during [setup step](#-download-database-and-dependencies). For each viral group, it needs four folowing files:

- model

  > random forest classifier model for identifying viral sequences

- customized.hmm (optional)

  > a collection of viral HMMs for gene annotation; if not specified, the one in `db/hmm/viral/combined.hmm` is used.

- hallmark-gene.list (optional)

  > names of hallmark gene hmm in the above viral hmm database file; These hallmark gene hmms can be collected by literature search or identified by comparing hallmark gene sequences (protein) against HMMs database above with `hmmsearch`; if not specified, no hallmark genes are counted in feature table

- rbs-prodigal-train.db (optional)

  > prodigal RBS (ribosomal binding site) motif training model; this can be produced with `-t` option in prodigal; This is useful feature for NCLDV due to large genome size for training; For other viral groups, it's OK to skip this file.

In this tutorial, I will show how to make `model` for *autolykiviridae*. 


```
# download sequences
wget https://github.com/jiarong/small-dataset/raw/master/vibrio_autolyki.fna.gz -O autolyki.fna.gz
# train feature file
virsorter train-feature --seqfile autolyki.fna.gz --hmm db/hmm/viral/combined.hmm --frags-per-genome 5 --jobs 4 -w autolyki-feature.out 
# check output
ls autolyki-feature.out
```

In the output directory (`autolyki-feature.out`), `all.pdg.ftr` is the feature file needed for next step, and `feature-importances.tsv` shows the importance of each feature used.  

To make the classifier model, we also need a feature file from cellular organisms. This can be done by collecting genomes from cellular organisms and repeat the above step. Note number of cellular genomes are very large (>200K). Here I will re-use the feature file I have prepared before. 

```
# fetch feature file for cellular organisms
wget https://zenodo.org/record/3823805/files/nonviral-common-random-fragments.ftr.gz?download=1 -O nonviral.ftr
# train the classifier model
virsorter train-model --viral-ftrfile autolyki-feature.out/all.pdg.ftr --nonviral-ftrfile nonviral.ftr --balanced --jobs 4 -w autolyki-model.out
```

In `autolyki-model.out`, `model` is the classifier model we need. Then put it in database directory

```
mkdir db/group/autolykiviridae
cp autolyki-model.out/model db/group/autolykiviridae
# reuse hallmark gene list from dsDNAphage due to their similarity
cp db/group/dsDNAphage/hallmark-gene.list db/group/autolykiviridae/
```

Now you can try this new classifier on the testing dataset:

```
virsorter run -w test.out -i test.fa --include-groups "autolykiviridae" -j 4 --min-score 0.8 all
```


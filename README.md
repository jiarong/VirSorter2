    #####################################################################
    ####  __   __(_) _ __  ___   ___   _ __ | |_  ___  _ __  |___ \  ####
    ####  \ \ / /| || '__|/ __| / _ \ | '__|| __|/ _ \| '__|   __) | ####
    ####   \ V / | || |   \__ \| (_) || |   | |_|  __/| |     / __/  ####
    ####    \_/  |_||_|   |___/ \___/ |_|    \__|\___||_|    |_____| ####
    #####################################################################


# VirSorter 2 

[![Version](https://anaconda.org/bioconda/virsorter/badges/version.svg)](https://anaconda.org/bioconda/virsorter)
[![Build Status](https://travis-ci.org/jiarong/VirSorter2.svg?branch=master)](https://travis-ci.org/jiarong/VirSorter2)
<!--
[![Bioconda](https://img.shields.io/conda/dn/bioconda/virsorter.svg?label=Bioconda )](https://anaconda.org/bioconda/virsorter)
-->


VirSorter2 applies a multi-classifier, expert-guided approach to detect diverse DNA and RNA virus genomes. It has made major updates to its [previous version](https://github.com/simroux/VirSorter):

- work with more viral groups including dsDNA phages, ssDNA viruses, RNA viruses, NCLDV (Nucleocytoviricota), *lavidaviridae* (virophages);
- apply machine learning to estimate viralness using genomic features including structural/functional/taxonomic annotation and viral hallmark genes;
- train with high quality virus genomes from metagenomes or other sources.


# Installation (tested on CentOS linux; should work in all linux; MacOS is not supported at the moment)

## Option 1 (bioconda: virsorter version 2.0.beta)

Conda is the easiest way to install VirSorter2. If you do not have conda installed, it can be installed following [this link](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).

```bash
conda install -c bioconda virsorter=2
```

## Option 2

The development version is most updated and recommended. To install the development version:

```bash
conda create -n vs2 -c bioconda -c conda-forge "python>=3.6" scikit-learn=0.22.1 imbalanced-learn pandas seaborn hmmer prodigal screed ruamel.yaml "snakemake>=5.16,<=5.26" click
conda activate vs2
git clone https://github.com/jiarong/VirSorter2.git
cd VirSorter2
pip install -e .
```

# Download database and dependencies

Before running VirSorter2, users must download all databases and install dependencies (takes 10+ mins, but this only need to be done once). The following command line downloads databases and dependencies to `db` directory, and its location is recorded in the tool configuration as a default, so you do not need to type `--db-dir` for other `virsorter` subcommands.

Note that this step takes > 10 mins. If you feel impatient and cancel the process, make sure to **remove the diretory specified by `-d/--db-dir` (`db` in this case) before running again**.

```bash
# just in case there is a failed attemp before; 
#   remove the whole diretory specified by -d
rm -rf db
# run setup
virsorter setup -d db -j 4
```

# Quick run

To run viral sequence identification on a test dataset:

```bash
# fetch testing data
wget -O test.fa https://raw.githubusercontent.com/jiarong/VirSorter2/master/test/8seq.fa
# run classification with 4 threads (-j) and test-out as output diretory (-w)
virsorter run -w test.out -i test.fa -j 4 all
ls test.out
```

Due to the large HMM database that VirSorter2 uses, this small dataset takes a few mins to finish. In the output directory (test.out), three files are useful:

- `final-viral-combined.fa`:  identified viral sequences
- `final-viral-score.tsv`:    table with score of each viral sequences across groups
- `final-viral-boundary.tsv`: table with boundary information

More details about each of these output files can be found [here](#detailed-description-on-output-files).

---
**NOTE**

Note that suffix `||full`, `||lt2gene` and `||{i}_partial` (`{i}` can be numbers starting from 0 to max number of viral fragments found in that contig) have been added to original sequence names to differentiate sub-sequences in case of multiple viral subsequences found in one contig.

---

# More options  

## choosing viral groups (`--include-groups`)

VirSorter2 finds all viral groups currently included (ssDNAphage, NCLDV , RNA, ssDNA virus, and *lavidaviridae*) by default. You can use `--include-groups` to select only specific groups. For those only interested in phage:

```bash
rm -rf test.out
virsorter run -w test.out -i test.fa --include-groups "dsDNAphage,ssDNA" -j 4 all
```

## re-run with different score cutoff (`--min-score` and `--classify`)

VirSorter2 takes one positional argument, `all` or `classify`. The default is `all`, which means running the whole pipeline, including 1) preprocessing, 2) annotation (feature extraction), and 3) classification. The main computational bottleneck is the annotation step, taking about 95% of CPU time. In case you just want to re-run with different score cutoff (--min-score), `classify` argument can skip the annotation steps, and only re-run only the classify step.

```bash
virsorter run -w test.out -i test.fa --include-groups "dsDNAphage,ssDNA" -j 4 --min-score 0.8 classify
```

The above overwrites the previous final output files. If you want to keep previous results, you can use `--label` to add a prefix to the new final output files.

```bash
virsorter run -w test.out -i test.fa --include-groups "dsDNAphage,ssDNA" -j 4 --min-score 0.9 --label rerun classify
```

## speed up a run (`--provirus-off`) 

In case you need to have some results quickly, there are two options: 1) turn off provirus step with `--provirus-off`; this reduces sensitivity on sequences that are only partially viral; 2) subsample ORFs from each sequence with `--max-orf-per-seq`; This option subsamples ORFs if a sequence has more ORFs than the number provided. Note that this option is only availale when `--provirus-off` is used. 

```bash
rm -rf test.out
virsorter run -w test.out -i test.fa --provirus-off --max-orf-per-seq 20 all
```

## Other options

You can run `virsorter run -h` to see all options. VirSorter2 is a wrapper around [snakemake](https://snakemake.readthedocs.io/en/stable/), a great pipeline management tool designed for reproducibility, and running on computer clusters. All snakemake options still work with VirSorter2, and users can simply append those snakemake option to virsorter options (after `all` or `classify`). For example, the `--forceall` snakemake option can be used to re-run the pipeline.

```bash
virsorter run -w test.out -i test.fa --provirus-off --max-orf-per-seq 20 all --forceall
```

When you re-run any VirSorter2 command, it will pick up at the step (rule in snakemake term) where it stopped last time. It will do nothing if it succeeded last time. The `--forceall` option can be used to enforce the re-run.

## DRAMv compatibility

[DRAMv](https://github.com/shafferm/DRAM) is a tool for annotating viral contigs identified by VirSorter. It needs two input files from VirSorter: 1) viral contigs, 2) `affi-contigs.tab` that have info on viral/nonviral and hallmark genes along contigs. In VirSorter2, these files can be generated by `--prep-for-dramv` flag.

```bash
rm -rf test.out
virsorter run --prep-for-dramv -w test.out -i test.fa -j 4 all
ls test.out/for-dramv
```

# Detailed description on output files

- final-viral-combined.fa

  > identified viral sequences, including two types:
  > - full sequences identified as viral (identified with suffix `||full`);
  > - partial sequences identified as viral (identified with suffix `||{i}_partial`); here `{i}` can be numbers starting from 0 to max number of viral fragments found in that contig;
  > - short (less than two genes) sequences with hallmark genes identified as viral (identified with suffix `||lt2gene`);

 
- final-viral-score.tsv

  > This table can be used for further screening of results. It includes the following columns:
  >   - sequence name
  >   - score of each viral sequences across groups (multiple columns)
  >   - max score across groups
  >   - max score group
  >   - contig length 
  >   - hallmark gene count
  >   - viral gene %
  >   - nonviral gene %

---
**NOTE**

Note that classifiers of different viral groups are not exclusive from each other, and may have overlap in their target viral sequence space, which means this information should not be used or considered as reliable taxonomic classification. We limit the purpose of VirSorter2 to viral idenfication only.

---

- final-viral-boundary.tsv

  > only some of the columns in this file might be useful:
  >   - seqname: original sequence name
  >   - trim\_orf\_index\_start, trim\_orf\_index\_end:  start and end ORF index on orignal sequence of identified viral sequence
  >   - trim\_bp\_start, trim\_bp\_end:  start and end position on orignal sequence of identified viral sequence
  >   - trim\_pr: score of final trimmed viral sequence
  >   - partial:  full sequence as viral or partial sequence as viral; this is defined when a full sequence has score > score cutoff, it is full (0), or else any viral sequence extracted within it is partial (1) 
  >   - pr\_full:  score of the original sequence
  >   - hallmark\_cnt:  hallmark gene count
  >   - group: the classifier of viral group that gives high score; this should **NOT** be used as reliable classification

---
**NOTE**

VirSorter2 tends to sometimes overestimate the size of viral sequence during provirus extraction procedure in order to achieve better sensitity. We recommend cleaning these provirus predictions to remove potential host genes on the edge of the predicted viral region, e.g. using a tool like CheckV (https://bitbucket.org/berkeleylab/checkv).

---

# Training customized classifiers

VirSorter2 currently has classifiers of five viral groups (dsDNAphage, NCLDV, RNA, ssNA virus, and *lavidaviridae*). It's designed for easy addition of more classifiers. The information of classifiers are store in the database (`-d`) specified during [setup step](#download-database-and-dependencies). For each viral group, it needs four files below:

- model

  > random forest classifier model for identifying viral sequences

- customized.hmm (optional)

  > a collection of viral HMMs for gene annotation; if not specified, the one in `db/hmm/viral/combined.hmm` is used.

- hallmark-gene.list (optional)

  > names of hallmark gene hmm in the above viral hmm database file; These hallmark gene hmms can be collected by literature search or identified by comparing hallmark gene sequences (protein) against HMMs database above with `hmmsearch`; if not specified, no hallmark genes are counted in feature table

- rbs-prodigal-train.db (optional)

  > prodigal RBS (ribosomal binding site) motif training model; this can be produced with `-t` option in prodigal; This is useful feature for NCLDV due to large genome size for training; For other viral groups, it's OK to skip this file.

In this tutorial, I will show how to make `model` for the *autolykiviridae* family.

First, prepare the dataset needed: 1) high quality viral genomes 2) protein sequence of hallmark gene; and install two more dependecies.

```bash
# download genome sequences
wget https://github.com/jiarong/small-dataset/raw/master/autolyki/vibrio_autolyki.fna.gz -O autolyki.fna.gz
# download hallmark gene seqs
wget https://raw.githubusercontent.com/jiarong/small-dataset/master/autolyki/DJR.fa -O DJR.fa
# download source code
git clone https://github.com/jiarong/VirSorter2.git
# install two more dependencies
conda install -c bioconda -y screed hmmer
```

Then identify hallmark gene HMMs by protein sequences of hallmark genes 

Note that we will need the VirSorter2 database here. If you skip the tutorial above, you can download the database by `virsorter setup -d db -j 4`. This will take 10+ mins.

```bash
# compare all HMMs and protein sequences of hallmark gene
# this will take 10+ mins due to large hmm database file
hmmsearch -T 50 --tblout DJR.hmmtbl --cpu 4 -o /dev/null db/hmm/viral/combined.hmm DJR.fa
# get HMMs names that are signicant hits with protein sequences of hallmark gene
python VirSorter2/virsorter/scripts/prepdb-train-get-seq-from-hmm-domtbl.py 50 DJR.hmmtbl > hallmark-gene.list
```

With `hallmark-gene.list` and the high quality genomes `autolyki.fna.gz`, you can train the features that are used for the classifier model.

```bash
# train feature file
virsorter train-feature --seqfile autolyki.fna.gz --hallmark hallmark-gene.list --hmm db/hmm/viral/combined.hmm --frags-per-genome 5 --jobs 4 -w autolyki-feature.out 
# check output
ls autolyki-feature.out
```

In the output directory (`autolyki-feature.out`), `all.pdg.ftr` is the feature file needed for next step.  

To make the classifier model, we also need a feature file from cellular organisms. This can be done by collecting genomes from cellular organisms and repeat the above step. Note number of cellular genomes are very large (>200K). Here I will re-use the feature file I have prepared before. 

```bash
# fetch feature file for cellular organisms
wget https://zenodo.org/record/3823805/files/nonviral-common-random-fragments.ftr.gz?download=1 -O nonviral.ftr.gz
gzip -d nonviral.ftr.gz
# train the classifier model
virsorter train-model --viral-ftrfile autolyki-feature.out/all.pdg.ftr --nonviral-ftrfile nonviral.ftr --balanced --jobs 4 -w autolyki-model.out
```

In `autolyki-model.out`, `feature-importances.tsv` shows the importance of each feature used. `model` is the classifier model we need. Then put the `model` and `hallmark-gene.list` in database directory as the existing viral groups. Note that **only letters** are allowed for group directory under `db/group/`.

```bash
### attention: only letters (both upper and lower case) are allowed in group names
mkdir db/group/autolykiviridae
cp autolyki-model.out/model db/group/autolykiviridae
cp hallmark-gene.list db/group/autolykiviridae/
```

Now you can try this new classifier on the testing dataset, and compare with `dsDNAphage` classifier:

```bash
# download the testing dataset
wget -O test.fa https://raw.githubusercontent.com/jiarong/VirSorter2/master/test/8seq.fa
# identify viral sequences in testing dataset; it takes 10+ mins;
virsorter run -w autolyki-model-test.out -i test.fa --include-groups "dsDNAphage,autolykiviridae" -j 4 --min-score 0.8 all
# check the scores in two classifiers
cat autolyki-model-test.out/final-viral-score.tsv
```


# Acknowledgement
VirSorter 2 is jointly developed by the Sullivan Lab at Ohio State University (https://u.osu.edu/viruslab/) and the Viral Genomics Group at the DOE Joint Genome Institute (https://jgi.doe.gov/our-science/scientists-jgi/viral-genomics/).  Funding was provided by NSF (#OCE1829831, #ABI1758974), the U.S. Department of Energy (#DE-SC0020173), and the Gordon and Betty Moore Foundation (#3790). The work conducted by the U.S. Department of Energy Joint Genome Institute is supported by the Office of Science of the U.S. Department of Energy under contract no. DE-AC02-05CH11231.

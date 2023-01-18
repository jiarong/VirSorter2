# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),

## [Unreleased]
- Make genbank file with pfam annotation: desc field has format "tax:viral||anno:XXXX" 
- Check "suspicious" words in each contig: methyltransferase, epimerase, endonuclease
- Mark short jobs as local_rules to optimaize for `virsorter run` --cluster mode
- Add "all" in --include-groups
- Change initial slide window to gene number based (4 genes) from bp size; apply the minimal prophage size at the end.

## [2.2.4] - 2022-01-18
- Fix issues with numpy v1.24 by limiting numpy version <1.24
- Change "singularity" to "apptainer" in `README`
- Add "PROVIRUS_MIN_HALLMARK_CNT" in `template-config-original.yaml` for tuning prophage search in `provirus.py`
- Add new cols: "final_max_score_group" and "final_max_score" in `final-viral-boundary.tsv` to be consistent with "max_score_group" and "max_score" in `final-viral-score.tsv`
- Remove records in `final-viral-boundary.tsv` but not in `final-viral-score.tsv`

## [2.2.3] - 2021-04-20
- Fix `--seqname-suffix-off` regex issues
- Fix issues with pandas v1.3.0
- Add email option (virsorter2 at gmail) for bug report

## [2.2.2] - 2021-04-20
- Fix bug with --seqname-suffix-off causing error with DRAMv

## [2.2.1] - 2021-03-31
- Fix minor bugs

## [2.2] - 2021-03-30
- Add --viral-gene-enrich-off option
- Add --seqname-suffix-off option
- Change default --include-groups from all groups to dsDNAphage and ssDNA
- Add FAQ, viral identification SOP link in README

## [2.1] - 2021-01-15
### Added
- Add --keep-original-seq option
- Add --high-confidence-only option
- Add --exclude-lt2gene option

### Changed
- Remove more itermediate files by default
- Change HMMER3.3.1 to HMMER3.3; HMMER3.3.1 throws seg fault with low complexity
sequences


## [2.0] - 2020-12-01
### Added
- Add --viral-gene-required
- Add --hallmark-required
- Add --hallmark-required-on-short
- Add --prep-for-dramv
- Add --label in virsorter run to allow adding prefix to output files
- Add shape to boundry file

### Changed
- Apply length prefilter to provirus.py; skip all seqs length < min(3000, MIN_GENOME_SIZE)
- Limit max file splits to 1000

## [2.0.beta] - 2020-06-19
### Added
- Make score table for trimmed viral seqs after provirus extraction
- Add length to score table
- Add subcommand `virsorter train-feature` and `virsorter train-model` for training customized classifiers of new viral groups, `virsorter config` for change configurations

## [2.0.alpha] - 2020-06-07
### Added
- First release 

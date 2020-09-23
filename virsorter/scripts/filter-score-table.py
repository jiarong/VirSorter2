#!/usr/bin/env python

import sys
import os
import screed
import numpy as np
import pandas as pd
import click

from ruamel.yaml import YAML


script_dir = os.path.dirname(os.path.abspath(__file__))
snakefile_dir = os.path.dirname(script_dir)
pkg_dir = os.path.dirname(snakefile_dir)
sys.path.append(pkg_dir)
from virsorter.config import get_default_config, set_logger

DEFAULT_CONFIG = get_default_config()
D = DEFAULT_CONFIG['GROUP_INFO']
DEFAULT_MIN_SIZE_ALLOWED_WO_HALLMARK_GENE = \
        DEFAULT_CONFIG['DEFAULT_MIN_SIZE_ALLOWED_WO_HALLMARK_GENE']

CONTEXT_SETTINGS = {'help_option_names': ['-h', '--help']}
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--hallmark-required', is_flag=True, default=False,
        help='require hallmark gene')
@click.option('--hallmark-required-on-short', is_flag=True, default=False,
        help='require hallmark gene on short seqs')
@click.option('--viral-gene-required', is_flag=True, default=False,
        help='require viral gene')
@click.argument('config', type=click.Path())
@click.argument('intable', type=click.Path())
@click.argument('inseqfile', type=click.Path())
@click.argument('outtable', type=click.Path()) 
@click.argument('outseqfile', type=click.Path()) 
def main(config, intable, inseqfile, outtable, outseqfile, hallmark_required,
        hallmark_required_on_short, viral_gene_required):
    '''Filter based on cols in score table

    Example:
        python filter-score-table.py <config.yaml> <score.tsv> <viral.fa> \
            <score-filtered.tsv> <viral-filtered.fa>

        <score.tsv>: score table with extra added info 
            (length, hallmark, viral, cellular)
        <viral.fa>: viral contig file with seqname matched to <score.tsv>
        <score-filtered.tsv>: filtered <score.tsv>
        <viral-filtered.fa>: filtered <viral.fa>

    '''
    config_f = config
    score_f = intable
    seqfile = inseqfile
    score_out_f = outtable
    seq_out_f = outseqfile

    config = YAML().load(open(config_f))
    df = pd.read_csv(score_f, sep='\t', header=0)
    sel = [True,] * len(df)
    if config['HALLMARK_REQUIRED']:
        sel = sel & (df['hallmark'] > 0)
    elif config['HALLMARK_REQUIRED_ON_SHORT']:
        length_cutoff_ser = df['max_score_group'].map(
                lambda x: D[x]['MIN_SIZE_ALLOWED_WO_HALLMARK_GENE'] \
                        if x in D \
                        else DEFAULT_MIN_SIZE_ALLOWED_WO_HALLMARK_GENE
        )
        sel_neg = (df['hallmark'] < 1) & (df['length'] < length_cutoff_ser)
        sel = sel & (~sel_neg)

    if config['VIRAL_GENE_REQUIRED']:
        sel = sel & (df['viral'] > 0)

    df = df.loc[sel,:]

    st = set(df['seqname'])
    # filter seqs 
    with screed.open(seqfile) as sp, open(seq_out_f, 'w') as fw:
        for rec in sp:
            header = rec.name
            #name, desc = header.split(None ,1)
            name = header.split(None ,1)[0]
            if not name in st:
                continue
            mes = f'>{name}\n{rec.sequence}\n'
            fw.write(mes)

    df.to_csv(score_out_f, sep='\t', na_rep='nan', 
            index=False, float_format='%.3f')

if __name__ == '__main__':
    main()

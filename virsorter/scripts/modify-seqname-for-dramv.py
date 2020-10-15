#!/usr/bin/env python

import sys
import os
import screed
import numpy as np
import pandas as pd
import click

script_dir = os.path.dirname(os.path.abspath(__file__))
snakefile_dir = os.path.dirname(script_dir)
pkg_dir = os.path.dirname(snakefile_dir)
sys.path.append(pkg_dir)
from virsorter.config import get_default_config, set_logger

DEFAULT_CONFIG = get_default_config()

CONTEXT_SETTINGS = {'help_option_names': ['-h', '--help']}
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--output', '-o', type=click.File('w'), default=sys.stdout,
        help='output seqfile with modified name for dramv')
@click.argument('inseqfile', type=click.Path())
@click.argument('intable', type=click.Path())
def main(intable, inseqfile, output):
    '''Filter based on cols in score table

    Example:
        python filter-score-table.py <viral.fa> <score.tsv> \
            -o <viral-name-modified-for-dramv.fa>

        <viral.fa>: viral contig file with seqname matched to <score.tsv>
        <score.tsv>: score table with extra added info 
            (length, hallmark, viral, cellular)
        <viral-name-modified-for-dramv.fa>: <viral.fa> with name modified 
            dramv (replace | with _ and add -cat[1-6] suffix

    '''
    score_f = intable
    seqfile = inseqfile
    fw = output

    df = pd.read_csv(score_f, sep='\t', header=0)
    df.index = df['seqname']

    with screed.open(seqfile) as sp:
        for rec in sp:
            header = rec.name
            #name, desc = header.split(None ,1)
            name = header.split(None ,1)[0]
            oriname, ty = name.rsplit('||', 1)
            provirus = False
            if ty.endswith('partial'):
                provirus = True

            has_hallmark = True if df.at[name, 'hallmark'] > 0 else False
            viral = df.at[name, 'viral']
            cellular = df.at[name, 'cellular']
            viral_enrich = True if viral > cellular else False

            if not provirus:
                if viral_enrich and has_hallmark:
                    contig_cat = 1
                elif viral_enrich or has_hallmark:
                    contig_cat = 2
                else:
                    contig_cat = 3
            else:
                if viral_enrich and has_hallmark:
                    contig_cat = 4
                elif viral_enrich or has_hallmark:
                    contig_cat = 5
                else:
                    contig_cat = 6

            name = name.replace('|', '_')
            mes = f'>{name}-cat_{contig_cat}\n{rec.sequence}\n'
            fw.write(mes)

if __name__ == '__main__':
    main()

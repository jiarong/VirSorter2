#!/usr/bin/env python

import sys
import os
import screed
import numpy as np
import pandas as pd
import click

@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('--suffix', help='suffix added to sequence name')
@click.argument('score-f', type=click.Path(exists=True))
@click.argument('seqfile', type=click.Path(exists=True))
@click.argument('outfile', type=click.Path())
def main(suffix, score_f, seqfile, outfile):
    '''Add sequence length to table

    \b
    <SCORE-F>: <all-fullseq-proba.tsv>, score table for all contigs >= 2 genes
    <SEQFILE>: <viral-fullseq.fa>, file with fullseq as viral seqname matched 
        to <all-fullseq-proba.tsv>
    <OUTFILE>: <viral-combined-proba.tsv>, score table for viral seqs with 
        suffix added to seqname
    '''
    st = set()
    with screed.open(seqfile) as sp:
        for rec in sp:
            header = rec.name
            name = header.split(None ,1)[0]
            st.add(name)

    df = pd.read_csv(score_f, sep='\t', header=0)
    df = df.loc[df['seqname'].isin(st), :]
    if suffix != None:
        ser = df['seqname'].map(lambda x: f'{x}{suffix}')
        df['seqname'] = ser

    df.to_csv(outfile, sep='\t', na_rep='NaN', 
            index=False, float_format='%.3f')

if __name__ == '__main__':
    main()

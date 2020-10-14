#!/usr/bin/env python

import sys
import os
import screed
import pandas as pd
import click
import numpy as np

script_dir = os.path.dirname(os.path.abspath(__file__))
snakefile_dir = os.path.dirname(script_dir)
pkg_dir = os.path.dirname(snakefile_dir)
sys.path.append(pkg_dir)
from virsorter.utils import TAXFILE_COLS, parse_hallmark_hmm

@click.command()
@click.option(
        '--hallmark',
        default=None,
        type=click.Path(),
        help='hallmark gene list file (hallmark-gene.list)'
)
@click.argument('taxfile', type=click.Path(exists=True))
@click.argument('outfile', type=click.Path())
def main(taxfile, outfile, hallmark):
    '''Add column showing if best hit hmm is a hallmark

    Example:
        python add-hallmark-to-taxfile.py --hallmark <hallmark>
            <taxfile> 
            <outfile>

        <TAXFILE>: all.pdg.hmm.tax which contain info on best hit hmm,
            its taxonomy association, bitscore
        <OUTFILE>: add column showing whether a best hit hmm is hallmark gene
    '''

    if hallmark != None:
        d_hallmark = parse_hallmark_hmm(hallmark)
    else:
        d_hallmark = {}

    # no header in all.pdg.hmm.tax
    df = pd.read_csv(taxfile, sep='\t', names=TAXFILE_COLS)
    cutoff_ser = df.loc[:,'hmm'].map(
            lambda x: d_hallmark.get(x, ['NA', np.inf])[1]
    )
    is_hallmark_ser = (df['score'] >= cutoff_ser)
    is_hallmark_ser = [1 if i else 0 for i in is_hallmark_ser]
    df['hallmark'] = is_hallmark_ser
    df.to_csv(outfile, sep='\t', na_rep='nan', index=False, header=False)

if __name__ == '__main__':
    main()

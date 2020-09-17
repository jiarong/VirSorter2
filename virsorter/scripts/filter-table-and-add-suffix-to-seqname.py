#!/usr/bin/env python

import sys
import os
import screed
import numpy as np
import pandas as pd


def main():
    '''Add sequence length to table

    Example:
        python filter-table-and-add-suffix-to-seqname.py \
            "||full" \
            <all-fullseq-proba.tsv> \
            <viral-fullseq.fa> \
            <viral-combined-proba.tsv>

        ||full: suffix to add to all items in seqname column
        <all-fullseq-proba.tsv>: score table for all contigs >= 2 genes
        <viral-fullseq.fa>: file with fullseq as viral seqname matched 
            to <all-fullseq-proba.tsv>
        <viral-combined-proba.tsv>: score table for viral seqs with 
            suffix added to seqname

    '''
    if len(sys.argv) != 5:
        mes = ('python {} "||full" <all-fullseq-proba.tsv> '
                '<viral-fullseq.fa> <viral-combined-proba.tsv>')
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    suffix = sys.argv[1]
    score_f = sys.argv[2]
    seqfile = sys.argv[3]
    outfile = sys.argv[4]

    st = set()
    with screed.open(seqfile) as sp:
        for rec in sp:
            header = rec.name
            name = header.split(None ,1)[0]
            st.add(name)

    df = pd.read_csv(score_f, sep='\t', header=0)
    df = df.loc[df['seqname'].isin(st), :]
    ser = df['seqname'].map(lambda x: f'{x}{suffix}')
    df['seqname'] = ser

    df.to_csv(outfile, sep='\t', na_rep='NaN', 
            index=False, float_format='%.3f')

if __name__ == '__main__':
    main()

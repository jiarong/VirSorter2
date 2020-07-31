#!/usr/bin/env python

import sys
import os
import screed
import pandas as pd


def main():
    '''Add sequence length to table

    Example:
        python add-length-to-table.py.py <infile.tsv> <viral-combined.fa> 
            <out.tsv>

        <infile.tsv>: merged classifcation results showing probability of being
            viral with different classifers
        <viral-combined.fa>: viral contig file with seqname matched 
            to <infile.tsv>

    '''
    if len(sys.argv) < 3:
        mes = '*** Usage: python {} <infile.tsv> <viral-combined.fa>\n'
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    infile = sys.argv[1]
    seqfile = sys.argv[2]
    outfile = sys.argv[3]
    d_len = {}
    with screed.open(seqfile) as sp:
        for rec in sp:
            l = rec.name.split(None, 1)
            name = l[0]
            length = len(rec.sequence)
            d_len[name] = length

    df = pd.read_csv(infile, sep='\t', header=0)
    old_cols = df.columns.tolist()
    _l = df['seqname'].map(d_len)
    df['length'] = _l
    cur_cols = ['seqname', 'length'] + old_cols[1:]
    df = df[cur_cols]
    df.to_csv(outfile, sep='\t', index=False, float_format='%.3f')

if __name__ == '__main__':
    main()

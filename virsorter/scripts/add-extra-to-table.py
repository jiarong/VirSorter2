#!/usr/bin/env python

import sys
import os
import screed
import numpy as np
import pandas as pd


def main():
    '''Add sequence length to table

    Example:
        python add-extra-to-table.py.py <score.tsv> <viral.fa> <out.tsv>

        <score.tsv>: merged classifcation results showing probability of being
            viral with different classifers
        <viral.fa>: viral contig file with seqname matched to <score.tsv>

    '''
    if len(sys.argv) != 4:
        mes = 'python {} <score.tsv> <viral.fa> <out.tsv>\n'
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    score_f = sys.argv[1]
    seqfile = sys.argv[2]
    outfile = sys.argv[3]

    d_info = {}
    with screed.open(seqfile) as sp:
        for rec in sp:
            header = rec.name
            name, desc = header.split(None ,1)
            _d = dict(i.split(':') for i in desc.split('||'))
            hallmark = _d.get('hallmark', np.nan)
            viral = _d.get('viral', np.nan)
            cellular = _d.get('cellular', np.nan)
            group = _d.get('group', 'NA')
            length = len(rec.sequence)
            d_info[name] = {'max_score_group': group, 'length':length, 
                    'hallmark':hallmark, 'viral':viral, 'cellular':cellular}

    df_info = pd.DataFrame.from_dict(d_info, orient='index')
    # seqname is index 
    df_info['seqname'] = df_info.index

    df = pd.read_csv(score_f, sep='\t', header=0)
    if len(df) == 0:
        cols = df.columns.values.tolist() + ['max_score', 'max_score_group', 
                'length', 'hallmark', 'viral', 'cellular']
        with open(outfile, 'w') as fw:
            fw.write('{}\n'.format('\t'.join(cols)))
            sys.exit(0)

    _df = df.drop('seqname', axis=1)
    _df.index = df.loc[:,'seqname']
    max_ser = _df.max(axis=1)
    df['max_score'] = max_ser.values

    max_score_group_ser = _df.idxmax(axis=1)
        
    for i in max_score_group_ser.index:
        try:
            if max_score_group_ser.loc[i] != df_info.loc[i, 'max_score_group']:
                df_info.loc[i, 'max_score_group'] = max_score_group_ser.loc[i]
        except ValueError as e:
            mes = ('*** Duplicate seqnames are found in input contig '
                   'sequence file; Please fix them and rerun..\n')
            sys.stderr.write(mes)
            sys.exit(1)

    df_merge = pd.merge(df, df_info, on=['seqname'], how='right')
    df_merge.to_csv(outfile, sep='\t', na_rep='NaN', 
            index=False, float_format='%.3f')

if __name__ == '__main__':
    main()

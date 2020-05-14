#!/usr/bin/env python

import sys
import os
import screed
import pandas as pd
from virsorter.config import DEFAULT_CONFIG

D = DEFAULT_CONFIG['GROUP_INFO']
DEFAULT_LEN_CUTOFF = 3000

def main():
    '''Remove seqs are shorter than a cutoff and without hallmark genes 
    
    Example:
        python remove-short-seq-wo-hallmark.py <file.hmkcnt> <contig.fa>

        <file.hmkcnt>: file with hallmark gene count of each contig in each viral group 
        <contig.fa>: contig sequence file

    '''
    if len(sys.argv) != 3:
        mes = '*** Usage: python {} <file.hmkcnt> <contig.fa>\n'
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    hmk_f = sys.argv[1]
    fa_f = sys.argv[2]

    df = pd.read_csv(hmk_f, sep='\t', header=0, index_col='seqname')
    group_ser = df.idxmax(axis=1)
    max_cnt_ser = df.max(axis=1)
    df_max = pd.concat([group_ser, max_cnt_ser], keys=['group', 'max_cnt'], axis=1)
    #print(df_max)

    st = set(df_max.index.unique())
    with screed.open(fa_f) as sp:
        for rec in sp:
            header = rec.name
            lis = header.split(None, 1)
            name = lis[0]
            seq = rec.sequence
            if not name in st:
                try:
                    _s = lis[1]
                    if ':' in _s:
                        _d = dict(i.split(':') for i in _s.split('||'))
                        group = _d['group']
                    else:
                        group = _s.strip()
                    cutoff = D[group]['MIN_SIZE_ALLOWED_WO_HALLMARK_GENE']
                except IndexError as e:
                    cutoff = DEFAULT_LEN_CUTOFF
                if len(seq) < cutoff:
                    continue
            elif df_max.at[name, 'max_cnt'] == 0:
                group = df_max.at[name, 'group']
                if len(seq) < D[group]['MIN_SIZE_ALLOWED_WO_HALLMARK_GENE']:
                    continue

            sys.stdout.write('>{}\n{}\n'.format(header, seq))

if __name__ == '__main__':
    main()

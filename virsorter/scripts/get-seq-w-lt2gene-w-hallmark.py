#!/usr/bin/env python

import sys
import os
import screed
import pandas as pd

script_dir = os.path.dirname(os.path.abspath(__file__))
snakefile_dir = os.path.dirname(script_dir)
pkg_dir = os.path.dirname(snakefile_dir)
sys.path.append(pkg_dir)
from virsorter.config import get_default_config

DEFAULT_CONFIG = get_default_config()

D = DEFAULT_CONFIG['GROUP_INFO']
DEFAULT_LEN_CUTOFF = 3000

def main():
    '''Get short contigs (with less than 2 genes) but has hallmark genes

    Example:
        python get-seq-w-lt2genes-w-hallmark.py \
                <file.hmkcnt> <all.proba.tsv> <contig.fa>

        <file.hmkcnt>: file with hallmark gene count for each contig from each
            virla groups (output from get-hallmark-cnt-for-each-seq.py)
        <all.proba.tsv>: file with proba being viral for each contig with more
            than 2 genes from each viral groups
        <contig.fa>: config sequence file

    '''
    if len(sys.argv) != 4:
        mes = '*** Usage: python {} <file.hmkcnt> <all.proba.tsv> <contig.fa>\n'
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    hmk_f = sys.argv[1]
    proba_f = sys.argv[2]
    fa_f = sys.argv[3]

    df = pd.read_csv(hmk_f, sep='\t', header=0, index_col='seqname')
    group_ser = df.idxmax(axis=1)
    max_cnt_ser = df.max(axis=1)
    df_max = pd.concat([group_ser, max_cnt_ser], 
            keys=['group', 'max_cnt'], axis=1)
    #print(df_max)
    sel = df_max['max_cnt'] >= 1
    df_max = df_max.loc[sel, :]
    seqnames_w_hm = df_max.index.unique()

    df_proba = pd.read_csv(proba_f, sep='\t', header=0, index_col='seqname')
    seqnames_w_mt2genes = df_proba.index.unique()

    st = set(seqnames_w_hm).difference(set(seqnames_w_mt2genes))
    with screed.open(fa_f) as sp:
        for rec in sp:
            header = rec.name
            lis = header.split(None, 1)
            name = lis[0]
            name = name.rsplit('||', 1)[0]
            seq = rec.sequence
            if not name in st:
                continue
            try:
                shape = lis[1]
            except IndexError as e:
                shape = 'NA'
            cnt = df_max.at[name, 'max_cnt']
            group = df_max.at[name, 'group']
            mes = '>{}  shape:{}||hm_cnt:{}||hm_group:{}\n{}\n'
            sys.stdout.write(mes.format(name, shape, cnt, group, seq))

if __name__ == '__main__':
    main()

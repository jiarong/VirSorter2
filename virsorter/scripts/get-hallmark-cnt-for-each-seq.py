#!/usr/bin/env python

import sys
import os
import screed
from collections import Counter
import click
import pandas as pd

def main():
    '''Get hallmark gene cnt for each contig

    Example:
        python get-hallmark-cnt-for-each-seq.py <outfile> \
                "A,B.." "<A.hallmark>,<B.hallmark>.." "<A.tax>,<B.tax>.."

        <outfile>: output file with hallmark gene counts for each contig from
            each viral groups
        "A,B..": list of viral group labels (comma separated)
        "<A.hallmark>,<B.hallmark>..": list of hallmark genes for viral groups
            (comma separated)
        "<A.tax>,<B.tax>..": list of .tax file from
            extract-feature-from-hmmout.py for viral groups (comma separated)
        
    '''
    if len(sys.argv) != 5:
        mes = ('*** Usage: python {} <outfile> "A,B.." '
                '"<A.hallmark>,<B.hallmark>.." '
                '"<A.tax>,<B.tax>.."\n')
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    outfile = sys.argv[1]
    groups = [ g.strip() for g in sys.argv[2].split(',') ]
    hallmark_fs = [ f.strip() for f in sys.argv[3].split(',') ]
    tax_fs = [ f.strip() for f in sys.argv[4].split(',') ]

    d_hallmark_cnt = {}
    for group, hallmark_f, tax_f in zip(groups, hallmark_fs, tax_fs):
        d_hallmark_hmm = {}
        with open(hallmark_f) as fp:
            for line in fp:
                if line.startswith('#'):
                    continue
                line = line.rstrip()
                lis = line.split('\t')
                hmm_name = lis[0]
                gene = lis[1]
                cutoff = float(lis[2])
                d_hallmark_hmm[hmm_name] = gene, cutoff


        d_hallmark_cnt[group] = {}
        with open(tax_f) as fp2:
            for line in fp2:
                if line.startswith('#'):
                    continue
                line = line.rstrip()
                lis = line.split('\t')
                name = lis[0]
                tax = lis[1]
                hmm = lis[2]
                try:
                    score = float(lis[3])
                except ValueError as e:
                    score = -10000000

                contig_w_rbs = name.rsplit('_', 1)[0]
                contig = contig_w_rbs.rsplit('||', 1)[0]

                cur_cnt = d_hallmark_cnt[group].setdefault(contig, 0)
                if hmm in d_hallmark_hmm:
                    if score >= d_hallmark_hmm[hmm][1]:
                        d_hallmark_cnt[group][contig] = cur_cnt + 1



    df = pd.DataFrame.from_dict(d_hallmark_cnt)
    df.to_csv(outfile, sep='\t', index=True, header=True, 
            index_label='seqname', float_format='%.0f') 

    #group_ser = df.idxmax(axis=1)
    #max_cnt_ser = df.max(axis=1)
    #df_max = pd.concat([group_ser, max_cnt_ser], keys=['group', 'max_cnt'], axis=1)
    #print(df_max)


if __name__ == '__main__':
    main()

#!/usr/bin/env python

import sys
import os
import pandas as pd

def main():
    if len(sys.argv) != 3:
        mes = '[INFO] Usage: python {} <final-viral-boundary.tsv.tmp> <final-viral-score.tsv>\n'
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    bdy_f = sys.argv[1]
    score_f = sys.argv[2]

    df1 = pd.read_csv(bdy_f, sep='\t', header=0)
    df2 = pd.read_csv(score_f, sep='\t', header=0)

    _d = {'seqname': 'seqname_new', 'max_score': 'final_max_score', 'max_score_group': 'final_max_score_group'}
    df2_sel = df2.loc[:, list(_d.keys())]
    df2_sel = df2_sel.rename(_d, axis=1)

    df_new = df1.merge(df2_sel, how='inner', on='seqname_new')

    df_new.to_csv('/dev/stdout', sep='\t', header=True, index=False, na_rep='NA')

if __name__ == '__main__':
    main()

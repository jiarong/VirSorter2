#!/usr/bin/env python

import sys
import os
import numpy as np
import pandas as pd

script_dir = os.path.dirname(os.path.abspath(__file__))
snakefile_dir = os.path.dirname(script_dir)
pkg_dir = os.path.dirname(snakefile_dir)
sys.path.append(pkg_dir)
from virsorter.config import get_default_config

DEFAULT_CONFIG = get_default_config()

PROVIRUS_MIN_PEAK_PROBA = DEFAULT_CONFIG['PROVIRUS_MIN_PEAK_PROBA']

def main():
    if len(sys.argv) != 5:
        mes = ('*** Usage: python {} '
                '"<A/all.pdg.prv>,<B/all.pdg.prv>" '
                '"A,B" '
                '<viral-partseq.tsv> '
                '<viral-fullseq.tsv>\n')

        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    fs = [f.strip() for f in sys.argv[1].split(',')]
    groups = [g.strip() for g in sys.argv[2].split(',')]
    f_partial = sys.argv[3]
    f_full = sys.argv[4]

    df_lis = []
    for group, f in zip(groups, fs):
        df = pd.read_csv(f, sep='\t', header=0)
        # remove "||rbs:common" or "||rbs:NCLDV" suffix from seqname
        seqnames = [name.rsplit('||', 1)[0] for name in df['seqname']]
        df['seqname'] = seqnames
        df['group'] = group
        df['size'] = df['trim_orf_index_end'] - df['trim_orf_index_start']
        df_lis.append(df)

    df_merged = pd.concat(df_lis)
    sel = df_merged['partial'] == 0
    df_merged_full = df_merged.loc[sel,:]
    df_merged_partial = df_merged.loc[~sel,:]

    # full
    # a copy is made in background with sorting
    seqnames_full = set(df_merged_full['seqname'].unique())
    _df = df_merged_full.sort_values(
            by=['seqname', 'trim_pr', 'size'], 
            ascending=False,
    )
    _df = _df.drop_duplicates(['seqname'])
    _df = _df.drop(['size'], axis=1)
    _df.to_csv(f_full, sep='\t', index=False, 
            na_rep='nan', float_format='%.3f')


    # partial
    with open(f_partial, 'w') as fw_partial:
        # ignore those already processed in full (whole seq as viral)
        sel_seqname = df_merged_partial['seqname'].isin(seqnames_full)
        df_merged_partial = df_merged_partial.loc[~sel_seqname,:]
	### moved the filtering to provirus.py
        #sel_peak_proba = \
        #        df_merged_partial['prox_pr_max'] >= PROVIRUS_MIN_PEAK_PROBA
        #df_merged_partial = df_merged_partial.loc[sel_peak_proba,:]

        # a copy is made in background with sorting
        df_merged_partial = df_merged_partial.sort_values(
                by=['seqname', 'size', 'prox_pr_max'], 
                ascending=False,
        )

	### moved the filtering to provirus.py
        #min_genome_sizes = df_merged_partial['group'].map(
        #        lambda x: DEFAULT_CONFIG['GROUP_INFO'][x]['MIN_GENOME_SIZE'])
        #sel_size = (
        #        (df_merged_partial['trim_bp_end'] - 
        #            df_merged_partial['trim_bp_start']) >= min_genome_sizes
        #)
        #df_merged_partial = df_merged_partial.loc[sel_size,:]
        df_merged_partial = df_merged_partial.drop(['size'], axis=1)
        # write the headers
        string_lis = df_merged_partial.columns
        fw_partial.write('{}\n'.format('\t'.join(string_lis)))

        seqnames = df_merged_partial['seqname']
        ind_starts = df_merged_partial['trim_orf_index_start']
        ind_ends = df_merged_partial['trim_orf_index_end']
        last_seqname = None
        for i in range(len(seqnames)):
            seqname = seqnames.iloc[i]

            cur_st = set(range(ind_starts.iloc[i], ind_ends.iloc[i]))
            if seqname != last_seqname:
                ser = df_merged_partial.iloc[i,:]
                string_lis = [
                        '{:.3f}'.format(i) if isinstance(i, float) \
                                else '{}'.format(i) for i in ser
                ]
                fw_partial.write('{}\n'.format('\t'.join(string_lis)))
                # reset for next new seqname
                existing_st = cur_st

            else:
                if len(cur_st.intersection(existing_st)) != 0:
                    continue
                ser = df_merged_partial.iloc[i,:]
                string_lis = [
                        '{:.3f}'.format(i) if isinstance(i, float) \
                                else '{}'.format(i) for i in ser
                ]
                fw_partial.write('{}\n'.format('\t'.join(string_lis)))
                existing_st.update(cur_st)

            last_seqname = seqname

if __name__ == '__main__':
    main()

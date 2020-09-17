#!/usr/bin/env python

import sys
import os
from collections import Counter
import screed
import pandas as pd
import numpy as np

script_dir = os.path.dirname(os.path.abspath(__file__))
snakefile_dir = os.path.dirname(script_dir)
pkg_dir = os.path.dirname(snakefile_dir)
sys.path.append(pkg_dir)
from virsorter.config import get_default_config, set_logger

DEFAULT_CONFIG = get_default_config()
TAX_FEATURE_LIST = DEFAULT_CONFIG['TAX_FEATURE_LIST']
FASTA_DESC_FORMAT_TEMPLATE = DEFAULT_CONFIG['FASTA_DESC_FORMAT_TEMPLATE']


def main():
    '''Add sequence length to table

    Example:
        python add-extra-to-fullseq-fasta-header.py \
                <viral-fullseq.fa> \
                "A/<all.pdg.ftr>,B/<all.pdg.ftr>" \
                "A,B"

        <viral-fullseq.fa>: viral contigs with score >= PROBA_CUTOFF
        <all.pdg.ftr>: feature table 
        A: viral group name
    '''

    if len(sys.argv) != 4:
        mes = ('python {} viral-fullseq.fa '
                '"A/<all.pdg.ftr>,B/<all.pdg.ftr>" '
                '"A,B"\n')
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    seqfile = sys.argv[1]
    ftr_fs = sys.argv[2]
    groups = sys.argv[3]

    d_group2name = {}
    with screed.open(seqfile) as sp:
        for rec in sp:
            header = rec.name
            name, desc = header.split(None ,1)
            _d = dict(i.split(':') for i in desc.split('||'))
            group = _d['group']
            st = d_group2name.setdefault(group, set())
            st.add(name)

    ftr_fs = [f.strip() for f in ftr_fs.split(',')]
    groups = [group.strip() for group in groups.split(',')]

    d_name2info = {}
    for ftr_f, group in zip(ftr_fs, groups):
        df_ftr = pd.read_csv(ftr_f, sep='\t', header=0)
        seqname_lis = df_ftr['seqname']
        seqname_ori_lis = [i.rsplit('||')[0] for i in seqname_lis]
        df_ftr.index = seqname_ori_lis

        name_st = d_group2name.get(group, set())
        for seqname_ori in name_st:
            hallmark = df_ftr.loc[seqname_ori, 'hallmark']
            viral = df_ftr.loc[seqname_ori, 'vir']
            cellular = sum(df_ftr.loc[seqname_ori, ['arc', 'bac', 'euk']])
            d_name2info[seqname_ori] = 1, np.nan, viral, cellular, hallmark


    with screed.open(seqfile) as sp:
        for rec in sp:
            header = rec.name
            name, desc = header.split(None ,1)
            d_desc = dict(i.split(':') for i in desc.split('||'))
            start_ind, end_ind, viral, cellular, hallmark = d_name2info[name]
            seq = rec.sequence
            length = len(seq)
            d_desc['start'] = 1
            d_desc['end'] = length
            d_desc['start_ind'] = start_ind
            d_desc['end_ind'] = end_ind
            d_desc['viral'] = viral
            d_desc['cellular'] = cellular

            try:
                d_desc['score']= float(d_desc['score'])
            except IndexError as e:
                d_desc['score']= np.nan
            d_desc['shape'] = d_desc.get('shape',  'NA')
            d_desc['hallmark'] = int(hallmark)

            desc = FASTA_DESC_FORMAT_TEMPLATE.format(**d_desc)
            mes = f'>{name}  {desc}\n{seq}\n'
            sys.stdout.write(mes)

if __name__ == '__main__':
    main()

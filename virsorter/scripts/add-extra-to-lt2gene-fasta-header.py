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
from virsorter.utils import (
        load_rbs_category, df_tax_per_config, parse_gff, 
        extract_feature_tax
)

DEFAULT_CONFIG = get_default_config()
TAX_FEATURE_LIST = DEFAULT_CONFIG['TAX_FEATURE_LIST']
FASTA_DESC_FORMAT_TEMPLATE = DEFAULT_CONFIG['FASTA_DESC_FORMAT_TEMPLATE']


def main():
    '''Add sequence length to table

    Example:
        python add-extra-to-lt2gene-fasta-header.py \
                <viral-lt2gene-w-hallmark.fa> \
                "A/<all.pdg.gff>,B/<all.pdg.gff>" \
                "A/<all.pdg.hmm.tax>,B/<all.pdg.hmm.ftr>" \
                "A,B"


        <viral-lt2gene-w-hallmark.fa>: viral contigs with less than \
                two genes and hallmark gene
        <all.pdg.gff>: gff file from prodigal
        <all.pdg.hmm.tax>: table with best hit of each gene, bit score, 
                and taxonomy
        A: viral group name


    '''
    if len(sys.argv) != 5:
        mes = ('python {} viral-lt2gene-w-hallmark.fa '
                '"A/<all.pdg.gff>,B/<all.pdg.gff>" '
                '"A/<all.pdg.hmm.tax>,B/<all.pdg.hmm.tax>" '
                '"A,B"\n')
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    seqfile = sys.argv[1]
    gff_fs = sys.argv[2]
    tax_fs = sys.argv[3]
    groups = sys.argv[4]

    d_group2name = {}
    with screed.open(seqfile) as sp:
        for rec in sp:
            header = rec.name
            name, desc = header.split(None ,1)
            _d = dict(i.split(':') for i in desc.split('||'))
            group = _d['group']
            st = d_group2name.setdefault(group, set())
            st.add(name)

    gff_fs = [f.strip() for f in gff_fs.split(',')]
    tax_fs = [f.strip() for f in tax_fs.split(',')]
    groups = [group.strip() for group in groups.split(',')]

    d_name2info = {}
    for gff_f, tax_f, group in zip(gff_fs, tax_fs, groups):
        gen_gff = parse_gff(gff_f) 
        seqname_lis = [l[0] for l in gen_gff]
        seqname_ori_lis = [i.rsplit('||')[0] for i in seqname_lis]
        d_name2cnt=Counter(seqname_ori_lis)

        name_st = d_group2name.get(group, set())
        for seqname_ori in name_st:
            total_gene_cnt = d_name2cnt[seqname_ori]
            ind = seqname_ori_lis.index(seqname_ori)
            seqname = seqname_lis[ind]
            df_tax_sel = df_tax_per_config(tax_f, seqname)
            sel_index_w_hallmark = [] # donot need hallmark cnt here 
            l_tax = extract_feature_tax(df_tax_sel, 
                    sel_index_w_hallmark=sel_index_w_hallmark, 
                    total_gene_cnt=total_gene_cnt)

            vir_ind = TAX_FEATURE_LIST.index('vir')
            arc_ind = TAX_FEATURE_LIST.index('arc')
            bac_ind = TAX_FEATURE_LIST.index('bac')
            euk_ind = TAX_FEATURE_LIST.index('euk')
            viral = l_tax[vir_ind] 
            cellular = l_tax[arc_ind] + l_tax[bac_ind] + l_tax[euk_ind]
            d_name2info[seqname_ori] = 1, total_gene_cnt, viral, cellular


    with screed.open(seqfile) as sp:
        for rec in sp:
            header = rec.name
            name, desc = header.split(None ,1)
            d_desc = dict(i.split(':') for i in desc.split('||'))
            start_ind, end_ind, viral, cellular = d_name2info[name]
            seq = rec.sequence
            length = len(seq)
            d_desc['start'] = 1
            d_desc['end'] = length
            d_desc['start_ind'] = start_ind
            d_desc['end_ind'] = end_ind
            d_desc['viral'] = viral
            d_desc['cellular'] = cellular
            d_desc['score'] = np.nan
            d_desc['hallmark'] = int(d_desc['hallmark'])

            desc = FASTA_DESC_FORMAT_TEMPLATE.format(**d_desc)
            mes = f'>{name}  {desc}\n{seq}\n'
            sys.stdout.write(mes)

if __name__ == '__main__':
    main()

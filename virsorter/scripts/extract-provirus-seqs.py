#!/usr/bin/env python

import sys
import os
import screed
import logging
import pandas as pd

script_dir = os.path.dirname(os.path.abspath(__file__))
snakefile_dir = os.path.dirname(script_dir)
pkg_dir = os.path.dirname(snakefile_dir)
sys.path.append(pkg_dir)
from virsorter.config import get_default_config, set_logger
from virsorter.utils import FASTA_DESC_FORMAT_TEMPLATE

DEFAULT_CONFIG = get_default_config()

set_logger()

def get_fasta_desc_d(ser, shape):
    d_desc = {}
    #trim_start_ind = ser.iloc[0]
    #trim_end_ind = ser.iloc[1]
    #trim_start_bp = ser.iloc[2]
    #trim_end_bp = ser.iloc[3]
    #full_start_ind = ser.iloc[-6]
    #full_end_ind = ser.iloc[-5]
    #group = ser.iloc[-1]
    #hallmark = ser.iloc[-2]
    #score = ser.iloc[4]
    trim_start_ind = ser.loc['trim_orf_index_start']
    trim_end_ind = ser.loc['trim_orf_index_end']
    trim_start_bp = ser.loc['trim_bp_start']
    trim_end_bp = ser.loc['trim_bp_end']
    full_start_ind = ser.loc['full_orf_index_start']
    full_end_ind = ser.loc['full_orf_index_end']
    group = ser.loc['group']
    viral = ser.loc['vir']
    cellular = sum(ser.loc[['arc', 'bac', 'euk']])
    hallmark = ser.loc['hallmark_cnt']
    score = ser.loc['trim_pr']

    d_desc['shape'] = shape
    d_desc['start'] = trim_start_bp
    d_desc['end'] = trim_end_bp
    d_desc['start_ind'] = trim_start_ind
    d_desc['end_ind'] = trim_end_ind
    d_desc['viral'] = viral
    d_desc['cellular'] = cellular
    d_desc['group'] = group
    d_desc['score'] = score
    d_desc['hallmark'] = hallmark

    return d_desc


def main():
    '''Extract substring sequences based in boundries

    Example:
        python extract-provirus-seqs.py \
                <contig.fa> <fullseq.tsv> <partseq.tsv> \
                <fullsesq-trim.fa> <partseq.fa>

        <config.fa>: contig file
        <fullseq.tsv>: full seq boudary info
        <partseq.tsv>: partial seq bourdary info
        <fullseq-trim.fa>: full seq trimmed according to boudary
        <partseq.fa>: partial seq extracted according to boudary

    '''
    if len(sys.argv) != 6:
        mes = ('*** python {} <input.contig> <input.full> <input.partial> '
                '<output.lytic> <output.lysogenic>\n')
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    seqfile = sys.argv[1]
    full_f = sys.argv[2]
    part_f = sys.argv[3]
    lytic_f = sys.argv[4]
    lyso_f = sys.argv[5]

    df_full = pd.read_csv(full_f, header=0, sep='\t', index_col='seqname')
    df_part = pd.read_csv(part_f, header=0, sep='\t', index_col='seqname')

    with screed.open(seqfile) as sp, \
            open(lytic_f, 'w') as fw_lytic, \
            open(lyso_f, 'w') as fw_lyso:
        st_full = set(df_full.index)
        st_part = set(df_part.index)
        for rec in sp:
            lis = rec.name.split(None, 1)
            seqname = lis[0]
            try:
                shape = lis[1]
            except IndexError as e:
                mes = '{} in {} does not have shape (linear or circular) info'
                logging.error(mes.format(seqname, seqfile))
                shape = 'NA'
            seqname = seqname.rsplit('||', 1)[0] # remove ||rbs:common
            if seqname in st_full:
                ser = df_full.loc[seqname, :]
                full_start_ind = ser.loc['full_orf_index_start']
                full_end_ind = ser.loc['full_orf_index_end']
                full_start_bp = ser.loc['full_bp_start']
                full_end_bp = ser.loc['full_bp_end']

                d_desc = get_fasta_desc_d(ser, shape)
                trim_start_bp = d_desc['start']
                trim_end_bp = d_desc['end'] 
                trim_start_ind = d_desc['start_ind']
                trim_end_ind = d_desc['end_ind']

                desc = FASTA_DESC_FORMAT_TEMPLATE.format(**d_desc)

                if shape == 'circular':
                    trim_start_bp_adj = trim_start_bp - full_start_bp + 1
                    trim_end_bp_adj = trim_end_bp - full_start_bp + 1
                else:
                    trim_start_bp_adj = trim_start_bp
                    trim_end_bp_adj = trim_end_bp
                
                if (trim_start_ind == full_start_ind and 
                        trim_end_ind == full_end_ind):
                    # not trimmed; save to lytic
                    seq = rec.sequence[(trim_start_bp_adj-1):trim_end_bp]
                else:
                    # trimmed; save to lytic (fullseq); decided not to 
                    #   interpret lytic or lyso here
                    seq = rec.sequence[(trim_start_bp_adj-1):trim_end_bp]

                mes = f'>{seqname}||full  {desc}\n{seq}\n'
                fw_lytic.write(mes)

            elif seqname in st_part:
                _df = df_part.loc[df_part.index == seqname, :]
                for i in range(len(_df)):
                    ser = _df.iloc[i, :]
                    d_desc = get_fasta_desc_d(ser, shape)
                    trim_start_bp = d_desc['start']
                    trim_end_bp = d_desc['end'] 
                    trim_start_ind = d_desc['start_ind']
                    trim_end_ind = d_desc['end_ind']
                    desc = FASTA_DESC_FORMAT_TEMPLATE.format(**d_desc)
                    seq = rec.sequence[(trim_start_bp-1):trim_end_bp]

                    #save to lyso
                    mes = f'>{seqname}||{i}_partial  {desc}\n{seq}\n'
                    fw_lyso.write(mes)

if __name__ == '__main__':
    main()

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
from virsorter.config import set_logger

set_logger()

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
            seqname = seqname.rsplit('||', 1)[0] # remove ||rbs:common
            if seqname in st_full:
                ser = df_full.loc[seqname, :]
                trim_start_ind = ser.iloc[0]
                trim_end_ind = ser.iloc[1]
                trim_start_bp = ser.iloc[2]
                trim_end_bp = ser.iloc[3]
                full_start_ind = ser.iloc[-6]
                full_end_ind = ser.iloc[-5]
                group = ser.iloc[-1]
                hallmark = ser.iloc[-2]
                score = ser.iloc[4]
                
                if (trim_start_ind == full_start_ind and 
                        trim_end_ind == full_end_ind):
                    # not trimmed; save to lytic
                    mes = ('>{}||full  '
                        'shape:{}||start:{}||end:{}||'
                        'group:{}||score:{}||hallmark:{}\n{}\n')
                    fw_lytic.write(
                            mes.format(seqname, shape, trim_start_bp, 
                                trim_end_bp, group, score, hallmark, 
                                rec.sequence)
                    )
                else:
                    # trimmed; save to lytic; decided not to 
                    #   interpret lytic or lyso here
                    mes = ('>{}||full  '
                        'shape:{}||start:{}||end:{}||'
                        'group:{}||score:{}||hallmark:{}\n{}\n')
                    fw_lytic.write(
                            mes.format(seqname, shape, trim_start_bp, 
                                trim_end_bp, group, score, hallmark, 
                                rec.sequence[(trim_start_bp-1):trim_end_bp])
                    )

            elif seqname in st_part:
                _df = df_part.loc[df_part.index == seqname, :]
                for i in range(len(_df)):
                    ser = _df.iloc[i, :]
                    trim_start_ind = ser.iloc[0]
                    trim_end_ind = ser.iloc[1]
                    trim_start_bp = ser.iloc[2]
                    trim_end_bp = ser.iloc[3]
                    full_start_ind = ser.iloc[-6]
                    full_end_ind = ser.iloc[-5]
                    group = ser.iloc[-1]
                    hallmark = ser.iloc[-2]
                    score = ser.iloc[4]
                    #save to lyso
                    mes = ('>{}||{}index_partial  '
                            'shape:{}||start:{}||end:{}||'
                            'group:{}||score:{}||hallmark:{}\n{}\n')
                    fw_lyso.write(
                            mes.format(seqname, i, shape, trim_start_bp, 
                                trim_end_bp, group, score, hallmark,
                                rec.sequence[(trim_start_bp-1):trim_end_bp])
                    )

if __name__ == '__main__':
    main()

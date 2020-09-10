#!/usr/bin/env python
import sys
import os
import warnings
from collections import OrderedDict, Counter
import logging

import click
import numpy as np
import pandas as pd

script_dir = os.path.dirname(os.path.abspath(__file__))
snakefile_dir = os.path.dirname(script_dir)
pkg_dir = os.path.dirname(snakefile_dir)
sys.path.append(pkg_dir)
from virsorter.config import get_default_config, set_logger
from virsorter.utils import (
        load_rbs_category, df_tax_per_config, parse_hallmark_hmm,
        parse_gff, extract_feature_gff, get_feature
)

DEFAULT_CONFIG = get_default_config() 

TOTAL_FEATURE_LIST = DEFAULT_CONFIG['TOTAL_FEATURE_LIST']
SELECT_FEATURE_LIST = DEFAULT_CONFIG['SELECT_FEATURE_LIST']

set_logger()

class gff_ftr:
    '''The class to process provirus extraction and trimming full seqs
    '''
    def __init__(self, gff_f, rbs_cat_f, outfile):
        '''Initialize object

        Args:
            gff_f (str): gff file from prodigal
            rbs_cat_f (str): file mapping rbs sites and categories (rbs directory of database)
            outfile (str): output file

        '''
        self.gff_f = gff_f
        self.rbs_cat_f = rbs_cat_f
        self.outfile = outfile

    def load_data(self):
        self.gff_gen = parse_gff(self.gff_f)
        self.rbs_cat_d = load_rbs_category(self.rbs_cat_f)
        self.gff_mat_colnames = ('orf_index', 'start', 'end', 'strand', 
                'partial', 'start_type', 'gc_cont', 'rbs_motif')

    def process_one_contig(self, seqname, mat):
        '''Process a contig

        Args:
            seqname (str): sequence name
            mat (list of list): list of list of info related to genomic
                features

        '''
        df_gff = pd.DataFrame(mat, columns=self.gff_mat_colnames)

        l = extract_feature_gff(df_gff, self.rbs_cat_d)
        # meet 1) >= 2 genes; 2) at least 1 full gene
        if len(l) != 0:
            if self.header_written == False:
                header_lis = TOTAL_FEATURE_LIST[:len(l)]
                self.fw.write(
                        'seqname\t{}\n'.format('\t'.join(header_lis))
                )
                self.header_written = True

            l = [str(i) for i in l]
            self.fw.write('{}\t{}\n'.format(seqname, '\t'.join(l)))

    def extract_gff_ftr(self):
        '''Iterate through each contigs in gff and boudary list in name2loc_d for each contig to get score
        '''
        mat = []
        last_seqname = None
        self.header_written = False
        with open(self.outfile, 'w') as self.fw:

            for lis in self.gff_gen:
                seqname = lis[0]
                # process a contig
                if last_seqname != None and last_seqname != seqname:
                    #logging.info('Processing {}'.format(last_seqname))
                    self.process_one_contig(last_seqname, mat)

                    # reset mat and last_seqname for next iter
                    mat = []
                    last_seqname = None

                # items in lis:
                # seqname, seqlen, orf_index, start, end, strand, partial, 
                # start_type, gc_cont, rbs_motif
                mat.append(lis[2:])
                last_seqname = seqname

            if len(mat) != 0:
                #logging.info('Processing {}'.format(last_seqname))
                self.process_one_contig(last_seqname, mat)


CONTEXT_SETTINGS = {'help_option_names': ['-h', '--help']}
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('rbs', type=click.Path()) 
@click.argument('gff', type=click.Path())
@click.argument('outfile', type=click.Path())
def main(gff, rbs, outfile):
    '''
    this function do the following:
        apply a classifier and assign classification score to each seq in combined seqfile of full and partial seqs
    
    \f
    :param rbs: rbs cattory file (rbs-catetory.tsv)
    :param gff: .gff file from prodigal
    :param outfile: feature from gff 
    :return: 0 if success else 1
    :rtype: int
    '''

    gff_ftr_ins = gff_ftr(gff, rbs, outfile)
    gff_ftr_ins.load_data()
    gff_ftr_ins.extract_gff_ftr()

if __name__ == '__main__':
    main()

#!/usr/bin/env python
import sys
import os
import warnings
from collections import OrderedDict, Counter
import logging

import screed
import joblib
import click
import sklearn
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

GENE_OVERLAP_MIN = DEFAULT_CONFIG['GENE_OVERLAP_MIN']
CLASSIFY_THREADS = DEFAULT_CONFIG['CLASSIFY_THREADS']
TOTAL_FEATURE_LIST = list(DEFAULT_CONFIG['TOTAL_FEATURE_LIST'])
SELECT_FEATURE_LIST = list(DEFAULT_CONFIG['SELECT_FEATURE_LIST'])

set_logger()


class clf_trim:
    '''The class to process provirus extraction and trimming full seqs
    '''
    def __init__(self, gff_f, tax_f, rbs_cat_f, model_f, 
            full_and_part_seqfile, outfile, hallmark_f=None, group=None):
        '''Initialize object

        Args:

            rbs_cat_d (dict): mapping hmm name to gene name and hmmsearch bit score
            cutoff from rbs_cat_d
            gff_f (str): gff file from prodigal
            tax_f (str): info loaded from .tax file (output of extract-feature-from-hmmout.py
            rbs_cat_f (str): file mapping rbs sites and categories (rbs directory of database)
            model_f (str): classifier model file
            full_and_part_seqfile (str): sequence file
            outfile (str): output file

            hallmark_f (str): hallmark gene list file
            group: viral group name

        '''
        self.gff_f = gff_f
        self.tax_f = tax_f
        self.rbs_cat_f = rbs_cat_f
        self.model_f = model_f
        self.outfile = outfile
        self.hallmark_f = hallmark_f
        self.full_and_part_seqfile = full_and_part_seqfile
        self.group = group

    def load_data(self):
        self.gff_gen = parse_gff(self.gff_f)

        df_tax_all = pd.read_csv(self.tax_f, sep='\t', header=None, 
                names=['orfname', 'tax', 'hmm', 'score'])
        try:
            seqnames, indice = zip(
                    *[orfname.rsplit('_', 1) 
                        for orfname in df_tax_all['orfname']]
            )
        except ValueError as e:
            seqnames = []
            indice = []

        df_tax_all['seqname'] = seqnames
        # convert to int to match orf_index in df_gff
        indice = [ int(i) for i in indice ]
        df_tax_all['orf_index'] = indice
        self.df_tax_all = df_tax_all

        self.rbs_cat_d = load_rbs_category(self.rbs_cat_f)

        if self.hallmark_f != None:
            self.d_hallmark_hmm = parse_hallmark_hmm(self.hallmark_f)
        else:
            self.d_hallmark_hmm = None

        self.name2loc_d = {}
        with screed.open(self.full_and_part_seqfile) as sp:
            for rec in sp:
                header = rec.name
                name, desc = header.split(None ,1)
                seqname_ori = name.rsplit('||', 1)[0]
                _d = dict(i.split(':') for i in desc.split('||'))
                seq_len = len(rec.sequence)

                self.name2loc_d.setdefault(seqname_ori, {})
                self.name2loc_d[seqname_ori].setdefault('seqname', [])
                self.name2loc_d[seqname_ori].setdefault('loc', [])

                self.name2loc_d[seqname_ori]['seqname'].append(name)
                self.name2loc_d[seqname_ori]['loc'].append(
                        (int(_d.get('start', 0)), 
                            int(_d.get('end', seq_len)))
                )

        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=DeprecationWarning)
            warnings.filterwarnings('ignore', category=FutureWarning)
            model = joblib.load(self.model_f)
            try:
                model.named_steps.gs.best_estimator_.set_params(
                        n_jobs=CLASSIFY_THREADS)
            except AttributeError as e:
                model.named_steps.rf.set_params(n_jobs=CLASSIFY_THREADS)

            self.model = model

        self.gff_mat_colnames = ('orf_index', 'start', 'end', 'strand', 
                'partial', 'start_type', 'gc_cont', 'rbs_motif')


    def process_one_contig(self, seqname, mat, start, end, name):
        '''Process a contig

        Args:
            seqname (str): sequence name
            mat (list of list): list of list of info related to genomic
                features

        '''
        df_gff = pd.DataFrame(mat, columns=self.gff_mat_colnames)

        sel_orf1 = df_gff['start'] >= start
        sel_orf2 = df_gff['end'] <= end

        df_gff_sel = df_gff.loc[sel_orf1 & sel_orf2,:]

        #print('start:', start)
        #print('end:', end)
        #print('start:', df_gff_sel['start'].iloc[0])
        #print('end:', df_gff_sel['end'].iloc[-1])
        #print('orf_index_start:', df_gff_sel['orf_index'].iloc[0])
        #print('orf_index_end:', df_gff_sel['orf_index'].iloc[-1])
        ##if seqname == 'Caudo-linear2||rbs:common':
        #if seqname == 'Caudo-provirus||rbs:common':
        #    print(df_gff.head())
        #    print(df_gff_sel.head())

        sel = (self.df_tax_all['seqname'] == seqname)
        df_tax = self.df_tax_all.loc[sel,:]
        sel = df_tax['orf_index'].isin(set(df_gff_sel['orf_index']))
        df_tax_sel = df_tax.loc[sel,:]

        if self.d_hallmark_hmm != None:
            hmms = df_tax['hmm']
            indice = df_tax['orf_index']
            score_cutoffs = hmms.map(
                    lambda x: self.d_hallmark_hmm.get(x, 
                        [np.nan, np.inf])[1],
                    na_action='ignore',
            )
            sel_score = df_tax['score'] > score_cutoffs
            sel_index_w_hallmark = indice.loc[sel_score].tolist()
        else:
            sel_index_w_hallmark = []

        # redo prediction here
        l = get_feature(df_gff_sel, df_tax_sel, self.rbs_cat_d, 
                sel_index_w_hallmark)
        if len(l) == 0:
            # does not meet 1) >= 2 genes; 2) at least 1 full gene
            return

        # classify

        res_lis = self.model.predict_proba([l,]) # [[0.84 0.16]]
        res = res_lis[0]  #[0.84 0.16]
        pr_full = res[1]       # 0.16
        # print(self.model.classes_)  # [0, 1] 1 is viral
        seqname_ori = seqname.rsplit('||', 1)[0]
        self.fw.write('{}\t{}\n'.format(name, pr_full))


    def classify_trimmed(self):
        '''Iterate through each contigs in gff and boudary list in name2loc_d for each contig to get score
        '''
        mat = []
        last_seqname = None
        with open(self.outfile, 'w') as self.fw:
            header_lis = [ 'seqname', self.group ]
            self.fw.write('{}\n'.format('\t'.join(header_lis)))

            for lis in self.gff_gen:
                seqname = lis[0]
                seqname_ori = seqname.rsplit('||', 1)[0]
                if not seqname_ori in self.name2loc_d:
                    continue
                # process a contig
                if last_seqname != None and last_seqname != seqname:
                    #logging.info('Processing {}'.format(last_seqname))
                    last_seqname_ori = last_seqname.rsplit('||', 1)[0]
                    name_lis = self.name2loc_d[last_seqname_ori]['seqname']
                    loc_lis = self.name2loc_d[last_seqname_ori]['loc']
                    for name, _lis in zip(name_lis, loc_lis) :
                        start, end = _lis
                        #print(len(mat), name, start, end)
                        self.process_one_contig(
                                last_seqname, mat, 
                                start, end, name,
                        )

                    # reset mat and last_seqname for next iter
                    mat = []
                    last_seqname = None

                mat.append(lis[2:])
                last_seqname = seqname

            if len(mat) != 0:
                #logging.info('Processing {}'.format(last_seqname))
                last_seqname_ori = last_seqname.rsplit('||', 1)[0]
                if last_seqname_ori in self.name2loc_d:
                    name_lis = self.name2loc_d[last_seqname_ori]['seqname']
                    loc_lis = self.name2loc_d[last_seqname_ori]['loc']
                    for name, _lis in zip(name_lis, loc_lis) :
                        start, end = _lis
                        self.process_one_contig(
                                last_seqname, mat, 
                                start, end, name,
                        )


CONTEXT_SETTINGS = {'help_option_names': ['-h', '--help']}
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--hallmark', default=None, type=click.Path(), 
        help='hallmark gene list file (hallmark-gene.list)')
@click.option('--group', default='viral', 
    help=('viral group of classifier'))
@click.argument('gff', type=click.Path())
@click.argument('tax', type=click.Path())
@click.argument('rbs', type=click.Path()) 
@click.argument('model', type=click.Path())
@click.argument('full_and_part_seqfile', type=click.Path())
@click.argument('outfile', type=click.Path())
def main(gff, tax, rbs, model, full_and_part_seqfile, outfile, hallmark, group):
    '''
    this function do the following:
        apply a classifier and assign classification score to each seq in combined seqfile of full and partial seqs
    
    \f
    :param hallmark: hallmark gene list file (hallmark-gene.list)
    :param full_and_part_seqfile: full and partial viral seqfile concatenated  with start and end info in header description
    :param group: viral group for classification 
    :param gff: .gff file from prodigal
    :param tax: .tax file from extract-feature-from-hmmout.py
    :param rbs: rbs cattory file (rbs-catetory.tsv)
    :param model: classifier model file
    :param outfile: output file for boundary info
    :return: 0 if success else 1
    :rtype: int
    '''

    clf_trim_ins = clf_trim(gff, tax, rbs, model, 
            full_and_part_seqfile, outfile, hallmark_f=hallmark, group=group)
    clf_trim_ins.load_data()
    clf_trim_ins.classify_trimmed()

if __name__ == '__main__':
    main()

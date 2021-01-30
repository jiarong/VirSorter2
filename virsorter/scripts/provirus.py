#!/usr/bin/env python
import sys
import os
import warnings
from collections import OrderedDict, Counter
import logging
import io

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
        parse_gff, extract_feature_gff, get_feature, GFF_PARSER_COLS
)


DEFAULT_CONFIG = get_default_config()

GROUP_DICT = DEFAULT_CONFIG['GROUP_INFO']
GENE_OVERLAP_MIN = DEFAULT_CONFIG['GENE_OVERLAP_MIN']
CLASSIFY_THREADS = DEFAULT_CONFIG['CLASSIFY_THREADS']
MIN_FRAC_OF_MAX_SCORE = DEFAULT_CONFIG['MIN_FRAC_OF_MAX_SCORE']
MAX_RETRY_TIMES = DEFAULT_CONFIG['MAX_RETRY_TIMES']
TOTAL_FEATURE_LIST = DEFAULT_CONFIG['TOTAL_FEATURE_LIST']
SELECT_FEATURE_LIST = DEFAULT_CONFIG['SELECT_FEATURE_LIST']
DEFAULT_MIN_GENOME_SIZE = DEFAULT_CONFIG['DEFAULT_MIN_GENOME_SIZE']
END_TRIM_OFF = DEFAULT_CONFIG['END_TRIM_OFF']
PROVIRUS_CHECK_MAX_FULLSEQ_PROBA = \
        DEFAULT_CONFIG['PROVIRUS_CHECK_MAX_FULLSEQ_PROBA']


set_logger()



class provirus(object):
    '''The class to process provirus extraction and trimming full seqs
    '''
    def __init__(self, gff_f, tax_f, rbs_cat_f, model_f, outfile, ftrfile,
            hallmark_f=None, fullseq_clf_f=None, group=None, proba=0.5):
        '''Initialize object

        Args:

            rbs_cat_d (dict): mapping hmm name to gene name and hmmsearch bit score
            cutoff from rbs_cat_d
            gff_f (str): gff file from prodigal
            tax_f (str): info loaded from .tax file (output of extract-feature-from-hmmout.py
            rbs_cat_f (str): file mapping rbs sites and categories (rbs directory of database)
            model_f (str): classifier model file
            outfile (str): output file
            ftr_f (str): feature file (genomic, taxonomic, and hallmark gene cnt)
            hallmark_f (str): hallmark gene list file
            fullseq_clf_f: file with probability of each contig being viral with different classifers
            group: viral group name
            proba: probabilty cutoff of being viral

        '''
        self.gff_f = gff_f
        self.tax_f = tax_f
        self.rbs_cat_f = rbs_cat_f
        self.model_f = model_f
        self.outfile = outfile
        self.ftr_f = ftrfile
        self.hallmark_f = hallmark_f
        self.fullseq_clf_f = fullseq_clf_f
        self.group = group
        self.proba = proba

    def load_data(self):
        self.gff_gen = parse_gff(self.gff_f)

        self.rbs_cat_d = load_rbs_category(self.rbs_cat_f)

        if self.hallmark_f != None:
            self.d_hallmark_hmm = parse_hallmark_hmm(self.hallmark_f)
        else:
            self.d_hallmark_hmm = None

        if self.fullseq_clf_f != None:
            df = pd.read_csv(self.fullseq_clf_f, sep='\t', header=0)
            # force seqname col to be str dtype in case seqname are 
            #   numbers only
            df = df.astype({'seqname': 'str'})
            decoy_lis = [i for i in df.columns if i.startswith('decoy')]
            df = df.drop(decoy_lis, axis=1)
            self.df_fullseq_clf = df

        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=DeprecationWarning)
            warnings.filterwarnings('ignore', category=FutureWarning)
            model = joblib.load(self.model_f)
            try:
                # pipe on grid search; legacy
                model.named_steps.gs.best_estimator_.set_params(
                        n_jobs=CLASSIFY_THREADS)
            except AttributeError as e:
                # grid search on pipe; new
                # steps = [('scaler', MinMaxScaler()), ('rf', clf_to_train)]
                model.named_steps.rf.set_params(n_jobs=CLASSIFY_THREADS)

            self.model = model

        self.gff_mat_colnames = GFF_PARSER_COLS[2:]
        self.hallmark_ftr_ind = SELECT_FEATURE_LIST.index('hallmark')
        self.arc_ind = SELECT_FEATURE_LIST.index('arc')
        self.bac_ind = SELECT_FEATURE_LIST.index('bac')
        self.euk_ind = SELECT_FEATURE_LIST.index('euk')
        self.vir_ind = SELECT_FEATURE_LIST.index('vir')
        self.mix_ind = SELECT_FEATURE_LIST.index('mix')
        self.unaligned_ind = SELECT_FEATURE_LIST.index('unaligned')


    def locate_ends(self, df_gff, df_tax, sel_index_w_hallmark, combs):
        '''Location ends of viral candidate region

        The longest regions with score within MIN_FRAC_OF_MAX_SCORE (95%) of
        best core is chosen

        Args:
            df_gff (pandas DataFrame): info from parse_gff
            df_tax (pandas DataFrame): info loaded from .tax file (output of
                extract-feature-from-hmmout.py
            sel_index_w_hallmark (pandas Series): index of hallmark genes in
                df_tax
            combs (list of list): list of start and end pairs

        Returns:
            list: list including final_start, final_end, final_ind_start,
                final_ind_end, pr, pr_max, ftr_lis

        '''
        indice = df_gff['orf_index']
        recs = []
        for start, end in combs:
            if end == 0:
                _end = None
            else:
                _end = end
            df_gff_sel = df_gff.iloc[start:_end,:]
            if len(df_gff_sel) <= 1:
                continue
            else:
                sel = df_tax['orf_index'].isin(
                        set(df_gff_sel['orf_index']))
                df_tax_sel = df_tax.loc[sel,:]
                l = get_feature(
                        df_gff_sel, df_tax_sel, 
                        self.rbs_cat_d, sel_index_w_hallmark
                )
                if len(l) == 0:
                    continue
                ind_start = indice.iloc[start]
                # end is exclusive for slicing;
                # ind_end is inclusive to match orf index
                ind_end = indice.iloc[end-1]
                recs.append((start, end, ind_start, ind_end, l))

        if len(recs) == 0:
            return []
        new_starts, new_ends, ind_starts, ind_ends, X = zip(*recs)
        sizes = np.array(ind_ends) - np.array(ind_starts) + 1

        res_lis = self.model.predict_proba(X)
        prs = list(zip(*res_lis))[1]
        _df = pd.DataFrame({'size':sizes, 'pr':prs})
        df_ori = _df
        pr_max = max(prs)
        cutoff = max(pr_max * MIN_FRAC_OF_MAX_SCORE, self.proba)
        #cutoff = pr_max * MIN_FRAC_OF_MAX_SCORE
        _df = _df.loc[_df['pr'] >= cutoff,:]
        if len(_df) == 0:
            return []
            
        # original index is still kept after selection
        ind = _df['size'].idxmax() 
        final_ind_start = ind_starts[ind]
        final_ind_end = ind_ends[ind]
        final_start = new_starts[ind]
        final_end = new_ends[ind]
        pr = _df['pr'].loc[ind]
        ftr_lis = X[ind]

        return (final_start, final_end, 
                final_ind_start, final_ind_end, pr, pr_max, ftr_lis)

    def trim_ends(self, df_gff, df_tax, sel_index_w_hallmark, seqname, 
            prox_pr, prox_pr_max, partial, 
            full_orf_index_start, full_orf_index_end,
            full_bp_start, full_bp_end, pr_full, 
            arc, bac, euk, vir, mix, unaligned, hallmark_cnt):
        '''Trim the ends of viral candidate regions

        The function writes the boudary info to output file, including 
            seqname, final_ind_start, final_ind_end,
            final_bp_start, final_bp_end, pr, pr_max,
            prox_ind_start, prox_ind_end,
            prox_bp_start, prox_bp_end,
            prox_pr, prox_pr_max,
            partial, full_orf_index_start, full_orf_index_end,
            full_bp_start, full_bp_end, pr_full, 
            arc, bac, euk, vir, mix, unaligned, hallmark_cnt

        '''
        # trim ends 5 or 10% of total genes
        indice = df_gff['orf_index']
        prox_ind_start = indice.iloc[0]
        prox_ind_end = indice.iloc[-1]
        prox_bp_start = df_gff['start'].iloc[0]
        prox_bp_end = df_gff['end'].iloc[-1]

        size = int(0.1* len(df_gff))
        #if size < 5:
        #    size = 5
        if END_TRIM_OFF:  #"true" in yaml ==> True (boolean) in python
            size = 0
        if size == 0:
            combs = [(0, 0),]
        else:
            starts = np.array(range(size))
            ends = -1*np.array(range(size))
            combs = [ (x, y) for x in starts for y in ends ]
        if size <= 10:
            # combinations
            lis = self.locate_ends(df_gff, df_tax, 
                    sel_index_w_hallmark, combs)
            try:
                _, _, final_ind_start, final_ind_end, pr, pr_max, ftr_lis = lis
            except ValueError as e:
                #assert len(lis) == 0:
                return
        else:
            # try start and end separately
            start_combs = [ (x, -(size-1)) for x in starts ]
            lis = self.locate_ends(df_gff, df_tax, 
                    sel_index_w_hallmark, start_combs)
            try:
                final_start, _, final_ind_start, _, _, _, _ = lis
            except ValueError as e:
                # assert len(lis) == 0
                return

            end_combs = [ (final_start, y) for y in ends ]
            lis = self.locate_ends(df_gff, df_tax, 
                            sel_index_w_hallmark, end_combs)
            _, final_end, _, final_ind_end, pr, pr_max, ftr_lis = lis

        final_bp_start = \
                df_gff['start'].loc[indice == final_ind_start].iloc[0]
        final_bp_end = df_gff['end'].loc[indice == final_ind_end].iloc[0]

        mes = ('{}\t{}\t{}\t{}\t{}\t{:.3f}\t{:.3f}'
                '\t{}\t{}\t{}\t{}\t{:.3f}\t{:.3f}'
                '\t{}\t{}\t{}\t{}\t{}\t{:.3f}'
                '\t{:.1f}\t{:.1f}\t{:.1f}\t{:.1f}\t{:.1f}\t{:.1f}\t{}\n')
        self.fw.write(
                mes.format(
                    seqname, final_ind_start, final_ind_end,
                    final_bp_start, final_bp_end, pr, pr_max,
                    prox_ind_start, prox_ind_end,
                    prox_bp_start, prox_bp_end,
                    prox_pr, prox_pr_max,
                    partial, full_orf_index_start, full_orf_index_end,
                    full_bp_start, full_bp_end, pr_full, 
                    arc, bac, euk, vir, mix, unaligned, int(hallmark_cnt)
                )
        )

        #print(ftr_lis)
        ftr_lis = [
                str(i) if isinstance(i, int) else '{:.3f}'.format(i) \
                        for i in ftr_lis
        ]
        mes = '{}\t{}\t{}\t{}\n'
        self.fw_ftr.write(
                mes.format(seqname, final_ind_start, 
                    final_ind_end, '\t'.join(ftr_lis))
        )

    def process_one_contig(self, seqname, seqlen, mat):
        '''Process a contig

        Args:
            seqname (str): sequence name
            mat (list of list): list of list of info related to genomic
                features

        '''
        try:
            MIN_GENOME_SIZE = GROUP_DICT[self.group]['MIN_GENOME_SIZE']
        except KeyError:
            MIN_GENOME_SIZE = DEFAULT_MIN_GENOME_SIZE

        # if proba table to get proba of whole seq is provided
        if self.fullseq_clf_f != None:
            seqname_ori = seqname.rsplit('||',1)[0]
            seqname_st = set(self.df_fullseq_clf['seqname'].unique())
            if not seqname_ori in seqname_st:
                return
            _df = self.df_fullseq_clf
            _df.index = _df['seqname']
            _df = _df.drop(['seqname'], axis=1)
            pr_full = _df.at[seqname_ori, self.group]
            pr_full_max_groups = max(_df.loc[seqname_ori,:])
            #if pr_full < self.proba and pr_full_max_groups >= self.proba:
            if pr_full < PROVIRUS_CHECK_MAX_FULLSEQ_PROBA and \
                    pr_full_max_groups >= max(self.proba,
                            PROVIRUS_CHECK_MAX_FULLSEQ_PROBA):
                # skip when another group is significant with full seq,
                #  so do not need go through sliding window, save computation
                #  in merge-provirus-from-groups.py step, provirus is selected
                #  based on longer length, 
                #  and full seq is always longer than partial, thus preferred
                return


            # prefilter out short contigs with low proba
            #  - too short for provirus extraction
            #  - unlikely to have proba > cutoff after trimming ends
            if seqlen < min(3000, MIN_GENOME_SIZE) and pr_full < self.proba:
                return

        df_gff = pd.DataFrame(mat, columns=self.gff_mat_colnames)
        df_tax = df_tax_per_config(self.tax_f, seqname)
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

        l = get_feature(df_gff, df_tax, self.rbs_cat_d, 
                sel_index_w_hallmark)
        if len(l) == 0:
            # does not meet 1) >= 2 genes; 2) at least 1 full gene
            return

        hallmark_cnt = l[self.hallmark_ftr_ind]
        arc = l[self.arc_ind]
        bac = l[self.bac_ind]
        euk = l[self.euk_ind]
        vir = l[self.vir_ind]
        mix = l[self.mix_ind]
        unaligned = l[self.unaligned_ind]

        if self.fullseq_clf_f == None:
            # redo prediction here
            # classify
            res_lis = self.model.predict_proba([l,]) # [[0.84 0.16]]
            res = res_lis[0]  #[0.84 0.16]
            pr_full = res[1]       # 0.16
            # print(self.model.classes_)  # [0, 1] 1 is viral

            # prefilter out short contigs with low proba
            #  - too short for provirus extraction
            #  - unlikely to have proba > cutoff after trimming ends
            if seqlen < min(3000, MIN_GENOME_SIZE) and pr_full < self.proba:
                return

        starts = df_gff['start']
        ends = df_gff['end']
        full_orf_index_start = df_gff['orf_index'].iloc[0]
        full_orf_index_end = df_gff['orf_index'].iloc[-1]
        full_bp_start = df_gff['start'].iloc[0]
        full_bp_end = df_gff['end'].iloc[-1]

        # fullseq
        #if pr_full >= self.proba:
        if pr_full >= PROVIRUS_CHECK_MAX_FULLSEQ_PROBA:
            # those with pr_full < self.proba will be filtered 
            #  through trim_ends()
            partial = 0
            if pr_full < self.proba:
                # those with pr_full < self.proba could be filtered 
                #  through trim_ends(), but filter here to save time 
                #  also keep consistent with provirus that only segments
                #  with > self.proba gets passed to trim_ends()
                return

            self.trim_ends(df_gff, df_tax, 
                    sel_index_w_hallmark, seqname, 
                    np.nan, np.nan,
                    partial, full_orf_index_start, full_orf_index_end, 
                    full_bp_start, full_bp_end, pr_full, 
                    arc, bac, euk, vir, mix, unaligned, hallmark_cnt)
        # partial
        else:
            # sliding windows
            ind_end = len(ends.loc[ends < MIN_GENOME_SIZE]) + 1
            if ind_end >= len(df_gff):
                # too short to be provirus; just trim_ends() as fullseq
                partial = 0
                self.trim_ends(df_gff, df_tax, 
                        sel_index_w_hallmark, seqname, 
                        np.nan, np.nan,
                        partial, full_orf_index_start, full_orf_index_end, 
                        full_bp_start, full_bp_end, pr_full, 
                        arc, bac, euk, vir, mix, unaligned, hallmark_cnt)
                return

            # first window index 0 -> ind
            ind_start = 0
            trigger = False
            trigger_cnt = 0
            provirus_cnt = 0
            retry_cnt = 0
            pr_max = -10000
            pr_last_valid = -10000
            while (ind_end < len(df_gff)):
                df_gff_sel = df_gff.iloc[ind_start:ind_end+1]
                sel = df_tax['orf_index'].isin(
                        set(df_gff_sel['orf_index']))
                df_tax_sel = df_tax.loc[sel,:]
                l = get_feature(
                        df_gff_sel, df_tax_sel, 
                        self.rbs_cat_d, sel_index_w_hallmark
                )

                hallmark_cnt = 0
                if len(l) == 0:
                    # do not meet 1) > 2 genes and 2) 1 full gene
                    #  this can happen if there are large genes
                    pr = 0
                else:
                    res_lis = self.model.predict_proba([l,])
                    res = res_lis[0]
                    pr = res[1]
                    hallmark_cnt = l[self.hallmark_ftr_ind]
                    arc = l[self.arc_ind]
                    bac = l[self.bac_ind]
                    euk = l[self.euk_ind]
                    vir = l[self.vir_ind]
                    mix = l[self.mix_ind]
                    unaligned = l[self.unaligned_ind]

                if trigger == False and pr < self.proba:
                    ind_start += 1
                    ind_end += 1
                elif trigger == False and pr >= self.proba:
                    trigger = True
                    trigger_cnt += 1
                    ind_end += 1
                    pr_last_valid = pr
                    if pr > pr_max:
                        pr_max = pr
                #elif trigger == True and pr >= self.proba:
                elif trigger == True and pr >= pr_max * 0.95:
                    ind_end += 1
                    retry_cnt = 0
                    pr_last_valid = pr
                    if pr > pr_max:
                        pr_max = pr
                #elif trigger == True and pr < self.proba:
                elif trigger == True and pr < pr_max * 0.95:
                    if retry_cnt <= MAX_RETRY_TIMES:
                        retry_cnt += 1
                        ind_end += 1
                    else:
                        # walk back retry_cnt extension
                        # and walk back 1 to previous ind_end; 
                        ind_end = ind_end - 1 - retry_cnt
                        df_gff_sel = df_gff.iloc[ind_start:ind_end+1]
                        sel = df_tax['orf_index'].isin(
                                set(df_gff_sel['orf_index']))
                        df_tax_sel = df_tax.loc[sel,:]
                        partial = 1
                        # require hallmark gene for provirus
                        if hallmark_cnt > 0:
                            provirus_cnt += 1
                            self.trim_ends(df_gff_sel, df_tax_sel, 
                                    sel_index_w_hallmark, seqname,
                                    pr_last_valid, pr_max,
                                    partial, 
                                    full_orf_index_start, full_orf_index_end, 
                                    full_bp_start, full_bp_end, 
                                    pr_full, 
                                    arc, bac, euk, vir, mix, unaligned,
                                    hallmark_cnt)


                        # set up for next provirus in the same contig
                        ind_start = ind_end + 1
                        _start = starts.iloc[ind_start]
                        ind_end = len(
                                ends.loc[ends < _start + MIN_GENOME_SIZE]) + 1
                        trigger = False
                        retry_cnt = 0
                        pr_max = -10000
                        pr_last_valid = -10000

            if trigger == True:
                # last segment is still provirus
                # walk back 1 and retry_cnt extension
                ind_end = ind_end - 1 - retry_cnt
                df_gff_sel = df_gff.iloc[ind_start:ind_end+1]
                sel = df_tax['orf_index'].isin(
                        set(df_gff_sel['orf_index']))
                df_tax_sel = df_tax.loc[sel,:]
                partial = 1
                # require hallmark gene for provirus
                if hallmark_cnt > 0:
                    provirus_cnt += 1
                    self.trim_ends(df_gff_sel, df_tax_sel,
                            sel_index_w_hallmark, seqname,
                            pr_last_valid, pr_max,
                            partial,
                            full_orf_index_start, full_orf_index_end, 
                            full_bp_start, full_bp_end, pr_full, 
                            arc, bac, euk, vir, mix, unaligned,
                            hallmark_cnt)

            if provirus_cnt == 0 and pr_full >= self.proba:
                # trigger_cnt > 0 <==> pr_full >= self.proba
                # add this condition to accommodate 
                # 1) provirus proba cutoff high, 
                #      triggering requiring hallmark
                # 2) seq has no hallmark 
                # 3) seq pr_full > proba cutoff;
                # 
                partial = 0
                self.trim_ends(df_gff, df_tax, 
                    sel_index_w_hallmark, seqname, 
                    np.nan, np.nan,
                    partial, full_orf_index_start, 
                    full_orf_index_end, full_bp_start, 
                    full_bp_end, pr_full, arc, bac, euk, 
                    vir, mix, unaligned, hallmark_cnt)

    def find_boundary(self):
        '''Iterate through each contigs in gff and get provirus boundary info
        '''
        mat = []
        last_seqname = None
        with open(self.outfile, 'w') as self.fw, \
                open(self.ftr_f, 'w') as self.fw_ftr:
            header_lis = [
                    'seqname', 'trim_orf_index_start', 'trim_orf_index_end',
                    'trim_bp_start', 'trim_bp_end', 'trim_pr', 'trim_pr_max',
                    'prox_orf_index_start', 'prox_orf_index_end',
                    'prox_bp_start', 'prox_bp_end',
                    'prox_pr', 'prox_pr_max',
                    'partial', 'full_orf_index_start', 'full_orf_index_end',
                    'full_bp_start', 'full_bp_end', 'pr_full', 
                    'arc', 'bac', 'euk', 'vir', 'mix', 'unaligned', 
                    'hallmark_cnt'
            ]
            self.fw.write('{}\n'.format('\t'.join(header_lis)))

            header_lis = [
                    'seqname', 'trim_orf_index_start', 'trim_orf_index_end', 
            ] + SELECT_FEATURE_LIST
            self.fw_ftr.write('{}\n'.format('\t'.join(header_lis)))

            for lis in self.gff_gen:
                seqname = lis[0]
                seqlen = lis[1]
                # process a contig
                if last_seqname != None and last_seqname != seqname:
                    #logging.info('Processing {}'.format(last_seqname))
                    self.process_one_contig(last_seqname, last_seqlen, mat)
                    # reset mat and last_seqname for next iter
                    mat = []
                    last_seqname = None
                    last_seqlen = None

                # do not need first two items: seqname, seqlen
                mat.append(lis[2:])
                last_seqname = seqname
                last_seqlen = seqlen

            if len(mat) != 0:
                #logging.info('Processing {}'.format(last_seqname))
                self.process_one_contig(last_seqname, last_seqlen, mat)


CONTEXT_SETTINGS = {'help_option_names': ['-h', '--help']}
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--hallmark', default=None, type=click.Path(), 
        help='hallmark gene list file (hallmark-gene.list)')
@click.option('--fullseq-clf', default=None, type=click.Path(), 
        help=('classification prob table of all groups for '
            'whole sequence (viral-fullseq-proba.tsv)'))
@click.option('--group', default=None, 
    help=('viral group of classifier, only needed when '
        '--fullseq-clf is provided'))
@click.option('--proba', default=0.5, type=click.FloatRange(min=0, max=1),
    help=('minimal proba required for provirus'))
@click.argument('gff', type=click.Path())
@click.argument('tax', type=click.Path())
@click.argument('rbs', type=click.Path()) 
@click.argument('model', type=click.Path())
@click.argument('outfile', type=click.Path())
@click.argument('ftrfile', type=click.Path()) 
def main(gff, tax, rbs, model, outfile, ftrfile, hallmark, fullseq_clf, group, proba):
    '''
    this function do the following:
        1) if a contig is classified as viral, it trims the ends to find a segment with optimal classifier score
        2) if a contig is classified as nonviral, it first do a slide window to find a segment with p > 0.5; and then extend till p fall below 0.5; then trim ends as in step 1)
    
    \f
    :param hallmark: hallmark gene list file (hallmark-gene.list)
    :param fullseq_clf: classification prob table for full seq (viral-fullseq-proba.tsv)
    :param group: viral group for classification (only needed when --fullseq-clf is provided)
    :param proba: minimal proba required for provirus
    :param gff: .gff file from prodigal
    :param tax: .tax file from extract-feature-from-hmmout.py
    :param rbs: rbs cattory file (rbs-catetory.tsv)
    :param model: classifier model file
    :param outfile: output file for boundary info
    :param ftr: file with features of selected bourdaries
    :return: 0 if success else 1
    :rtype: int
    '''
    if fullseq_clf != None and  group == None:
        mes = '***--group is required when --fullseq_clf is provided\n'
        sys.stderr.write(mes)
        sys.exit(1)

    provirus_ins = provirus(gff, tax, rbs, model, outfile, ftrfile,
        hallmark_f=hallmark, fullseq_clf_f=fullseq_clf, 
        group=group, proba=proba)
    provirus_ins.load_data()
    provirus_ins.find_boundary()

if __name__ == '__main__':
    main()

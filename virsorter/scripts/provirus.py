#!/usr/bin/env python
import sys
import os
import warnings
from collections import OrderedDict, Counter
import logging

import joblib
import click
import sklearn
import numpy as np
import pandas as pd

from virsorter.config import DEFAULT_CONFIG, set_logger

GROUP_DICT = DEFAULT_CONFIG['GROUP_INFO']
GENE_OVERLAP_MIN = DEFAULT_CONFIG['GENE_OVERLAP_MIN']
CLASSIFY_THREADS = DEFAULT_CONFIG['CLASSIFY_THREADS']
TAXON_LIST = DEFAULT_CONFIG['TAXON_LIST']
MIN_FRAC_OF_MAX_SCORE = DEFAULT_CONFIG['MIN_FRAC_OF_MAX_SCORE']
MAX_RETRY_TIMES = DEFAULT_CONFIG['MAX_RETRY_TIMES']
TOTAL_FEATURE_LIST = DEFAULT_CONFIG['TOTAL_FEATURE_LIST']
SELECT_FEATURE_LIST = DEFAULT_CONFIG['SELECT_FEATURE_LIST']

set_logger()

def load_rbs_category(f):
    '''load into rbs site and rbs catetory mapping into a dict
    '''
    d = {}
    with open(f) as fp:
        for line in fp:
            if line.startswith('#'):
                continue
            line = line.rstrip()
            site, cat = line.split('\t', 1)
            d[site] = cat

    return d

def parse_hallmark_hmm(hallmark_f):
    '''map hmm name to gene name and hmmsearch bit score cutoff
    '''
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

    return d_hallmark_hmm

def parse_gff(gff):
    ''' parsing gff file from prodigal to info related genomic features.
    Args:
        gff (str): gff file from prodial output

    Returns:
        list of tuple: list of genomic features from each contig including
            seqname, seqlen, orf_index, start, end, strand, partial,
            start_type, gc_cont, rbs_motif

    '''
    with open(gff) as fp:
        for line in fp:
            if line.startswith('# Sequence Data:'):
                """
                The fields in this header are as follows:

                seqnum: An ordinal ID for this sequence, beginning at 1.
                seqlen: Number of bases in the sequence.
                seqhdr: The entire FASTA header line.
                version: Version of Prodigal used to analyze this sequence.
                run_type: "Ab initio" for normal mode, "Anonymous" for anonymous mode.
                model (Anonymous mode only): Information about the preset training file used to analyze the sequence.
                gc_cont: % GC content of the sequence.
                transl_table: The genetic code used to analyze the sequence.
                uses_sd: Set to 1 if Prodigal used its default RBS finder, 0 if it scanned for other motifs.
                """
                line = line.rstrip()
                s = line.split('Sequence Data:')[1].strip()
                try:
                    seq_data = OrderedDict(
                            i.strip().split('=') for i in s.split(';')
                    )
                except ValueError as e:
                    seq_data = OrderedDict()
                    # in case there are ";" in the seq description
                    for i in s.split(';'):
                        if not '=' in i:
                            continue
                        m, n = i.strip().split('=',1)
                        seq_data[m] = n

                seqlen = int(seq_data['seqlen'])
                seqhdr = seq_data['seqhdr'].strip('"')
                seqname = seqhdr.split(None, 1)[0]
                continue

            elif line.startswith('#'):
                continue

            line = line.rstrip()
            """
            [0]Seqname\t[1]Prodigal_v2.6.3\t[2]CDS\t[3]start\t[4]end\t[5]score\t[6]strand\t[7]frame\t[8]semicolon-delimited-string\n   
            --------------
            fields in semicolon-delimited-string:
            [0]ID: A unique identifier for each gene, consisting of the ordinal ID of the sequence and an ordinal ID of that gene within the sequence (separated by an underscore). For example, "4_1023" indicates the 1023rd gene in the 4th sequence in the file.
            [1]partial: An indicator of if a gene runs off the edge of a sequence or into a gap. A "0" indicates the gene has a true boundary (a start or a stop), whereas a "1" indicates the gene is "unfinished" at that edge (i.e. a partial gene). For example, "01" means a gene is partial at the right boundary, "11" indicates both edges are incomplete, and "00" indicates a complete gene with a start and stop codon.
            [2]start_type: The sequence of the start codon (usually ATG, GTG, or TTG). If the gene has no start codon, this field will be labeled "Edge".
            [3]rbs_motif: The RBS motif found by Prodigal (e.g. "AGGA" or "GGA", etc.)
            [4]rbs_spacer: The number of bases between the start codon and the observed motif.
            [5]gc_cont: The GC content of the gene sequence.
            [6]conf: A confidence score for this gene, representing the probability that this gene is real, i.e. 78.3% means Prodigal believes that gene is real 78.3% of the time and a false positive 21.7% of the time.
            [7]score: The total score for this gene.
            [8]cscore: The hexamer coding portion of the score, i.e. how much this gene looks like a true protein.
            [9]sscore: A score for the translation initiation site for this gene; it is the sum of the following three fields.
            [10]rscore: A score for the RBS motif of this gene.
            [11]uscore: A score for the sequence surrounding the start codon.
            [12]tscore: A score for the start codon type (ATG vs. GTG vs. TTG vs. Nonstandard).
            """
            items = line.split('\t')

            start = int(items[3])
            end = int(items[4])
            strand = items[6]

            last = items[8] # ; delimited string
            sub_items = OrderedDict(
                    i.strip().split('=') for i in last.rstrip(';').split(';')
            )
            ind = sub_items['ID']
            #"4_1023" is the 1023rd gene in the 4th sequence in the file
            orf_index = int(ind.split('_')[1])
            partial =  sub_items['partial']
            if partial != '00':
                partial = 1
            else:
                partial = 0

            rbs_motif = sub_items['rbs_motif']
            #print(rbs_motif)
            start_type = sub_items['start_type']
            # gc_cont=0.482; convert to perc
            gc_cont = 100*float(sub_items['gc_cont']) 

            yield (seqname, seqlen, orf_index, start, end, strand, partial, 
                    start_type, gc_cont, rbs_motif)

def extract_feature_gff(df, rbs_cat_d):
    '''Extract genomic features from info from parse_gff

    Args:
        df (pandas DataFrame): info from parse_gff
        rbs_cat_d (dict): mapping hmm name to gene name and hmmsearch bit score
            cutoff from rbs_cat_d

    Returns:
        list: list of genomic features or empty if less than two orfs

    '''
    # require at least two genes and one needs to be full
    partial_lis = df['partial']
    cnt_partial_per_contig = sum(partial_lis)
    cnt_full_per_contig = len(partial_lis) - cnt_partial_per_contig
    if cnt_full_per_contig >= 1 \
            and cnt_partial_per_contig + cnt_full_per_contig >= 2:
        ends = df['end']
        starts = df['start']
        gene_sizes = abs(ends - starts)
        seqlen = abs(ends.iloc[-1] - starts.iloc[0])
        # select only full gene for mean size
        sel = (partial_lis == 0)
        gene_sizes = gene_sizes.loc[sel]
        size = np.mean(gene_sizes)

        # pd.Series arithmatic operation automatically match index
        #   need to convert to np.array first
        gene_spacings = np.array(starts.iloc[1:]) - np.array(ends.iloc[:-1])
        spacing = np.mean(gene_spacings)

        total_cnt = cnt_partial_per_contig + cnt_full_per_contig
        # gene cnt per 1kb
        density = 1000*total_cnt/seqlen

        # gene overlap perc
        strands = df['strand']
        cnt_strand_switch = 0
        cnt_gene_overlap = 0
        for i in range(len(strands)):
            if i == 0:
                continue
            if strands.iloc[i] != strands.iloc[i-1]:
                cnt_strand_switch += 1
            else:
                spa = gene_spacings[i-1]   # has been converted to np.array
                if spa < -1*GENE_OVERLAP_MIN:
                    cnt_gene_overlap += 1

        gene_overlap_perc = 100.0*cnt_gene_overlap/(total_cnt-1)
        strand_switch_perc = 100.0*cnt_strand_switch/(total_cnt-1)

        rbs_list = df['rbs_motif']
        rbs_cnt_d = Counter(rbs_list)
        cat_cnt_d = {}
        rbs_cat_list = sorted(set(rbs_cat_d.values()))
        seen_cat_st = set()
        for key in rbs_cnt_d:
            cnt = rbs_cnt_d[key]
            try:
                cat = rbs_cat_d[key]
            except KeyError as e:
                #mes = '{} is not in {}\n'
                #if not key in seen_cat_st:
                #    logging.info(
                #            mes.format(key, 
                #                os.path.basename(rbs_cat_f))
                #    )
                seen_cat_st.add(key)
                cat = 'Other'
            current = cat_cnt_d.get(cat, 0)
            cat_cnt_d[cat] = current + cnt

        cat_cnt_list = [cat_cnt_d.get(cat, 0) 
                            for cat in rbs_cat_list]
        _sum = sum(cat_cnt_list)
        cat_perc_list = [100.0*i/_sum for i in cat_cnt_list]

        start_codon_list = df['start_type']
        start_codon_d = Counter(start_codon_list)
        start_codon_num = len(start_codon_list)
        atg_perc = 100.0*start_codon_d['ATG']/start_codon_num
        gtg_perc = 100.0*start_codon_d['GTG']/start_codon_num
        ttg_perc = 100.0*start_codon_d['TTG']/start_codon_num

        gc_list = df['gc_cont']
        gc_sd = np.std(gc_list)
        gc_mean = np.mean(gc_list)
        
        _l = [seqlen, size, gene_overlap_perc, density, 
                strand_switch_perc,
                atg_perc, gtg_perc, ttg_perc, gc_mean, gc_sd, ]

        _l.extend(cat_perc_list)
        return _l
    else:
        return []

def get_feature(df_gff, df_tax, rbs_cat_d, sel_index_w_hallmark):
    '''Get all features (genomic and taxonomic and hallmark gene cnt)

    Args:
        df_gff (pandas DataFrame): info from parse_gff
        df_tax (pandas DataFrame): info loaded from .tax file (output of
            extract-feature-from-hmmout.py
        rbs_cat_d (dict): mapping hmm name to gene name and hmmsearch bit score
            cutoff from rbs_cat_d
        sel_index_w_hallmark: index of hallmark genes in df_tax

    '''
    df_gff_sel = df_gff
    df_tax_sel = df_tax
    l = extract_feature_gff(df_gff_sel, rbs_cat_d)
    if len(l) == 0:
        # does not meet 1) >= 1 full gene 2) >= 2 total genes
        return l

    if len(df_tax_sel) == 0:
        # all unaligned
        l_tax = (0, 0, 0, 0, 0, 100, 0)
    else:
        if len(sel_index_w_hallmark) != 0:
            indice = df_tax_sel['orf_index']
            hallmark_st = set(indice).intersection(
                    set(sel_index_w_hallmark))
            hallmark_cnt = len(hallmark_st)
        else:
            hallmark_cnt = 0

        tax_lis = df_tax_sel['tax']
        aligned = len(tax_lis)
        _d = Counter(tax_lis)
        total = len(df_gff_sel)
        l_tax = [ 100.0*_d.get(key, 0)/total for key in TAXON_LIST ]
        unaligned = total - aligned
        unaligned_perc = 100.0*unaligned/total
        l_tax.append(unaligned_perc)
        l_tax.append(hallmark_cnt)

    l.extend(l_tax)
    ser = pd.Series(l, index=TOTAL_FEATURE_LIST)
    l = ser.loc[ser.index.isin(SELECT_FEATURE_LIST)].tolist()
    return l

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

        df_tax_all = pd.read_csv(self.tax_f, sep='\t', header=None, 
                names=['orfname', 'tax', 'hmm', 'score'])
        seqnames, indice = zip(
                *[orfname.rsplit('_', 1) for orfname in df_tax_all['orfname']]
        )
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

        if self.fullseq_clf_f != None:
            df = pd.read_csv(self.fullseq_clf_f, sep='\t', header=0)
            decoy_lis = [i for i in df.columns if i.startswith('decoy')]
            df = df.drop(decoy_lis, axis=1)
            self.df_fullseq_clf = df

        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=DeprecationWarning)
            warnings.filterwarnings('ignore', category=FutureWarning)
            model = joblib.load(self.model_f)
            model.named_steps.gs.best_estimator_.set_params(
                    n_jobs=CLASSIFY_THREADS)
            self.model = model

        self.gff_mat_colnames = ('orf_index', 'start', 'end', 'strand', 
                'partial', 'start_type', 'gc_cont', 'rbs_motif')


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
                df_gff_sel = df_gff.iloc[start:,:]
            else:
                df_gff_sel = df_gff.iloc[start:end,:]
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
        #print(_df)
        pr_max = max(prs)
        cutoff = pr_max * MIN_FRAC_OF_MAX_SCORE
        _df = _df.loc[_df['pr'] >= cutoff,:]
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
            full_bp_start, full_bp_end, pr_full):
        '''Trim the ends of viral candidate regions

        The function writes the boudary info to output file, including 
            seqname, final_ind_start, final_ind_end,
            final_bp_start, final_bp_end, pr, pr_max,
            prox_ind_start, prox_ind_end,
            prox_bp_start, prox_bp_end,
            prox_pr, prox_pr_max,
            partial, full_orf_index_start, full_orf_index_end,
            full_bp_start, full_bp_end, pr_full

        '''
        # trim ends 5 or 10% of total genes
        indice = df_gff['orf_index']
        prox_ind_start = indice.iloc[0]
        prox_ind_end = indice.iloc[-1]
        prox_bp_start = df_gff['start'].iloc[0]
        prox_bp_end = df_gff['end'].iloc[-1]

        size = int(0.1* len(df_gff))
        if size < 5:
            size = 5
        starts = np.array(range(size))
        ends = -1*np.array(range(size))
        combs = [ (x, y) for x in starts for y in ends ]
        if size <= 10:
            # combinations
            lis = self.locate_ends(df_gff, df_tax, 
                    sel_index_w_hallmark, combs)
            if len(lis) == 0:
                return
            _, _, final_ind_start, final_ind_end, pr, pr_max, ftr_lis = \
                    lis
        else:
            # try start and end separately
            start_combs = [ (x, -(size-1)) for x in starts ]
            lis = self.locate_ends(df_gff, df_tax, 
                    sel_index_w_hallmark, start_combs)
            final_start, _, final_ind_start, _, _, _, _ = lis

            end_combs = [ (final_start, y) for y in ends ]
            lis = self.locate_ends(df_gff, df_tax, 
                            sel_index_w_hallmark, end_combs)
            _, final_end, _, final_ind_end, pr, pr_max, ftr_lis = lis

        final_bp_start = \
                df_gff['start'].loc[indice == final_ind_start].iloc[0]
        final_bp_end = df_gff['end'].loc[indice == final_ind_end].iloc[0]
        mes = ('{}\t{}\t{}\t{}\t{}\t{:.3f}\t{:.3f}'
                '\t{}\t{}\t{}\t{}\t{:.3f}\t{:.3f}'
                '\t{}\t{}\t{}\t{}\t{}\t{:.3f}\n')
        self.fw.write(
                mes.format(
                    seqname, final_ind_start, final_ind_end,
                    final_bp_start, final_bp_end, pr, pr_max,
                    prox_ind_start, prox_ind_end,
                    prox_bp_start, prox_bp_end,
                    prox_pr, prox_pr_max,
                    partial, full_orf_index_start, full_orf_index_end,
                    full_bp_start, full_bp_end, pr_full
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

    def process_one_contig(self, seqname, mat):
        '''Process a contig

        Args:
            seqname (str): sequence name
            mat (list of list): list of list of info related to genomic
                features

        '''
        df_gff = pd.DataFrame(mat, columns=self.gff_mat_colnames)
        sel = (self.df_tax_all['seqname'] == seqname)
        df_tax = self.df_tax_all.loc[sel,:]
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
            if pr_full < self.proba and pr_full_max_groups >= self.proba:
                # skip when another group is significant with full seq,
                #  so do not need go through sliding window, save computation
                #  in merge-provirus-from-groups.py step, provirus is selected
                #  based on longer length, 
                #  and full seq is always longer than partial, thus preferred
                return

        else:
            # redo prediction here
            l = get_feature(df_gff, df_tax, self.rbs_cat_d, 
                    sel_index_w_hallmark)
            if len(l) == 0:
                # does not meet 1) >= 2 genes; 2) at least 1 full gene
                return

            # classify

            res_lis = self.model.predict_proba([l,]) # [[0.84 0.16]]
            res = res_lis[0]  #[0.84 0.16]
            pr_full = res[1]       # 0.16
            # print(self.model.classes_)  # [0, 1] 1 is viral

        starts = df_gff['start']
        ends = df_gff['end']
        full_orf_index_start = df_gff['orf_index'].iloc[0]
        full_orf_index_end = df_gff['orf_index'].iloc[-1]
        full_bp_start = df_gff['start'].iloc[0]
        full_bp_end = df_gff['end'].iloc[-1]

        if pr_full >= self.proba:
            partial = 0
            self.trim_ends(df_gff, df_tax, 
                    sel_index_w_hallmark, seqname, 
                    np.nan, np.nan,
                    partial, full_orf_index_start, full_orf_index_end, 
                    full_bp_start, full_bp_end, pr_full)
        else:
            # sliding windows
            MIN_GENOME_SIZE = GROUP_DICT[self.group]['MIN_GENOME_SIZE']
            ind_end = len(ends.loc[ends < MIN_GENOME_SIZE]) + 1
            if ind_end >= len(df_gff):
                # too short to be provirus
                return
            # first window index 0 -> ind
            ind_start = 0
            trigger = False
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

                hallmark_exist= False
                if len(l) == 0:
                    # do not meet 1) > 2 genes and 2) 1 full gene
                    #  this can happen if there are large genes
                    pr = 0
                else:
                    res_lis = self.model.predict_proba([l,])
                    res = res_lis[0]
                    pr = res[1]
                    if l[-1] > 0:
                        hallmark_exist = True

                if trigger == False and pr < self.proba:
                    ind_start += 1
                    ind_end += 1
                elif trigger == False and pr >= self.proba:
                    trigger = True
                    ind_end += 1
                    pr_last_valid = pr
                    if pr > pr_max:
                        pr_max = pr
                elif trigger == True and pr >= self.proba:
                    ind_end += 1
                    retry_cnt = 0
                    pr_last_valid = pr
                    if pr > pr_max:
                        pr_max = pr
                elif trigger == True and pr < self.proba:
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
                        if hallmark_exist:
                            self.trim_ends(df_gff_sel, df_tax_sel, 
                                    sel_index_w_hallmark, seqname,
                                    pr_last_valid, pr_max,
                                    partial, 
                                    full_orf_index_start, full_orf_index_end, 
                                    full_bp_start, full_bp_end, pr_full)

                        # set up for next provirus in the same contig
                        ind_start = ind_end + 1
                        _start = starts.iloc[ind_start]
                        ind_end = len(
                                ends.loc[ends < _start + MIN_GENOME_SIZE]) + 1
                        trigger = False
                        retry_cnt = 0
                        pr_max = -10000
                        pr_last_valid = -10000

            if trigger == False:
                # not provirus
                return

            # walk back 1 and retry_cnt extension
            ind_end = ind_end - 1 - retry_cnt
            df_gff_sel = df_gff.iloc[ind_start:ind_end+1]
            sel = df_tax['orf_index'].isin(
                    set(df_gff_sel['orf_index']))
            df_tax_sel = df_tax.loc[sel,:]
            partial = 1
            # require hallmark gene for provirus
            if hallmark_exist:
                self.trim_ends(df_gff_sel, df_tax_sel,
                        sel_index_w_hallmark, seqname,
                        pr_last_valid, pr_max,
                        partial,
                        full_orf_index_start, full_orf_index_end, 
                        full_bp_start, full_bp_end, pr_full)


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
            ]
            self.fw.write('{}\n'.format('\t'.join(header_lis)))

            header_lis = [
                    'seqname', 'trim_orf_index_start', 'trim_orf_index_end', 
            ] + SELECT_FEATURE_LIST
            self.fw_ftr.write('{}\n'.format('\t'.join(header_lis)))

            for lis in self.gff_gen:
                seqname = lis[0]
                # process a contig
                if last_seqname != None and last_seqname != seqname:
                    #logging.info('Processing {}'.format(last_seqname))
                    self.process_one_contig(last_seqname, mat)
                    # reset mat and last_seqname for next iter
                    mat = []
                    last_seqname = None

                mat.append(lis[2:])
                last_seqname = seqname

            if len(mat) != 0:
                #logging.info('Processing {}'.format(last_seqname))
                self.process_one_contig(last_seqname, mat)


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

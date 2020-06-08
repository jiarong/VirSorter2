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
from virsorter.config import DEFAULT_CONFIG, set_logger

GENE_OVERLAP_MIN = DEFAULT_CONFIG['GENE_OVERLAP_MIN']
CLASSIFY_THREADS = DEFAULT_CONFIG['CLASSIFY_THREADS']
TAXON_LIST = DEFAULT_CONFIG['TAXON_LIST']
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

        self.name2loc_d = {}
        with screed.open(self.full_and_part_seqfile) as sp:
            for rec in sp:
                header = rec.name
                name, desc = header.split(None ,1)
                seqname_ori = name.rsplit('||', 1)[0]
                _d = dict(i.split(':') for i in desc.split('||'))

                self.name2loc_d.setdefault(seqname_ori, {})
                self.name2loc_d[seqname_ori].setdefault('seqname', [])
                self.name2loc_d[seqname_ori].setdefault('loc', [])

                self.name2loc_d[seqname_ori]['seqname'].append(name)
                self.name2loc_d[seqname_ori]['loc'].append(
                        (int(_d.get('start')), int(_d.get('end')))
                )

        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=DeprecationWarning)
            warnings.filterwarnings('ignore', category=FutureWarning)
            model = joblib.load(self.model_f)
            model.named_steps.gs.best_estimator_.set_params(
                    n_jobs=CLASSIFY_THREADS)
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

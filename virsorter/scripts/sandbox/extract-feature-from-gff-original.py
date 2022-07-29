#!/usr/bin/env python
# Usage: python <thisfile> <rbs-cat.tsv> <in.gff> <feature-table.tsv>

import sys
import os
import screed
from collections import OrderedDict, Counter
import logging

script_dir = os.path.dirname(os.path.abspath(__file__))
snakefile_dir = os.path.dirname(script_dir)
pkg_dir = os.path.dirname(snakefile_dir)
sys.path.append(pkg_dir)
from virsorter.config import set_logger

GENE_OVERLAP_MIN=5

set_logger()

def mean(data):
    '''
    Return the sample arithmetic mean of data.
    '''
    n = len(data)
    if n < 1:
        raise ValueError('mean requires at least one data point')
    return 1.0*sum(data)/n

def _ss(data):
    '''Return sum of square deviations of sequence data.
    '''
    c = mean(data)
    ss = sum((x-c)**2 for x in data)
    return ss

def stddev(data, ddof=0):
    '''Calculates the standard deviation

    Calculates population by default; specify ddof=1 to compute the sample
    standard deviation.
    '''
    n = len(data)
    if n < 2:
        raise ValueError('variance requires at least two data points')
    ss = _ss(data)
    pvar = 1.0*ss/(n-ddof)
    return pvar**0.5

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

def main():
    '''Extract genomic features from gff file

    Example:
        python extract-feature-from-gff.py <rbs-category.tsv> \
                <in.gff> <feature-table.tsv>
        
        <rbs-category.tsv>: rsb site and category mapping file in rbs directory
            in database
        <in.gff>: gff file
        <feature-table.tsv>: tab delimited file with genomic features
        
    '''
    if len(sys.argv) != 4:
        mes = (
                '*** Usage: python {} '
                '<rbs-category.tsv> ' 
                '<in.gff> <feature-table.tsv>\n'
                )
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    rbs_cat_f = sys.argv[1]
    in_gff = sys.argv[2]
    out_gff = sys.argv[3]
    if out_gff == '-':
        out_gff = '/dev/stdout'

    rbs_cat_d = load_rbs_category(rbs_cat_f)
    rbs_cat_list = sorted(set(rbs_cat_d.values()))
    seen_cat_st = set()
    with open(in_gff) as fp, open(out_gff, 'w') as fw:
        _l = ['contig_size', 'gene_size', 'gene_overlap_perc', 
                'density', 'strand_switch_perc',
                'atg_perc', 'gtg_perc', 'ttg_perc', 'gc_mean', 'gc_sd']
        _l += ['rbs_{}'.format(i) for i in rbs_cat_list]
        _s = '\t'.join(['{}'.format(i) for i in _l])
        fw.write('{}\t{}\n'.format('seqname', _s))

        # initialize
        cnt1 = 0    # cnt partial genes in whole file
        cnt2 = 0    # cnt totol contigs
        cnt3 = 0    # cnt contigs with < 2 genes in whole file

        seqname = None
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
                cnt2 += 1
                if seqname == None:
                    pass
                # require at least two genes and one needs to be full
                elif cnt_full_per_contig >= 1 \
                        and cnt_partial_per_contig + cnt_full_per_contig >= 2:
                    size = mean(gene_sizes)

                    spacing = mean(gene_spacings)

                    total_cnt = cnt_partial_per_contig + cnt_full_per_contig
                    # gene cnt per 1kb
                    density = 1000*total_cnt/seqlen

                    # gene overlap perc
                    gene_overlap_perc = 100.0*cnt_gene_overlap/(total_cnt-1)

                    strand_switch_perc = 100.0*cnt_strand_switch/(total_cnt-1)

                    rbs_cnt_d = Counter(rbs_list)
                    cat_cnt_d = {}
                    for key in rbs_cnt_d:
                        cnt = rbs_cnt_d[key]
                        try:
                            cat = rbs_cat_d[key]
                        except KeyError as e:
                            mes = '{} is not in {}\n'
                            if not key in seen_cat_st:
                                logging.info(
                                        mes.format(key, 
                                            os.path.basename(rbs_cat_f))
                                )
                            seen_cat_st.add(key)
                            cat = 'Other'
                        current = cat_cnt_d.get(cat, 0)
                        cat_cnt_d[cat] = current + cnt

                    cat_cnt_list = [cat_cnt_d.get(cat, 0) 
                                        for cat in rbs_cat_list]
                    _sum = sum(cat_cnt_list)
                    cat_perc_list = [100.0*i/_sum for i in cat_cnt_list]

                    start_codon_d = Counter(start_codon_list)
                    start_codon_num = len(start_codon_list)
                    atg_perc = 100.0*start_codon_d['ATG']/start_codon_num
                    gtg_perc = 100.0*start_codon_d['GTG']/start_codon_num
                    ttg_perc = 100.0*start_codon_d['TTG']/start_codon_num

                    gc_sd = stddev(gc_list)
                    gc_mean = mean(gc_list)

                    print('====>bbbb')
                    print(cat_cnt_d)
                    print(len(gc_list))
                    print(start_codon_num)
                    print(total_cnt)
                    
                    _l = [seqlen, size, gene_overlap_perc, density, 
                            strand_switch_perc,
                            atg_perc, gtg_perc, ttg_perc, gc_mean, gc_sd, ]
                    _l += cat_perc_list
                    _s = '\t'.join(['{:.1f}'.format(i) for i in _l])
                    fw.write('{}\t{}\n'.format(seqname, _s))
                else:
                    mes = ('Contig removed due to < 2 total genes or '
                            'has no full gene: {}')
                    logging.info(mes.format(seqname))
                    cnt3 += 1

                # reintialize after output
                gene_sizes = []
                gene_spacings = []
                rbs_list = []
                start_codon_list = []
                gc_list = []
                cnt_partial_per_contig = 0
                cnt_full_per_contig = 0
                cnt_strand_switch = 0
                cnt_gene_overlap = 0
                first_gene = True

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

            if not first_gene:
                prev_end = end
                prev_strand = strand
            else:
                pass

            start = int(items[3])
            end = int(items[4])
            strand = items[6]
            if not first_gene:
                gene_spacings.append(start - prev_end)
                if strand != prev_strand:
                    cnt_strand_switch += 1
                else:
                    if start - prev_end <= -1*GENE_OVERLAP_MIN:
                        # gene overlap cnt on the same stran
                        # more common in virus
                        cnt_gene_overlap += 1
            else:
                first_gene = False

            last = items[8] # ; delimited string
            sub_items = OrderedDict(
                    i.strip().split('=') for i in last.rstrip(';').split(';')
            )
            partial =  sub_items['partial']
            rbs_motif = sub_items['rbs_motif']
            #print(rbs_motif)
            start_type = sub_items['start_type']
            # gc_cont=0.482; convert to perc
            gc_cont = 100*float(sub_items['gc_cont']) 

            gc_list.append(gc_cont)

            if partial != '00':
                # discard partial gene for average gene size
                cnt1 += 1
                cnt_partial_per_contig += 1
            else:
                cnt_full_per_contig += 1
                # only add full gene to average gene size calc
                gene_sizes.append(abs(end - start))
                rbs_list.append(rbs_motif)
                start_codon_list.append(start_type)

        # save results right after loop
        if seqname == None:
            pass
        # require at least two genes and one needs to be full
        elif cnt_full_per_contig >= 1 \
                and cnt_partial_per_contig + cnt_full_per_contig >= 2:
            size = mean(gene_sizes)
            spacing = mean(gene_spacings)
            total_cnt = cnt_partial_per_contig + cnt_full_per_contig
            # gene cnt per 1kb
            density = 1000*total_cnt/seqlen

            # gene overlap perc
            gene_overlap_perc = 100.0*cnt_gene_overlap/(total_cnt-1)
            strand_switch_perc = 100.0*cnt_strand_switch/(total_cnt-1)

            rbs_cnt_d = Counter(rbs_list)
            cat_cnt_d = {}
            for key in rbs_cnt_d:
                cnt = rbs_cnt_d[key]
                try:
                    cat = rbs_cat_d[key]
                except KeyError as e:
                    mes = '{} is not in {}'
                    if not key in seen_cat_st:
                        logging.info(
                                mes.format(key, os.path.basename(rbs_cat_f))
                        )
                    seen_cat_st.add(key)
                    cat = 'Other'
                current = cat_cnt_d.get(cat, 0)
                cat_cnt_d[cat] = current + cnt

            cat_cnt_list = [cat_cnt_d.get(cat, 0) 
                                for cat in rbs_cat_list]
            _sum = sum(cat_cnt_list)
            cat_perc_list = [100.0*i/_sum for i in cat_cnt_list]

            start_codon_d = Counter(start_codon_list)
            start_codon_num = len(start_codon_list)
            atg_perc = 100.0*start_codon_d['ATG']/start_codon_num
            gtg_perc = 100.0*start_codon_d['GTG']/start_codon_num
            ttg_perc = 100.0*start_codon_d['TTG']/start_codon_num

            gc_sd = stddev(gc_list)
            gc_mean = mean(gc_list)

            _l = [seqlen, size, gene_overlap_perc, density, strand_switch_perc,
                    atg_perc, gtg_perc, ttg_perc, gc_mean, gc_sd]
            _l += cat_perc_list
            _s = '\t'.join(['{:.1f}'.format(i) for i in _l])
            fw.write('{}\t{}\n'.format(seqname, _s))
        else:
            mes = ('Contig removed due to < 2 total genes '
                    'or has no full gene: {}')
            logging.info(mes.format(seqname))
            cnt3 += 1

    mes = '# of partial genes removed from all contigs: {}'
    logging.info(mes.format(cnt1))
    mes = '# of total contigs: {}'
    logging.info(mes.format(cnt2))
    mes = ('# of contigs removed due to having < 2 genes or has ' 
            'no full gene: {}')
    logging.info(mes.format(cnt3))

if __name__ == '__main__':
    main()

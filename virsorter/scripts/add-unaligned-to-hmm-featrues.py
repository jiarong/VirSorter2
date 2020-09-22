#!/usr/bin/env python
# by gjr

import sys
import os
import screed
from collections import Counter
import click

script_dir = os.path.dirname(os.path.abspath(__file__))
snakefile_dir = os.path.dirname(script_dir)
pkg_dir = os.path.dirname(snakefile_dir)
sys.path.append(pkg_dir)
from virsorter.utils import TAXON_LIST

#TAXON_LIST = ['arc', 'bac', 'euk', 'vir', 'mixed']

@click.command()
@click.option('--hallmark',default=None,
        type=click.Path(),
        help='hallmark gene list file (hallmark-gene.list)')
@click.argument('faa', type=click.Path(exists=True)) 
@click.argument('tax', type=click.Path(exists=True))
def main(faa, tax, hallmark):
    '''
    This script combines info from FAA (prodigal-out.faa) and 
    TAX (prodigal-out.hmm.tax) and HALLMARK (hallmark-gene.list, if exists)
    and produce a taxonomic feature file

    \f
    :param faa: faa file output from prodigal
    :param tax: tax file output from extract-feature-from-hmmout.py
    :param hallmark: hallmark-gene.list for a viral group
    :return: 0 if success, 1 if fail
    :rtype: int
    '''
    faa = faa
    hmm = tax
    hallmark_f = hallmark
    if hmm == '-':
        hmm = '/dev/stdin'

    sys.stdout.write(
            '#seqname\t{}\tunaligned\thallmark\n'.format('\t'.join(TAXON_LIST)))

    d_cnt = {}
    with screed.open(faa) as fp:
        for rec in fp:
            name = rec.name.split(None, 1)[0]
            contig = name.rsplit('_', 1)[0]
            cnt = d_cnt.get(contig, 0)
            d_cnt[contig] = cnt + 1

    d_hallmark_hmm = {}
    if hallmark != None:
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

    d_tax = {}
    d_hallmark_cnt = {}
    with open(hmm) as fp2:
        for line in fp2:
            if line.startswith('#'):
                continue
            line = line.rstrip()
            lis = line.split('\t')
            name = lis[0]
            tax = lis[1]
            hmm = lis[2]
            try:
                score = float(lis[3])
            except ValueError as e:
                score = -10000000

            contig = name.rsplit('_', 1)[0]

            lis2 = d_tax.get(contig, [])
            d_tax[contig] = lis2 + [tax,]

            cur_cnt = d_hallmark_cnt.setdefault(contig, 0)
            if hmm in d_hallmark_hmm:
                if score >= d_hallmark_hmm[hmm][1]:
                    d_hallmark_cnt[contig] = cur_cnt + 1

    # contig with no gene aligned are not in d_tax
    for contig in d_cnt:
        if contig in d_tax:
            tax_lis = d_tax[contig]
            aligned = len(tax_lis)
            _d = Counter(tax_lis)

            total = d_cnt[contig]
            _l = [ 100.0*_d.get(key, 0)/total for key in TAXON_LIST ]
            unaligned = total - aligned
            unaligned_perc = 100.0*unaligned/total
            _l.append(unaligned_perc)

            _s = '\t'.join(['{:.1f}'.format(i) for i in _l])
            sys.stdout.write(
                    '{}\t{}\t{}\n'.format(contig, _s, d_hallmark_cnt[contig])
                    )
        else:
            mes = '{}\t0.0\t0.0\t0.0\t0.0\t0.0\t100.0\t0\n'
            sys.stdout.write(mes.format(contig))


if __name__ == '__main__':
    main()

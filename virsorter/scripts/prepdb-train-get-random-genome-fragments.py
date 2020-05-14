#!/usr/bin/env python
# Usage: python <thisfile> 
#               <id2taxid2taxon.tab> 
#               <refeq-genomes.dir> 
#               "taxon_level"
#               number-of-contigs-per-taxon-level
#
# This script is used when all refseq genome has been download locally

import sys
import os
import glob
import screed
import random

TAXON_LEVELS = ['superkingdom','phylum','class', 'order','family',
                    'genus','species']

N_CONTIGS_PER_TAXON = 5   #random contig among genomes in a genus
N_STARTS_PER_CONTIG = 2   #random starting postion on each contig
N_LENGTHS_PER_START = 2   #random length at each starting position

MIN_SIZE = 1000

def main():
    if len(sys.argv) != 4:
        mes = (
                '*** Usage: python {} <id2taxid2taxon.tab> '
                '<refeq-genomes.dir> "taxon_level"\n'
        )
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    tab_f = sys.argv[1]
    genome_dir = sys.argv[2].rstrip('/')
    taxon_level = sys.argv[3]

    assert taxon_level in TAXON_LEVELS, \
            '*** {} must be one of {}'.format(
                    taxon_level, 
                    ','.join(TAXON_LEVELS)
            )
    ind = TAXON_LEVELS.index(taxon_level)
    ind += 1   # assembly id and taxon id cols are before taxons

    d = {}
    with open(tab_f) as fp:
        for line in fp:
            if line.startswith('#'):
                continue
            line = line.rstrip()
            li = line.split()
            acc = li[0]
            genus = li[ind]
            st = d.setdefault(genus, set())
            f = glob.glob('{}/{}*_genomic.fna.gz'.format(genome_dir, acc))
            assert len(f) == 1, \
                '*** {} has > 1 genomes matched:\n {}'.format(
                        acc, '\n'.join(f)
                )
            f = f[0]
            d[genus].add(f)

    cnt_written = 0
    length_li = []
    for genus in d:
        fs = d[genus]
        cnt = 0
        for f in fs:
            with screed.open(f) as fp:
                for rec in fp:
                    cnt += 1

        rand_contig_li = random.sample(range(cnt), N_CONTIGS_PER_TAXON)
        cnt2 = 0
        for f in fs:
            with screed.open(f) as fp:
                for rec in fp:
                    seq = rec.sequence
                    if cnt2 in set(rand_contig_li):
                        length = len(seq)
                        rand_start_li = []
                        while 1:
                            li = random.sample(range(length), 1)
                            i = li[0]
                            if length - i < MIN_SIZE:
                                continue
                            rand_start_li.append(i)
                            if len(rand_start_li) == N_STARTS_PER_CONTIG:
                                break

                        for i in rand_start_li:
                            rand_length_li = random.sample(
                                    range(i, length - 1000), 
                                    N_LENGTHS_PER_START
                            )
                            for j in rand_length_li:
                                sys.stdout.write(
                                        '>{}_{}_{}\n{}\n'.format(
                                            rec.name, i, j, seq[i, j]
                                        )
                                )
                                cnt_written += 1
                                length_li.append((j - i))
                    cnt2 += 1


    sys.stderr.write(
            '*** total sequence randomly sampled: {}'.format(cnt_written)
    )
    mini = min(length_li)
    maxi = max(length_li)
    mean = sum(length_li)*1.0/len(length_li)
    sys.stderr.write((
            '*** Among randomly sampled fragments:\n'
            '***   Max: {}\n'
            '***   Min: {}\n'
            '***    Mean: {:.1f}').format(maxi, mini, mean)
    )

if __name__ == '__main__':
    main()

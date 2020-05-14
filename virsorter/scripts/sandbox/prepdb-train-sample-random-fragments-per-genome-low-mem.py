#!/usr/bin/env python

import sys
import os
import screed
import random

N_FRAGMENTS_PER_GENOME = 4

MIN_SIZE = 1000

def main():
    if len(sys.argv) < 3:
        mes = (
                '*** Usage: python {} '
                '<nonviral-sampled-accession-with-cnt.list> '
                '<genome1.fna.gz> <genome2.fna.gz>..\n'
        )
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    list_f = sys.argv[1]
    genome_fs = sys.argv[2:]

    d = {}
    with open(list_f) as fp:
        for line in fp:
            if line.startswith('#'):
                continue
            line = line.rstrip()
            li = line.split()
            acc = li[0]
            cnt = int(li[1])
            d[acc] = cnt

    cnt_written = 0
    length_li = []
    for f in genome_fs:
        bname = os.path.basename(f)
        #GCF_900097015.1_PADLG01_genomic.fna.gz
        acc = '_'.join(bname.split('_')[:2])
        if acc not in d:
            continue

        contig_length_li = []
        with screed.open(f) as fp:
            for n, rec in enumerate(fp):
                length = len(rec.sequence)
                contig_length_li.append(length)

        rand_contig_li = []
        rand_start_li = []
        rand_end_li = []
        while 1:
            # pick a contig
            ind = random.choice(range(len(contig_length_li)))
            length = contig_length_li[ind]
            if length <= MIN_SIZE:
                continue
            rand_contig_li.append(ind)
            start = random.choice(range(length - MIN_SIZE))
            rand_start_li.append(start)
            end = random.choice(range((start + MIN_SIZE), length))
            rand_end_li.append(end)
            if len(rand_start_li) == \
                    N_FRAGMENTS_PER_GENOME * d[acc]:
                break

        with screed.open(f) as fp:
            for n, rec in enumerate(fp):
                seq = rec.sequence
                if not n in set(rand_contig_li):
                    continue

                for i,j in enumerate(rand_contig_li):
                    if j != n:
                        continue
                    start = rand_start_li[i]
                    end = rand_end_li[i]

                    sys.stdout.write(
                            '>{}_{}_{}\n{}\n'.format(
                                rec.name, start, end, seq[start:end]
                            )
                    )
                    cnt_written += 1
                    length_li.append((end - start))

    assert cnt_written == sum(
            [d[acc] * N_FRAGMENTS_PER_GENOME for acc in d]
    )
    sys.stderr.write(
            '*** Total sequence randomly sampled: {}'.format(cnt_written)
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

#!/usr/bin/env python

import sys
import os
import screed
import random
from collections import Counter


MIN_SIZE = 1000
# stop if > MAX # of trial on a genme
# indicating bad genome quality e.g. < MIN_SIZE bp or too many "N" in genome
MAX_FIALED_TRIAL = 20

def main():
    '''Sample random DNA fragments from refseq genomes
    
    argv[1]: fragments per genome 
    argv[2]: tab delimited file with genome as col 1, weight as col 2 
    argv[3:]: genome file paths; space separated 

    Each file is treated as a genome; This script only works with refseq
    genome downloaded by ncbi-genomes-download due to specific file name
    conventions related to assembly accession # used in argv[2]
    '''
    if len(sys.argv) < 3:
        mes = (
                '*** Usage: python {} cnt-per-genome '
                '<nonviral-sampled-accession-with-cnt.list> '
                '<genome1.fna.gz> <genome2.fna.gz>..\n'
        )
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    N_FRAGMENTS_PER_GENOME = int(sys.argv[1])
    list_f = sys.argv[2]
    genome_fs = sys.argv[3:]


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
    expected = 0
    length_li = []
    f_cnt = 0
    for f in genome_fs:
        bname = os.path.basename(f)
        #sys.stderr.write('*** start processing {}..\n'.format(bname))
        #GCF_900097015.1_PADLG01_genomic.fna.gz
        acc = '_'.join(bname.split('_')[:2])
        if acc not in d:
            continue

        expected += d[acc] * N_FRAGMENTS_PER_GENOME
        d_seq = {}
        with screed.open(f) as fp:
            for n, rec in enumerate(fp):
                d_seq[rec.name] = rec.sequence

        rand_start_li = []
        rand_end_li = []
        cnt_per_genome = 0
        cnt_failed_trial = 0
        while 1:
            # pick a contig
            if cnt_failed_trial > MAX_FIALED_TRIAL:
                # stop if > MAX # of trial on a genme
                #    indicating bad genome quality e.g. < MIN_SIZE bp
                #    or too many "N" in genome
                mes = (
                        '*** Warning: Only {} of {} sampled from {}, '
                        'check genome quality, e.g. contig length, '
                        '"N" in sequence\n'
                )
                sys.stderr.write(
                        mes.format(
                            cnt_per_genome, 
                            N_FRAGMENTS_PER_GENOME, 
                            bname
                        )
                )
                break

            name = random.choice(list(d_seq.keys()))
            seq = d_seq[name]
            length = len(seq)
            if length <= MIN_SIZE:
                cnt_failed_trial += 1
                continue
            start = random.choice(range(length - MIN_SIZE))
            rand_start_li.append(start)
            end = random.choice(range((start + MIN_SIZE), length))
            rand_end_li.append(end)
            subseq = seq[start:end]
            d_cnt_letter = Counter(list(subseq.upper()))
            cnt_n = d_cnt_letter.get('N', 0)
            if cnt_n*1.0/len(subseq) > 0.01:
                cnt_failed_trial += 1
                continue
            seqid, desc = name.split(None, 1)
            sys.stdout.write(
                    '>{}||s:{}||e:{}||index:{}||accession:{}  {}\n{}\n'.format(
                        seqid, start, end, cnt_written, acc, desc, subseq
                    )
            )
            cnt_per_genome += 1
            cnt_written += 1
            length_li.append((end - start))

            if cnt_per_genome == N_FRAGMENTS_PER_GENOME * d[acc]:
                break

    mes = '*** Total written ({}) vs. Total expected fragments ({})\n'
    sys.stderr.write(mes.format(cnt_written, expected))

    sys.stderr.write(
            '*** Total sequence randomly sampled: {}\n'.format(cnt_written)
    )

    mini = min(length_li)
    maxi = max(length_li)
    mean = sum(length_li)*1.0/len(length_li)
    sys.stderr.write((
            '*** Among randomly sampled fragments:\n'
            '***   Max: {}\n'
            '***   Min: {}\n'
            '***   Mean: {:.1f}\n').format(maxi, mini, mean)
    )

if __name__ == '__main__':
    main()

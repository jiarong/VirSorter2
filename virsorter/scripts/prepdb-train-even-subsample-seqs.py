#!/usr/bin/env python
# Usage: python <thisfile> num <seqfile1.fa> <seqfile2.fa> ..

import sys
import os
import screed
import random

def main():
    if len(sys.argv) < 3:
        mes = (
                '*** Usage: python {} num <seqfile1> <seqfile2> ..\n'
                '***    Use "max" for num if subsample to max '
                '(the seq # of seqfile with least seq #)\n\n'
        )
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    num = sys.argv[1]
    fs = sys.argv[2:]

    length_li = []

    for f in fs:
        with screed.open(f) as fp:
            cnt = 0
            for rec in fp:
                cnt += 1

            length_li.append(cnt)

    min_size = min(length_li)
    if num.lower() == 'max':
        size = min_size
    elif int(num) > min_size:
        size = min_size
    else:
        size = int(num)

    sys.stderr.write('*** Evenly subsample to {} seqs\n'.format(size))
    for n, f in enumerate(fs):
        total = length_li[n]
        l = random.sample(range(total), size)
        st = set(l) 
        bname = os.path.basename(f)
        with screed.open(f) as fp, open('{}.sub'.format(bname), 'w') as fw:
            for m, rec in enumerate(fp):
                if m in st:
                    fw.write('>{}\n{}\n'.format(rec.name, rec.sequence))

if __name__ == '__main__':
    main()

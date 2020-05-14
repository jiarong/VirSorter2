#!/usr/bin/env python

import sys
import os
import numpy as np

def main():
    if len(sys.argv) != 2:
        mes = '*** Usage: python {} <genomesize.list>\n'
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)


    input = sys.argv[1]
    if input == '-':
        input = '/dev/stdin'

    lis = []
    with open(input) as fp:
        for line in fp:
            if line.startswith('#'):
                continue
            size = int(line.rstrip())
            lis.append(size)

    mes = 'min\t10pct\t25pct\tmean\tmedian\t75pct\t90pct\tmax\n'

    sys.stdout.write(mes)
    mes = '{:.0f}\t{:.0f}\t{:.0f}\t{:.0f}\t{:.0f}\t{:.0f}\t{:.0f}\t{:.0f}\n'
    sys.stdout.write(
            mes.format(
                min(lis), np.percentile(lis, 10), np.percentile(lis, 25), np.mean(lis), np.median(lis), np.percentile(lis, 75), np.percentile(lis, 90), max(lis))
    )

if __name__ == '__main__':
    main()

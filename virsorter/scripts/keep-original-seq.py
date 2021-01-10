#!/usr/bin/env python

import sys
import os
import screed
import logging


script_dir = os.path.dirname(os.path.abspath(__file__))
snakefile_dir = os.path.dirname(script_dir)
pkg_dir = os.path.dirname(snakefile_dir)
sys.path.append(pkg_dir)
from virsorter.config import set_logger

set_logger()

def main():
    if len(sys.argv) != 3:
        mes = ('*** Usage: python {} <trimmed-seqfile.fa> '
                                        '<original-seqfile.fa>\n')
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    clf_f = sys.argv[1]
    faa_f = sys.argv[2]

    if clf_f == '-':
        clf_f = '/dev/stdin'

    with screed.open(clf_f) as fp, screed.open(faa_f) as sp:
        d = {}
        for rec in fp:
            header = rec.name
            lis = header.split(None, 1)
            seqname = lis[0]
            if seqname.endswith('_partial'):
                mes = f'>{seqname}\n{rec.sequence}\n'
                sys.stdout.write(mes)
                continue

            seqname_ori = seqname.rsplit('||', 1)[0]
            d[seqname_ori] = seqname

        for rec in sp:
            header = rec.name
            lis = header.split(None, 1)
            seqname_ori = lis[0]
            if not seqname_ori in d:
                continue
            seqname = d[seqname_ori]
            mes = f'>{seqname}\n{rec.sequence}\n'
            sys.stdout.write(mes)

if __name__ == '__main__':
    main()

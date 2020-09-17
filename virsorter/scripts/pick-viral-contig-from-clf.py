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
    if len(sys.argv) != 4:
        mes = '*** Usage: python {} proba_cutoff <file.clf> <contig.fa>\n'
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    cutoff = float(sys.argv[1])
    clf_f = sys.argv[2]
    faa_f = sys.argv[3]

    with open(clf_f) as fp, screed.open(faa_f) as sp:
        d = {}
        for n, line in enumerate(fp):
            line = line.rstrip()
            lis = line.split('\t')
            if n == 0:
                headers = lis
                continue

            seqname = lis[0]
            proba_lis = [float(i) for i in lis[1:]]
            maxi = max(proba_lis)
            if maxi >= cutoff:
                ind = proba_lis.index(maxi)
                # ind is for lis[1:]
                group = headers[ind+1]
                d[seqname] = group, maxi

        for rec in sp:
            header = rec.name
            lis = header.split(None, 1)
            name = lis[0]
            try:
                desc = lis[1]
            except IndexError as e:
                mes = '{} in {} does not have shape (linear or circular) info'
                logging.error(mes.format(seqname, seqfile))
                sys.exit(1)
            name = name.rsplit('||', 1)[0]
            if not name in d:
                continue
            group, score = d[name]
            mes = (f'>{name}  shape:{desc}||group:{group}||score:{score}\n'
                    f'{rec.sequence}\n')
            sys.stdout.write(mes)

if __name__ == '__main__':
    main()

#!/usr/bin/env python

import sys
import os
import screed


def main():
    '''Add suffix to sequence name

    Example: 
        python add-suffix-seqname.py <config.fa> suffix

    argv[1]: contig sequence file
    argv[2]: suffix
    output to stdout

    '''
    if len(sys.argv) != 3:
        mes = ('*** Usage: python {} <contig.fa> suffix\n')
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    infile = sys.argv[1]
    suffix = sys.argv[2]
    with screed.open(infile) as sp:
        for rec in sp:
            header = rec.name
            name = header.split(None, 1)[0]
            seq = rec.sequence
            sys.stdout.write('>{}{}\n{}\n'.format(name, suffix, seq))

if __name__ == '__main__':
    main()

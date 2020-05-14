#!/usr/bin/env python

import sys
import os
import screed


def main():
    '''Get sequence names from sequence file (fasta or fastq)

    Example:
        python get-list-from-seqfile.py <seqfile>

        <seqfile>: sequence file

    '''

    if len(sys.argv) != 2:
        mes = ('*** usage: python {} <seqfile>\n')
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    infile = sys.argv[1]
    with screed.open(infile) as sp:
        for rec in sp:
            header = rec.name
            name = header.split(None, 1)[0]
            seq = rec.sequence
            sys.stdout.write('{}\n'.format(name))

if __name__ == '__main__':
    main()

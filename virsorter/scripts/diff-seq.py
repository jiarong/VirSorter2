#!/usr/bin/env python

import sys
import os
import screed


def main():
    '''Print the seqs with unique the first seqfile.

    Example:
        python diff-seq.py <seqfile1.fa> <seqfile2.fa>

    '''
    if len(sys.argv) != 3:
        mes = ('*** Usage: python {} <seqfile1.fa> <seqfile2.fa>\n')
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    infile1 = sys.argv[1]
    infile2 = sys.argv[2]
    with screed.open(infile1) as sp1, screed.open(infile2) as sp2:
        st_from_infile2 = set([rec.name.split(None, 1)[0] for rec in sp2])
        for rec in sp1:
            name = rec.name.split(None, 1)[0]
            if name in st_from_infile2:
                continue
            seq = rec.sequence
            sys.stdout.write('>{}\n{}\n'.format(rec.name, seq))

if __name__ == '__main__':
    main()

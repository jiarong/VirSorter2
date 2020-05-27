#!/usr/bin/env python

import sys
import os
import screed


def main():
    if len(sys.argv) != 2:
        mes = '*** Usage: python {} <seqfile>\n'
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    seqfile = sys.argv[1]
    with screed.open(seqfile) as sp:
        for rec in sp:
            header =  rec.name
            sys.stdout.write('{}\t{}\n'.format(header, len(rec.sequence)))

if __name__ == '__main__':
    main()


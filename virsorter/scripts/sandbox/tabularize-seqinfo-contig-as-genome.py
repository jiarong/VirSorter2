#!/usr/bin/env python


import sys
import os
import screed

GROUPS = ['caudovirales', 'NCLDV', 'RNA', 'ssDNA', 'lavidaviridae']

def main():
    if len(sys.argv) < 2:
        mes = '*** Usage: python {} <seqfile1> <seqfile2> ..\n'
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    for f in sys.argv[1:]:
        _l = f.split('/')
        group = _l[-3]
        if not group in GROUPS:
            sys.stderr.write(
                    '*** {} not found in {}; skipping..'.format(group, Groups)
            )
            continue
        source = _l[-2]
        with screed.open(f) as sp:
            for rec in sp:
                name = rec.name.split()[0]
                mes = '{}\t{}\t{}\t{}\n'.format(
                        group, source, name, len(rec.sequence)
                )
                sys.stdout.write(mes)


if __name__ == '__main__':
    main()

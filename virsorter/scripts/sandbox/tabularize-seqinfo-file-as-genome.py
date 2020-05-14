#!/usr/bin/env python


import sys
import os
import screed

GROUPS = ['caudovirales', 'NCLDV', 'RNA', 'ssDNA', 'lavidaviridae']

def main():
    if len(sys.argv) < 4:
        mes = '*** Usage: python {} group source <seqfile1> <seqfile2> ..\n'
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    group = sys.argv[1]
    source = sys.argv[2]

    for f in sys.argv[3:]:

        #_l = os.path.basename(f).split('_')
        #genome_id = '_'.join(_l[:2])
        genome_id = os.path.basename(f).split('.fna.gz')[0]
        size = sum([len(rec.sequence) for rec in screed.open(f)])
        mes = '{}\t{}\t{}\t{}\n'.format(
                group, source, genome_id, size
        )
        sys.stdout.write(mes)


if __name__ == '__main__':
    main()

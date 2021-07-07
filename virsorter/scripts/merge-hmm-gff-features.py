#!/usr/bin/env python
# by gjr

import sys
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
snakefile_dir = os.path.dirname(script_dir)
pkg_dir = os.path.dirname(snakefile_dir)
sys.path.append(pkg_dir)
from virsorter.config import get_default_config

DEFAULT_CONFIG = get_default_config()

TOTAL_FEATURE_LIST = list(DEFAULT_CONFIG['TOTAL_FEATURE_LIST'])
SELECT_FEATURE_LIST = list(DEFAULT_CONFIG['SELECT_FEATURE_LIST'])

def main():
    if len(sys.argv) != 3:
        mes = '*** python {} <file.gff.features> <file.hmm.features>\n'
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    gff = sys.argv[1]
    hmm = sys.argv[2]

    if hmm == '-':
        hmm = '/dev/stdin'
    
    d = {}
    with open(hmm) as fp:
        for line in fp:
            if line.startswith('#'):
                continue
            line = line.rstrip()
            lis = line.split('\t', 1)
            name = lis[0]
            tax = lis[1]
            d[name] = tax

    with open(gff) as fp2:
        for line in fp2:
            if line.startswith('#'):
                continue
            line = line.rstrip()
            if line.startswith('seqname'):
                sys.stdout.write(
                        'seqname\t{}\n'.format(
                            '\t'.join(TOTAL_FEATURE_LIST)
                        )
                )
                continue

            name, other = line.split('\t', 1)
            tax = d[name]
            tax_item_size = len(tax.split('\t'))
            # if tax does not have hallmark cnt, add 0
            if tax_item_size == 6:
                tax = '{}\t0'.format(tax)
            sys.stdout.write('{}\t{}\n'.format(line, tax))


if __name__ == '__main__':
    main()

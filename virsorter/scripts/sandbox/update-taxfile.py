#!/usr/bin/env python

import sys
import os

def main():
    if len(sys.argv) != 3:
        mes = '*** Usage: python {} <new-viral-hits.list> <all.pdg.hmm.tax>\n'
        sys.stdout.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    list_f = sys.argv[1]
    tax_f = sys.argv[2]

    st = set()
    with open(tax_f) as fp:
        for line in fp:
            sys.stdout.write(line)
            if line.startswith('#'):
                continue
            line = line.rstrip()
            lis = line.split()
            st.add(lis[0])

    with open(list_f) as fp1:
        for line in fp1:
            if line.startswith('#'):
                continue
            item = line.strip()
            if item in st:
                continue
            sys.stdout.write('{}\tvir\t-\t-\n'.format(item))

if __name__ == '__main__':
    main()

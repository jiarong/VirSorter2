#!/usr/bin/env python

import sys
import os
import random
import collections

TAXON_LEVELS = ['superkingdom','phylum','class', 'order','family',
                    'genus','species']

def main():
    if len(sys.argv) != 4:
        mes = (
                '*** Usage: python {} <id2taxid2taxon.tab> cnt '
                '"taxon_level"\n'
                '***   taxon_level options:\n'
                '***   {}\n'
        )
        sys.stderr.write(
                mes.format(
                    os.path.basename(sys.argv[0]), 
                    repr(TAXON_LEVELS)
                )
        )
        sys.exit(1)

    tab_f = sys.argv[1]
    n = int(sys.argv[2])
    taxon_level = sys.argv[3]

    assert taxon_level in TAXON_LEVELS, \
            '*** {} must be one of {}'.format(
                    taxon_level, 
                    ','.join(TAXON_LEVELS)
            )
    ind = TAXON_LEVELS.index(taxon_level)
    ind += 2   # assembly id and taxon id cols are before taxons

    d = {}
    with open(tab_f) as fp:
        for line in fp:
            if line.startswith('#'):
                continue
            line = line.rstrip()
            li = line.split()
            acc = li[0]
            genus = li[ind]
            st = d.setdefault(genus, set())
            d[genus].add(acc)

    for genus in d:
        accs = d[genus]
        #print('===>' + genus)
        #print('===>' + repr(accs))
        if len(accs) > n:
            rand_li = random.sample(accs, n)
        else:
            accs_li = list(accs)
            rand_li = [random.choice(accs_li) for i in range(n)]

        d_cnt = collections.Counter(rand_li)
        for acc in rand_li:
            sys.stdout.write('{}\t{}\n'.format(acc, d_cnt[acc]))

if __name__ == '__main__':
    main()

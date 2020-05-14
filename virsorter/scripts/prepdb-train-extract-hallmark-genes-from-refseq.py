#!/usr/bin/env python

import sys
import os
import screed
import re
import errno

def main():
    if len(sys.argv) != 4:
        mes = '*** python {} <refseq.faa> <gene-anno.tsv> <outdir>\n'
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    seqfile = sys.argv[1]
    annofile = sys.argv[2]
    outdir = sys.argv[3]

    d = {}
    st = set()
    with open(annofile) as fp1:
        for line in fp1:
            if line.startswith('#'):
                continue
            line = line.rstrip()
            if not line:
                continue
            items = line.split('\t')
            name = items[0]
            name = '-'.join(name.split())
            d[name] = []
            for i in items[1:]:
                for j in i.split('*'):
                    j = j.strip()
                    if j == '':
                        continue
                    st.add(j)
                    if j == 'structur':
                        st.add('structure')
                        st.add('structural')

                pat = i.replace('*', '.*') 
                d[name].append(pat)

    try:
        try:
            os.makedirs(outdir)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

        d_w = {}
        for key in d:
            d_w[key] = open('{}/{}.fa'.format(outdir, key), 'w')

        for n, rec in enumerate(screed.open(seqfile)):
            try:
                name, desc = rec.name.split(None, 1)
            except IndexError as e:
                sys.stderr.write('*** No annotation for {}\n'.format(rec.name))
                continue

            words = set(desc.split())
            ol = st.intersection(words)
            if len(ol) == 0:
                continue

            trigger = False
            for key in d:
                pats = d[key]
                for pat in pats:
                    res = re.search(pat, desc)
                    if res:
                        d_w[key].write(
                                '>{}\n{}\n'.format(rec.name, rec.sequence)
                        )
                        trigger = True
                        break

                if trigger:
                    break

    finally:
        for key in d_w:
            d_w[key].close()

if __name__ == '__main__':
    main()

#!/usr/bin/env python

import sys
import os
import screed


def main():
    '''
    pick seqs that are in list file from many seqfiles
    '''
    if len(sys.argv) < 3:
        mes = ('*** Usage: python {} <seqname.list> '
                    '<seqfile1.fa> <seqfile2.fa> ..\n')
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    list_f = sys.argv[1]
    seqfiles = sys.argv[2:]
    st = set()
    mes = (
        '*** Only take first part of seqname; e.g. '
        'NZ_JH600071.1||s:3689||e:6051||index:6682||accession:GCF_000245015.1' 
        ' ==> NZ_JH600071.1\n')
    sys.stderr.write(mes)
    with open(list_f) as fp:
        for line in fp:
            if line.startswith('#'):
                continue
            line = line.rstrip()

            name = line.split('||', 1)[0]
            st.add(name)

    for seqfile in seqfiles:
        with screed.open(seqfile) as sp:
            for rec in sp:
                name = rec.name.split(None, 1)[0]
                if not name in st:
                    continue
                seq = rec.sequence
                sys.stdout.write('>{}\n{}\n'.format(rec.name, seq))

if __name__ == '__main__':
    main()

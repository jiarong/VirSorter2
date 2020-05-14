#!/usr/bin/env python
# remove gaps - or . in aligned fasta file
# by gjr

import sys
import os
import screed

def main():
    '''
    remove gaps in seqs
    '''
    if len(sys.argv) != 3:
        mes = 'Usage: python {} <gaped.file> <nogap.file>\n'
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    infile = sys.argv[1]
    outfile = sys.argv[2]
    
    if infile == '-':
        infile = '/dev/stdin'
    if outfile == '-':
        outfile = '/dev/stdout'

    with open(outfile, 'w') as fw, screed.open(infile) as sp:
        for record in sp:
            _l = record['name'].split(None, 1)
            name = _l[0]
            try:
                desc = _l[1]
            except IndexError as e:
                desc = ''
            seq = record['sequence']

            seq1 = seq.replace('-','').replace('.','')
            lis = seq1.split()
            seq1 = ''.join(lis)
            fw.write('>{}  {}\n{}\n'.format(name, desc, seq1))

if __name__ == '__main__':
    main()

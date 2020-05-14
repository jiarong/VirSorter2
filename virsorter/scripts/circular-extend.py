#!/usr/bin/env python

import sys
import os
import screed


def main():
    ''' Extend contigs by stitching 5' frag to 3'

    Example:
        python circular-extend.py <circular.fa> <circular.extend.fa>
    
    argv[1]: circular sequence file where sequences are shown as linear
        fasta, and the shared subtring between 5' and 3' is clipped at
        one end.  
    argv[2]: file with two duplicated sequences stitching one's 5'
        end stitched to the other's 3' end. Thus the gene where the
        circular genome break can be recovered.

    '''
    if len(sys.argv) != 3:
        mes = ('*** Usage: python {} <circular.fa> <circular.extend.fa>')
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    infile = sys.argv[1]
    outfile = sys.argv[2]
    with open(outfile, 'w') as fw:
        for rec in screed.open(infile):
            seq = rec.sequence
            fw.write('>{}\n{}{}\n'.format(rec.name, seq, seq))

if __name__ == '__main__':
    main()

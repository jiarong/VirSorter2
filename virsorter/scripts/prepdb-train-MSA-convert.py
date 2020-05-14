#!/usr/bin/env python
# parse tabular output from HMMER, do identity filter and MSA convert
# by gjr; Jan 25, 12

# Parse tabular output from HMMER, do length filter (30 bp minimum)
# Convert .sto to .afa

import sys
import os
import screed

N = 10 # e-value cutoff
M = None # bit score cutoff

def main():
    #read seqs into a dict, not memory efficient
    if len(sys.argv) != 3:
        mes = 'python {} <file.sto> <file.afa>"\n'
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    infile = sys.argv[1]
    outfile = sys.argv[2]

    if infile == '-':
        infile = '/dev/stdin'
    if outfile == '-':
        outfile = '/dev/stdout'

    with open(infile) as fp, open(outfile, 'w') as fw: 
        d_seq = {}
        for line in fp:
            if line.startswith('#'):
                continue
            if line.startswith('//'):
                #end of file
                continue
            line = line.rstrip()
            if not line:
                continue
            # assume no space in names
            name, subseq = line.split()

            if name in d_seq:
                curseq = d_seq[name]
                # combined parts of interleaved MSA
                d_seq[name] = curseq + subseq
            else:
                d_seq[name] = subseq

        for key in d_seq:
            seq = d_seq[key]
            fw.write('>{}\n{}\n'.format(key, seq))

if __name__ == '__main__':
    main()


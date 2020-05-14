#!/usr/bin/env python
# subsample genes from each contig from prodigal faa output
# gene name has the format: Name.1 Name.2 Name.3

import sys
import os
import screed
import random

def main():
    '''Subsample n sequences

    Example::

        $ python {} n <in.faa>
        
        n: number of subsampled seqs
        <in.faa>: input sequence file 
        
    '''
    if len(sys.argv) != 3:
        mes = '*** python {} n <in.faa>\n'
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    subsample_n = int(sys.argv[1])
    infile = sys.argv[2]

    if infile == '-':
        infile = '/dev/stdin'

    random.seed(99)

    d = {}
    with screed.open(infile) as fp:
        for n, rec in enumerate(fp):
            name = rec.name.split(None, 1)[0]
            contig = name.rsplit('_', 1)[0]
            d.setdefault(contig, [])
            d[contig].append(n)

    d_subsample = {}
    for contig in d:
        l = d[contig]
        if len(l) <= subsample_n:
            d_subsample[contig] = l
            continue
        subsample_l = random.sample(l, subsample_n)
        d_subsample[contig] = subsample_l

    l_flatten = [ i for l in d_subsample.values() for i in l ]
    st = set(l_flatten)
    
    with screed.open(infile) as fp:
        for n, rec in enumerate(fp):
            if n not in st:
                continue
            mes = '>{}\n{}\n'
            sys.stdout.write(mes.format(rec.name, rec.sequence))

if __name__ == '__main__':
    main()

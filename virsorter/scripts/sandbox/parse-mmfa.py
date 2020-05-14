#!/usr/bin/env python
# Usage: python <thisfile> <mmfa> <clustered.dir> <unclustered.fa>

import os
import sys
import screed
from collections import Counter

def get_clu_from_mmfa(f):
    with open(f) as fp:
        pre = None
        rep = None
        d = {}
        for n, line in enumerate(fp):
            line = line.rstrip()
            if not line:
                continue
            if not line.startswith('>'):
                continue
            name = line.split(None, 1)[0] # take str before 1st space
            name = name[1:] # remove >
            if pre == name:
                rep = pre
                continue
            if rep:
                d[pre] = rep
            pre = name

        if rep:
            d[pre] = rep

    return d

def main():
    if len(sys.argv) != 4:
        mes = '*** Usage: python {} <mmfa> <clustered.dir> <unclustered.fa>\n'
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    mmfa = sys.argv[1]
    clu_outdir = sys.argv[2]
    unclu_f = sys.argv[3]

    os.makedirs(clu_outdir)

    # d: dict with seqname as key, representative seqname as val
    d = get_clu_from_mmfa(mmfa) 
    d_cnt = Counter(d.values())

    # make a dict with repseq as key, fw as val
    # might be io intensive
    d_fw = {}
    for n, key in enumerate(d_cnt):
        if d_cnt[key] < 2:
            continue
        d_fw[key] = open('{}/{}.fa'.format(clu_outdir, key), 'w')

    try:
        with screed.open(mmfa) as fp, open(unclu_f, 'w') as fw:
            for rec in fp:
                if len(rec.sequence) == 0:
                    continue
                name = rec.name
                rep = d[name] 
                if rep in d_fw:
                    d_fw[rep].write('>{}\n{}\n'.format(name, rec.sequence))
                else:
                    fw.write('>{}\n{}\n'.format(name, rec.sequence))
        
    finally:
        for key in d_fw:
            d_fw[key].close()
    
if __name__ == '__main__':
    main()

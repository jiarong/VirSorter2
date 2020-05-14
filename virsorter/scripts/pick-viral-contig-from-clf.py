#!/usr/bin/env python

import sys
import os
import screed


def main():
    if len(sys.argv) != 4:
        mes = '*** Usage: python {} proba_cutoff <file.clf> <contig.fa>\n'
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    cutoff = float(sys.argv[1])
    clf_f = sys.argv[2]
    faa_f = sys.argv[3]

    with open(clf_f) as fp, screed.open(faa_f) as sp:
        d = {}
        for n, line in enumerate(fp):
            line = line.rstrip()
            lis = line.split('\t')
            if n == 0:
                headers = lis
                continue

            seqname = lis[0]
            proba_lis = [float(i) for i in lis[1:]]
            maxi = max(proba_lis)
            if maxi >= cutoff:
                ind = proba_lis.index(maxi)
                # ind is for lis[1:]
                group = headers[ind+1]
                d[seqname] = group

        for rec in sp:
            header = rec.name
            name = header.split(None, 1)[0]
            name = name.rsplit('||', 1)[0]
            if not name in d:
                continue
            sys.stdout.write('>{}  group:{}\n{}\n'.format(name, d[name], rec.sequence))

if __name__ == '__main__':
    main()

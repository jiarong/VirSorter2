#!/usr/bin/env python

import sys
import os
import screed


def main():
    if len(sys.argv) != 6:
        mes = '*** Usage: python {} proba_cutoff <file.clf> "A,B" "<A.pdg.faa>,<B.pdg.faa" "<A.pdg.vir.faa>,<B.pdg.vir.faa>"\n'
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    cutoff = float(sys.argv[1])
    clf_f = sys.argv[2]
    tags = [i.strip() for i in sys.argv[3].split(',')]
    faa_fs = [i.strip() for i in sys.argv[4].split(',')]
    viral_faa_fs = [i.strip() for i in sys.argv[5].split(',')]

    assert len(tags) == len(faa_fs) == len(viral_faa_fs)

    with open(clf_f) as fp:
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

    for tag, faa_f, viral_faa_f in zip(tags, faa_fs, viral_faa_fs):
        with screed.open(faa_f) as sp, open(viral_faa_f, 'w') as fw:
            for rec in sp:
                header = rec.name
                name = header.split(None, 1)[0]
                # remove "_i" from prodigal
                contig_name_w_suffix = name.rsplit('_', 1)[0]
                # remove "||rbs:common" or "||rbs:NCLDV"
                contig_name = name.rsplit('||', 1)[0]
                if not contig_name in d:
                    continue
                if d[contig_name] != tag:
                    continue
                fw.write('>{}  {}\n{}\n'.format(name, tag, rec.sequence))

if __name__ == '__main__':
    main()

#!/usr/bin/env python

import sys
import os


def open_gff_by_contig(f):
    '''Get gff records by contig.

    Args:
        f (str): gff file from prodigal.

    Returns:
        generator: generator of gff records by contig.
    
    '''
    lines = []
    first = True
    with open(f) as fp:
        for line in fp:
            if line.startswith('# Sequence Data:'):
                if not first:
                    yield ''.join(lines)
                    lines = []
                else:
                    first = False

            lines.append(line)

        yield ''.join(lines)


def main():
    '''Split gff by num contig records per split

    Example:
        python  <in.gff> <outdir> num

    '''
    if len(sys.argv) != 4:
        mes = '*** Usage: python {} <in.gff> <outdir> num\n'
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    gff_f = sys.argv[1]
    outdir = sys.argv[2]
    n = int(sys.argv[3])

    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)

    try:
        fw_lis = []
        split_idx = 0
        bname = os.path.basename(gff_f)
        fw = open('{}/{}.{}.split'.format(outdir, bname, split_idx), 'w')
        fw_lis.append(fw)
        gen = open_gff_by_contig(gff_f)
        for i, rec in enumerate(gen):
            if i < (split_idx + 1) * n:
                fw.write(rec)
            else:
                split_idx += 1
                fw = open(
                        '{}/{}.{}.split'.format(outdir, bname, split_idx),
                        'w',
                )
                fw_lis.append(fw)
                fw.write(rec)

    finally:
        _l = [fw.close() for fw in fw_lis]


if __name__ == '__main__':
    main()
                
            

            
                    

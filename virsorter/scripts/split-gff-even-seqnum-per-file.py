#!/usr/bin/env python

import sys
import os
import logging

script_dir = os.path.dirname(os.path.abspath(__file__))
snakefile_dir = os.path.dirname(script_dir)
pkg_dir = os.path.dirname(snakefile_dir)
sys.path.append(pkg_dir)
from virsorter.config import get_default_config, set_logger

DEFAULT_CONFIG = get_default_config()
MAX_SPLIT = DEFAULT_CONFIG['MAX_SPLIT']
GFF_SEQNUM_PER_SPLIT = DEFAULT_CONFIG['GFF_SEQNUM_PER_SPLIT']

set_logger()

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
    if n > MAX_SPLIT:
        mes = (
            f'Too many ({n}) splits requested on GFF file for provirus '
            'extraction that may deteriate file system; reducing it '
            f'to {MAX_SPLIT}..\n; If you are running in cluster mode '
            '(virsorter run --cluster), this causes the run time for '
            '"rule provirus" on each split to increase and the default '
            'walltime might may become not enough..\n'
        )
        logging.warning(mes)
        n = MAX_SPLIT

    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)

    try:
        split_idx = 0
        bname = os.path.basename(gff_f)
        fw = open('{}/{}.{}.split'.format(outdir, bname, split_idx), 'w')
        gen = open_gff_by_contig(gff_f)
        for i, rec in enumerate(gen):
            if i < (split_idx + 1) * n:
                fw.write(rec)
            else:
                fw.close()
                split_idx += 1
                fw = open(
                        '{}/{}.{}.split'.format(outdir, bname, split_idx),
                        'w',
                )
                fw.write(rec)

    finally:
        fw.close()


if __name__ == '__main__':
    main()

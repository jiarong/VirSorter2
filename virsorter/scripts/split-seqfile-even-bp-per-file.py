#!/usr/bin/env python
# split file in to even pieces, splited base on file size, save time
#   on iterating through big file

import sys
import screed 
import errno
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

def split_by_bp_per_group(f, n, outdir='.'):
    ''' evenly split with n bp seq in each (except last one).

    Args:
        param1 (str): seqfile path.
        param2 (int): bp per split.
        param3 (str): output directory path.

    Returns:
        int: 0 for sucess, other for failure.

    '''
    outdir = outdir.rstrip('/')
    logging.info('%s bp per group' %n)
    group = 0
    try:
        os.mkdir(outdir)
    except OSError as e:
        if e.errno == errno.EEXIST and os.path.isdir(outdir):
            pass
        else:
            raise

    try:
        fw = open('%s/%s.%d.split' %(outdir, os.path.basename(f), group), 'w')
        total = 0
        for m, record in enumerate(screed.open(f)):
            name = record.name
            seq = record.sequence
            total = total+len(seq)
            fw.write('>%s\n%s\n' %(name, seq))
            if total >= n:
                logging.info(
                        '%s/%s.%d.split written..' \
                                %(outdir, os.path.basename(f), group)
                )
                fw.close()
                group += 1
                total = 0
                fw = open(
                        '%s/%s.%d.split' %(
                            outdir, os.path.basename(f), group), 
                        'w')
        if total == 0:
            # remove if empty file
            os.remove('%s/%s.%d.split' %(outdir, 
                            os.path.basename(f), group))
            group -= 1
        else:
            logging.info('%s/%s.%d.split written' %(
                outdir, os.path.basename(f), group))

        logging.info('%s is splitted into %d' %(
            os.path.basename(f), group+1))
    finally:
        fw.close()


def split_by_seq_num_per_group(f, n, outdir='.'):
    ''' evenly split with n seqs in each (except last one).

    Args:
        param1 (str): seqfile path.
        param2 (int): number of seqs per split.
        param3 (str): output directory path.

    Returns:
        int: 0 for sucess, other for failure.

    '''
    outdir = outdir.rstrip('/')
    sys.stderr.write('*** %s seqs per group\n' %n)
    group = 0
    try:
        os.mkdir(outdir)
    except OSError as e:
        if e.errno == errno.EEXIST and os.path.isdir(outdir):
            pass
        else:
            raise

    try:
        fw = open('%s/%s.%d.split' %(outdir, os.path.basename(f), group), 'w')
        total = 0
        for m, record in enumerate(screed.open(f)):
            name = record.name
            seq = record.sequence
            total += 1
            fw.write('>%s\n%s\n' %(name, seq))
            if total >= n:
                sys.stderr.write(
                        '*** %s/%s.%d.split written..\n' %(
                            outdir, os.path.basename(f), group))
                fw.close()
                group += 1
                total = 0
                fw = open(
                        '%s/%s.%d.split' %(
                            outdir, os.path.basename(f), group
                            ), 'w')

        if total == 0:
            # remove if empty file
            os.remove('%s/%s.%d.split' %(outdir, 
                            os.path.basename(f), group))
        else:
            sys.stderr.write(
                    '*** %s/%s.%d.split written\n' %(
                        outdir, os.path.basename(f), group)
                    )

    finally:
        fw.close()

def main():
    '''evenly split a seqfile with n seqs in each

    Example::

        $ python split-seqfile-even-bp-per-file.py \
                <seqfile> <outdir> num

        <seqfiel>: sequence file to split
        <outdir>: output direcotry
        num: bp of sequences per file

    '''
    if len(sys.argv) != 4:
        mes = '*** Usage: python {} <seqfile> <outdir> num\n'
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    seqfile = sys.argv[1]
    outdir = sys.argv[2]
    n = int(sys.argv[3])
    if n > MAX_SPLIT:
        seqfile_bname = os.path.basename(seqfile)
        if seqfile_bname.endswith('.faa'):
            mes = (
                f'Too many ({n}) splits requested on {seqfile_bname} that '
                'may deteriate file system; reducing it '
                f'to {MAX_SPLIT}..\n; If you are running in cluster mode '
                '(virsorter run --cluster), this causes the run time '
                'of "rule hmmsearch" for each split to increase and the '
                'default walltime might may become not enough..\n'
            )
        else:
            mes = (
                f'Too many ({n}) splits requested on {seqfile_bname} that '
                'may deteriate file system; reducing it '
                f'to {MAX_SPLIT}..\n; If you are running in cluster mode '
                '(virsorter run --cluster), this causes the run time '
                'of "rule gene_call" for each split to increase and the '
                'default walltime might may become not enough..\n'
            )

        logging.warning(mes)
        n = MAX_SPLIT

    split_by_bp_per_group(seqfile, n, outdir)

if __name__ == '__main__':
    main()


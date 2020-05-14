#!/usr/bin/env python
# split file in to even pieces, splited base on file size, save time
#   on iterating through big file

import sys
import screed 
import errno
import os
import logging
from virsorter.config import set_logger

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

    split_by_bp_per_group(seqfile, n, outdir)

if __name__ == '__main__':
    main()


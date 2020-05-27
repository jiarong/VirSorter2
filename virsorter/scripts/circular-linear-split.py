#!/usr/bin/env python

import sys
import os
import logging
import screed

from virsorter.config import set_logger
set_logger()


def main():
    ''' separate circular and linear contigs.

    Circular and linear sequences needs to be handeled differently. This script
    split linear and circular sequences into two different files, and also
    remove the overlapping substring (duplicated) between 5' and 3' at 3' end
    in circular sequence.

    Example:
        python circular-linear-split.py \
                <infile.fa> <circular.fa> <linear.fa> \
                <seqname-length.map> tag min-length

        <infile.fa>: input contig sequence file.
        <circular.fa>: output circular sequence file.
        <circular.fa>: output linear sequence file.
        <seqname-length.map>: tab separated output file with seqname as first
            col, sequence length as second col.
        tag: suffix added to sequence names
        min_length: minimal seq length required (>=0, int)

    '''
    if len(sys.argv) != 7:
        mes = ('*** Usage: python {} '
                              '<infile.fa> <circular.fa> '
                              '<linear.fa> <seqname-length.map> '
                              'tag min-length\n')
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    infile = sys.argv[1]
    circular = sys.argv[2]
    linear = sys.argv[3]
    id_len_map = sys.argv[4]
    tag = sys.argv[5]
    min_length = int(sys.argv[6])
    with open(circular, 'w') as fw1, \
           open(linear, 'w') as fw2, \
           open(id_len_map, 'w') as fw3:
        cnt1 = 0
        cnt2 = 0
        cnt3 = 0
        for rec in screed.open(infile):
            header = rec.name
            name = header.split(None, 1)[0]
            seq = rec.sequence
            if len(seq) < min_length:
                cnt3 += 1
                continue
            prefix = seq[:10]
            i = seq.rfind(prefix, 10)
            # found a match at another location 
            if i != -1:
                suffix = seq[i:]
                # should allow mismatch in future
                if suffix == seq[:len(suffix)]:
                    seq = seq[:i]
                    fw1.write('>{}{}  circular\n{}\n'.format(name, tag, seq))
                    cnt1 += 1
                else:
                    fw2.write('>{}{}  linear\n{}\n'.format(name, tag, seq))

            else:
                fw2.write('>{}{}  linear\n{}\n'.format(name, tag, seq))

            fw3.write('{}\t{}\n'.format(name, len(seq)))
            cnt2 += 1

        logging.info('# of seqs < {} bp and removed: {}'.format(
            min_length, cnt3))
        logging.info('# of circular seqs: {}'.format(cnt1))
        logging.info('# of linear seqs  : {}'.format(cnt2 - cnt1))

if __name__ == '__main__':
    main()

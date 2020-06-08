#!/usr/bin/env python
# by gjr

import sys
import os
import screed
import logging

script_dir = os.path.dirname(os.path.abspath(__file__))
snakefile_dir = os.path.dirname(script_dir)
pkg_dir = os.path.dirname(snakefile_dir)
sys.path.append(pkg_dir)
from virsorter.config import set_logger

N = 10 # e-value cutoff
M = None # bit score cutoff

set_logger()

def get_dict_whole_seq(fs, tags):
    """Parse .tblout file from hmmsearch

    Args:
        fs (str, list of str): tblout files 
        tags (str, list of str): labels (str) matching fs

    Returns:
        dict: a dictionary with sequence name as key (str)
            and e-value, bit score, tag as value (tuple)

    """
    d = {}
    st = set()
    for tag, f in zip(tags, fs):
        with open(f) as fp:
            for line in fp:
                if line.startswith('#'):
                    continue
                line = line.strip()
                lis = line.split()
                try:
                    name = lis[0]
                    hmm_name = lis[2]
                    st.add(name)
                    e_val = float(lis[4])
                    bit = float(lis[5])
                except IndexError as e:
                    mes = ('There is issue with format of {}; '
                            'should be domain tabular output format; ')
                    logging.error(mes.format(f))
                if M != None:
                    if bit < M:
                        continue  
                else:
                    if e_val > N:
                        continue
                old_bit = d.get(name, [1000, -1000, '', ''])[1]

                if bit > old_bit:
                    d[name] = e_val, bit, hmm_name, tag
                else:
                    continue

    total_hits = len(st)
    filtered_hits = len(d)
    mes = '{} out of {} hits are kept after bit score filter ({:.1f})'
    logging.info(mes.format(filtered_hits, total_hits, M))
    return d


def main():
    '''Extract taxonomic feature from hmmsearch outputs

    The script assigns each sequence to the taxonomic group (tag) where its
    best hit belongs.

    Example:
        python extract-feature-from-hmmout.py \
                bit_score_cutoff "<A.tblout>,<B.tblout>.." "A,B.."

        bit_score_cutoff: hmmsearch score cutoff to be considered as a hit
        "<A.tblout>,<B.tblout>..": list of hmmsearch tabular output (--tblout),
            comma separated, no space in between
        "A,B..": list of labels matching hmmsearch tabular output, comma
            separated, no space in between

    '''
    if len(sys.argv) != 4:
        mes = ('python {} bit-score-cutoff '
                '"<arc.tblout>,<bac.tblout>.." "arc.tag,bac.tag.."\n')
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    cutoff = float(sys.argv[1])
    global M
    M = cutoff # change default bit score cutoff
    
    tblout_lis = [ i.strip() for i in sys.argv[2].split(',') ]
    tag_lis = [ i.strip() for i in sys.argv[3].split(',') ]

    d_tblout = get_dict_whole_seq(tblout_lis, tag_lis)  
    #dict of hmm hits, may take big memory if too many hits

    for key in d_tblout:
        tag = d_tblout[key][-1]
        hmm_name = d_tblout[key][-2]
        score = d_tblout[key][1]
        hallmark_cnt = 0

        outline = '{}\t{}\t{}\t{}\n'
        sys.stdout.write(outline.format(key, tag, hmm_name, score))

if __name__ == '__main__':
    main()


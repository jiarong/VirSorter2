#!/usr/bin/env python
# parse tabular output from HMMER, do identity filter and MSA convert
# by gjr; Jan 25, 12

# Parse tabular output from HMMER, do length filter (30 bp minimum)
# Convert .sto to .fa
# Usage: python get-seq-from-hmmout.py <hmmdomtblout> <hmmseqout.sto> <outfile.fa>

import sys
import os
import screed

N = 1000000 # e-value cutoff
#N = 1 # e-value cutoff
M = 50 # bit score cutoff
#M = 30 # bit score cutoff

def get_dict_whole_seq(fs):
    """
    Parse .tblout file from hmmsearch

    Parameters:
    -----------
    f : files
        .tblout file list

    Returns:
    --------
    dict:
        a dictionary with sequence name as key (str)
        and e-value, bit score as value (tuple)

    """

    d = {}
    with open(f) as fp:
        for line in fp:
            if line.startswith('#'):
                continue
            line = line.strip()
            lis = line.split()
            name = lis[0]
            e_val = float(lis[4])
            bit = float(lis[5])
            if M != None:
                if bit < M:
                    continue  
            else:
                if e_val > N:
                    continue
            d[name] = e_val, bit

        return d

def get_hmm_name_st_from_tbl(f):
    """
    Parse .tblout file from hmmsearch

    Parameters:
    -----------
    f : file
        .tblout file

    Returns:
    --------
    set:
        hmm names

    """
    st = set()
    with open(f) as fp:
        for line in fp:
            if line.startswith('#'):
                continue
            line = line.strip()
            lis = line.split()
            seq_name = lis[0]
            hmm_name = lis[2]
            e_val = float(lis[4])
            bit = float(lis[5])
            if M != None:
                if bit < M:
                    continue  
            else:
                if e_val > N:
                    continue

            st.add(hmm_name)

    return st

def get_dict_domain(fs, d_len, cutoff):
    """
    Parse .hmmtblout file from hmmsearch

    Parameters:
    -----------
    fs : files
         .hmmtblout file list
    d_len : dict
         group name as key, average length as value

    Returns:
    --------
    dict:
        a dictionary with sequence name as key (str)
        and e-value, bit score, read length as value (tuple)

    """
    d = {}
    for f in fs:
        #sys.stderr.write('parsing {} ..\n'.format(os.path.basename(f)))
        bname = os.path.basename(f)
        tag = bname.rsplit('.', 1)[0]
        with open(f) as fp:
            for line in fp:
                if line.startswith('#'):
                    continue
                line = line.strip()
                lis = line.split()
                qname = lis[0]  
                # MATOU-v1_62589762_gc0_2
                # for each contig, only the combo of hmm, genetic code, 
                #   and frame with best score will be recorded
                name = qname.split('_gc', 1)[0]
                qlen = int(lis[2])
                e_val = float(lis[6])
                bit = float(lis[7])
                dom_bit = float(lis[13])
                qs = int(lis[17])
                qe = int(lis[18])
                ali = (qe - qs + 1)

                if M != None:
                    if bit < M:
                        continue  
                else:
                    if e_val > N:
                        continue

                old_bit = d.get(name, 
                        [1000, -1000, -1000, -1, -1, '', ''])[1]
                old_dom_bit = d.get(name, 
                        [1000, -1000, -1000, -1, -1, '', ''])[2]

                sto_name = '{}/{}-{}'.format(qname, qs, qe)
                if bit > old_bit:
                    d[name] = e_val, bit, dom_bit, qlen, ali, sto_name, tag
                elif bit == old_bit and dom_bit > old_dom_bit:
                    d[name] = e_val, bit, dom_bit, qlen, ali, sto_name, tag


    #print(d['Picorna-Calici_spider134013_Hubei_picorna-like_virus_55_len9795'])
    to_remove = set()
    for name in d:
        tag = d[name][-1]
        ali = d[name][4]
        if ali < d_len[tag] * cutoff:
            to_remove.add(name)

    for name in to_remove:
        _t = d.pop(name)

    return d

def get_set_domain_no_filter(fs):
    """
    Parse .hmmtblout file from hmmsearch

    Parameters:
    -----------
    fs : files
         .hmmtblout file list

    Returns:
    --------
    set:
        a set of sequence name

    """
    st = set()
    for f in fs:
        #sys.stderr.write('parsing {} ..\n'.format(os.path.basename(f)))
        bname = os.path.basename(f)
        tag = bname.rsplit('.', 1)[0]
        with open(f) as fp:
            for line in fp:
                if line.startswith('#'):
                    continue
                line = line.strip()
                lis = line.split()
                qname = lis[0]  
                # MATOU-v1_62589762_gc0_2
                # for each contig, only the combo of hmm, genetic code, 
                #   and frame with best score will be recorded
                name = qname.split('_gc', 1)[0]
                qlen = int(lis[2])
                e_val = float(lis[6])
                bit = float(lis[7])
                dom_bit = float(lis[13])
                qs = int(lis[17])
                qe = int(lis[18])
                ali = (qe - qs + 1)

                if M != None:
                    if bit < M:
                        continue  
                else:
                    if e_val > N:
                        continue

                st.add(qname)

    return st

def main():
    #read seqs into a dict, not memory efficient
    if len(sys.argv) != 3:
        mes = 'python {} scorecutoff <gene.domtblout> \n'
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    cutoff = int(sys.argv[1])
    domtblout_f = sys.argv[2]
    fname = os.path.basename(domtblout_f)
    tag = fname.rsplit('.', 1)[0]

    global M
    M = cutoff

    st = get_hmm_name_st_from_tbl(domtblout_f)

    for hmm_name in st:
        sys.stdout.write('{}\t{}\t{}\n'.format(hmm_name, tag, M))

if __name__ == '__main__':
    main()


#!/usr/bin/env python
# Usage: python <thisfile> <Pfam-A.hmm.gz> <Pfam-A.regions.uniprot.tsv.gz> \
#                 <uniprot_sprot_archaea.dat.gz> ..

import sys
import os
import gzip
import bz2
from collections import Counter

magic_dict = {
    b'\x1f\x8b\x08': 'gz',
    b'\x42\x5a\x68': 'bz2',
    b'\x50\x4b\x03\x04': 'zip'
    }

max_len = max(len(x) for x in magic_dict)

def file_type(filename):
    # have to open as binary here so whole file as bytes
    #   coz magic number is not valid in unicode
    with open(filename, 'rb') as f:
        file_start = f.read(max_len)
    for magic, filetype in magic_dict.items():
        if file_start.startswith(magic):
            return filetype
    return "no match"

def parse_pfam2uniprot(f):
    d_pfam2uniprot = {}
    st_uniprot_acc_in_pfam = set()
    try:
        if file_type(f) == 'gz':
            # default is 'rb' (binary mode), which creates byte string
            # 'rt' (text mode) creates str (unicode)
            fp = gzip.open(f, 'rt')
        elif file_type(f) == 'bz2':
            fp = bz2.open(f, 'rt')
        else:
            fp = open(f)

        for line in fp:
            if line.startswith('#'):
                continue
            if line.startswith('uniprot_acc'):
                continue
            line = line.rstrip()
            li = line.split()
            uniprot_acc = li[0]
            pfam_acc = li[4]
            st = d_pfam2uniprot.setdefault(pfam_acc, set())
            d_pfam2uniprot[pfam_acc].add(uniprot_acc)
            st_uniprot_acc_in_pfam.add(uniprot_acc)

    finally:
        fp.close()

    return st_uniprot_acc_in_pfam, d_pfam2uniprot

def main():
    if len(sys.argv) != 4:
        mes = ('*** Usage: python {} <Pfam-A.hmm.gz> '
                        '<Pfam-A.regions.uniprot.tsv.gz> '
                        '<uniprot-acc2taxon.list> \n')
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    pfam_hmm_f = sys.argv[1]
    pfam2uniprot_map_f = sys.argv[2]
    uniprot_acc2tax_f = sys.argv[3]

    st_uniprot_acc_in_pfam, d_pfam2uniprot = \
            parse_pfam2uniprot(pfam2uniprot_map_f)
    
    
    d_uniprot2tax = {}
    kindoms = set()
    try:
        if uniprot_acc2tax_f == '-':
            f = '/dev/stdin'
        else:
            f = uniprot_acc2tax_f

        if file_type(f) == 'gz':
            fp = gzip.open(f, 'rt')
        elif file_type(f) == 'bz2':
            fp = bz2.open(f, 'rt')
        else:
            fp = open(f)

        for line in fp:
            if line.startswith('#'):
                continue
            line = line.rstrip()
            acc, taxon = line.split('\t')
            if not acc in st_uniprot_acc_in_pfam:
                continue
            taxon_levels = taxon.split(';')
            kindom = taxon_levels[0]

            kindoms.add(kindom)
            d_uniprot2tax[acc] = kindom

    finally:
        fp.close()

    d_pfam2tax = {}
    cnt_total = 0
    cnt_missing = 0
    for key in d_pfam2uniprot:

        li = []
        for uniprot in d_pfam2uniprot[key]:
            cnt_total += 1
            try:
                li += [d_uniprot2tax[uniprot],]
            except KeyError as e:
                cnt_missing += 1
                #mes= (
                #        '*** Warning: {} in {} NOT found in'
                #        '<uniprot-acc2taxon.list>, ignored..\n'
                #)
                #sys.stderr.write(
                #        mes.format(uniprot, 
                #            os.path.basename(pfam2uniprot_map_f)
                #        )
                #)
        total = len(li)
        d_kindom_cnt = Counter(li)
        d_kindom_frac = dict(
                (kindom, d_kindom_cnt[kindom]*1.0/total) 
                for kindom in d_kindom_cnt
        )
        d_pfam2tax[key] = d_kindom_frac


    mes= (
        '*** Warning: {} out of {} in {} are NOT found in'
        '<uniprot-acc2taxon.list>\n'
    )
    sys.stderr.write(
            mes.format(cnt_missing, 
                cnt_total, 
                os.path.basename(pfam2uniprot_map_f)
            )
    ) 

    try:
        d_kindom_fp = {}
        kindoms.add("Mixed")
        for kindom in kindoms:
            kindom = kindom.split()[0]
            if top_kindom in set(['unclassified', 'other']):
                # no unclassified or other > 0.9 in Pfam-A 32.0
                #   downloaded on Aug 2, 2019
                # put this conditional here just in case
                top_kindom = 'Mixed'

            d_kindom_fp[kindom] = open('Pfam-A-{}.hmm'.format(kindom), 'w') 

        if file_type(pfam_hmm_f) == 'gz':
            fp = gzip.open(pfam_hmm_f, 'rt')
        elif file_type(pfam_hmm_f) == 'bz2':
            fp = bz2.open(pfam_hmm_f, 'rt')
        else:
            fp = open(pfam_hmm_f)

        s = ''
        for line in fp:
            if line.startswith('ACC'):
                line1 = line.rstrip()
                acc_with_version = line1.split()[1]
                acc = acc_with_version.rsplit('.', 1)[0]

            if line == '//\n':
                d_kindom_frac = d_pfam2tax[acc]
                sorted_items = sorted(
                        d_kindom_frac.items(), 
                        key=lambda x: x[1]
                        )
                top_kindom, relabund = sorted_items[-1]
                if relabund < 0.9:
                    top_kindom = 'Mixed'
                elif top_kindom in set(['unclassified', 'other']):
                    # no unclassified or other > 0.9 in Pfam-A 32.0
                    #   downloaded on Aug 2, 2019
                    # put this conditional here just in case
                    top_kindom = 'Mixed'
                #if top_kindom == 'Viruses' and relabund >= 0.9:
                #    continue
                d_kindom_fp[top_kindom].write('{}//\n'.format(s))
                s = ''
                continue

            s = '{}{}'.format(s,line)

    finally:
        fp.close()
        for kindom in d_kindom_fp:
            d_kindom_fp[kindom].close()

if __name__ == '__main__':
    main()

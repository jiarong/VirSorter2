#!/usr/bin/env python

import sys
import os
import gzip
import bz2


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

def main():
    if len(sys.argv) < 2:
        mes = '*** Usage: python {} <uniprot_sprot_archaea.dat.gz> ..\n'
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    uniprot_fs = sys.argv[1:]

    accs = []
    taxon_levels = []
    trigger = False
    for f in uniprot_fs:
        sys.stderr.write('*** Processing: {}\n'.format(os.path.basename(f)))
        try:
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

                if line.startswith('AC'):
                    # AC   Q9V2L2; G8ZFP4;
                    s = line.split(None, 1)[1]
                    s = s.rstrip(';')
                    li = [i.strip() for i in s.split(';')]
                    accs = li
                    continue

                if line.startswith('OC'):
                    if not trigger:
                        trigger = True
                    s = line.split(None, 1)[1]
                    s = s.rstrip(';')
                    s = s.rstrip('.')
                    li = [i.strip() for i in s.split(';')]
                    taxon_levels += li

                elif trigger:
                    for i in accs:
                        ###
                        # only take the domain level here
                        ###
                        sys.stdout.write('{}\t{}\n'.format(i, taxon_levels[0]))

                    taxon_levels = []
                    trigger = False

        finally:
            fp.close()

if __name__ == '__main__':
    main()

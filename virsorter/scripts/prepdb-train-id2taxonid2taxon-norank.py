#!/usr/bin/env python
# make table of id, taxonid, lineage
# by gjr; 061715

import sys
import os
import ete3
from ete3 import NCBITaxa



def main():
    if len(sys.argv) != 3:
        mes = 'Usage: python {} <id2taxonid.table> <id2taxonid2taxon.table>\n'
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    
    infile = sys.argv[1]
    outfile = sys.argv[2]

    if infile == '-':
        infile = '/dev/stdin'

    if outfile == '-':
        outfile = '/dev/stdout'

    ncbi = NCBITaxa()

    with open(infile) as fp, open(outfile, 'w') as fw:
        cnt = 0
        for line in fp:
            line = line.rstrip()
            id, taxonid = line.split()
            try:
                lineage_taxid = ncbi.get_lineage(taxonid)
            except ValueError as e:
                mes = '*** Invalid TaxonID (not in NCBI Taxonomy db): {}\n'
                sys.stderr.write(mes.format(taxonid))
                continue

            lineage_tax_d = ncbi.get_taxid_translator(lineage_taxid)
            lineage_rank_d = ncbi.get_rank(lineage_taxid)

            lineage_tax_lis = [
                    lineage_tax_d[taxid] for taxid in lineage_taxid 
                    ]
            lineage_rank_lis = [
                    lineage_rank_d[taxid].replace(' ', '') \
                            for taxid in lineage_taxid
                    ]

            # skip 'root'
            lis = zip(lineage_rank_lis[1:], lineage_tax_lis[1:])
            lis2 = ['{}:{}'.format(x, y) for x, y in lis]

            fw.write('{}\t{}\t{}\n'.format(id, taxonid, 
                                            '\t'.join(lis2)))


if __name__ == '__main__':
    main()

#!/usr/bin/env python
# make table of id, taxonid, lineage
# by gjr; 061715

import sys
import os
import ete3
from ete3 import NCBITaxa



def main():
    if len(sys.argv) != 4:
        mes = 'Usage: python {} "phylum,genus" <id2taxonid.table> <id2taxonid2taxon.table>\n'
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    
    REQUIRED_RANK = [i.strip() for i in sys.argv[1].strip().split(',')]
    infile = sys.argv[2]
    outfile = sys.argv[3]

    if infile == '-':
        infile = '/dev/stdin'

    if outfile == '-':
        outfile = '/dev/stdout'

    ncbi = NCBITaxa()
    rank_lis_needed = ['superkingdom','phylum','class',
                            'order','family','genus','species']


    with open(infile) as fp, open(outfile, 'w') as fw:
        cnt = 0
        for line in fp:
            line = line.rstrip()
            id, taxonid = line.split()
            try:
                lineage_num = ncbi.get_lineage(taxonid)
            except ValueError as e:
                mes = '*** Invalid TaxonID (not in NCBI Taxonomy db): {}\n'
                sys.stderr.write(mes.format(taxonid))
                continue
            lineage_tax = ncbi.translate_to_names(lineage_num)


            rank_dict = ncbi.get_rank(lineage_num)
            st_rank = set(rank_dict.values())
        
            if not set(REQUIRED_RANK).issubset(st_rank):
                mes = (
                        '*** Irregular taxon format with TaxonID ({})'
                        '(skipped): \n{}\n{}\n'
                )
                sys.stderr.write(mes.format(taxonid, repr(lineage_tax), repr(st_rank)))
                cnt += 1
                continue


            n_temp = -1
            index_temp = -1
            lis_lineage_taxon = []
            for n, idd in enumerate(lineage_num):
                rank = rank_dict[idd]
                if rank in set(rank_lis_needed):
                    index = rank_lis_needed.index(rank)
                    if index != index_temp + 1:
                        skipped  = index -1 - index_temp
                        # add skipped levels
                        while 1:
                            lis_lineage_taxon.append('Other')
                            skipped -= 1
                            if skipped == 0:
                                break
                            n_temp += 1


                    #assert n_temp < n, 'n_temp: {}; n: {}'.format(n_temp, n)
                    lis_lineage_taxon.append(lineage_tax[n])
                    n_temp = n
                    index_temp = index
    
            l1 = len(lis_lineage_taxon)
            l2 = len(rank_lis_needed)

            if l1 != l2:
                lis_lineage_taxon.extend(['Other']*(l2 - l1))

            assert len(lis_lineage_taxon) == len(rank_lis_needed)
            fw.write('{}\t{}\t{}\n'.format(id, taxonid, 
                                            '\t'.join(lis_lineage_taxon)))

        mes = (
                '*** Number of taxonid with irregular taxonomy '
                '(without {} info): {}'
        )
        sys.stderr.write(mes.format(', '.join(REQUIRED_RANK), cnt))


if __name__ == '__main__':
    main()

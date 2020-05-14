#!/usr/bin/env python
# Usage: python <thisfile> <in.gff> <out.gff>

import sys
import os
import screed
from collections import OrderedDict
import logging
from virsorter.config import set_logger

set_logger()

def main():
    ''' remove partial gene in gff file from prodigal.

    There are partial genes resulted from break point of circular
    genomes. This script remove these partial genes from the 5' and 3'
    of extended circular sequences output from circular-extend.py 

    Example:
        python circular-remove-partial-gene.py <in.gff> <out.gff>

        <in.gff>: gff file of extended circular sequence file output 
                from circular-extend.py
        <out.gff>: gff file with partial genes removed from 5' and 3'
            ends
        
    '''
    if len(sys.argv) != 3:
        mes = '*** Usage: python {} <in.gff> <out.gff>\n'
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    in_gff = sys.argv[1]
    out_gff = sys.argv[2]
    if out_gff == '-':
        out_gff = '/dev/stdout'
    with open(in_gff) as fp, open(out_gff, 'w') as fw:
        cnt1 = 0
        cnt2 = 0
        cnt3 = 0
        lines = []
        cnt_gene = 0
        seqname = None
        for line in fp:
            if line.startswith('# Sequence Data:'):
                """
                The fields in this header are as follows:

                seqnum: An ordinal ID for this sequence, beginning at 1.
                seqlen: Number of bases in the sequence.
                seqhdr: The entire FASTA header line.
                version: Version of Prodigal used to analyze this sequence.
                run_type: "Ab initio" for normal mode, "Anonymous" for anonymous mode.
                model (Anonymous mode only): Information about the preset training file used to analyze the sequence.
                gc_cont: % GC content of the sequence.
                transl_table: The genetic code used to analyze the sequence.
                uses_sd: Set to 1 if Prodigal used its default RBS finder, 0 if it scanned for other motifs.
                """
                if seqname == None:
                    pass
                elif cnt_gene >= 2:
                    fw.write(''.join(lines)) 
                else:
                    mes = 'Conitg removed due to < 2 genes: {}'
                    logging.info(mes.format(seqname))
                    cnt3 += 1

                lines = []
                cnt_gene = 0
                line = line.rstrip()
                s = line.split('Sequence Data:')[1].strip()
                seq_data = OrderedDict(
                        i.strip().split('=') for i in s.split(';')
                        )
                seqlen = seq_data['seqlen']
                seqhdr = seq_data['seqhdr'].strip('"')
                seqname = seqhdr.split(None, 1)[0]
                # circular was duplicated after removing overlap
                seqlen = int(int(seqlen)/2.0)
                seq_data['seqlen'] = str(seqlen)

                l = ['='.join((i, seq_data[i])) for i in seq_data]
                line = '# Sequence Data: {}\n'.format(';'.join(l))
                lines += [line,]
                continue

            elif line.startswith('#'):
                fw.write(line)
                lines += [line,]
                continue

            line = line.rstrip()
            """
            [0]Seqname\t[1]Prodigal_v2.6.3\t[2]CDS\t[3]start\t[4]end\t[5]score\t[6]strand\t[7]frame\t[8]semicolon-delimited-string\n   
            --------------
            fields in semicolon-delimited-string:
            [0]ID: A unique identifier for each gene, consisting of the ordinal ID of the sequence and an ordinal ID of that gene within the sequence (separated by an underscore). For example, "4_1023" indicates the 1023rd gene in the 4th sequence in the file.
            [1]partial: An indicator of if a gene runs off the edge of a sequence or into a gap. A "0" indicates the gene has a true boundary (a start or a stop), whereas a "1" indicates the gene is "unfinished" at that edge (i.e. a partial gene). For example, "01" means a gene is partial at the right boundary, "11" indicates both edges are incomplete, and "00" indicates a complete gene with a start and stop codon.
            [2]start_type: The sequence of the start codon (usually ATG, GTG, or TTG). If the gene has no start codon, this field will be labeled "Edge".
            [3]rbs_motif: The RBS motif found by Prodigal (e.g. "AGGA" or "GGA", etc.)
            [4]rbs_spacer: The number of bases between the start codon and the observed motif.
            [5]gc_cont: The GC content of the gene sequence.
            [6]conf: A confidence score for this gene, representing the probability that this gene is real, i.e. 78.3% means Prodigal believes that gene is real 78.3% of the time and a false positive 21.7% of the time.
            [7]score: The total score for this gene.
            [8]cscore: The hexamer coding portion of the score, i.e. how much this gene looks like a true protein.
            [9]sscore: A score for the translation initiation site for this gene; it is the sum of the following three fields.
            [10]rscore: A score for the RBS motif of this gene.
            [11]uscore: A score for the sequence surrounding the start codon.
            [12]tscore: A score for the start codon type (ATG vs. GTG vs. TTG vs. Nonstandard).
            """
            items = line.split('\t')
            last = items[8] # ; delimited string

            sub_items = OrderedDict(
                    i.strip().split('=') for i in last.rstrip(';').split(';')
                    )
            partial =  sub_items['partial']
            if partial != '00':
                # discard partial gene in circular contigs
                cnt1 += 1
                continue
            start = int(items[3])
            if start > seqlen:
                # discard genes on extended part
                continue
            end = int(items[4])
            frame = int(items[7])
            ###
            # decide NOT to move the end of last gene pass boudry
            #   to the begin of seq, coz it cause issue with gene size
            #   estimation
            ###
            if start <= seqlen and end > seqlen:
                # a gene extends from end (3') to start (5') 
                # change alignment coordinate
                cnt2 += 1
            #    end = end - seqlen
            #    items[4] = str(end)
            #    line = '\t'.join(items)

            line = '{}\n'.format(line)
            lines += [line,]
            cnt_gene += 1

        if seqname == None:
            pass
        if cnt_gene >= 2:
            fw.write(''.join(lines)) 
        else:
            mes = 'Conitg removed due to < 2 genes: {}'
            logging.info(mes.format(seqname))
            cnt3 += 1

    mes = '# of partial genes removed from circular contigs: {}'
    logging.info(mes.format(cnt1))
    mes = '# of genes cover both 5\' and 3\' from circular contigs: {}'
    logging.info(mes.format(cnt2))
    mes = '# of contigs removed due to having < 2 genes: {}'
    logging.info(mes.format(cnt3))

if __name__ == '__main__':
    main()

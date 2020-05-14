#!/usr/bin/env python
# Usage: python <thisfile> <file.gff> <in.fa> <out.fa>

import sys
import os
import screed
from collections import OrderedDict

def gff_parser(f):
    '''parse gff from prodigal

    args:
        f (str): gff file from prodigal

    returns:
        list of str: a list seqname of orf in faa file
            ({contig-seqname}_{index} from prodigal

    '''
    with open(f) as fp:
        for line in fp:
            if line.startswith('#'):
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

            seqname = items[0]
            index = sub_items['ID'].split('_')[1]
            seqname_in_fasta = '{}_{}'.format(seqname, index) 
            yield seqname_in_fasta

def main():
    '''filter seqnames in faa based on records in gff
    
    Example:
        python filter-seqs-by-gff.py <gff> <in.faa> <out.faa>

        <gff>: gff file from prodigal
        <in.faa>: faa file from prodigal
        <out.faa>: faa file with orfs not in gff removed

    '''
    if len(sys.argv) != 4:
        mes = '*** Usage: python {} <gff> <in.fa> <out.fa>\n'
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    
    gff = sys.argv[1]
    in_fa = sys.argv[2]
    out_fa = sys.argv[3]
    if out_fa == '-':
        out_fa = '/dev/stdout'

    st = set(gff_parser(gff))
    with screed.open(in_fa) as fp, open(out_fa, 'w') as fw:
        for rec in fp:
            header = rec.name
            name = header.split(None, 1)[0]
            if not name in st:
                continue
            fw.write('>{}\n{}\n'.format(name, rec.sequence))

if __name__ == '__main__':
    main()

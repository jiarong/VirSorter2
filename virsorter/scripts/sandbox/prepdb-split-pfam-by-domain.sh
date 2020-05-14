#!/bin/bash
#SBATCH --nodes=1 --time=48:00:00 --mem=5G


cd ~/Documents/dev/VirSorter2/prepdb/prep-pfam

### download pfam
#wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/Pfam-A.regions.uniprot.tsv.gz
#wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/Pfam-A.hmm.dat.gz
#wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/Pfam-A.hmm.gz

### download unprot mapping data (both swissprot and tremble, all taxa)
wget -r -np -nH --cut-dirs=6 --accept "uniprot*.dat.gz" --reject uniprot_trembl_bacteria.dat.gz ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions

### download tremble bacteria separately in another script due to large size
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_bacteria.dat.gz

time  python ../../scripts/prepdb-parse-uniprot.py uniprot_*.dat.gz > uniprot-acc2taxon.list
time python ../../scripts/prepdb-split-pfam-by-domain.py Pfam-A.hmm.gz Pfam-A.regions.uniprot.tsv.gz uniprot-acc2taxon.list

scontrol show job ${SLURM_JOBID}

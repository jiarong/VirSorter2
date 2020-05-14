cd /mnt/home/guojiaro/Documents/dev/VirSorter2/test

Outdir=test-out
rm -rf $Outdir
#/usr/bin/time -v virsorter run -w $Outdir -d ../db -i 4seq.fa -j 8 --hallmark-required-on-short all
/usr/bin/time -v virsorter run -w $Outdir -d ../db -i 4seq.fa -j 8 --hallmark-required-on-short --provirus all

#Outdir=test
#python virsorter/scripts/provirus.py $Outdir/iter-0/caudovirales/all.pdg.gff $Outdir/iter-0/caudovirales/all.pdg.hmm.tax db/rbs/rbs-catetory.tsv db/group/caudovirales/model zzz zzz.ftr --hallmark db/group/caudovirales/hallmark-gene.list
#python virsorter/scripts/provirus.py $Outdir/iter-0/caudovirales/all.pdg.gff $Outdir/iter-0/caudovirales/all.pdg.hmm.tax db/rbs/rbs-catetory.tsv db/group/caudovirales/model zzz zzz.ftr --hallmark db/group/caudovirales/hallmark-gene.list --fullseq-clf test/iter-0/viral-contig-proba.tsv --group Caudo
#python virsorter/get-hallmark-cnt-for-each-seq.py zz.hallmark.tsv  "caudovirales,NCLDV" "db/group/caudovirales/hallmark-gene.list,db/group/NCLDV/hallmark-gene.list" "test2/iter-0/caudovirales/all.pdg.hmm.tax,test2/iter-0/NCLDV/all.pdg.hmm.tax"
#python virsorter/scripts/remove-short-seq-wo-hallmark.py zz.hallmark.tsv test2/final-viral-fullseq-contig.fa

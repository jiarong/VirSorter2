cd /mnt/home/guojiaro/Documents/dev/VirSorter2/test

Outdir=test-out
rm -rf $Outdir
#/usr/bin/time -v virsorter run -w $Outdir -d ../db -i 4seq.fa -j 8 --hallmark-required-on-short all
/usr/bin/time -v virsorter run -w $Outdir -d ../db -i 4seq.fa -j 8 --hallmark-required-on-short --provirus all

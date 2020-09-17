
Outdir=test-out
rm -rf $Outdir
/usr/bin/time -v virsorter run  -w $Outdir -i 8seq.fa -j 4 all
#/usr/bin/time -v virsorter run  -w $Outdir -i 8seq.fa -j 4 classify
#/usr/bin/time -v virsorter run  -w $Outdir -i 8seq.fa -j 4 --hallmark-required-on-short all

#!/bin/bash

for sample in `ls -d NA*/ | cut -d'/' -f 1`
do
if [ -s $sample/background.coverage_results.bed ] && [ -f $sample/cal.target.done ] ; then
  echo $sample
  back_med=`cat $sample/background.coverage_results.bed |\
    awk '{print $4}' | sort -n |\
	  awk '{ a[i++]=$1; } END { if (i%2==0) { print ( a[int(i/2)] + a[int((i-1)/2)])/2; } else { print a[int(i/2)]; } }'`

  sex=$(awk '{print $NF}' $sample/sex.determine.txt)

  dj=$(cut -f 4 $sample/DJ.coverage_results.bed)
  djcopy=`echo "scale=2; 2*$dj/$back_med" | bc`

  fiveS=$(cut -f 4 $sample/5S.coverage_results.bed)
  fiveScopy=`echo "scale=2; 2*$fiveS/$back_med" | bc`

  rdna=$(cut -f 5 $sample/rDNA.coverage_results.bed)
  rdnacopy=`echo "scale=2; 2*$rdna/$back_med" | bc`

  rdna18s=$(cut -f 5 $sample/rDNA_18S.coverage_results.bed)
  rdna18scopy=`echo "scale=2; 2*$rdna18s/$back_med" | bc`

  chr13_PHR=$(cut -f 5 $sample/PHR_keep.coverage_results.bed)
  chr13_PHRcopy=`echo "scale=2; 2*$chr13_PHR/$back_med" | bc`

  chr13_PHR_arm1=$(fgrep "chr13_PHR_arm1" $sample/PHR_forCov.coverage_results.bed | cut -f 5)
  chr13_PHR_arm1_copy=`echo "scale=2; 2*$chr13_PHR_arm1/$back_med" | bc`

  chr13_PHR_arm2=$(fgrep "chr13_PHR_arm2" $sample/PHR_forCov.coverage_results.bed | cut -f 5)
  chr13_PHR_arm2_copy=`echo "scale=2; 2*$chr13_PHR_arm2/$back_med" | bc`

  echo -e "$sample\t$sex\t$djcopy\t$fiveScopy\t$rdnacopy\t$rdna18scopy\t$chr13_PHRcopy\t$chr13_PHR_arm1_copy\t$chr13_PHR_arm2_copy" > $sample/$sample.copynum_dip.txt
fi
done
out=copynum_dip.txt

cat */*.$out | sed -e '1i Sample\tSex\tDJ\t5S\trDNA\trDNA18S\tChr13_PHR\tChr13_PHR_arm1\tChr13_PHR_arm2' > $out
awk 'NR==1 || $NF>3.2 && $NF<4.2' $out > $out.PHR
awk 'NR==1 || $3>8.5 && $3<9.5' $out > $out.DJ_9
awk 'NR==1 || $3<8.5' $out > $out.ROB

cat $out

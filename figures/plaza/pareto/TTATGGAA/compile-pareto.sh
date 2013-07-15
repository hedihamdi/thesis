#!/bin/bash
cd /home/santolin/these/programs/scangen/scangen-2.0/build/;make;cd -
ln -sf /home/santolin/these/programs/scangen/scangen-2.0/build/src/scangen

for k in enhancers #targets
do
   echo $k
   if [ "$k" == "enhancers" ]; then
      pos='/home/santolin/these/files/droso/plaza/coords/201012/coords-pos.dat'
      neg='/home/santolin/these/files/droso/plaza/coords/201203-full-neg/neg-dm2.dat'
      ftarg='roc-enhancers-info.dat'
   elif [ "$k" == "targets" ];then
      pos='/home/santolin/these/files/droso/plaza/target-genes/align/regs-pos-align-no-training-set.dat'
      neg='/home/santolin/these/files/droso/plaza/target-genes/align/regs-neg-align.dat'
      ftarg='roc-target-info.dat'
   fi
   echo pos: $pos
   echo neg: $neg

   [ -e $ftarg ] && rm $ftarg


   sw=1000

   n=1

   for tmot in 5 6 7 8 9 10 11 12
   do
      mot="/home/santolin/these/files/droso/plaza/TTATGGAA/t$tmot/TTATGGAA-wname.dat"
      for j in $tmot-5 $tmot-4 $tmot-3 $tmot-2 $tmot-1 $tmot
      do
         t=`echo $j | bc -l`
         echo $mot, thres $t
         if [ "$k" == "enhancers" ]; then
            ./scangen --roc --info --poscoords $pos --negcoords $neg -s 1 --mots-w-names $mot -t $t -n $n
         elif [ "$k" == "targets" ];then
            ./scangen --roc --info --posalign $pos --negalign $neg -s 1 --mots-w-names $mot -t $t -n $n
         fi
         awk 'NR>1{print evol,tmot,tdet,$0}' evol="TTATGGAA" tmot=$tmot tdet=$t score/roc-infos.dat >> $ftarg
         if [ "$k" == "enhancers" ]; then
            ./scangen --roc --info --poscoords $pos --negcoords $neg -s 1 --mots-w-names $mot -t $t -n $n --wocons
         elif [ "$k" == "targets" ];then
            ./scangen --roc --info --posalign $pos --negalign $neg -s 1 --mots-w-names $mot -t $t -n $n  --wocons
         fi
         awk 'NR>1{print evol,tmot,tdet,$0}' evol="TTATGGAA_no_cons" tmot=$tmot tdet=$t score/roc-infos.dat >> $ftarg
      done
   done

done


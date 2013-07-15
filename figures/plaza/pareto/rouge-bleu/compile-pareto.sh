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


   mot="bestmots.dat"
   n=2
   for t in 5 6 7 8 9
   do
      awk 'NR==1' bestmots.dat > foo
      echo foo, thres $t
      if [ "$k" == "enhancers" ]; then
         ./scangen --roc --info --poscoords $pos --negcoords $neg -s 1 --mots-w-names foo -t $t -n $n
      elif [ "$k" == "targets" ];then
         ./scangen --roc --info --posalign $pos --negalign $neg -s 1 --mots-w-names foo -t $t -n $n
      fi
      awk 'NR>1{print evol,tmot,tdet,$0}' evol="rouge" tmot=7 tdet=$t score/roc-infos.dat >> $ftarg
      if [ "$k" == "enhancers" ]; then
         ./scangen --roc --info --poscoords $pos --negcoords $neg -s 1 --mots-w-names foo -t $t -n $n --wocons
      elif [ "$k" == "targets" ];then
         ./scangen --roc --info --posalign $pos --negalign $neg -s 1 --mots-w-names foo -t $t -n $n  --wocons
      fi
      awk 'NR>1{print evol,tmot,tdet,$0}' evol="rouge_nocons" tmot=7 tdet=$t score/roc-infos.dat >> $ftarg
      
      awk 'NR==2' bestmots.dat > foo
      echo foo, thres $t
      if [ "$k" == "enhancers" ]; then
         ./scangen --roc --info --poscoords $pos --negcoords $neg -s 1 --mots-w-names foo -t $t -n $n
      elif [ "$k" == "targets" ];then
         ./scangen --roc --info --posalign $pos --negalign $neg -s 1 --mots-w-names foo -t $t -n $n
      fi
      awk 'NR>1{print evol,tmot,tdet,$0}' evol="bleu" tmot=7 tdet=$t score/roc-infos.dat >> $ftarg
      if [ "$k" == "enhancers" ]; then
         ./scangen --roc --info --poscoords $pos --negcoords $neg -s 1 --mots-w-names foo -t $t -n $n --wocons
      elif [ "$k" == "targets" ];then
         ./scangen --roc --info --posalign $pos --negalign $neg -s 1 --mots-w-names foo -t $t -n $n  --wocons
      fi
      awk 'NR>1{print evol,tmot,tdet,$0}' evol="bleu_nocons" tmot=7 tdet=$t score/roc-infos.dat >> $ftarg
   done
   echo "combi 7_7 6_6 1 1 0.5 3 5" >> $ftarg
   echo "combi_nocons 7_7 6_6 1 1 0.5 4 5" >> $ftarg

done


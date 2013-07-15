#!/bin/bash

cd /home/santolin/these/programs/scangen/scangen-2.0/build/;make;cd -
ln -sf /home/santolin/these/programs/scangen/scangen-2.0/build/src/scangen

mot='bestmots.dat'
opt='bestmots-optthr.dat'
n=2
poscoords='/home/santolin/these/files/droso/plaza/coords/201012/coords-pos.dat'
negcoords='/home/santolin/these/files/droso/plaza/coords/201012/coords-neg.dat'

[ -d score/dispscore ] && rm -r score/dispscore
./scangen --dispscore --weightmots --poscoords $poscoords --negcoords $negcoords -s 1 --mots-w-name $mot --opt-thr-file $opt
./scangen --dispscore --weightmots --poscoords $poscoords --negcoords $negcoords -s 1 --mots-w-name $mot --opt-thr-file $opt --wocons
cd score/dispscore/
   for i in *.tex
   do
      pdflatex $i
   done
cd ../..

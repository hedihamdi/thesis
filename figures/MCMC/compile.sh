#!/bin/bash
#cd /home/santolin/these/programs/imogene/build/;make install;cd -
imogene test --MCMC -a align -s "eutherian" --bsinit 'GGTGCTGAGT' -t 10 -e 1 #--nbmotsmin 50

#./plotSampling.R

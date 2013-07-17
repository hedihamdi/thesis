#!/bin/bash

awk '{if ($1 ~ ">"){printf("%s\t",substr($1,2))} else{printf("%s\n",$1)}}' fasta/mef3-pos-coordsavailable.dat > mef3-pos-coordsavailable.dat
awk '{if ($1 ~ ">"){printf("%s\t",substr($1,2))} else{printf("%s\n",$1)}}' fasta/mef3-pos-outoflearning.dat > mef3-pos-outoflearning.dat

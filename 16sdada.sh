#! /usr/bin/env bash

#BSUB -J 16sdadasp
#BSUB -n 8
#BSUB -R "span[hosts=1]"
#BSUB -W 24:00
#BSUB -o dada16sp_%J.out
#BSUB -e dada16sp_%J.err

Rscript 16S_dadasp.R

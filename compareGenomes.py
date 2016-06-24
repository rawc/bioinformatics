#!/usr/bin/env python3
# Name: Louis Chatta (lchatta)
# Group Members: (lchatta, cjavery)

from genomeAnalyzer import genomeCompare
# By creating the genomeCompare class with an init 
# We can easily add functionality to multiple files
print ('This program will compare the following two files:')
genomeCompare(filenames=['testGenome.fa','haloVolc1_1-genes.fa'],compare=True)


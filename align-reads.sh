#!/bin/bash
# Align reads in reads.fasta and put them in msa.fasta

mafft --globalpair --retree 2 --maxiterate 2 --anysymbol reads.fasta | \
  tr '\n' ' ' | sed -E 's/$/\n/;s/>(\w+) /\n>\1\n/g;s/ //g' | tail -n+2 > msa.fasta

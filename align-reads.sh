#!/bin/bash
# Align reads in reads.fasta and put them in msa.fasta

mafft --globalpair --retree 2 --maxiterate 2 --anysymbol reads.fasta > msa.fasta

#!/bin/bash

set -euxo pipefail

# Convert .wav files to a single file -> reads.fasta
./voice2fasta.sh

# Align reads.fasta using mafft -> msa.fasta
./align-reads.sh

# Call the consensus on msa.fasta
./consensus.py msa.fasta | tee consensus.fasta

# Re-render slides
quarto render sequence-alignment.qmd

# Send notification when done
notify-send -u normal "Pipeline is complete!"

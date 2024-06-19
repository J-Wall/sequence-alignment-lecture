#!/bin/bash
set -euxo pipefail

rm -f reads.fasta

for n in $(basename -a -s .wav reads/*.wav)
do
  echo '>'$n | tee -a reads.fasta
  whisper reads/$n.wav \
      --model tiny.en -f srt --output_dir .srt --condition_on_previous_text False --fp16 False \
      | cut -f5- -d' ' | sed -E 's/([A-Z])/\L\1/g;s/ /_/g;s/[^a-z_]//g;s/__+/_/g' | tee -a reads.fasta
done

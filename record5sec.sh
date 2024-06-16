#!/bin/bash
# Usage: record5sec.sh NAME

mkdir -p reads
rm -f reads/$1.wav
ffmpeg -f pulse -i 0 -t 5 reads/$1.wav

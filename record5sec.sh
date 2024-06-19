#!/bin/bash
# Usage: record5sec.sh NAME
# Need to edit this script to use the right audio device.
# run `arecord -l` to see what's available

mkdir -p reads
rm -f reads/$1.wav
ffmpeg -f pulse -i 0 -t 5 reads/$1.wav

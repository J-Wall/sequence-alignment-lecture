#!/usr/bin/env python
# Usage: consensus.py msa.fasta

from collections import Counter
import sys


def parse_fasta(fp):
    with open(fp, "r") as f:
        next(f)  # skip first header

        acc: list[str] = []
        for l in f:
            if l.startswith(">"):
                yield "".join(acc)
                acc = []
                continue

            acc.append(l.strip())

    yield "".join(acc)


def main():
    seqs = parse_fasta(sys.argv[1])
    seq = next(seqs)
    counts_per_pos = [Counter(c) for c in seq]

    for seq in seqs:
        for counter, c in zip(counts_per_pos, seq):
            counter[c] += 1

    consensus = "".join(counter.most_common(1)[0][0] for counter in counts_per_pos)
    print(">consensus", consensus, ">without_gaps", consensus.replace("-", ""), sep="\n")


if __name__ == "__main__":
    main()

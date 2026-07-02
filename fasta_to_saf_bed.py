#!/usr/bin/env python3

import sys
from pathlib import Path

if len(sys.argv) < 2:
    print("Usage: python fasta_to_saf_bed.py <transcript.fa>")
    sys.exit(1)

fasta_file = sys.argv[1]
base = Path(fasta_file).stem

saf_out = base + ".saf"
bed_out = base + ".bed"

def parse_fasta_headers(fasta_path):
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                header = line[1:].strip()
                parts = header.split()
                gene = parts[0]  # e.g. 4933401J01Rik
                # Extract length field (format: length:1070)
                length_field = [x for x in parts if x.startswith("length:")]
                if length_field:
                    length = int(length_field[0].split(":")[1])
                else:
                    length = None
                yield gene, length

with open(saf_out, "w") as saf, open(bed_out, "w") as bed:
    for gene, length in parse_fasta_headers(fasta_file):
        if length is None:
            continue  # skip if no length info
        # SAF: GeneID=1:length  GeneID  1  length  .
        saf.write(f"{gene}=1:{length}\t{gene}\t1\t{length}\t.\n")
        # BED: GeneID  0  length  GeneID=0:length  .  +
        bed.write(f"{gene}\t0\t{length}\t{gene}=0:{length}\t.\t+\n")

print(f"✅ Generated {saf_out} and {bed_out}")

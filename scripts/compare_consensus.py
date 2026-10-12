#!/usr/bin/env python3
import csv
import sys

BASES = set("ACGT")

def read_fasta(path):
    records = {}
    name, seq = None, []

    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(seq).upper()
                name = line[1:].split()[0]  # first part of FASTA header = sample ID
                seq = []
            else:
                seq.append(line)

    if name is not None:
        records[name] = "".join(seq).upper()

    return records

if len(sys.argv) != 3:
    sys.exit("Usage: python compare_consensus.py pipeline_A.fasta pipeline_B.fasta")

a = read_fasta(sys.argv[1])
b = read_fasta(sys.argv[2])

with open("consensus_comparison.csv", "w", newline="") as out, \
     open("snp_positions.csv", "w", newline="") as positions:

    writer = csv.writer(out)
    pos_writer = csv.writer(positions)

    writer.writerow([
        "sample", "length_A", "length_B",
        "Ns_A", "Ns_B", "delta_Ns_B_minus_A",
        "unambiguous_SNPs"
    ])
    pos_writer.writerow(["sample", "position_1_based", "base_A", "base_B"])

    for sample in sorted(set(a) & set(b)):
        seq_a, seq_b = a[sample], b[sample]

        if len(seq_a) != len(seq_b):
            writer.writerow([
                sample, len(seq_a), len(seq_b),
                seq_a.count("N"), seq_b.count("N"),
                seq_b.count("N") - seq_a.count("N"),
                "length_mismatch"
            ])
            continue

        snps = []
        for i, (base_a, base_b) in enumerate(zip(seq_a, seq_b), start=1):
            if base_a in BASES and base_b in BASES and base_a != base_b:
                snps.append((i, base_a, base_b))

        writer.writerow([
            sample,
            len(seq_a),
            len(seq_b),
            seq_a.count("N"),
            seq_b.count("N"),
            seq_b.count("N") - seq_a.count("N"),
            len(snps)
        ])

        for position, base_a, base_b in snps:
            pos_writer.writerow([sample, position, base_a, base_b])

    missing_a = sorted(set(b) - set(a))
    missing_b = sorted(set(a) - set(b))

    if missing_a:
        print("Only in pipeline B:", ", ".join(missing_a))
    if missing_b:
        print("Only in pipeline A:", ", ".join(missing_b))
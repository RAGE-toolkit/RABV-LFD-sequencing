#!/usr/bin/env python3

import csv
import re
import sys
from collections import defaultdict

BASES = set("ACGT")


def read_fasta(path):
    records = {}
    header = None
    sequence = []

    with open(path) as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    records[header] = "".join(sequence).upper()

                # Use the first word in the FASTA header as the sequence ID
                header = line[1:].split()[0]
                sequence = []

            else:
                sequence.append(line)

    if header is not None:
        records[header] = "".join(sequence).upper()

    return records


def count_ns(sequence):
    """Return leading, trailing, and internal N counts."""

    stripped = sequence.strip("N")

    if not stripped:
        return len(sequence), 0, 0

    leading = len(sequence) - len(sequence.lstrip("N"))
    trailing = len(sequence) - len(sequence.rstrip("N"))
    internal = stripped.count("N")

    return leading, trailing, internal


def parse_header(header):
    """
    Accepted formats:

    S4-D0
        -> sample S4, Brain, no temperature, day 0

    S4-FTA-RT-D1
        -> sample S4, FTA, RT, day 1

    S4-FTA-37-D1
        -> sample S4, FTA, 37, day 1

    S4-RDT-RT-D1
        -> sample S4, RDT, RT, day 1

    S4-RDT-37-D1
        -> sample S4, RDT, 37, day 1

    RADDL4B-2026-016
        -> Brain sample, no temperature/day

    RADDL4B-2026-016-RDT
        -> RDT specimen, no temperature/day
    """

    # FTA/RDT + RT/37 + day
    specimen_with_temp_and_day = re.match(
        r"^(.+)[-_](FTA|RDT)[-_](RT|37)[-_]D(\d+)$",
        header
    )

    if specimen_with_temp_and_day:
        return (
            specimen_with_temp_and_day.group(1),
            specimen_with_temp_and_day.group(2),
            specimen_with_temp_and_day.group(3),
            int(specimen_with_temp_and_day.group(4)),
        )

    # FTA/RDT + day, without temperature
    specimen_with_day = re.match(
        r"^(.+)[-_](FTA|RDT)[-_]D(\d+)$",
        header
    )

    if specimen_with_day:
        return (
            specimen_with_day.group(1),
            specimen_with_day.group(2),
            "",
            int(specimen_with_day.group(3)),
        )

    # Brain baseline + D0
    baseline_with_day = re.match(
        r"^(.+)[-_]D0$",
        header
    )

    if baseline_with_day:
        return (
            baseline_with_day.group(1),
            "Brain",
            "",
            0,
        )

    # Field RDT with no temperature/day
    specimen_no_day = re.match(
        r"^(.+)[-_](FTA|RDT)$",
        header
    )

    if specimen_no_day:
        return (
            specimen_no_day.group(1),
            specimen_no_day.group(2),
            "",
            "",
        )

    # Anything else is treated as a Brain sample
    return header, "Brain", "", ""


if len(sys.argv) != 2:
    sys.exit(
        "Usage: python compare_sample_specimens.py consensus_sequences.fasta"
    )


records = read_fasta(sys.argv[1])
samples = defaultdict(list)


for header, sequence in records.items():

    sample_id, medium, temperature, day = parse_header(header)

    samples[sample_id].append({
        "header": header,
        "sequence": sequence,
        "medium": medium,
        "temperature": temperature,
        "day": day,
    })


with open("comparison_to_baseline.csv", "w", newline="") as out, \
     open("snp_positions_to_baseline.csv", "w", newline="") as snp_out:

    writer = csv.writer(out)
    snp_writer = csv.writer(snp_out)

    writer.writerow([
        "sample_id",
        "baseline",
        "specimen",
        "medium",
        "temperature",
        "day",
        "baseline_leading_Ns",
        "baseline_trailing_Ns",
        "baseline_internal_Ns",
        "specimen_leading_Ns",
        "specimen_trailing_Ns",
        "specimen_internal_Ns",
        "delta_internal_Ns",
        "unambiguous_SNPs",
    ])

    snp_writer.writerow([
        "sample_id",
        "baseline",
        "specimen",
        "medium",
        "temperature",
        "day",
        "position_1_based",
        "baseline_base",
        "specimen_base",
    ])

    for sample_id, specimens in sorted(samples.items()):

        baselines = [
            x for x in specimens
            if x["medium"] == "Brain"
        ]

        if len(baselines) != 1:
            print(
                f"Skipping {sample_id}: expected one brain baseline, "
                f"found {len(baselines)}."
            )
            continue

        baseline = baselines[0]
        baseline_seq = baseline["sequence"]

        for specimen in specimens:

            if specimen["medium"] == "Brain":
                continue

            specimen_seq = specimen["sequence"]

            # For now, require identical sequence lengths.
            if len(baseline_seq) != len(specimen_seq):
                print(
                    f"Skipping {specimen['header']}: sequence length differs "
                    f"from baseline {baseline['header']}."
                )
                continue

            snps = []

            for position, (base0, base1) in enumerate(
                zip(baseline_seq, specimen_seq),
                start=1,
            ):

                # Count only A/C/G/T differences as SNPs.
                # Differences involving N are ignored.
                if (
                    base0 in BASES
                    and base1 in BASES
                    and base0 != base1
                ):
                    snps.append(
                        (position, base0, base1)
                    )

            baseline_lead, baseline_trail, baseline_internal = count_ns(
                baseline_seq
            )

            specimen_lead, specimen_trail, specimen_internal = count_ns(
                specimen_seq
            )

            writer.writerow([
                sample_id,
                baseline["header"],
                specimen["header"],
                specimen["medium"],
                specimen["temperature"],
                specimen["day"],
                baseline_lead,
                baseline_trail,
                baseline_internal,
                specimen_lead,
                specimen_trail,
                specimen_internal,
                specimen_internal - baseline_internal,
                len(snps),
            ])

            for position, base0, base1 in snps:

                snp_writer.writerow([
                    sample_id,
                    baseline["header"],
                    specimen["header"],
                    specimen["medium"],
                    specimen["temperature"],
                    specimen["day"],
                    position,
                    base0,
                    base1,
                ])
from Bio import SeqIO

for record in SeqIO.parse(
    "/Volumes/KB/PHL_September_15_2026/outputs/run1_complete/consensus.fa",
    "fasta"
):
    if record.id.startswith("S4"):
        print(record.id, len(record.seq))
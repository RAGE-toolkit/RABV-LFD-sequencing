from pathlib import Path
import pandas as pd

files = sorted(Path(".").rglob("*.amplicon_depths.tsv"))

tables = []

for file in files:
    sample = file.name.removesuffix(".amplicon_depths.tsv")

    # Skip NTC and NEC samples
    if sample.startswith(("NTC", "NEC")):
        continue

    df = pd.read_csv(file, sep="\t")
    df = df.rename(columns={"mean_depth": sample})
    tables.append(df)

combined = tables[0]

for df in tables[1:]:
    combined = combined.merge(df, on=["chrom", "amplicon"], how="outer")

combined = combined.sort_values(["chrom", "amplicon"])
combined.to_csv("combined_amplicon_depths.tsv", sep="\t", index=False)

print(f"Combined {len(tables)} samples, excluding NTC and NEC.")
print("Saved combined_amplicon_depths.tsv")
from pathlib import Path
import pandas as pd

# Find all post-normalisation coverage files recursively
files = sorted(
    Path(".").rglob("*.post-normalisation.coverage.tsv")
)

if not files:
    raise FileNotFoundError(
        "No *.post-normalisation.coverage.tsv files found"
    )

combined = None
excluded = 0
included = 0

for file in files:

    # Sample ID = filename without suffix
    sample_id = file.name.replace(
        ".post-normalisation.coverage.tsv",
        ""
    )

    # Skip NTC and NEC samples
    if sample_id.startswith(("NTC", "NEC")):
        excluded += 1
        continue

    # Read coverage file
    df = pd.read_csv(file, sep="\t")

    # Keep common coordinates and rename depth column
    df = df[["chrom", "pos", "depth"]].rename(
        columns={"depth": sample_id}
    )

    # Merge samples using chrom + pos
    if combined is None:
        combined = df
    else:
        combined = combined.merge(
            df,
            on=["chrom", "pos"],
            how="outer"
        )

    included += 1

if combined is None:
    raise ValueError("No samples remaining after excluding NTC and NEC.")

# Sort by genomic position
combined = combined.sort_values(["chrom", "pos"])

# Write combined output
output = Path("all_seq_combined_post_normalisation_coverage.tsv")

combined.to_csv(
    output,
    sep="\t",
    index=False
)

print(f"Total files found: {len(files)}")
print(f"Samples included: {included}")
print(f"NTC/NEC samples excluded: {excluded}")
print(f"Output: {output}")
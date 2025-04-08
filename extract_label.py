from Bio import SeqIO
import pandas as pd
import os
from pathlib import Path


# Set your base data folder
BASE_DIR = Path("/nfs_share/students/jinhyun/covid")
variant_dirs = ["wildtype", "alpha", "beta", "delta", "omicron"]
standard_dir = BASE_DIR / "standard"


annotations = []

for variant in variant_dirs:
    variant_path = BASE_DIR / variant

    if not variant_path.is_dir():
        print(f"⚠️ Skipping non-existent folder: {variant_path}")
        continue

    fasta_files = list(variant_path.glob("*.fasta")) + list(variant_path.glob("*.fa"))
    print(f"🔍 Found {len(fasta_files)} FASTA files in {variant}/")

    
    for fasta_file in fasta_files:
        count = 0
        full_path = os.path.join(variant_path, fasta_file)
        for record in SeqIO.parse(full_path, "fasta"):
            sample_id = record.id
            annotations.append({
                "sample_id": sample_id,
                "target": variant,
                "fasta_path": full_path
            })
            count += 1
        print(f"✅ Parsed {count} sequences from {fasta_file.name}")


# ✅ Manually map standard files to variant labels
standard_file_mapping = {
    "OL689430.1.fasta": "alpha",   # B.1.1.7
    "MZ433432.1.fasta": "beta",    # B.1.351
    "OK091006.1.fasta": "delta",   # B.1.617.2
    "OW998817.1.fasta": "omicron"  # BA.1
}

# Annotate standard files using mapping
for filename, label in standard_file_mapping.items():
    full_path = os.path.join(standard_dir, filename)
    for record in SeqIO.parse(full_path, "fasta"):
        sample_id = record.id
        annotations.append({
            "sample_id": sample_id,
            "target": label,
            "fasta_path": full_path
        })

# Convert to DataFrame
df = pd.DataFrame(annotations)

# Check for duplicates
duplicates = df[df.duplicated("sample_id", keep=False)]
if not duplicates.empty:
    print("⚠️ Duplicate sample IDs found:")
    print(duplicates.sort_values("sample_id"))
else:
    print("✅ No duplicate sample IDs found.")

# Save to Excel
df.to_excel("results/sample_annotations.xlsx", index=False)
print("✅ Annotation file saved to 'results/sample_annotations.xlsx'")

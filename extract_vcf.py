from Bio import AlignIO
import pandas as pd

alignment_ref = AlignIO.read("/nfs_share/students/jinhyun/covid/reference/NC_045512.2.fasta", "fasta")
alignment = AlignIO.read("results/aligned.fasta", "fasta")
reference = str(alignment[0].seq)  # Wuhan reference
variant_dict = {}

assert len(alignment_ref[0].seq) == len(reference)
print("Reference length:", len(reference))
print("Number of samples:", len(alignment)-1)
for record in alignment[1:]:
    sample_id = record.id
    sample_seq = str(record.seq)
    sample_variants = []
    for pos, (ref_base, sample_base) in enumerate(zip(reference, sample_seq)):
        ref_base = ref_base.upper()
        sample_base = sample_base.upper()

        if sample_base == 'N':
            continue

        if ref_base not in {'A', 'T', 'C', 'G', '-'}:
            print(f"Unknown base found at ref: {ref_base}")
            raise
        if sample_base not in {'A', 'T', 'C', 'G', 'Y', 'W', 'R', 'M', 'K', 'S', '-'}:
            print(f"Unknown base found at sample: {sample_base}")
            raise

        if ref_base != sample_base:
            if sample_base == '-':
                variant = f"{pos+1}_{ref_base}>DEL"
            elif ref_base == '-':
                variant = f"{pos+1}_INS"
            else:
                variant = f"{pos+1}_{ref_base}>{sample_base}"
            sample_variants.append(variant)

    variant_dict[sample_id] = sample_variants

def extract_pos(variant_str):
    return int(variant_str.split('_')[0])
all_variants = sorted(set(var for variant_list in variant_dict.values() for var in variant_list), key=extract_pos)

variant_matrix = pd.DataFrame(0, index=variant_dict.keys(), columns=all_variants)

for sample_id, variants in variant_dict.items():
    variant_matrix.loc[sample_id, variants] = 1

variant_matrix = variant_matrix.reset_index()
variant_matrix = variant_matrix.rename(columns={'index': 'sample_id'})

variant_matrix.to_excel("results/variant_matrix.xlsx", index=False)
print("✅ 'results/variant_matrix.xlsx' saved successfully.")

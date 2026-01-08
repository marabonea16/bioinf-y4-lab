"""
Exercise 10.2 — Identify top SNP–Gene correlations

TODO:
- încărcați matricea integrată multi-omics
- împărțiți rândurile în SNPs vs gene (după indice sau după nume)
- calculați corelații între fiecare SNP și fiecare genă
- filtrați |r| > 0.5
- exportați snp_gene_pairs_<handle>.csv
"""

from pathlib import Path
import pandas as pd

HANDLE = "marabonea16"
JOINT_CSV = Path(f"labs/10_integrative/submissions/{HANDLE}/multiomics_concat_{HANDLE}.csv")

OUT_CSV = Path(f"labs/10_integrative/submissions/{HANDLE}/snp_gene_pairs_{HANDLE}.csv")

# Load the integrated multi-omics matrix
df_joint = pd.read_csv(JOINT_CSV, index_col=0)

print(f"Loaded joint matrix with shape: {df_joint.shape}")
print(f"Features: {df_joint.index.tolist()}")
print(f"Samples: {df_joint.columns.tolist()}")

# Split rows into SNPs vs Genes based on feature names
snp_rows = [idx for idx in df_joint.index if idx.startswith('SNP')]
gene_rows = [idx for idx in df_joint.index if idx.startswith('GENE')]

print(f"\nFound {len(snp_rows)} SNPs and {len(gene_rows)} genes")

df_snp = df_joint.loc[snp_rows]
df_gene = df_joint.loc[gene_rows]

# Compute correlation matrix between all features
# Transpose so samples are rows, then compute correlation matrix
corr_matrix = df_joint.T.corr()

# Extract SNP-gene correlations (exclude SNP-SNP and gene-gene correlations)
snp_gene_corr = corr_matrix.loc[snp_rows, gene_rows]

print(f"\nSNP-Gene correlation matrix shape: {snp_gene_corr.shape}")
print("\nCorrelation matrix:")
print(snp_gene_corr)

# Filter pairs with |r| > 0.5
pairs = []
for snp in snp_rows:
    for gene in gene_rows:
        r = snp_gene_corr.loc[snp, gene]
        if abs(r) > 0.5:
            pairs.append({
                'SNP': snp,
                'Gene': gene,
                'Correlation': r,
                'Abs_Correlation': abs(r)
            })

# Convert to DataFrame and sort by absolute correlation
df_pairs = pd.DataFrame(pairs)
df_pairs = df_pairs.sort_values('Abs_Correlation', ascending=False)

# Remove the helper column before saving
df_pairs = df_pairs[['SNP', 'Gene', 'Correlation']]

# Export results
df_pairs.to_csv(OUT_CSV, index=False)

print(f"\n{'='*60}")
print(f"Found {len(df_pairs)} SNP-gene pairs with |r| > 0.5")
print(f"{'='*60}")
print("\nTop SNP-Gene pairs:")
print(df_pairs.to_string(index=False))
print(f"\nResults saved to: {OUT_CSV}")

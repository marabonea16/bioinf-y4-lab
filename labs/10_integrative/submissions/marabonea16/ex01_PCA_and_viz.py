"""
Exercise 10 — PCA Single-Omics vs Joint

TODO:
- încărcați SNP și Expression
- normalizați fiecare strat (z-score)
- rulați PCA pe:
    1) strat SNP
    2) strat Expression
    3) strat Joint (concat)
- generați 3 figuri PNG
- comparați vizual distribuția probelor
"""

from pathlib import Path
import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

HANDLE = "marabonea16"

SNP_CSV = Path(f"data/work/{HANDLE}/lab10/snp_matrix_{HANDLE}.csv")
EXP_CSV = Path(f"data/work/{HANDLE}/lab10/expression_matrix_{HANDLE}.csv")

OUT_DIR = Path(f"labs/10_integrative/submissions/{HANDLE}")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Load data
df_snp = pd.read_csv(SNP_CSV, index_col=0)
df_exp = pd.read_csv(EXP_CSV, index_col=0)

# Both datasets already have samples as columns, no transpose needed

# Find common samples
common_samples = df_snp.columns.intersection(df_exp.columns)
print(f"Common samples: {len(common_samples)}")

df_snp = df_snp[common_samples]
df_exp = df_exp[common_samples]

# Z-score normalization (per feature, across samples)
df_snp_norm = df_snp.sub(df_snp.mean(axis=1), axis=0).div(df_snp.std(axis=1), axis=0)
df_exp_norm = df_exp.sub(df_exp.mean(axis=1), axis=0).div(df_exp.std(axis=1), axis=0)

# Concatenate features
df_joint = pd.concat([df_snp_norm, df_exp_norm], axis=0)

# Save concatenated multiomics data
df_joint.to_csv(OUT_DIR / f"multiomics_concat_{HANDLE}.csv")
print(f"Saved: {OUT_DIR / f'multiomics_concat_{HANDLE}.csv'}")

# =============================================================================
# PCA Analysis
# =============================================================================

# 1. PCA on SNP layer only
pca_snp = PCA(n_components=2)
proj_snp = pca_snp.fit_transform(df_snp_norm.T)

plt.figure(figsize=(8, 6))
plt.scatter(proj_snp[:, 0], proj_snp[:, 1], c='red', s=100, alpha=0.6)
for i, sample in enumerate(common_samples):
    plt.annotate(sample, (proj_snp[i, 0], proj_snp[i, 1]), fontsize=8)
plt.title(f"PCA on SNP Layer\nExplained variance: {pca_snp.explained_variance_ratio_}")
plt.xlabel(f"PC1 ({pca_snp.explained_variance_ratio_[0]:.2%})")
plt.ylabel(f"PC2 ({pca_snp.explained_variance_ratio_[1]:.2%})")
plt.tight_layout()
plt.savefig(OUT_DIR / f"pca_snp_{HANDLE}.png", dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved: {OUT_DIR / f'pca_snp_{HANDLE}.png'}")

# 2. PCA on Expression layer only
pca_exp = PCA(n_components=2)
proj_exp = pca_exp.fit_transform(df_exp_norm.T)

plt.figure(figsize=(8, 6))
plt.scatter(proj_exp[:, 0], proj_exp[:, 1], c='green', s=100, alpha=0.6)
for i, sample in enumerate(common_samples):
    plt.annotate(sample, (proj_exp[i, 0], proj_exp[i, 1]), fontsize=8)
plt.title(f"PCA on Expression Layer\nExplained variance: {pca_exp.explained_variance_ratio_}")
plt.xlabel(f"PC1 ({pca_exp.explained_variance_ratio_[0]:.2%})")
plt.ylabel(f"PC2 ({pca_exp.explained_variance_ratio_[1]:.2%})")
plt.tight_layout()
plt.savefig(OUT_DIR / f"pca_expr_{HANDLE}.png", dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved: {OUT_DIR / f'pca_expr_{HANDLE}.png'}")

# 3. PCA on Joint (integrated) data
pca_joint = PCA(n_components=2)
proj_joint = pca_joint.fit_transform(df_joint.T)

plt.figure(figsize=(8, 6))
plt.scatter(proj_joint[:, 0], proj_joint[:, 1], c='blue', s=100, alpha=0.6)
for i, sample in enumerate(common_samples):
    plt.annotate(sample, (proj_joint[i, 0], proj_joint[i, 1]), fontsize=8)
plt.title(f"PCA on Joint Multi-Omics\nExplained variance: {pca_joint.explained_variance_ratio_}")
plt.xlabel(f"PC1 ({pca_joint.explained_variance_ratio_[0]:.2%})")
plt.ylabel(f"PC2 ({pca_joint.explained_variance_ratio_[1]:.2%})")
plt.tight_layout()
plt.savefig(OUT_DIR / f"pca_joint_{HANDLE}.png", dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved: {OUT_DIR / f'pca_joint_{HANDLE}.png'}")

# =============================================================================
# Cross-Omics Correlation
# =============================================================================

# Compute correlation between each SNP and each gene across samples
corr_matrix = pd.concat([df_snp_norm, df_exp_norm], axis=0).T.corr()

# Extract SNP-gene correlations (exclude SNP-SNP and gene-gene)
snp_rows = [idx for idx in df_joint.index if idx.startswith('SNP')]
gene_rows = [idx for idx in df_joint.index if idx.startswith('GENE')]

snp_gene_corr = corr_matrix.loc[snp_rows, gene_rows]

# Filter pairs with |r| > 0.5
pairs = []
for snp in snp_rows:
    for gene in gene_rows:
        r = snp_gene_corr.loc[snp, gene]
        if abs(r) > 0.5:
            pairs.append({'SNP': snp, 'Gene': gene, 'Correlation': r})

df_pairs = pd.DataFrame(pairs)
df_pairs = df_pairs.sort_values('Correlation', key=abs, ascending=False)
df_pairs.to_csv(OUT_DIR / f"snp_gene_pairs_{HANDLE}.csv", index=False)
print(f"\nSaved: {OUT_DIR / f'snp_gene_pairs_{HANDLE}.csv'}")
print(f"\nFound {len(df_pairs)} SNP-gene pairs with |r| > 0.5:")
print(df_pairs)

print("\n=== Exercise completed successfully! ===")

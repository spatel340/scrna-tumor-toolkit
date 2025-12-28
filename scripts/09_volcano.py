import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

DE = Path("results/tables/de_malignant_vs_immune_cluster_pseudobulk.tsv")
OUTDIR = Path("figures/de")
OUTDIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(DE, sep="\t").dropna(subset=["padj", "logFC_malignant_minus_immune"])

# -log10(padj) with floor to avoid inf
df["neglog10_padj"] = -np.log10(np.clip(df["padj"].values, 1e-300, 1.0))

# thresholds
padj_thr = 0.05
lfc_thr = 0.5

df["sig"] = (df["padj"] < padj_thr) & (np.abs(df["logFC_malignant_minus_immune"]) > lfc_thr)

# label top genes
top_up = df[(df["padj"] < padj_thr) & (df["logFC_malignant_minus_immune"] > lfc_thr)].head(8)
top_dn = df[(df["padj"] < padj_thr) & (df["logFC_malignant_minus_immune"] < -lfc_thr)].head(8)
label_df = pd.concat([top_up, top_dn], axis=0)

plt.figure(figsize=(7.5, 6.5))
# background
plt.scatter(df["logFC_malignant_minus_immune"], df["neglog10_padj"], s=10, alpha=0.35, color="#9aa0a6", linewidths=0)
# significant
sig = df[df["sig"]]
plt.scatter(sig["logFC_malignant_minus_immune"], sig["neglog10_padj"], s=12, alpha=0.65, color="#d93025", linewidths=0)

# threshold lines
plt.axhline(-np.log10(padj_thr), color="black", lw=1, ls="--", alpha=0.6)
plt.axvline(lfc_thr, color="black", lw=1, ls="--", alpha=0.6)
plt.axvline(-lfc_thr, color="black", lw=1, ls="--", alpha=0.6)

# labels
for _, r in label_df.iterrows():
    plt.text(r["logFC_malignant_minus_immune"], r["neglog10_padj"], str(r["gene"]),
             fontsize=8, ha="left", va="bottom")

plt.title("Malignant vs Immune (cluster pseudobulk DE)")
plt.xlabel("logFC (malignant − immune)")
plt.ylabel("−log10(FDR)")
plt.tight_layout()

out_png = OUTDIR / "volcano_malignant_vs_immune.png"
plt.savefig(out_png, dpi=220)
print("Wrote:", out_png)

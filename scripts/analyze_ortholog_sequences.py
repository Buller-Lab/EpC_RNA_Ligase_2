import argparse
import os
from Bio import SeqIO
from Bio.Align import PairwiseAligner
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from scipy import stats
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform
# Publication style: Arial, palette #E4E8F0 -> #4B6087, #262626 for text and lines.
PALETTE = ["#E4E8F0", "#BFC9DB", "#8094B7", "#657DA7", "#4B6087"]
INK = "#262626"
STYLE_RC = {
    "font.family": "Arial",
    "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
    "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": INK,
    "xtick.color": INK, "ytick.color": INK, "lines.color": INK, "patch.edgecolor": INK,
}
HEATMAP_CMAP = LinearSegmentedColormap.from_list("publication_blue", PALETTE)
parser = argparse.ArgumentParser()
parser.add_argument("fasta", help="Input FASTA file")
parser.add_argument("--outdir", default=".", help="Output directory (default: current directory)")
args = parser.parse_args()
os.makedirs(args.outdir, exist_ok=True)

def out(name):
    return os.path.join(args.outdir, name)
records = list(SeqIO.parse(args.fasta, "fasta"))
ids = [r.id for r in records]
# Friendly figure labels for special entries; the identity matrix keeps raw IDs.
DISPLAY_OVERRIDES = {"Rnl2_27": "EpC (Rnl2_27)"}
display_ids = [DISPLAY_OVERRIDES.get(i, i) for i in ids]
seqs = [str(r.seq) for r in records]
n = len(seqs)
aligner = PairwiseAligner()
aligner.mode = "global"
matrix = np.zeros((n, n))
for i in range(n):
    for j in range(i, n):
        if i == j:
            matrix[i, j] = 100.0
            continue
        alignment = aligner.align(seqs[i], seqs[j])[0]
        alen = alignment.shape[1]
        identities = alignment.counts().identities
        pid = identities / alen * 100 if alen > 0 else 0.0
        matrix[i, j] = pid
        matrix[j, i] = pid
df = pd.DataFrame(matrix, index=ids, columns=ids)
df.to_excel(out("identity_matrix.xlsx"))
sns.set_theme(style="white", rc=STYLE_RC)
figsize = (max(10, n * 0.6), max(8, n * 0.6))
plt.figure(figsize=figsize)
ax = sns.heatmap(df, cmap=HEATMAP_CMAP, annot=True, fmt=".1f", square=True, linewidths=0.5, cbar_kws={"label": "Identity (%)", "shrink": 0.8}, vmin=0, vmax=100, annot_kws={"size": 9}, xticklabels=display_ids, yticklabels=display_ids)
plt.title("Pairwise Protein Sequence Identity (%)", fontsize=14, pad=20)
plt.xticks(rotation=90, ha="right", fontsize=9)
plt.yticks(rotation=0, fontsize=9)
plt.tight_layout()
plt.savefig(out("sequence_identity_heatmap.png"), dpi=300, bbox_inches="tight")
plt.savefig(out("sequence_identity_heatmap.svg"), bbox_inches="tight")
plt.close()
lengths = np.array([len(s) for s in seqs])
id_values = matrix[np.triu_indices(n, k=1)] if n > 1 else np.array([100.0])
total_pairs = len(id_values) if n > 1 else 0
length_min = np.min(lengths)
length_max = np.max(lengths)
length_mean = np.mean(lengths)
length_median = np.median(lengths)
length_std = np.std(lengths)
length_skew = stats.skew(lengths)
length_kurt = stats.kurtosis(lengths)
if n > 1:
    identity_min = np.min(id_values)
    identity_max = np.max(id_values)
    identity_mean = np.mean(id_values)
    identity_median = np.median(id_values)
    identity_std = np.std(id_values)
    identity_skew = stats.skew(id_values)
    identity_kurt = stats.kurtosis(id_values)
    above_30 = int(np.sum(id_values > 30))
    above_50 = int(np.sum(id_values > 50))
    above_70 = int(np.sum(id_values > 70))
    above_90 = int(np.sum(id_values > 90))
else:
    identity_min = identity_max = identity_mean = identity_median = identity_std = identity_skew = identity_kurt = 100.0
    above_30 = above_50 = above_70 = above_90 = 0
plt.figure(figsize=(10, 6))
sns.histplot(lengths, bins=min(30, max(5, n//2)), kde=True, color=PALETTE[2], edgecolor=INK, line_kws={"color": INK})
plt.title("Sequence Length Distribution")
plt.xlabel("Length (amino acids)")
plt.ylabel("Frequency")
plt.tight_layout()
plt.savefig(out("sequence_length_distribution.png"), dpi=300, bbox_inches="tight")
plt.savefig(out("sequence_length_distribution.svg"), bbox_inches="tight")
plt.close()
if n > 1:
    distance_matrix = 100 - matrix
    condensed_dist = squareform(distance_matrix)
    linkage_mat = linkage(condensed_dist, method="average")
    plt.figure(figsize=(max(12, n * 0.35), max(8, n * 0.25)))
    dendrogram(linkage_mat, labels=display_ids, orientation="left", leaf_font_size=max(5, 9 - n//8),
               color_threshold=0, above_threshold_color=INK)
    plt.title("UPGMA Phylogenetic Dendrogram (Distance = 100 - % Identity)")
    plt.xlabel("Distance")
    plt.tight_layout()
    plt.savefig(out("phylogenetic_dendrogram.png"), dpi=300, bbox_inches="tight")
    plt.savefig(out("phylogenetic_dendrogram.svg"), bbox_inches="tight")
    plt.close()
else:
    plt.figure(figsize=(6, 4))
    plt.text(0.5, 0.5, "Single sequence - no tree", ha="center", va="center", fontsize=14)
    plt.axis("off")
    plt.savefig(out("phylogenetic_dendrogram.png"), dpi=300, bbox_inches="tight")
    plt.savefig(out("phylogenetic_dendrogram.svg"), bbox_inches="tight")
    plt.close()
total_residues = int(np.sum(lengths))
pct_above_50 = (above_50 / total_pairs * 100) if total_pairs > 0 else 0
report = f"""Protein Sequence Pairwise Identity Analysis Report
Generated on: {__import__('datetime').datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

Input FASTA file: {args.fasta}
Number of protein sequences: {n}
Total amino acid residues: {total_residues}

SEQUENCE LENGTH STATISTICS
Min length: {length_min} aa
Max length: {length_max} aa
Mean length: {length_mean:.2f} aa
Median length: {length_median:.2f} aa
Standard deviation: {length_std:.2f} aa
Skewness: {length_skew:.4f}
Kurtosis: {length_kurt:.4f}

PAIRWISE IDENTITY STATISTICS (off-diagonal pairs: {total_pairs})
Min identity: {identity_min:.2f}%
Max identity: {identity_max:.2f}%
Mean identity: {identity_mean:.2f}%
Median identity: {identity_median:.2f}%
Standard deviation: {identity_std:.2f}%
Skewness: {identity_skew:.4f}
Kurtosis: {identity_kurt:.4f}

Threshold analysis:
Pairs >30% identity: {above_30} ({(above_30/total_pairs*100) if total_pairs>0 else 0:.1f}%)
Pairs >50% identity: {above_50} ({pct_above_50:.1f}%)
Pairs >70% identity: {above_70} ({(above_70/total_pairs*100) if total_pairs>0 else 0:.1f}%)
Pairs >90% identity: {above_90} ({(above_90/total_pairs*100) if total_pairs>0 else 0:.1f}%)

OUTPUT FILES GENERATED
- identity_matrix.xlsx          : Complete pairwise % identity matrix (Excel)
- sequence_identity_heatmap.png/.svg : Heatmap with publication blue colormap
- sequence_length_distribution.png/.svg : Histogram + KDE of sequence lengths
- phylogenetic_dendrogram.png/.svg   : UPGMA hierarchical clustering tree (distance-based)
- description.txt               : This full analysis report

All computations used global Needleman-Wunsch alignment with match=1, mismatch=0.
Identity % = (identical residues / alignment length) x 100.
Phylogenetic tree constructed via average-linkage (UPGMA) on (100 - identity) distances.
"""
with open(out("description.txt"), "w") as f:
    f.write(report)
print(f"Analysis complete. Files written to {args.outdir}: identity_matrix.xlsx, sequence_identity_heatmap.png/.svg, sequence_length_distribution.png/.svg, phylogenetic_dendrogram.png/.svg, description.txt")
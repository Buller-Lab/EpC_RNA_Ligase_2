import argparse
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

# Publication style: Arial, palette #E4E8F0 -> #4B6087, #262626 for text and lines.
PALETTE = ["#E4E8F0", "#BFC9DB", "#8094B7", "#657DA7", "#4B6087"]
INK = "#262626"
STYLE_RC = {
    "font.family": "Arial",
    "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
    "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": INK,
    "xtick.color": INK, "ytick.color": INK, "lines.color": INK, "patch.edgecolor": INK,
}
mpl.rcParams.update(STYLE_RC)

parser = argparse.ArgumentParser(
    description="Plot per-residue RMSF (mean +/- SD over replicate rows) for a single dataset."
)
parser.add_argument("-i", "--input", required=True,
                    help="CSV with one 'replicate N' row per replicate and one column per residue number")
parser.add_argument("-o", "--output", default=None,
                    help="Output base path without extension (default: next to the input CSV)")
parser.add_argument("--ylabel", default="RMSF (nm)", help="Y-axis label (default: 'RMSF (nm)')")
args = parser.parse_args()

df = pd.read_csv(args.input, index_col=0)
# Keep only the replicate rows; precomputed summary rows (average, stdev) are ignored.
df = df[df.index.str.startswith("replicate")]
residues = df.columns.astype(int)
mean = df.mean(axis=0).to_numpy()
std = df.std(axis=0, ddof=0).to_numpy()  # population SD, matches the 'stdev' row of the input

output_path = args.output or os.path.splitext(args.input)[0] + "_plot"

fig, ax = plt.subplots(figsize=(10, 4), facecolor="white")
ax.fill_between(residues, mean - std, mean + std, color=PALETTE[1], linewidth=0, label="± SD")
ax.plot(residues, mean, color=PALETTE[4], linewidth=1.8, label=f"Mean (n = {len(df)})")

ax.set_xlim(residues.min(), residues.max())
ax.set_ylim(bottom=0)
ax.set_xlabel("Residue number", fontsize=12, labelpad=8, weight="bold")
ax.set_ylabel(args.ylabel, fontsize=12, labelpad=8, weight="bold")
for tick in ax.get_xticklabels() + ax.get_yticklabels():
    tick.set_fontweight("bold")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.legend(loc="upper left", frameon=False, fontsize=10)

plt.tight_layout()
plt.savefig(f"{output_path}.png", format="png", dpi=450, bbox_inches="tight", facecolor="white")
plt.savefig(f"{output_path}.svg", format="svg", bbox_inches="tight", facecolor="white")
plt.close(fig)
print(f"Saved RMSF plot as: {output_path}.png and {output_path}.svg")

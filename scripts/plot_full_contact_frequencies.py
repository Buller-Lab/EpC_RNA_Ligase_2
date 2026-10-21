import argparse
import os
import re
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch

# Publication style: Arial, palette #E4E8F0 -> #4B6087, #262626 for text and lines.
PALETTE = ["#E4E8F0", "#BFC9DB", "#8094B7", "#657DA7", "#4B6087"]
INK = "#262626"
STYLE_RC = {
    "font.family": "Arial",
    "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
    "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": INK,
    "xtick.color": INK, "ytick.color": INK, "lines.color": INK, "patch.edgecolor": INK,
}

parser = argparse.ArgumentParser(
    description="Plot the full per-contact frequency analysis for every simulation "
                "system (one landscape subplot per sheet, stacked vertically)."
)
parser.add_argument("-i", "--input", required=True, help="compact_vs_loose xlsx file")
parser.add_argument("-n", "--top", type=int, default=35,
                    help="keep only the top N residues (by mean frequency) per system "
                         "(0 = all; default 35)")
args = parser.parse_args()

aa_3_to_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
}

# Colors taken from panel b) of the publication figure (dummy.png) for the
# three shared categories; two distinct muted tones added for the categories
# that only appear in the full set.
COLOR_MAP = {
    "conserved": "#8094b7",   # blue   (panel b)
    "tested": "#dfb295",      # orange (panel b)
    "consensus": "#569f85",   # green  (panel b)
    "potential": "#a68cb7",   # purple (added)
    "TPW": "#d6a3a3",         # red    (added)
}

LEGEND_LABELS = {
    "conserved": "highly conserved (>95% occurence)",
    "tested": "already tested (previous rounds)",
    "consensus": "consensus (in compact & loose)",
    "potential": "potential (new candidate)",
    "TPW": "R2.1 (S57T_S62P_N65W)",
}

# Legend / stacking order; only categories actually present are drawn.
ANNOT_ORDER = ["conserved", "tested", "consensus", "potential", "TPW"]


def format_residue(part):
    """GLU205 -> E205, nucleotides such as U22 are kept as-is."""
    part = part.strip()
    match = re.match(r"^([A-Z]{3})(\d+)$", part, re.IGNORECASE)
    if match and match.group(1).upper() in aa_3_to_1:
        return aa_3_to_1[match.group(1).upper()] + match.group(2)
    return part


def format_contact(contact_str):
    """"[GLU205, U22]" -> "E205": keep only the amino-acid partner, drop the nucleotide."""
    cleaned = contact_str.replace("[", "").replace("]", "")
    for raw in cleaned.split(","):
        raw = raw.strip()
        match = re.match(r"^([A-Z]{3})(\d+)$", raw, re.IGNORECASE)
        if match and match.group(1).upper() in aa_3_to_1:
            return format_residue(raw)
    # Fallback: no amino acid found, keep the formatted pair.
    return "–".join(format_residue(p) for p in cleaned.split(",") if p.strip())


def load_sheet(xls, sheet_name):
    df = pd.read_excel(xls, sheet_name=sheet_name)
    df = df.dropna(subset=["Contact", "Average"]).copy()
    df["Label"] = df["Contact"].apply(format_contact)
    df["Stddev"] = df["Stddev"].fillna(0)

    # One bar per residue: residues contacting more than one nucleotide are
    # averaged (mean contact frequency / mean stddev), as in the condensed plot.
    # Sorting first lets "first" pick the annotation of the strongest contact.
    df = df.sort_values(by="Average", ascending=False)
    df = df.groupby("Label", as_index=False).agg(
        Average=("Average", "mean"),
        Stddev=("Stddev", "mean"),
        Annotation=("Annotation", "first"),
        N=("Average", "size"),
    )
    df = df.sort_values(by="Average", ascending=False).reset_index(drop=True)
    if args.top and args.top > 0:
        df = df.head(args.top).reset_index(drop=True)
    return df


file_path = args.input
xls = pd.ExcelFile(file_path)
systems = [(name, load_sheet(xls, name)) for name in xls.sheet_names]

max_bars = max(len(df) for _, df in systems)
present_annots = set()
for _, df in systems:
    present_annots.update(df["Annotation"].dropna().unique())

# Full-width landscape figure, one row per simulation system stacked above each other.
fig_width = max(14.0, max_bars * 0.30 + 2.5)
row_height = 4.6
fig_height = row_height * len(systems)

sns.set_theme(style="white", font_scale=1.0, rc=STYLE_RC)
fig, axes = plt.subplots(
    len(systems), 1, figsize=(fig_width, fig_height), facecolor="white"
)
if len(systems) == 1:
    axes = [axes]

for ax, (name, df) in zip(axes, systems):
    x_positions = np.arange(len(df))
    bar_colors = [COLOR_MAP.get(a, PALETTE[0]) for a in df["Annotation"]]

    ax.bar(
        x_positions,
        df["Average"],
        yerr=df["Stddev"],
        capsize=2,
        color=bar_colors,
        edgecolor=INK,
        linewidth=0.6,
        error_kw={"ecolor": INK, "linewidth": 0.6},
    )

    ax.set_ylim(0, 1.1)
    ax.set_ylabel("Contact Frequency", size=12, weight="bold", labelpad=10)
    ax.set_title(f"{name} system", size=13, weight="bold", pad=8, loc="left")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["bottom"].set_linewidth(0.8)

    ax.set_xlim(-0.8, max_bars - 0.2)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(df["Label"], rotation=90, ha="center",
                       fontsize=10, fontweight="bold")
    ax.tick_params(axis="x", pad=2)
    ax.tick_params(axis="y", pad=3, labelsize=10)
    for label in ax.get_yticklabels():
        label.set_weight("bold")

legend_elements = [
    Patch(facecolor=COLOR_MAP[a], edgecolor=INK, linewidth=0.8, label=LEGEND_LABELS[a])
    for a in ANNOT_ORDER if a in present_annots
]

plt.tight_layout(pad=1.6, h_pad=2.4, rect=(0, 0, 1, 0.95))
# Horizontal legend above the top subplot so it never overlaps the bars.
fig.legend(
    handles=legend_elements, loc="upper center", bbox_to_anchor=(0.5, 1.0),
    ncol=len(legend_elements), frameon=False, fontsize=11,
    columnspacing=1.8, handletextpad=0.6,
)

input_dir = os.path.dirname(os.path.abspath(file_path))
base_name = os.path.splitext(os.path.basename(file_path))[0] + "_full_barplot"
output_path = os.path.join(input_dir, base_name)

plt.savefig(f"{output_path}.png", format="png", dpi=450, bbox_inches="tight", facecolor="white")
plt.savefig(f"{output_path}.svg", format="svg", bbox_inches="tight", facecolor="white")

plt.show()

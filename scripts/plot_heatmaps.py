#!/usr/bin/env python3
"""Plot an activity heatmap of RNA ligase 2 orthologs across assay conditions,
with rows ordered by a sequence-based (phylogenetic) clustering.

The input CSV (semicolon-separated by default) holds one row per enzyme, an
``Enzyme`` identifier column holding publication IDs (e.g. ``T4_Rnl2``,
``Rnl2_02`` ... ``Rnl2_30``; the reference Erwinia phage Cronus ligase is
``Rnl2_27``), one numeric column per assay condition, and an (empty)
``sequence`` column. Sequences are looked up by matching the publication ID
directly against the orthologs FASTA (which uses the same IDs).

Example:
    python scripts/plot_heatmaps.py evolutionary_analysis/Top15_rLI2_orthologs/top15_rLI2_orthologs_diff_conditions.csv
"""

import argparse
import os
import re
import sys
import shutil
import tempfile
import subprocess

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator
from scipy.cluster.hierarchy import dendrogram, linkage

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_FASTA = os.path.join(HERE, "..", "evolutionary_analysis", "rLI2_orthologs.fasta")

# Friendly display labels for special entries (data files keep the publication ID).
DISPLAY_OVERRIDES = {"Rnl2_27": "EpC (Rnl2_27)"}

VALUE_LABEL = "Conversion (%)"

# Publication style: Arial, palette #E4E8F0 -> #4B6087, #262626 for text and lines.
PALETTE = ["#E4E8F0", "#BFC9DB", "#8094B7", "#657DA7", "#4B6087"]
INK = "#262626"
mpl.rcParams.update({
    "font.family": "Arial",
    "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
    "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": INK,
    "xtick.color": INK, "ytick.color": INK, "lines.color": INK, "patch.edgecolor": INK,
})
HEATMAP_CMAP = LinearSegmentedColormap.from_list("publication_blue", PALETTE)

# --- Layout ---
CELL_HEIGHT = 0.26
COL_WIDTH = 0.55
DENDRO_WIDTH = 0.90
CBAR_SPACER = 0.35   # gap between dendrogram and colorbar
CBAR_WIDTH = 0.18
LEFT_MARGIN = 1.8
RIGHT_MARGIN = 0.9
TOP_MARGIN = 0.9
BOTTOM_MARGIN = 1.6


def clean_label(s):
    """Normalize whitespace (incl. non-breaking spaces) in a column label."""
    return re.sub(r"\s+", " ", str(s).replace("\xa0", " ")).strip()


def format_condition_label(c):
    """Render condition labels with a uniform 1 decimal for both time and %,
    replacing '(v/v)' with 'lys' (v/v is stated in the figure legend),
    e.g. '1 h (1% (v/v))' -> '1.0 h (1.0% lys)', so all x labels align."""
    m = re.match(r"([\d.]+)\s*h\s*\(([\d.]+)\s*%\s*\(v/v\)\)", c)
    if not m:
        return c
    return f"{float(m.group(1)):.1f} h ({float(m.group(2)):.1f}% lys)"


def phylogenetic_order(records, aligner):
    """Align sequences and return (leaf_order, linkage_Z, ids) via average-linkage on BLOSUM62 distances."""
    if shutil.which(aligner) is None:
        sys.exit(f"Error: '{aligner}' not found in PATH.")

    fasta_in = tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False)
    for rec in records:
        SeqIO.write(rec, fasta_in, "fasta")
    fasta_in.close()
    aligned = fasta_in.name + ".aln"

    subprocess.run(
        [aligner, "-i", fasta_in.name, "-o", aligned,
         "--force", "--auto", "--outfmt=fasta"],
        check=True,
    )
    alignment = MultipleSeqAlignment(list(SeqIO.parse(aligned, "fasta")))
    os.unlink(fasta_in.name)
    os.unlink(aligned)

    dm = DistanceCalculator("blosum62").get_distance(alignment)
    ids = list(dm.names)
    n = len(ids)
    condensed = [dm[ids[i], ids[j]] for i in range(n) for j in range(i + 1, n)]
    Z = linkage(condensed, method="average")
    order = dendrogram(Z, labels=ids, no_plot=True)["ivl"]
    return order, Z, ids


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("csv", help="Input CSV of enzymes x conditions")
    p.add_argument("--fasta", default=DEFAULT_FASTA,
                   help="Orthologs FASTA for sequence lookup (default: evolutionary_analysis/rLI2_orthologs.fasta)")
    p.add_argument("--sep", default=";", help="CSV field separator (default: ';')")
    p.add_argument("--aligner", default="clustalo", help="MSA executable (default: clustalo)")
    p.add_argument("--vmin", type=float, default=0.0)
    p.add_argument("--vmax", type=float, default=100.0)
    p.add_argument("--output", default=None,
                   help="Output figure base path (without extension). Default: alongside the CSV.")
    a = p.parse_args()

    df = pd.read_csv(a.csv, sep=a.sep)
    df.columns = [clean_label(c) for c in df.columns]
    df["Enzyme"] = df["Enzyme"].astype(str).str.strip()

    value_cols = [c for c in df.columns if c not in ("Enzyme", "sequence")]
    if not value_cols:
        sys.exit("Error: no value columns found (only Enzyme/sequence present).")
    for c in value_cols:
        df[c] = pd.to_numeric(
            df[c].astype(str).str.replace("\xa0", " ").str.strip(), errors="coerce"
        )
    df = df.dropna(subset=value_cols, how="all").reset_index(drop=True)

    # Look up each enzyme's sequence by publication ID; assign unique internal ids
    fasta_lookup = {rec.id: rec for rec in SeqIO.parse(a.fasta, "fasta")}
    records, internal_ids, display_labels = [], [], []
    seen = {}
    for _, row in df.iterrows():
        enzyme = row["Enzyme"]
        if enzyme not in fasta_lookup:
            sys.exit(f"Error: '{enzyme}' not found in {a.fasta}")
        seen[enzyme] = seen.get(enzyme, 0) + 1
        uid = enzyme if seen[enzyme] == 1 else f"{enzyme}__{seen[enzyme]}"
        label = DISPLAY_OVERRIDES.get(enzyme, enzyme)
        if seen[enzyme] > 1:
            label = f"{label} ({seen[enzyme]})"
        records.append(SeqRecord(Seq(str(fasta_lookup[enzyme].seq).replace("*", "")),
                                 id=uid, description=""))
        internal_ids.append(uid)
        display_labels.append(label)

    df["_uid"] = internal_ids
    df["_label"] = display_labels
    print(f"Loaded {len(df)} enzymes, {len(value_cols)} conditions")

    order, Z, ids = phylogenetic_order(records, a.aligner)
    df_ordered = df.set_index("_uid").loc[order].reset_index()

    # --- Figure ---
    # Column order: heatmap | dendrogram (flush against the heatmap) | spacer | colorbar.
    n_rows = len(df_ordered)
    n_cols = len(value_cols)
    width_total = (LEFT_MARGIN + n_cols * COL_WIDTH + DENDRO_WIDTH
                   + CBAR_SPACER + CBAR_WIDTH + RIGHT_MARGIN)
    height_total = TOP_MARGIN + BOTTOM_MARGIN + n_rows * CELL_HEIGHT

    fig = plt.figure(figsize=(width_total, height_total))
    gs = fig.add_gridspec(
        1, 4, width_ratios=[n_cols * COL_WIDTH, DENDRO_WIDTH, CBAR_SPACER, CBAR_WIDTH],
        wspace=0.0,
        left=LEFT_MARGIN / width_total, right=1 - RIGHT_MARGIN / width_total,
        top=1 - TOP_MARGIN / height_total, bottom=BOTTOM_MARGIN / height_total,
    )
    ax_hm = fig.add_subplot(gs[0, 0])
    ax_dend = fig.add_subplot(gs[0, 1])
    ax_cbar = fig.add_subplot(gs[0, 3])  # gs[0, 2] is an empty spacer

    sns.heatmap(
        df_ordered[value_cols].values,
        cmap=HEATMAP_CMAP, vmin=a.vmin, vmax=a.vmax,
        yticklabels=df_ordered["_label"],
        xticklabels=[format_condition_label(c) for c in value_cols],
        linewidths=1.2, linecolor="white",
        cbar=True, cbar_ax=ax_cbar, cbar_kws={"label": VALUE_LABEL},
        ax=ax_hm,
    )
    ax_hm.tick_params(axis="both", labelsize=10)
    ax_hm.tick_params(axis="x", rotation=90)
    plt.setp(ax_hm.get_yticklabels(), rotation=0)
    # Highlight the EpC reference row in bold.
    for lbl in ax_hm.get_yticklabels():
        if "EpC" in lbl.get_text():
            lbl.set_fontweight("bold")
    for spine in ax_hm.spines.values():
        spine.set_visible(True)
        spine.set_edgecolor(INK)
        spine.set_linewidth(1.0)

    dendrogram(Z, labels=ids, ax=ax_dend, orientation="right",
               color_threshold=0, above_threshold_color=INK)
    for line in ax_dend.get_lines():
        line.set_color(INK)
        line.set_linewidth(1.0)
    ax_dend.invert_yaxis()
    ax_dend.axis("off")

    base = a.output
    if base is None:
        d = os.path.dirname(os.path.abspath(a.csv))
        base = os.path.join(d, os.path.splitext(os.path.basename(a.csv))[0]
                            + "_heatmap_phylogenetic_order")
    plt.savefig(f"{base}.svg", bbox_inches="tight")
    plt.savefig(f"{base}.png", bbox_inches="tight", dpi=400)
    plt.close(fig)
    print(f"Figure saved: {base}.svg and {base}.png")


if __name__ == "__main__":
    main()

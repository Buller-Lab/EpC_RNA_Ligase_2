#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot top 10 coevolving residues mapped to a reference FASTA sequence."
    )
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-f", "--fasta", required=True)
    parser.add_argument("-p", "--position", type=int, required=True)
    return parser.parse_args()

def main():
    args = parse_args()

    if not os.path.exists(args.input) or not os.path.exists(args.fasta):
        print("Error: Input matrix or FASTA file does not exist.", file=sys.stderr)
        sys.exit(1)

    try:
        records = list(SeqIO.parse(args.fasta, "fasta"))
        if not records:
            print("Error: FASTA file is empty.", file=sys.stderr)
            sys.exit(1)
        ref_seq = str(records[0].seq).upper()
    except Exception as e:
        print(f"Error reading FASTA: {e}", file=sys.stderr)
        sys.exit(1)

    clean_seq = ref_seq.replace("-", "")
    if args.position < 1 or args.position > len(clean_seq):
        print(f"Error: Position {args.position} out of bounds for sequence of length {len(clean_seq)}.", file=sys.stderr)
        sys.exit(1)

    target_aa = clean_seq[args.position - 1]
    
    aln_pos = 0
    ref_pos = 0
    target_aln_idx = None
    
    for i, char in enumerate(ref_seq):
        if char != "-":
            ref_pos += 1
            if ref_pos == args.position:
                target_aln_idx = aln_pos
                break
        if char != "-": 
            aln_pos += 1
        else:
            pass

    print(f"Reading matrix data from {args.input}...")
    try:
        df = pd.read_excel(args.input, index_col=0)
    except Exception as e:
        print(f"Error reading Excel file: {e}", file=sys.stderr)
        sys.exit(1)

    labels = df.index.astype(str).tolist()
    
    pos_to_label = {}
    for label in labels:
        num_str = "".join([c for c in label if c.isdigit()])
        if num_str:
            pos_to_label[int(num_str)] = label

    matrix_pos = target_aln_idx + 1
    if matrix_pos not in pos_to_label:
        print(f"Error: Calculated matrix position {matrix_pos} ({target_aa}{args.position}) not found in Excel labels.", file=sys.stderr)
        sys.exit(1)

    target_label = pos_to_label[matrix_pos]
    print(f"Mapping Target: Reference position {args.position} -> Matrix column {target_label}")

    series = df[target_label].drop(index=target_label)
    top_10 = series.nlargest(10).reset_index()
    top_10.columns = ["Residue", "Z_Score"]

    n_bars = len(top_10)
    cell_size = 0.48
    fig_width = max(8.5, min(12.0, n_bars * cell_size + 3.5))
    full_height = max(7.0, n_bars * cell_size + 2.2)
    fig_height = full_height / 1.25

    fig, ax = plt.subplots(figsize=(fig_width, fig_height), facecolor="white")
    sns.set_theme(style="white", font_scale=1.0)

    bar_color = "#89b3d0" 
    x_positions = np.arange(n_bars)

    bars = ax.bar(
        x_positions,
        top_10["Z_Score"],
        color=bar_color,
        edgecolor="#404040",
        linewidth=0.8
    )

    max_val = top_10["Z_Score"].max()
    ax.set_ylim(0, max_val * 1.15)
    ax.set_ylabel("Coevolution Z-score", size=12, weight="bold", labelpad=12)
    ax.set_title(f"Top 10 Coevolving Residues for {target_aa}{args.position} (Col {target_label})", size=13, weight="bold", pad=15)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["bottom"].set_linewidth(0.8)

    ax.set_xticks(x_positions)
    ax.set_xticklabels([f"$\\bf{{{label}}}$" for label in top_10["Residue"]], rotation=45, ha="right")

    ax.tick_params(axis='x', pad=5, labelsize=11)
    ax.tick_params(axis='y', pad=3, labelsize=10.5)

    for label in ax.get_yticklabels():
        label.set_weight("bold")

    plt.tight_layout(pad=2.2)

    input_dir = os.path.dirname(os.path.abspath(args.input))
    base_name = f"ref_pos_{args.position}_{target_aa}_top10_barplot"
    output_path = os.path.join(input_dir, base_name)

    plt.savefig(f"{output_path}.png", format="png", dpi=450, bbox_inches="tight", facecolor="white")
    plt.savefig(f"{output_path}.svg", format="svg", bbox_inches="tight", facecolor="white")
    
    print(f"Completed Successfully.\nAssets generated at:\n - {output_path}.png\n - {output_path}.svg")
    plt.close()

if __name__ == "__main__":
    main()
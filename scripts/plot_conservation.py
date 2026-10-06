import argparse
import os
import sys
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

custom_cmap = LinearSegmentedColormap.from_list(
    "white_blue",
    ["#FFFFFF", "#89B3D0", "#4D8BD0", "#1C3A5E"]
)


def style_and_save_heatmap(ax, fig, cax, base_output_path, suffix):
    ax.tick_params(axis='x', pad=6, labelsize=11)
    ax.tick_params(axis='y', pad=4, labelsize=11)
    
    for tick in ax.get_xticklabels():
        tick.set_fontweight("bold")
        tick.set_fontfamily("Arial")
    for tick in ax.get_yticklabels():
        tick.set_fontweight("bold")
        tick.set_fontfamily("Arial")
        
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=10)
    for label in cbar.ax.get_yticklabels():
        label.set_fontweight("bold")
        label.set_fontfamily("Arial")
        
    cbar.set_label("Frequency (%)", size=11, weight="bold", labelpad=12)

    plt.tight_layout(pad=2.2)
    full_path = f"{base_output_path}_{suffix}" if suffix else base_output_path
    
    plt.savefig(f"{full_path}.png", format="png", dpi=450, bbox_inches="tight", facecolor="white")
    plt.savefig(f"{full_path}.svg", format="svg", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved: {full_path}.png and .svg")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("xlsx_file")
    parser.add_argument("--positions", required=True, help="Comma-separated list, e.g. 45,46,47,120,121")
    parser.add_argument("--output", default="conservation_heatmap", help="Output base name or file path without extension")
    args = parser.parse_args()

    selected = [int(x.strip()) for x in args.positions.split(",")]

    try:
        df = pd.read_excel(args.xlsx_file, sheet_name="Conservation", header=None)
    except Exception as e:
        print(f"Error reading Excel sheet: {e}", file=sys.stderr)
        sys.exit(1)

    positions = [int(x) for x in df.iloc[0, 1:].tolist()]
    ref_aas = df.iloc[1, 1:].tolist()

    data_start_row = 2
    aa_labels = df.iloc[data_start_row:, 0].tolist()
    matrix = df.iloc[data_start_row:, 1:].astype(float).values

    col_idx = [positions.index(p) for p in selected if p in positions]
    if len(col_idx) != len(selected):
        missing = [p for p in selected if p not in positions]
        print(f"Warning: positions not found in matrix: {missing}")

    sub_matrix = matrix[:, col_idx].T
    y_labels = [f"({ref}){p}" for p, ref in zip(selected, [ref_aas[i] for i in col_idx])]

    n_pos = len(selected)
    n_aa = len(aa_labels)

    cell_size = 0.45
    fig_height = max(3.5, n_pos * cell_size + 2.0)

    fig_width1 = max(8.5, min(14.5, n_aa * cell_size + 3.0))
    sns.set_theme(style="white")
    
    fig1, ax1 = plt.subplots(figsize=(fig_width1, fig_height), facecolor="white")
    divider1 = make_axes_locatable(ax1)
    cax1 = divider1.append_axes("right", size="3%", pad="1.5%")

    sns.heatmap(
        sub_matrix,
        xticklabels=aa_labels,
        yticklabels=y_labels,
        cmap=custom_cmap,
        annot=True,
        fmt=".1f",
        square=True,
        linewidths=1.5,
        linecolor="white",
        cbar_ax=cax1,
        vmin=0,
        vmax=100,
        annot_kws={"size": 12, "weight": "bold", "family": "Arial"},
        ax=ax1
    )
    
    ax1.set_xticklabels(aa_labels, rotation=0, ha="center")
    base_output = os.path.splitext(args.output)[0]
    style_and_save_heatmap(ax1, fig1, cax1, base_output, "")

    all_top_aas = set()
    top3_per_row = []
    
    for i in range(n_pos):
        row = sub_matrix[i, :]
        valid_idx = [j for j in range(len(row)) if aa_labels[j] != "-"]
        if valid_idx:
            top3_idx = sorted(valid_idx, key=lambda j: row[j], reverse=True)[:3]
        else:
            top3_idx = []
        top_aas = [aa_labels[idx] for idx in top3_idx]
        top3_per_row.append((top3_idx, top_aas))
        all_top_aas.update(top_aas)
    
    all_top_aas = sorted(list(all_top_aas))
    n_top = len(all_top_aas)
    
    consistent_matrix = np.zeros((n_pos, n_top + 1))
    col_map = {aa: idx for idx, aa in enumerate(all_top_aas)}
    
    for i in range(n_pos):
        row = sub_matrix[i, :]
        top3_idx, _ = top3_per_row[i]
        for idx in top3_idx:
            aa = aa_labels[idx]
            if aa in col_map:
                consistent_matrix[i, col_map[aa]] = row[idx]
        consistent_matrix[i, -1] = max(0.0, 100.0 - np.sum(consistent_matrix[i, :-1]))

    x_labels_top3 = all_top_aas + ["Other"]
    fig_width2 = max(6.5, min(12.5, (n_top + 1) * cell_size + 3.0))

    fig2, ax2 = plt.subplots(figsize=(fig_width2, fig_height), facecolor="white")
    divider2 = make_axes_locatable(ax2)
    cax2 = divider2.append_axes("right", size="3%", pad="1.5%")

    sns.heatmap(
        consistent_matrix,
        xticklabels=x_labels_top3,
        yticklabels=y_labels,
        cmap=custom_cmap,
        annot=True,
        fmt=".1f",
        square=True,
        linewidths=1.5,
        linecolor="white",
        cbar_ax=cax2,
        vmin=0,
        vmax=100,
        annot_kws={"size": 12, "weight": "bold", "family": "Arial"},
        ax=ax2
    )

    ax2.set_xticklabels(x_labels_top3, rotation=0, ha="center")
    style_and_save_heatmap(ax2, fig2, cax2, base_output, "top3")


if __name__ == "__main__":
    main()
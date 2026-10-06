import argparse
import sys
import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from Bio import AlignIO
import os
import string

mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

custom_cmap = LinearSegmentedColormap.from_list(
    "white_blue",
    ["#FFFFFF", "#89B3D0", "#4D8BD0", "#1C3A5E"]
)

def get_alignment_column_indices(aln, ref_id):
    ref_row = None
    for i, rec in enumerate(aln):
        if ref_id in rec.id:
            ref_row = i
            break
            
    if ref_row is None:
        print(f"Error: Reference ID '{ref_id}' not found in alignment headers.", file=sys.stderr)
        print("Available headers:", file=sys.stderr)
        for rec in aln[:10]:
            print(f"  {rec.id}", file=sys.stderr)
        sys.exit(1)
            
    ref_seq = str(aln[ref_row].seq).upper()
    aln_idx_to_ref_pos = {}
    ref_pos_to_aln_idx = {}
    
    current_ref_pos = 0
    for aln_idx, aa in enumerate(ref_seq):
        if aa != "-":
            current_ref_pos += 1
            aln_idx_to_ref_pos[aln_idx] = current_ref_pos
            ref_pos_to_aln_idx[current_ref_pos] = aln_idx
            
    return ref_pos_to_aln_idx, aln_idx_to_ref_pos

def main():
    parser = argparse.ArgumentParser(description="Plot co-occurring mutation frequency heatmaps from FASTA alignments.")
    parser.add_argument("fasta_aln", help="Path to the FASTA alignment file")
    parser.add_argument("--ref_id", required=True, help="Reference sequence ID in alignment")
    parser.add_argument("--focus", type=int, required=True, help="Focus reference position index")
    parser.add_argument("--focus_muts", required=True, help="Comma-separated list of focus mutations (e.g. E,K)")
    parser.add_argument("--targets", required=True, help="Comma-separated list of target positions to test against")
    parser.add_argument("--output", default="coev", help="Output directory or prefix")
    args = parser.parse_args()

    focus_muts = [x.strip().upper() for x in args.focus_muts.split(",")]
    target_positions = [int(x.strip()) for x in args.targets.split(",")]

    if len(target_positions) < 1:
        print("Error: Please provide at least one target position.", file=sys.stderr)
        sys.exit(1)

    muts_str = "".join(focus_muts)

    output_path = args.output
    if output_path.endswith(('/', '\\')) or os.path.isdir(output_path) or '.' not in os.path.basename(output_path):
        output_dir = output_path
        base_name = f"Foc{args.focus}_{muts_str}_vs_Targets_{'_'.join(map(str, target_positions))}"
    else:
        output_dir = os.path.dirname(output_path) or "."
        base_name = os.path.basename(output_path)

    os.makedirs(output_dir, exist_ok=True)
    full_output_name = os.path.join(output_dir, base_name)

    aln = AlignIO.read(args.fasta_aln, "fasta")
    ref_pos_to_aln_idx, _ = get_alignment_column_indices(aln, args.ref_id)
    
    if args.focus not in ref_pos_to_aln_idx:
        raise ValueError(f"Focus position {args.focus} not found or is a gap.")
    for t in target_positions:
        if t not in ref_pos_to_aln_idx:
            raise ValueError(f"Target position {t} not found or is a gap.")

    focus_aln_col = ref_pos_to_aln_idx[args.focus]
    aa_order = list("ACDEFGHIKLMNPQRSTVWY-")

    def build_full_matrix(target_pos):
        target_aln_col = ref_pos_to_aln_idx[target_pos]
        matrix_counts = {m: {aa_t: 0 for aa_t in aa_order} for m in focus_muts}
        total_per_focus = {m: 0 for m in focus_muts}
        
        for rec in aln:
            f_aa = str(rec.seq[focus_aln_col]).upper()
            t_aa = str(rec.seq[target_aln_col]).upper()
            if f_aa in focus_muts and t_aa in aa_order:
                matrix_counts[f_aa][t_aa] += 1
                total_per_focus[f_aa] += 1

        matrix_pct = []
        for m in focus_muts:
            row_vals = []
            denom = total_per_focus[m]
            for aa_t in aa_order:
                pct = (matrix_counts[m][aa_t] / denom * 100) if denom > 0 else 0.0
                row_vals.append(pct)
            matrix_pct.append(row_vals)
            
        return np.array(matrix_pct)

    matrices = [build_full_matrix(t) for t in target_positions]

    y_labels = focus_muts
    num_targets = len(target_positions)
    num_aa = len(aa_order)

    sns.set_theme(style="white")

    # Layout with increased spacing
    ax_height_per = 2.9
    subplot_vgap = 0.95          # Increased vertical gap between subplots
    top_margin = 0.8
    bottom_margin = 1.3
    left_margin = 1.2
    cbar_width = 0.23
    cbar_gap = 0.20
    right_margin = cbar_gap + cbar_width + 0.7

    ax_width = num_aa * (ax_height_per / len(focus_muts) * 0.96)
    fig_width = left_margin + ax_width + right_margin
    fig_height = top_margin + num_targets * ax_height_per + (num_targets - 1) * subplot_vgap + bottom_margin

    fig = plt.figure(figsize=(fig_width, fig_height), facecolor="white")

    axes = []
    cbar_axes = []

    for i in range(num_targets):
        y_pos = bottom_margin + (num_targets - 1 - i) * (ax_height_per + subplot_vgap)
        
        # Main heatmap axis
        ax = fig.add_axes([
            left_margin / fig_width,
            y_pos / fig_height,
            ax_width / fig_width,
            ax_height_per / fig_height
        ])
        axes.append(ax)
        
        # Individual colorbar axis
        cax = fig.add_axes([
            (left_margin + ax_width + cbar_gap) / fig_width,
            y_pos / fig_height,
            cbar_width / fig_width,
            ax_height_per / fig_height
        ])
        cbar_axes.append(cax)

    # Plot each target
    letters = string.ascii_lowercase

    for i, (target_pos, matrix) in enumerate(zip(target_positions, matrices)):
        ax = axes[i]
        cax = cbar_axes[i]

        sns.heatmap(
            matrix,
            xticklabels=aa_order,
            yticklabels=y_labels,
            cmap=custom_cmap,
            annot=True,
            fmt=".1f",
            square=False,
            linewidths=1.5,
            linecolor="white",
            vmin=0,
            vmax=100,
            annot_kws={"size": 15.5, "weight": "bold", "family": "Arial"},
            cbar=True,
            cbar_ax=cax,
            ax=ax
        )

        # Axis labels
        ax.set_xlabel(f"Position {target_pos}", fontsize=12, labelpad=9, weight="bold")
        ax.set_ylabel(f"Position {args.focus}", fontsize=12, labelpad=9, weight="bold")

        ax.tick_params(axis='x', pad=6, labelsize=11, rotation=0)
        ax.tick_params(axis='y', pad=6, labelsize=11)

        # Bold ticks
        for tick in ax.get_xticklabels():
            tick.set_fontweight("bold")
            tick.set_fontfamily("Arial")
        for tick in ax.get_yticklabels():
            tick.set_fontweight("bold")
            tick.set_fontfamily("Arial")

        # Panel label (a, b, c, ...) in top-left corner
        ax.text(-0.12, 1.08, f"({letters[i]})", transform=ax.transAxes,
                fontsize=14, fontweight="bold", va='top', ha='right')

        # Colorbar
        cbar = ax.collections[0].colorbar
        cbar.set_label("Frequency (%)", size=11, weight="bold", labelpad=10)
        cbar.ax.tick_params(labelsize=10)
        for label in cbar.ax.get_yticklabels():
            label.set_fontweight("bold")
            label.set_fontfamily("Arial")

    plt.savefig(f"{full_output_name}.png", dpi=450, bbox_inches="tight", facecolor="white")
    plt.savefig(f"{full_output_name}.svg", bbox_inches="tight", facecolor="white")
    plt.close()

    print(f"Saved plot with {num_targets} subplots (each with own colorbar + panel label) as: {full_output_name}.png")

if __name__ == "__main__":
    main()
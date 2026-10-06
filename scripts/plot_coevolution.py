import argparse
import sys
import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from Bio import AlignIO
import os

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

    if len(target_positions) < 2:
        print("Error: Please provide at least two target positions.", file=sys.stderr)
        sys.exit(1)

    t1, t2 = target_positions[0], target_positions[1]
    muts_str = "".join(focus_muts)

    output_path = args.output
    if output_path.endswith(('/', '\\')) or os.path.isdir(output_path) or '.' not in os.path.basename(output_path):
        output_dir = output_path
        base_name = f"Foc{args.focus}_{muts_str}_vs_Targets_{t1}_{t2}"
    else:
        output_dir = os.path.dirname(output_path) or "."
        base_name = os.path.basename(output_path)

    os.makedirs(output_dir, exist_ok=True)
    full_output_name = os.path.join(output_dir, base_name)

    aln = AlignIO.read(args.fasta_aln, "fasta")
    ref_pos_to_aln_idx, _ = get_alignment_column_indices(aln, args.ref_id)
    
    if args.focus not in ref_pos_to_aln_idx or t1 not in ref_pos_to_aln_idx or t2 not in ref_pos_to_aln_idx:
        raise ValueError("Focus or target positions map to a gap or are out of index scope.")

    focus_aln_col = ref_pos_to_aln_idx[args.focus]
    aa_order = list("ACDEFGHIKLMNPQRSTVWY-")

    def build_matrix_with_other(target_pos):
        target_aln_col = ref_pos_to_aln_idx[target_pos]
        matrix_counts = {m: {aa_t: 0 for aa_t in aa_order} for m in focus_muts}
        total_per_focus = {m: 0 for m in focus_muts}
        
        for rec in aln:
            f_aa = str(rec.seq[focus_aln_col]).upper()
            t_aa = str(rec.seq[target_aln_col]).upper()
            if f_aa in focus_muts and t_aa in aa_order:
                matrix_counts[f_aa][t_aa] += 1
                total_per_focus[f_aa] += 1

        total_counts_per_aa = {aa_t: sum(matrix_counts[m][aa_t] for m in focus_muts) for aa_t in aa_order}
        sorted_aas = [aa for aa, count in sorted(total_counts_per_aa.items(), key=lambda item: item[1], reverse=True) if count > 0]
        
        top_2_aas = sorted_aas[:2]
        remaining_aas = sorted_aas[2:]
        
        matrix_pct = []
        for m in focus_muts:
            row_vals = []
            denom = total_per_focus[m]
            
            for aa_t in top_2_aas:
                pct = (matrix_counts[m][aa_t] / denom * 100) if denom > 0 else 0.0
                row_vals.append(pct)
                
            if remaining_aas:
                other_count = sum(matrix_counts[m][aa_t] for aa_t in remaining_aas)
                other_pct = (other_count / denom * 100) if denom > 0 else 0.0
                row_vals.append(other_pct)
            else:
                row_vals.append(0.0)
                
            matrix_pct.append(row_vals)
            
        display_cols = top_2_aas + ["Other"]
        return np.array(matrix_pct), display_cols

    sub_matrix1, display_cols1 = build_matrix_with_other(t1)
    sub_matrix2, display_cols2 = build_matrix_with_other(t2)
    
    y_labels = focus_muts

    sns.set_theme(style="white")

    ax_height = 2.5
    top_margin = 0.5
    bottom_margin = 1.0
    left_margin = 1.0
    subplot_gap = 0.4
    cbar_gap = 0.2
    cbar_width = 0.15
    right_margin = cbar_gap + cbar_width + 0.8
    
    fig_height = ax_height + top_margin + bottom_margin
    num_rows = len(focus_muts)
    cell_size = ax_height / num_rows
    
    ax_width1 = len(display_cols1) * cell_size
    ax_width2 = len(display_cols2) * cell_size
    fig_width = left_margin + ax_width1 + subplot_gap + ax_width2 + right_margin

    fig = plt.figure(figsize=(fig_width, fig_height), facecolor="white")

    x_pos1 = left_margin
    x_pos2 = left_margin + ax_width1 + subplot_gap
    x_pos_cbar = x_pos2 + ax_width2 + cbar_gap

    ax1 = fig.add_axes([x_pos1 / fig_width, bottom_margin / fig_height, ax_width1 / fig_width, ax_height / fig_height])
    ax2 = fig.add_axes([x_pos2 / fig_width, bottom_margin / fig_height, ax_width2 / fig_width, ax_height / fig_height])
    cbar_ax = fig.add_axes([x_pos_cbar / fig_width, bottom_margin / fig_height, cbar_width / fig_width, ax_height / fig_height])

    sns.heatmap(
        sub_matrix1,
        xticklabels=display_cols1,
        yticklabels=y_labels,
        cmap=custom_cmap,
        annot=True,
        fmt=".1f",
        square=False,
        linewidths=1.5,
        linecolor="white",
        vmin=0,
        vmax=100,
        annot_kws={"size": 17, "weight": "bold", "family": "Arial"},
        cbar=False,
        ax=ax1
    )

    sns.heatmap(
        sub_matrix2,
        xticklabels=display_cols2,
        yticklabels=False,
        cmap=custom_cmap,
        annot=True,
        fmt=".1f",
        square=False,
        linewidths=1.5,
        linecolor="white",
        vmin=0,
        vmax=100,
        annot_kws={"size": 17, "weight": "bold", "family": "Arial"},
        cbar_ax=cbar_ax,
        ax=ax2
    )

    ax1.set_xlabel(f"Position {t1}", fontsize=11, labelpad=8, weight="bold")
    ax1.set_ylabel(f"Position {args.focus}", fontsize=11, labelpad=8, weight="bold")
    ax1.tick_params(axis='x', pad=6, labelsize=11)
    ax1.tick_params(axis='y', pad=6, labelsize=11)
    ax1.set_xticklabels(display_cols1, rotation=0)
    ax1.set_yticklabels(y_labels, rotation=0)

    ax2.set_xlabel(f"Position {t2}", fontsize=11, labelpad=8, weight="bold")
    ax2.set_ylabel("")
    ax2.tick_params(axis='x', pad=6, labelsize=11)
    ax2.tick_params(axis='y', left=False)
    ax2.set_xticklabels(display_cols2, rotation=0)

    for ax in [ax1, ax2]:
        for tick in ax.get_xticklabels():
            tick.set_fontweight("bold")
            tick.set_fontfamily("Arial")
        for tick in ax.get_yticklabels():
            tick.set_fontweight("bold")
            tick.set_fontfamily("Arial")

    cbar = ax2.collections[0].colorbar
    cbar.set_label("Frequency (%)", size=11, weight="bold", labelpad=12)
    cbar.ax.tick_params(labelsize=10)

    for label in cbar.ax.get_yticklabels():
        label.set_fontweight("bold")
        label.set_fontfamily("Arial")

    plt.savefig(f"{full_output_name}.png", dpi=450, bbox_inches="tight", facecolor="white")
    plt.savefig(f"{full_output_name}.svg", bbox_inches="tight", facecolor="white")
    plt.close()

    print(f"Saved combined plot as: {full_output_name}.png")

if __name__ == "__main__":
    main()
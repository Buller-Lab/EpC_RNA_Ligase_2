import argparse
import os
import re
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", required=True)
args = parser.parse_args()

aa_3_to_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
}

VALID_ANNOTATIONS = ["conserved", "tested", "consensus"]

def extract_amino_acid(contact_str):
    if pd.isna(contact_str):
        return None
    cleaned = contact_str.replace("[", "").replace("]", "")
    parts = [p.strip() for p in cleaned.split(",")]
    for part in parts:
        match = re.match(r"^([A-Z]{3})(\d+)$", part, re.IGNORECASE)
        if match:
            aa_3 = match.group(1).upper()
            num = match.group(2)
            if aa_3 in aa_3_to_1:
                return f"{aa_3}{num}"
    return None

def to_one_letter(aa_str):
    match = re.match(r"^([A-Z]{3})(\d+)$", aa_str)
    if match:
        return aa_3_to_1[match.group(1)] + match.group(2)
    return aa_str

def process_sheet(df):
    df["AminoAcid"] = df["Contact"].apply(extract_amino_acid)
    df = df.dropna(subset=["AminoAcid"])
    df = df[df["Annotation"].isin(VALID_ANNOTATIONS)]
    sheet_grouped = df.groupby("AminoAcid").agg({"Average": "mean", "Annotation": "first"}).reset_index()
    return sheet_grouped

file_path = args.input
xls = pd.ExcelFile(file_path)

df_compact = process_sheet(pd.read_excel(xls, sheet_name="compact"))
df_loose = process_sheet(pd.read_excel(xls, sheet_name="loose"))

df_compact["Sheet_Count"] = 1
df_loose["Sheet_Count"] = 1

df_combined = pd.concat([df_compact, df_loose])

counts = df_combined.groupby("AminoAcid")["Sheet_Count"].sum()
shared_amino_acids = counts[counts == 2].index

df_combined = df_combined[df_combined["AminoAcid"].isin(shared_amino_acids)]

final_df = df_combined.groupby("AminoAcid").agg({"Average": ["mean", "std"], "Annotation": "first"})
final_df.columns = ["Final_Average", "Final_Stddev", "Annotation"]
final_df = final_df.reset_index()
final_df["Final_Stddev"] = final_df["Final_Stddev"].fillna(0)
final_df["Residue"] = final_df["AminoAcid"].apply(to_one_letter)
final_df = final_df.sort_values(by="Final_Average", ascending=False).reset_index(drop=True)

color_map = {
    "conserved": "#a68cb7",
    "tested": "#dfb295",
    "consensus": "#89b3d0"
}

legend_labels = {
    "conserved": "highly conserved (>95% occurence)",
    "tested": "already tested (previous rounds)",
    "consensus": "consensus (in compact & loose)"
}

n_bars = len(final_df)
cell_size = 0.48
fig_width = max(10.5, min(14.5, n_bars * cell_size + 3.5))
full_height = max(8.0, n_bars * cell_size + 2.2)
fig_height = full_height / 4.0

if n_bars <= 8:
    y_tick_size = 11
    label_size = 11
else:
    y_tick_size = 10.5
    label_size = 12

fig, ax = plt.subplots(figsize=(fig_width, fig_height), facecolor="white")
sns.set_theme(style="white", font_scale=1.0)

bar_colors = [color_map.get(annot, "#7f7f7f") for annot in final_df["Annotation"]]
x_positions = np.arange(len(final_df))

bars = ax.bar(
    x_positions,
    final_df["Final_Average"],
    yerr=final_df["Final_Stddev"],
    capsize=3,
    color=bar_colors,
    edgecolor="#404040",
    linewidth=0.8,
    error_kw={"ecolor": "#404040", "linewidth": 0.8}
)

consensus_indices = final_df[final_df["Annotation"] == "consensus"].index[:4]
for idx in consensus_indices:
    val = final_df.loc[idx, "Final_Average"]
    err = final_df.loc[idx, "Final_Stddev"]
    ax.plot(idx, val + err + 0.10, marker='*', color='#404040', markersize=7, zorder=5)

ax.set_ylim(0, 1.1)
ax.set_ylabel("Contact Frequency", size=label_size, weight="bold", labelpad=12)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_linewidth(0.8)
ax.spines["bottom"].set_linewidth(0.8)

ax.set_xticks(x_positions)
ax.set_xticklabels([f"$\\bf{{{label}}}$" for label in final_df["Residue"]], rotation=45, ha="right")

ax.tick_params(axis='x', pad=5, labelsize=11)
ax.tick_params(axis='y', pad=3, labelsize=y_tick_size)

for label in ax.get_yticklabels():
    label.set_weight("bold")

legend_elements = [
    Patch(facecolor=color_map[annot], edgecolor="#404040", linewidth=0.8, label=legend_labels[annot]) 
    for annot in VALID_ANNOTATIONS if annot in final_df["Annotation"].unique()
]
ax.legend(handles=legend_elements, loc="upper right", frameon=False, fontsize=10)

plt.tight_layout(pad=2.2)

input_dir = os.path.dirname(os.path.abspath(file_path))
base_name = os.path.splitext(os.path.basename(file_path))[0] + "_barplot"
output_path = os.path.join(input_dir, base_name)

plt.savefig(f"{output_path}.png", format="png", dpi=450, bbox_inches="tight", facecolor="white")
plt.savefig(f"{output_path}.svg", format="svg", bbox_inches="tight", facecolor="white")

plt.show()
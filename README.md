## This is the GitHub Repository for the publication:

# Evolved RNA ligase 2 facilitates therapeutic RNA oligonucleotide manufacturing

In this repository, we provide the code and data used for the computationally guided enzyme discovery and engineering of an RNA ligase 2 from *Erwinia* phage Cronus. The computational workflows provided here supported five rounds of evolution to engineer an RNA ligase capable of assembling a 100-nt single guide RNA from two SPOS-derived 50-mers, as well as handling fully modified small interfering RNAs carrying 2’ sugar modifications.

Specifically, this repository contains scripts for:
1. Co-evolution analysis utilizing an adapted GREMLIN_v2 approach [1].
2. Consensus of interacting residues derived from Molecular Dynamics (MD) simulations.
3. Conservation analysis of homologous sequences utilizing MMseqs2 [2].
4. Post-analysis of MD simulations, including binding energy calculations utilizing Rosetta [3].
5. Ortholog sequence analysis (pairwise identity, length distribution, and phylogenetic clustering).

[1] Kamisetty, Hetunandan, et al. "Assessing the utility of coevolution-based residue–residue contact predictions in a sequence-and structure-rich era." Proceedings of the National Academy of Sciences 110.40 (2013): 15674-15679. *(Update if a different GREMLIN reference is preferred)*

[2] Steinegger, Martin, and Johannes Söding. "MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets." Nature biotechnology 35.11 (2017): 1026-1028.

[3] Chaudhury, Sidhartha, et al. "PyRosetta: a script-based interface for implementing molecular modeling algorithms using Rosetta." Bioinformatics 26.5 (2010): 689-691.

[4] Özgen, Fatma Feyza, et al. "Biomass‐Derived Diformylxylose as a Renewable Solvent for Biocatalysis Applications." ChemSusChem 19.4 (2026): e202502273.  https://doi.org/10.1002/cssc.202502273

[5] Honda Malca, Sumire, Peter Stockinger, Miquel Estévez-Gay, and Rebecca Buller. "Discovery and engineering of polymerases and ligases for the synthesis of modified nucleic acids." Current Opinion in Chemical Biology 93 (2026): 102697. https://doi.org/10.1016/j.cbpa.2026.102697

# Installation

We recommend running this code on UNIX-based systems such as Ubuntu. This repository can be downloaded to your local machine via the command:
```bash
git clone https://github.com/Buller-Lab/EpC_RNA_Ligase_2
```
then navigate into the cloned repository with:
```bash
cd EpC_RNA_Ligase_2
```

this should only take a few seconds.

# System Requirements

## Hardware requirements

This code was developed and tested on the following hardware:

- CPU: AMD Ryzen Threadripper 3970X 32-Core Processor
- Memory: 128 GiB RAM
- GPU: 2x NVIDIA GeForce RTX 3090

## Software requirements
All dependencies are specified in `environment.yml` (standard data-science stack
plus BioPython, and the external tools MMseqs2, MAFFT and Clustal Omega, as well
as MDTraj/PyEMMA and PyRosetta for the MD post-analysis). Create the environment
with conda or mamba/micromamba:

```bash
# with conda (or mamba)
conda env create -f environment.yml
conda activate rnalig

# or with micromamba
micromamba create -y -n rnalig -f environment.yml
micromamba activate rnalig
```

Notes:
- PyRosetta is pulled from the RosettaCommons channel declared in `environment.yml`
  and requires accepting their license. If you do not need the MD binding-energy
  step, you may remove the `pyrosetta`, `mdtraj` and `pyemma` lines.
- For Workflow 3 (Conservation Analysis) you also need an MMseqs2 sequence
  database; follow the official MMseqs2 documentation to download it.

## Ortholog numbering

The RNA ligase 1 and 2 orthologs follow the publication IDs of Table S3
(`T4_Rnl1`/`T4_Rnl2`, `Rnl1_02`…`Rnl1_22`, `Rnl2_02`…`Rnl2_30`). The reference
enzyme from *Erwinia* phage Cronus (EpC) is `Rnl2_27` (shown as "EpC (Rnl2_27)"
in figures). Two sequences listed in Table S3 are marked *excluded* (Escherichia
phage T2 and *Methanothermobacter thermautotrophicus*) and are not part of the
panel. The residue-level coevolution and conservation analyses use UniProt
accessions and the `EpC_rLI2` reference, and are therefore unaffected by this
numbering.

# Instructions for use

## 1. Co-evolution analysis
This analysis is based on GREMLIN_v2 (slightly adapted code to provide z-scores). The input Multiple Sequence Alignment (MSA) was generated with the MAFFT online server (with default settings) and the alignment fasta file was manually adapted to remove linebreaks.

**Run the co-evolution analysis:**
```bash
python scripts/coevolution_analysis.py evolutionary_analysis/msa_for_coevolution.fasta \
    --ref evolutionary_analysis/EpC_rLI2.fasta \
    --outdir evolutionary_analysis/R4_coevolving_positions
```
This writes `coevolution.xlsx` to the output directory.

**Visualize co-evolution z-scores for a selected set of positions (e.g., with position 285):**
```bash
python scripts/plot_zscores.py -i evolutionary_analysis/R4_coevolving_positions/coevolution.xlsx \
    -f evolutionary_analysis/EpC_rLI2.fasta \
    -p 285 -s 322 319 262 280 237 277 270
```

**Visualize the top 10 co-evolving residues for a given position:**
```bash
python scripts/plot_top_zscores.py -i evolutionary_analysis/R4_coevolving_positions/coevolution.xlsx \
    -f evolutionary_analysis/EpC_rLI2.fasta \
    -p 285
```

**Exemplary command for co-evolving amino acids (focusing on variants that have E, D, and K in position 285):**
```bash
python scripts/plot_coevolution.py --ref_id EpC_rLI2 \
    --focus 285 \
    --focus_muts E,D,K \
    --targets 237,277 \
    --output evolutionary_analysis/R4_coevolving_positions evolutionary_analysis/msa_for_coevolution.fasta
```

**Plot co-occurring mutation frequencies across all amino acids for the target positions:**
```bash
python scripts/plot_coevolution_all.py --ref_id EpC_rLI2 \
    --focus 285 \
    --focus_muts E,D,K \
    --targets 322,319,262,280,237,277,270 \
    --output evolutionary_analysis/R4_coevolving_positions evolutionary_analysis/msa_for_coevolution.fasta
```

## 2. Consensus of Interacting Residues (derived from MD)
This section handles the analysis of conserved contacts (these residues will not be mutated).

**Identify conserved contacts:**
```bash
python scripts/conservation_analysis.py evolutionary_analysis/rLI2_homologs/hits.fasta \
    --ref evolutionary_analysis/EpC_rLI2.fasta \
    --outdir evolutionary_analysis/R5_conserved_contacts
```
This writes `04_conservation.xlsx` to the output directory.

**Plot the conservation heatmap for specific positions:**
```bash
python scripts/plot_conservation.py evolutionary_analysis/R5_conserved_contacts/04_conservation.xlsx \
    --output evolutionary_analysis/R5_conserved_contacts/heatmap \
    --positions 55,56,64,66,205,226,228,232,234
```
The `--output` value is a base path (without extension); the script saves both `.png` and `.svg` files.

**Get consensus contacts (compact vs. loose):**
```bash
python scripts/get_consensus_contacts.py -i molecular_analysis/MD_contact_frequencies/compact_vs_loose_contact_frequencies.xlsx
```

## 3. Conservation Analysis
*Note: The result of this analysis depends on the pool of homologous sequences. We realized considerable differences upon recent reproduction, mostly caused by new/removed UniProt entries.*

To reproduce the workflow with your own database (results may differ based on sequence entries):

**Step 1:** Create a conda environment for MMseqs2 and download/specify the database according to the MMseqs2 documentation.

**Step 2:** Extract homologs (acknowledgement: imidase & polymerase review [5]). The `get_homologs.py` script is part of the Buller-Lab enzyme discovery toolkit, available at https://buller-lab.github.io/DNA-RNA_Polymerases_Ligases.
```bash
python scripts/get_homologs.py -i evolutionary_analysis/rLI2_orthologs.fasta \
    -d /mnt/bkup/tools/UniProt \
    -o evolutionary_analysis/rLI2_homologs
```

**Step 3:** Run the conservation analysis to design the R5 combinatorial library (select mutations for MD-derived positions).
```bash
python scripts/conservation_analysis.py evolutionary_analysis/rLI2_homologs/hits.fasta \
    --ref evolutionary_analysis/EpC_rLI2.fasta \
    --outdir evolutionary_analysis/R5_combin_library
```

**Step 4:** Plot the conservation heatmap for the selected positions.
```bash
python scripts/plot_conservation.py evolutionary_analysis/R5_combin_library/04_conservation.xlsx \
    --output evolutionary_analysis/R5_combin_library/heatmap \
    --positions 67,230,233,238
```

## 4. MD Simulations (Post-analysis)
MD simulations were performed as described in the methods section, based on previously published code (see Reference [4] - DFX paper). The post-analysis calculates binding energies and maps the energy landscape. The docked input structures for the compact and loose systems are provided under `molecular_analysis/MD_input_structures/`; `scripts/change_numbering.py` is a helper to apply a constant offset to all residue numbers in a PDB file (writes `<input>_renumbered.pdb` next to the input):
```bash
python scripts/change_numbering.py molecular_analysis/MD_input_structures/28_docked.pdb 100
```

**Identify the top 100 frames:**
For all variants in the two systems (compact vs loose), the top 100 frames were identified based on the distance of the ligation site to the catalytic site.
```bash
python scripts/identify_top100_frames.py
```

**Calculate binding energies:**
The resulting conformational ensembles were used as input for binding energy calculation utilizing PyRosetta.
```bash
python scripts/calculate_binding_energy.py
```

**Plot the per-residue RMSF:**
The RMSF input (`molecular_analysis/MD_rmsf_analysis/EpC_RMSF_analysis.csv`) holds one `replicate N` row per replicate and one column per residue number. The plot shows the mean over replicates with a ±SD band.
```bash
python scripts/plot_rmsf.py -i molecular_analysis/MD_rmsf_analysis/EpC_RMSF_analysis.csv \
    -o molecular_analysis/MD_rmsf_analysis/rmsf_plot
```
The `-o` value is a base path (without extension); the script saves both `.png` and `.svg` files.

## 5. Ortholog Sequence Analysis
For both RNA ligase 1 (rLI1) and RNA ligase 2 (rLI2) orthologs, we computed pairwise sequence identities (global Needleman-Wunsch alignment), length distributions, and a UPGMA phylogenetic dendrogram. The precomputed results and the orthologous sequences are provided under `orthologs/rLI1/` and `orthologs/rLI2/`.

**Run the analysis:**
```bash
python scripts/analyze_ortholog_sequences.py orthologs/rLI2/rLI2_orthologs.fasta \
    --outdir orthologs/rLI2
```
This generates `identity_matrix.xlsx`, `sequence_identity_heatmap`, `sequence_length_distribution`, and `phylogenetic_dendrogram` (each as `.png` and `.svg`), and a `description.txt` summary report. The same command applies to rLI1:
```bash
python scripts/analyze_ortholog_sequences.py orthologs/rLI1/rLI1_orthologs.fasta \
    --outdir orthologs/rLI1
```

**Amino acid conservation across the rLI2 orthologs:**
```bash
python scripts/conservation_analysis.py evolutionary_analysis/rLI2_orthologs.fasta \
    --ref evolutionary_analysis/EpC_rLI2.fasta \
    --outdir evolutionary_analysis/rLI2_orthologs
```
This writes `04_conservation.xlsx` to the output directory. Plot the conservation heatmap for the selected positions:
```bash
python scripts/plot_conservation.py evolutionary_analysis/rLI2_orthologs/04_conservation.xlsx \
    --output evolutionary_analysis/rLI2_orthologs/heatmap \
    --positions 57,59,62,65,109,215,218,221,236,296
```

**Activity heatmap of top orthologs across assay conditions (phylogenetically ordered):**
This plots the measured activity of the top Rnl2 orthologs across the different assay conditions (time points and enzyme loadings), with rows ordered by a sequence-based clustering of the orthologs. Enzymes are referenced by their publication IDs (`Rnl2_<N>`, `T4_Rnl2`, and `Rnl2_27` for the EpC reference) and sequences are looked up directly from the orthologs FASTA. Alignment uses Clustal Omega (`clustalo` must be in `PATH`). The sequences are read from `evolutionary_analysis/rLI2_orthologs.fasta` by default (override with `--fasta`); the figure is written next to the input CSV (override with `--output`):
```bash
python scripts/plot_heatmaps.py evolutionary_analysis/Top15_rLI2_orthologs/top15_rLI2_orthologs_diff_conditions.csv
```

# References

If you utilize this code, please cite:
Sumire Honda Malca, Peter Stockinger, Daniela Milbredt, Nadine Duss, Michael Niklaus, David Patsch, Lisa Schelbert, Nicolas Imstepf, Kaja Stalder, Nicolò Zuin Fantoni, Irene Marzuoli, Steven Paul Hanlon, Stefan G. Koenig, Filippo Sladojevich, Hans Iding, Serena Bisagni, and Rebecca Buller


"Evolved RNA ligase 2 facilitates therapeutic RNA oligonucleotide manufacturing." (DOI to be added upon publication)

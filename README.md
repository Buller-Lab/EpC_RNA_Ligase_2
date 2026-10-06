## This is the GitHub Repository for the publication: 

# Evolved RNA ligase 2 facilitates therapeutic RNA oligonucleotide manufacturing
[![DOI](https://zenodo.org/badge/DOI/[INSERT_ZENODO_DOI].svg)](https://doi.org/NSERT_ZENODO_DOI])

In this repository, we provide the code and data used for the computationally guided zyme discovery and engineering of an RNA ligase 2 from *Erwinia* phage Cronus. The mputational workflows provided here supported five rounds of evolution to engineer an RNA gase capable of assembling a 100-nt single guide RNA from two SPOS-derived 50-mers, as ll as handling fully modified small interfering RNAs carrying 2’ sugar modifications. 

Specifically, this repository contains scripts for:
1. Co-evolution analysis utilizing an adapted GREMLIN_v2 approach [1].
2. Consensus of interacting residues derived from Molecular Dynamics (MD) simulations.
3. Conservation analysis of homologous sequences utilizing MMseqs2 [2].
4. Post-analysis of MD simulations, including binding energy calculations utilizing Rosetta [3].

[1] Kamisetty, Hetunandan, et al. "Assessing the utility of coevolution-based sidue–residue contact predictions in a sequence-and structure-rich era." Proceedings of e National Academy of Sciences 110.40 (2013): 15674-15679. *(Update if a different EMLIN reference is preferred)*

[2] Steinegger, Martin, and Johannes Söding. "MMseqs2 enables sensitive protein sequence arching for the analysis of massive data sets." Nature biotechnology 35.11 (2017): 26-1028.

[3] Chaudhury, Sidhartha, et al. "PyRosetta: a script-based interface for implementing lecular modeling algorithms using Rosetta." Bioinformatics 26.5 (2010): 689-691.

[4] [INSERT DFX PAPER REFERENCE CITED IN MD SIMULATIONS SECTION]

# Installation

We recommend running this code on UNIX-based systems such as Ubuntu. This repository can  downloaded to your local machine via the command:
```bash
git clone [https://github.com/Buller-Lab/](https://github.com/Buller-Lab/)NSERT_REPO_NAME]
```
then navigate into the cloned repository with:
```bash
cd [INSERT_REPO_NAME]
```

this should only take a few seconds.

# System Requirements

## Hardware requirements

This code was developed and tested on the following hardware:

- CPU: AMD Ryzen Threadripper 3970X 32-Core Processor
- Memory: 128 GiB RAM
- GPU: 2x NVIDIA GeForce RTX 3090

## Software requirements
To create conda environments with the necessary dependencies, run:

```bash
# [INSERT CONDA ENVIRONMENT CREATION COMMANDS HERE, e.g., environment.yml containing Rosetta and standard data science packages]
# Note: For Workflow 3 (Conservation Analysis), a separate environment with MMseqs2 is quired. Please follow the official MMseqs2 documentation to install and download the cessary databases.
```

# Instructions for use

## 1. Co-evolution analysis
This analysis is based on GREMLIN_v2 (slightly adapted code to provide z-scores). The put Multiple Sequence Alignment (MSA) was generated with the MAFFT online server (with fault settings) and the alignment fasta file was manually adapted to remove linebreaks.

**Run the co-evolution analysis:**
```bash
python scripts/coevolution_analysis.py evolutionary_analysis/msa_for_coevolution.fasta \
    --ref evolutionary_analysis/EpC_rLI2.fasta \
    --outdir evolutionary_analysis/R4_coevolving_positions
```

**Visualize top co-evolving residues (e.g., with position 285) based on z-scores:**
```bash
python scripts/plot_zscores.py -i evolutionary_analysis/R4_coevolving_positions/evolution.xlsx \
    -f evolutionary_analysis/EpC_rLI2.fasta \
    -p 285 -s 322 319 262 280 237 277 270
```

**Exemplary command for co-evolving amino acids (focusing on variants that have E, D, and in position 285):**
```bash
python scripts/plot_coevolution.py --ref_id EpC_rLI2 \
    --focus 285 \
    --focus_muts E,D,K \
    --targets 237,277 \
    --output evolutionary_analysis/R4_coevolving_positions evolutionary_analysis/a_for_coevolution.fasta
```

## 2. Consensus of Interacting Residues (derived from MD)
This section handles the analysis of conserved contacts (these residues will not be tated).

**Identify conserved contacts:**
```bash
python scripts/conservation_analysis.py evolutionary_analysis/rLI2_homologs/hits.fasta \
    --ref evolutionary_analysis/EpC_rLI2.fasta \
    --outdir evolutionary_analysis/R5_conserved_contacts
```

**Plot the conservation heatmap for specific positions:**
```bash
python scripts/plot_conservation.py evolutionary_analysis/R5_conserved_contacts/_conservation.xlsx \
    --output evolutionary_analysis/R5_conserved_contacts/heatmap.png \
    --positions 55,56,64,66,205,226,228,232,234
```

**Get consensus contacts (compact vs. loose):**
```bash
python scripts/get_consensus_contacts.py -i molecular_analysis/MD_contact_frequencies/mpact_vs_loose_contact_frequencies.xlsx 
```

## 3. Conservation Analysis
*Note: The result of this analysis depends on the pool of homologous sequences. We alized considerable differences upon recent reproduction, mostly caused by new/removed iProt entries.* 

To reproduce the workflow with your own database (results may differ based on sequence tries):

**Step 1:** Create a conda environment for MMseqs2 and download/specify the database cording to the MMseqs2 documentation.

**Step 2:** Extract homologs (acknowledgement: imidase & polymerase review).
```bash
python scripts/get_homologs.py -i evolutionary_analysis/rLI2_orthologs.fasta \
    -d /mnt/bkup/tools/UniProt \
    -o evolutionary_analysis/rLI2_homologs
```

**Step 3:** Run the conservation analysis to design the R5 combinatorial library (select tations for MD-derived positions).
```bash
python scripts/conservation_analysis.py evolutionary_analysis/rLI2_homologs/hits.fasta \
    --ref evolutionary_analysis/EpC_rLI2.fasta \
    --outdir evolutionary_analysis/R5_combin_library
```

**Step 4:** Plot the conservation heatmap for the selected positions.
```bash
python scripts/plot_conservation.py evolutionary_analysis/R5_combin_library/_conservation.xlsx \
    --output evolutionary_analysis/R5_combin_library/heatmap.png \
    --positions 67,230,233,238
```

## 4. MD Simulations (Post-analysis)
MD simulations were performed as described in the methods section, based on previously blished code (see Reference [4] - DFX paper). The post-analysis calculates binding ergies and maps the energy landscape.

**Identify the top 100 frames:**
For all variants in the two systems (compact vs loose), the top 100 frames were entified based on the distance of the ligation site to the catalytic site.
```bash
python scripts/identify_top_100_frames.py
```

**Calculate binding energies:**
The resulting conformational ensembles were used as input for binding energy calculation ilizing PyRosetta.
```bash
python scripts/calculate_binding_energy.py
```

# References

If you utilize this code, please cite:
Sumire Honda Malca, Peter Stockinger, Daniela Milbredt, Nadine Duss, Michael Niklaus, David Patsch, Lisa Schelbert, Nicolas Imstepf, Kaja Stalder, Irene Marzuoli, Steven Paul Hanlon, Stefan G. Koenig, Filippo Sladojevich, Hans Iding, Serena Bisagni, and Rebecca Buller. "Evolved RNA ligase 2 facilitates therapeutic RNA oligonucleotide manufacturing." [INSERT PUBLICATION VENUE/PREPRINT SERVER, YEAR, DOI].
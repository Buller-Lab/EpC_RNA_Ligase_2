#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ------------------------------------------------------------
# "THE BEERWARE LICENSE" (Revision 42):
# <so@g.harvard.edu> and <pkk382@g.harvard.edu> wrote this code.
# As long as you retain this notice, you can do whatever you want
# with this stuff. If we meet someday, and you think this stuff
# is worth it, you can buy us a beer in return.
# --Sergey Ovchinnikov and Peter Koo
#
# This implementation comes from Peter Stockinger (peter.stockinger@roche.com).
# The beerware licence now extends to him and he would love to buy and drink
# a beer for/with Sergey and Peter Koo.
# ------------------------------------------------------------

import argparse
import sys
import os
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize
from scipy.spatial.distance import pdist, squareform
from Bio import SeqIO, AlignIO
from openpyxl import Workbook

# --- Global Configurations ---
alphabet = "ARNDCQEGHILKMFPSTWYV-"
states = len(alphabet)
a2n = {a: n for n, a in enumerate(alphabet)}

def aa2num(aa):
    return a2n.get(aa, a2n['-'])

# --- Alignment Data Processing ---
def filt_gaps(msa, gap_cutoff=0.5):
    tmp = (msa == states - 1).astype(float)
    non_gaps = np.where(np.sum(tmp.T, -1).T / msa.shape[0] < gap_cutoff)[0]
    return msa[:, non_gaps], non_gaps

def get_eff(msa, eff_cutoff=0.8):
    msa_sm = 1.0 - squareform(pdist(msa, "hamming"))
    msa_w = (msa_sm >= eff_cutoff).astype(float)
    msa_w = 1 / np.sum(msa_w, -1)
    return msa_w

def mk_msa(seqs):
    msa_ori = np.array([[aa2num(aa) for aa in seq] for seq in seqs])
    msa, v_idx = filt_gaps(msa_ori, 0.5)
    msa_weights = get_eff(msa, 0.8)
    ncol = msa.shape[1]
    w_idx = v_idx[np.stack(np.triu_indices(ncol, 1), -1)]
    return {
        "msa_ori": msa_ori,
        "msa": msa,
        "weights": msa_weights,
        "neff": np.sum(msa_weights),
        "v_idx": v_idx,
        "w_idx": w_idx,
        "nrow": msa.shape[0],
        "ncol": ncol,
        "ncol_ori": msa_ori.shape[1]
    }

# --- Modernized GREMLIN Vector Engine (NumPy + SciPy Optimize) ---
def gremlin_numpy(msa_data, maxiter=100):
    ncol = msa_data["ncol"]
    nrow = msa_data["nrow"]
    neff = msa_data["neff"]
    msa = msa_data["msa"]
    weights = msa_data["weights"]

    # One-hot encoding of MSA sequences
    oh_msa = np.eye(states)[msa]  # Shape: (nrow, ncol, states)

    # Pre-calculate flat dimensions for packing vectors
    v_dim = ncol * states
    w_dim = ncol * states * ncol * states

    # Symmetrization wrapper
    def unpack_params(x):
        V = x[:v_dim].reshape((ncol, states))
        W_raw = x[v_dim:].reshape((ncol, states, ncol, states))
        # Enforce zero diagonal blocks and symmetry
        W = W_raw * (1.0 - np.eye(ncol)[:, None, :, None])
        W = 0.5 * (W + W.transpose((2, 3, 0, 1)))
        return V, W

    # Objective Loss & Gradient Evaluator
    def objective_func(x):
        V, W = unpack_params(x)

        # Vectorized Hamiltonian mapping: V + dot(OH_MSA, W)
        # Tensordot across axis (ncol, states)
        vw = V + np.tensordot(oh_msa, W, axes=((1, 2), (2, 3))) # Shape: (nrow, ncol, states)

        # Partition functions (local Z) and scores
        log_z = scipy_logsumexp(vw, axis=2) # Shape: (nrow, ncol)
        prob = np.exp(vw - log_z[:, :, None]) # Shape: (nrow, ncol, states)

        # Pseudo-Log-Likelihood
        h = np.sum(oh_msa * vw, axis=(1, 2))
        z = np.sum(log_z, axis=1)
        pll = np.sum((h - z) * weights)

        # L2 Regularization
        l2_v = 0.01 * np.sum(V**2)
        l2_w = 0.01 * np.sum(W**2) * 0.5 * (ncol - 1) * (states - 1)
        
        loss = -(pll / neff) + (l2_v + l2_w) / neff

        # Vectorized Gradient Calculation
        # Compute difference between observed counts and model expectations
        diff = (prob - oh_msa) * weights[:, None, None]
        
        grad_V = np.sum(diff, axis=0) / neff + (0.02 * V) / neff
        grad_W = np.tensordot(oh_msa.transpose((1, 2, 0)), diff, axes=((2,), (0,))) / neff
        grad_W = grad_W * (1.0 - np.eye(ncol)[:, None, :, None])
        grad_W = 0.5 * (grad_W + grad_W.transpose((2, 3, 0, 1)))
        grad_W = grad_W + (0.02 * W * 0.5 * (ncol - 1) * (states - 1)) / neff

        return loss, np.concatenate([grad_V.ravel(), grad_W.ravel()])

    # Helper function to avoid external dependencies
    def scipy_logsumexp(a, axis=None):
        a_max = np.amax(a, axis=axis, keepdims=True)
        out = np.log(np.sum(np.exp(a - a_max), axis=axis, keepdims=True))
        out += a_max
        return np.squeeze(out, axis=axis)

    # Initialize V via sequence profile frequencies
    pseudo_count = 0.01 * np.log(neff)
    counts = np.sum(oh_msa * weights[:, None, None], axis=0) + pseudo_count
    v_ini = np.log(counts)
    v_ini -= np.mean(v_ini, axis=-1, keepdims=True)

    x0 = np.concatenate([v_ini.ravel(), np.zeros(w_dim)])

    # Use High-Performance L-BFGS optimization framework
    res = minimize(objective_func, x0, method='L-BFGS-B', jac=True, options={'maxiter': maxiter})
    
    final_V, final_W = unpack_params(res.x)
    tri = np.triu_indices(ncol, 1)
    W_triur = final_W[tri[0], :, tri[1], :]

    return {"v": final_V, "w": W_triur, "v_idx": msa_data["v_idx"], "w_idx": msa_data["w_idx"]}

def normalize(x):
    x = stats.boxcox(x - np.amin(x) + 1.0)[0]
    return (x - np.mean(x)) / np.std(x)

def get_mtx(mrf):
    raw = np.sqrt(np.sum(mrf["w"][:, :-1, :-1]**2, axis=(1, 2)))
    raw_sq = squareform(raw)
    ap_sq = np.sum(raw_sq, axis=0, keepdims=True) * np.sum(raw_sq, axis=1, keepdims=True) / np.sum(raw_sq)
    apc = squareform(raw_sq - ap_sq, checks=False)
    return {"zscore": normalize(apc)}

# --- Runtime Entrypoint ---
def main():
    p = argparse.ArgumentParser(description="GREMLIN optimization from a pre-aligned FASTA file.")
    p.add_argument("alignment_input", help="Path to the pre-aligned input file (FASTA format).")
    p.add_argument("--ref", required=True, help="Path to the reference sequence file (FASTA format) to identify the target track.")
    p.add_argument("--outdir", required=True, help="Directory to save the resulting coevolution Excel sheet.")
    p.add_argument("--opt_iter", type=int, default=50, help="SciPy optimization iterations.")
    a = p.parse_args()

    os.makedirs(a.outdir, exist_ok=True)

    # Read the target/reference sequence structure
    refs = list(SeqIO.parse(a.ref, "fasta"))
    if not refs:
        print("Error: Reference file empty.", file=sys.stderr)
        sys.exit(1)
    ref = refs[0]
    rseq = str(ref.seq).upper()

    # Read the direct aligned input
    try:
        aln = AlignIO.read(a.alignment_input, "fasta")
    except Exception as e:
        print(f"Error reading alignment input: {str(e)}", file=sys.stderr)
        sys.exit(1)

    # Locate which index row in the input alignment matches our reference track
    ref_row = 0
    ref_found = False
    for i, rec in enumerate(aln):
        if rec.id == ref.id or str(rec.seq).upper().replace("-", "") == rseq:
            ref_row = i
            ref_found = True
            break
            
    if not ref_found:
        print("Warning: Direct reference matching sequence structural match was not explicitly pinpointed. Defaulting to index 0.", file=sys.stderr)

    # --- Modern NumPy Co-evolution Stage ---
    print("Executing native GREMLIN optimization engine (L-BFGS-B)...")
    seqs_str = [str(rec.seq).upper() for rec in aln]
    gremlin_msa_payload = mk_msa(seqs_str)
    
    mrf = gremlin_numpy(gremlin_msa_payload, maxiter=a.opt_iter)
    
    # Save raw W matrix
    np.save(os.path.join(a.outdir, "raw_W_matrix.npy"), mrf["w"])
    
    mtx = get_mtx(mrf)
    
    zscore_sq = squareform(mtx["zscore"])
    ref_aln = str(aln[ref_row].seq).upper()
    aln_col_to_ref_pos = {}
    ref_positions_labels = []
    pos = 0
    for c in range(len(ref_aln)):
        ch = ref_aln[c]
        if ch == "-":
            aln_col_to_ref_pos[c] = None
            continue
        pos += 1
        aln_col_to_ref_pos[c] = pos
        ref_positions_labels.append(f"{ch}{pos}")

    full_coev_matrix = np.zeros((pos, pos))
    v_idx = mrf["v_idx"]
    
    for i in range(len(v_idx)):
        for j in range(len(v_idx)):
            ref_i = aln_col_to_ref_pos[v_idx[i]]
            ref_j = aln_col_to_ref_pos[v_idx[j]]
            if ref_i is not None and ref_j is not None:
                full_coev_matrix[ref_i - 1, ref_j - 1] = zscore_sq[i, j]

    f5 = os.path.join(a.outdir, "coevolution.xlsx")
    wb_coev = Workbook()
    ws_coev = wb_coev.active
    ws_coev.title = "Coevolution_Zscores"
    
    ws_coev.append([""] + ref_positions_labels)
    for idx, label in enumerate(ref_positions_labels):
        row_vals = [round(val, 4) for val in full_coev_matrix[idx, :]]
        ws_coev.append([label] + row_vals)
    wb_coev.save(f5)

    print(f"Completed Successfully.\nCo-evolution Matrix written to: {f5}")

if __name__ == "__main__":
    main()
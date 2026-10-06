#!/usr/bin/env python3
import argparse
import sys
import os
import shutil
import tempfile
import subprocess
from collections import Counter
from Bio import SeqIO, AlignIO
import pandas as pd
from openpyxl import Workbook


def filter_length(records, ref_len, tol):
    if tol <= 0:
        return records
    min_l = int(ref_len * (1 - tol))
    max_l = int(ref_len * (1 + tol))
    return [r for r in records if min_l <= len(r.seq) <= max_l]


def cluster_mmseqs(in_fasta, prefix, min_id, cov, thr):
    tmp = tempfile.mkdtemp(prefix="mmseqs_")
    db = os.path.join(tmp, prefix)
    cmd = [
        "mmseqs", "easy-cluster", in_fasta, db, tmp,
        "--min-seq-id", str(min_id),
        "-c", str(cov),
        "--cov-mode", "1",
        "--threads", str(thr),
        "--remove-tmp-files", "1"
    ]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"MMseqs2 error:\n{e.stderr}", file=sys.stderr)
        shutil.rmtree(tmp, ignore_errors=True)
        sys.exit(1)
    rep = f"{db}_rep_seq.fasta"
    if not os.path.exists(rep):
        shutil.rmtree(tmp, ignore_errors=True)
        sys.exit(1)
    out = os.path.join(os.path.dirname(in_fasta), f"{prefix}_reps.fasta")
    shutil.move(rep, out)
    shutil.rmtree(tmp, ignore_errors=True)
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("input")
    p.add_argument("--ref", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--len_tol", type=float, default=0.10)
    p.add_argument("--mmseqs_id", type=float, default=0.95)
    p.add_argument("--mmseqs_cov", type=float, default=0.80)
    p.add_argument("--threads", type=int, default=8)
    a = p.parse_args()

    os.makedirs(a.outdir, exist_ok=True)

    hits = list(SeqIO.parse(a.input, "fasta"))
    refs = list(SeqIO.parse(a.ref, "fasta"))
    if not refs:
        sys.exit(1)
    ref = refs[0]
    rlen = len(ref.seq)

    seen = set()
    uniq = []
    for r in hits:
        s = str(r.seq).upper()
        if s not in seen:
            seen.add(s)
            uniq.append(r)

    flt = filter_length(uniq, rlen, a.len_tol)
    f1 = os.path.join(a.outdir, "01_filtered.fasta")
    SeqIO.write(flt, f1, "fasta")

    reps_path = cluster_mmseqs(f1, "02", a.mmseqs_id, a.mmseqs_cov, a.threads)
    reps = list(SeqIO.parse(reps_path, "fasta"))

    final = [ref]
    rseq = str(ref.seq).upper()
    for r in reps:
        if str(r.seq).upper() != rseq:
            final.append(r)

    f2 = os.path.join(a.outdir, "02_reps.fasta")
    SeqIO.write(final, f2, "fasta")

    f3 = os.path.join(a.outdir, "03_aln.fasta")
    cmd = ["mafft", "--auto", "--thread", str(a.threads), f2]
    try:
        res = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        with open(f3, "w") as fh:
            fh.write(res.stdout)
    except Exception as e:
        print(str(e), file=sys.stderr)
        sys.exit(1)

    aln = AlignIO.read(f3, "fasta")
    ref_row = 0
    for i, rec in enumerate(aln):
        if rec.id == ref.id or str(rec.seq).upper().replace("-", "") == rseq:
            ref_row = i
            break

    ref_aln = str(aln[ref_row].seq).upper()
    aa_list = list("ACDEFGHIKLMNPQRSTVWY") + ["-"]
    data = []
    pos = 0
    n = len(aln)
    for c in range(len(ref_aln)):
        ch = ref_aln[c]
        if ch == "-":
            continue
        pos += 1
        col = [str(rec.seq[c]).upper() for rec in aln]
        cnt = Counter(col)
        for aa in aa_list:
            pct = (cnt.get(aa, 0) / n * 100) if n > 0 else 0
            data.append({"pos": pos, "ref": ch, "aa": aa, "pct": round(pct, 1)})

    df = pd.DataFrame(data)
    wide = df.pivot(index="aa", columns="pos", values="pct")
    pos_list = list(range(1, pos + 1))
    ref_dict = df.drop_duplicates("pos").set_index("pos")["ref"].to_dict()
    ref_row_vals = [ref_dict.get(p, "-") for p in pos_list]

    f4 = os.path.join(a.outdir, "04_conservation.xlsx")
    wb = Workbook()
    ws = wb.active
    ws.title = "Conservation"
    ws.append([""] + pos_list)
    ws.append(["Ref"] + ref_row_vals)
    for aa in aa_list:
        vals = [round(wide.loc[aa, p], 1) if p in wide.columns else 0.0 for p in pos_list]
        ws.append([aa] + vals)
    wb.save(f4)

    print(f"Done. Alignment: {f3}")
    print(f"Conservation: {f4}")


if __name__ == "__main__":
    main()
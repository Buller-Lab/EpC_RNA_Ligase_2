#!/usr/bin/env python3
"""Rename RNA ligase ortholog FASTA headers from internal 'Excelzyme' IDs
(rLI1_NN / rLI2_NN) to the publication IDs defined in Table S3, dropping the
sequences marked 'excluded'.

The mapping is read from evolutionary_analysis/ortholog_id_mapping.csv
(columns: excelzyme, publication, number, species, accession).

Header numbers are matched by family + integer, tolerating the historical
suffixes used in the raw files, e.g. 'rLI2_8_mod_(+1)' -> rLI2_08 -> excluded,
'rLI2_28_(+1)' -> rLI2_28 -> Rnl2_27.

Usage:
    python scripts/apply_publication_ids.py <in.fasta> <out.fasta> [--dry-run]
"""
import argparse
import csv
import os
import re
import sys
from Bio import SeqIO

MAPPING_CSV = os.path.join(os.path.dirname(__file__), "..",
                           "evolutionary_analysis", "ortholog_id_mapping.csv")


def load_mapping(path):
    mp = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            mp[row["excelzyme"].strip()] = row["publication"].strip()
    return mp


def excelzyme_id(header):
    m = re.search(r"rLI([12])_(\d+)", header)
    if not m:
        return None
    return f"rLI{m.group(1)}_{int(m.group(2)):02d}"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("infile")
    ap.add_argument("outfile")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    mp = load_mapping(os.path.abspath(MAPPING_CSV))
    kept, dropped = [], []
    for rec in SeqIO.parse(a.infile, "fasta"):
        exc = excelzyme_id(rec.id)
        if exc is None or exc not in mp:
            sys.exit(f"Error: cannot map header '{rec.id}' (parsed '{exc}')")
        pub = mp[exc]
        if pub.lower() == "excluded":
            dropped.append((rec.id, exc))
            continue
        rec.id = pub
        rec.name = pub
        rec.description = ""
        kept.append(rec)

    print(f"{a.infile}: {len(kept)} kept, {len(dropped)} dropped")
    for old, exc in dropped:
        print(f"  dropped {old} ({exc}, excluded)")
    if not a.dry_run:
        SeqIO.write(kept, a.outfile, "fasta")
        print(f"  wrote {a.outfile}")


if __name__ == "__main__":
    main()

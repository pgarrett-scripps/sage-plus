#!/usr/bin/env python3
"""Species-entrapment FDP for Sage Plus peptide / protein / protein-group q-values.

A human-only sample is searched against human + yeast + E. coli. Any yeast or
E. coli identification is false. The combined estimator (Wen et al. 2025,
FDRBench) is FDP = N_e * (1 + 1/r) / (N_t + N_e), where r is the entrapment to
target database size ratio; N_e / (N_t + N_e) is a strict lower bound.

usage: entrapment.py FASTA results.sage.parquet [label]
"""

import json
import re
import sys

import duckdb

MASS = dict(G=57.02146, A=71.03711, S=87.03203, P=97.05276, V=99.06841, T=101.04768,
            C=160.03065, L=113.08406, I=113.08406, N=114.04293, D=115.02694, Q=128.05858,
            K=128.09496, E=129.04259, M=131.04049, H=137.05891, F=147.06841, R=156.10111,
            Y=163.06333, W=186.07931, U=150.95364, O=237.14773)
TARGET_SPECIES = ("human", "irt")


def read_fasta(path):
    name, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line)
    if name:
        yield name, "".join(seq)


def digest(seq, missed=1, lo=7, hi=50):
    sites = [0] + [i + 1 for i, a in enumerate(seq[:-1]) if a in "KR" and seq[i + 1] != "P"] + [len(seq)]
    for i in range(len(sites) - 1):
        for j in range(i + 1, min(i + 2 + missed, len(sites))):
            p = seq[sites[i]:sites[j]]
            if lo <= len(p) <= hi and all(a in MASS for a in p):
                m = sum(MASS[a] for a in p) + 18.01056
                if 500 <= m <= 5000:
                    yield p.replace("I", "L")


def species(accession):
    return accession.split(":", 1)[0]


def ratios(fasta):
    pep_species = {}
    for name, seq in read_fasta(fasta):
        sp = species(name) in TARGET_SPECIES
        for p in set(digest(seq)):
            pep_species[p] = pep_species.get(p, False) or sp
    n_t = sum(pep_species.values())
    n_e = len(pep_species) - n_t
    prots = [species(n) in TARGET_SPECIES for n, _ in read_fasta(fasta)]
    return {"peptide_r": n_e / n_t, "protein_r": (len(prots) - sum(prots)) / sum(prots)}


def fdp(n_t, n_e, r):
    n = n_t + n_e
    return {"targets": n_t, "entrapments": n_e, "total": n,
            "lower_bound": n_e / n if n else None,
            "combined": n_e * (1 + 1 / r) / n if n else None}


def is_target_list(proteins):
    return any(species(p) in TARGET_SPECIES for p in re.split("[;/]", proteins))


def main():
    fasta, parquet = sys.argv[1], sys.argv[2]
    r = ratios(fasta)
    con = duckdb.connect()
    con.execute(f"create view t as select * from read_parquet('{parquet}') where not is_decoy")
    out = {"ratios": r, "levels": {}}
    for q in (0.005, 0.01, 0.02, 0.05):
        level = {}
        rows = con.execute(f"select proteins from t where spectrum_q <= {q}").fetchall()
        t = sum(is_target_list(p) for (p,) in rows)
        level["psm"] = fdp(t, len(rows) - t, r["peptide_r"])
        rows = con.execute(
            f"select distinct peptide, proteins from t where peptide_q <= {q}").fetchall()
        t = sum(is_target_list(p) for _, p in rows)
        level["peptide"] = fdp(t, len(rows) - t, r["peptide_r"])
        rows = con.execute(
            f"select distinct proteins from t where num_proteins = 1 and protein_q <= {q}").fetchall()
        t = sum(is_target_list(p) for (p,) in rows)
        level["protein"] = fdp(t, len(rows) - t, r["protein_r"])
        rows = con.execute(
            f"select distinct protein_groups from t where num_protein_groups = 1 "
            f"and protein_groups is not null and protein_group_q <= {q}").fetchall()
        t = sum(is_target_list(p) for (p,) in rows)
        level["protein_group"] = fdp(t, len(rows) - t, r["protein_r"])
        out["levels"][str(q)] = level
    json.dump(out, sys.stdout, indent=1)
    print()


if __name__ == "__main__":
    main()

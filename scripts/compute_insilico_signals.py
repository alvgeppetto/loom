#!/usr/bin/env python3
"""Compute additional in silico characterization for all 52 novel CRISPR targets.

Signals computed:
  1. Homopolymer runs (max run length per nucleotide)
  2. RNAfold MFE (ViennaRNA) for all 52
  3. Melting temperature (nearest-neighbor method, SantaLucia 1998)

Output: data/crispr_guides/insilico_characterization_52.json
"""
import json
import math
import re
import subprocess
import sys
from pathlib import Path

DATA = Path("data/crispr_guides")

# ── Load targets ──────────────────────────────────────────────────────
targets = json.load(open(DATA / "novel_targets_52.json"))["regions"]
same_corpus = json.load(open(DATA / "novel_targets_same_corpus.json"))
gc_data = json.load(open(DATA / "gc_content_52.json"))

# Build lookup by position for same-corpus conservation
# same_corpus is a flat list with guide_id containing the position
sc_lookup = {}
for e in same_corpus:
    # Extract position from guide_id like "ORF1b (nsp12-16, RdRp)_16813"
    pos_str = e["guide_id"].rsplit("_", 1)[-1]
    sc_lookup[int(pos_str)] = e
gc_lookup = {e["position"]: e for e in gc_data["targets"]}


# ── 1. Homopolymer runs ──────────────────────────────────────────────
def max_homopolymer(seq: str) -> dict:
    """Return max run length for each nucleotide and overall."""
    result = {}
    for base in "ACGT":
        runs = re.findall(f"{base}+", seq)
        result[base] = max((len(r) for r in runs), default=0)
    result["max_run"] = max(result.values())
    result["max_base"] = max(result, key=result.get)
    return result


# ── 2. Melting temperature (nearest-neighbor, SantaLucia 1998) ────────
# ΔH in cal/mol, ΔS in cal/(mol·K) for DNA/DNA duplexes
NN_DH = {
    "AA": -7900, "TT": -7900, "AT": -7200, "TA": -7200,
    "CA": -8500, "TG": -8500, "GT": -8400, "AC": -8400,
    "CT": -7800, "AG": -7800, "GA": -8200, "TC": -8200,
    "CG": -10600, "GC": -9800, "GG": -8000, "CC": -8000,
}
NN_DS = {
    "AA": -22.2, "TT": -22.2, "AT": -20.4, "TA": -21.3,
    "CA": -22.7, "TG": -22.7, "GT": -22.4, "AC": -22.4,
    "CT": -21.0, "AG": -21.0, "GA": -22.2, "TC": -22.2,
    "CG": -27.2, "GC": -24.4, "GG": -19.9, "CC": -19.9,
}
# Initiation parameters
INIT_DH = {"G": 100, "C": 100, "A": 2300, "T": 2300}
INIT_DS = {"G": -2.8, "C": -2.8, "A": 4.1, "T": 4.1}


def calc_tm(seq: str, oligo_conc_nM: float = 250.0, na_conc_M: float = 0.05) -> float:
    """Nearest-neighbor Tm for a short oligonucleotide."""
    seq = seq.upper()
    dH = INIT_DH[seq[0]] + INIT_DH[seq[-1]]  # initiation
    dS = INIT_DS[seq[0]] + INIT_DS[seq[-1]]

    for i in range(len(seq) - 1):
        pair = seq[i:i + 2]
        dH += NN_DH.get(pair, -8000)  # fallback for ambiguous
        dS += NN_DS.get(pair, -22.0)

    # Salt correction (Owczarzy 2004 simplified)
    dS += 0.368 * (len(seq) - 1) * math.log(na_conc_M)

    # Self-complementary check
    complement = str.maketrans("ACGT", "TGCA")
    is_self_comp = seq == seq[::-1].translate(complement)

    R = 1.987  # cal/(mol·K)
    if is_self_comp:
        ct = oligo_conc_nM * 1e-9
    else:
        ct = oligo_conc_nM * 1e-9 / 4.0

    tm_K = dH / (dS + R * math.log(ct))
    return round(tm_K - 273.15, 1)


# ── 3. RNAfold MFE for all 52 ────────────────────────────────────────
def batch_rnafold(sequences: list[str]) -> list[float]:
    """Run RNAfold on all sequences using ViennaRNA Python API."""
    import RNA
    mfes = []
    for seq in sequences:
        (_structure, mfe) = RNA.fold(seq)
        mfes.append(round(mfe, 1))
    return mfes


# ── Run everything ────────────────────────────────────────────────────
sequences = [t["best_seq"] for t in targets]
positions = [t["best_pos"] for t in targets]

# Homopolymer analysis
print("Computing homopolymer runs...", flush=True)
homopolymers = [max_homopolymer(seq) for seq in sequences]

# Melting temperatures
print("Computing melting temperatures...", flush=True)
tms = [calc_tm(seq) for seq in sequences]

# RNAfold
print("Running RNAfold on all 52 targets...", flush=True)
mfes = batch_rnafold(sequences)
assert len(mfes) == 52, f"Expected 52 MFE values, got {len(mfes)}"

# ── Assemble results ─────────────────────────────────────────────────
results = []
for i, t in enumerate(targets):
    pos = t["best_pos"]
    sc = sc_lookup.get(pos, {})
    gc = gc_lookup.get(pos, {})

    entry = {
        "rank": i + 1,
        "position": pos,
        "sequence": t["best_seq"],
        "gene": t["gene"],
        "pam": t["pam_types"][0] if t["pam_types"] else "unknown",
        "conservation_pct": sc.get("conservation_pct", t.get("ncov_conservation_pct")),
        "gc_pct": gc.get("gc_pct"),
        "mfe_kcal": mfes[i],
        "tm_celsius": tms[i],
        "max_homopolymer_run": homopolymers[i]["max_run"],
        "max_homopolymer_base": homopolymers[i]["max_base"],
        "homopolymer_detail": {
            b: homopolymers[i][b] for b in "ACGT"
        },
    }
    results.append(entry)

# ── Summary statistics ────────────────────────────────────────────────
poly_flagged = [r for r in results if r["max_homopolymer_run"] >= 4]
mfe_flagged = [r for r in results if r["mfe_kcal"] <= -5.0]
tm_range = (min(r["tm_celsius"] for r in results), max(r["tm_celsius"] for r in results))

output = {
    "description": "In silico characterization of all 52 novel SARS-CoV-2 CRISPR targets",
    "date": "2026-03-08",
    "signals": ["homopolymer_runs", "rnafold_mfe", "melting_temperature"],
    "conditions": {
        "rnafold": "ViennaRNA RNAfold, default parameters, --noPS",
        "tm": "nearest-neighbor (SantaLucia 1998), 250 nM oligo, 50 mM Na+",
        "homopolymer": "max consecutive identical nucleotides"
    },
    "summary": {
        "total_targets": 52,
        "homopolymer_ge4": len(poly_flagged),
        "homopolymer_flagged_positions": [r["position"] for r in poly_flagged],
        "mfe_le_minus5": len(mfe_flagged),
        "tm_range_celsius": list(tm_range),
        "tm_mean_celsius": round(sum(r["tm_celsius"] for r in results) / len(results), 1),
    },
    "targets": results,
}

outpath = DATA / "insilico_characterization_52.json"
json.dump(output, open(outpath, "w"), indent=2)
print(f"\nWrote {outpath}")

# ── Print summary table ──────────────────────────────────────────────
print(f"\n{'Rank':>4} {'Pos':>6} {'Gene':<22} {'GC%':>4} {'MFE':>6} {'Tm':>5} {'MaxRun':>6} {'Base':>4}")
print("-" * 70)
for r in results[:12]:  # show top 12
    print(f"{r['rank']:4d} {r['position']:6d} {r['gene']:<22} "
          f"{r['gc_pct']:4.0f} {r['mfe_kcal']:6.1f} {r['tm_celsius']:5.1f} "
          f"{r['max_homopolymer_run']:6d} {r['max_homopolymer_base']:>4}")

print(f"\nHomopolymer ≥4: {len(poly_flagged)} targets")
for r in poly_flagged:
    print(f"  pos {r['position']} ({r['gene']}): {r['max_homopolymer_run']}× {r['max_homopolymer_base']} "
          f"in {r['sequence']}")
print(f"\nMFE ≤ -5.0: {len(mfe_flagged)} targets")
print(f"Tm range: {tm_range[0]}–{tm_range[1]}°C")

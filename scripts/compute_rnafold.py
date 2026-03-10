#!/usr/bin/env python3
"""Run RNAfold on top 8 CRISPR candidates."""
import json
import RNA

top8 = [
    (10456, "ACTATTAAGGGTTCATTCCT", "ORF1a", "Cas12a"),
    (17431, "CTGCTCAATTACCTGCACCA", "nsp13", "SpCas9"),
    (17561, "TGCTGAAATTGTTGACACTG", "nsp13", "SpCas9"),
    (19218, "ATTGTTTGTAGATTTGACAC", "nsp14", "SpCas9"),
    (16813, "GAGAGTACACCTTTGAAAAA", "nsp13", "SpCas9"),
    (10531, "TGTTACATGCACCATATGGA", "ORF1a", "Cas12a"),
    (17078, "TGGTACTGGTAAGAGTCATT", "nsp13", "SpCas9"),
    (12888, "ACCTTGTAGGTTTGTTACAG", "ORF1a", "SpCas9"),
]

results = []
print("=== RNAfold for Top 8 Candidates ===")
print(f"{'Pos':>6} {'Protein':<8} {'PAM':<7} {'GC%':>4} {'MFE':>7} {'Structure':<25}")
for pos, seq, protein, pam in top8:
    gc = (seq.count("G") + seq.count("C")) / len(seq) * 100
    (ss, mfe) = RNA.fold(seq)
    (ens_ss, ens_dG) = RNA.pf_fold(seq)
    print(f"{pos:>6} {protein:<8} {pam:<7} {gc:4.0f} {mfe:>7.2f} {ss:<25}")
    results.append({
        "position": pos,
        "sequence": seq,
        "protein": protein,
        "pam": pam,
        "gc_pct": round(gc, 1),
        "mfe_structure": ss,
        "mfe_kcal": round(mfe, 2),
        "ensemble_dG": round(ens_dG, 2),
        "accessible": mfe > -5.0,
    })

accessible = sum(1 for r in results if r["accessible"])
print(f"\n{accessible}/8 predicted accessible (MFE > -5.0 kcal/mol)")

with open("data/crispr_guides/rnafold_top8.json", "w") as f:
    json.dump({"description": "RNAfold MFE prediction for top 8 SARS-CoV-2 CRISPR targets",
               "tool": "ViennaRNA RNAfold",
               "results": results}, f, indent=2)
print("Wrote rnafold_top8.json")

#!/usr/bin/env python3
"""Audit all gap claims in the PubMed scan results."""
import json

with open("web/data/pubmed-scan-results.json") as f:
    data = json.load(f)

print("=== ALL GAP CLAIMS ===\n")

print("--- APPLICATION GAPS (20 claimed) ---")
for opp in data["opportunities"]:
    if opp["type"] == "application_gap":
        print(f"  {opp['pathogen']:30s} | {opp['area']}")

print("\n--- GENE GAPS (43 claimed) ---")
for opp in data["opportunities"]:
    if opp["type"] == "gene_gap":
        print(f"  {opp['pathogen']:30s} | {opp['gene']}")

print("\n--- PATHOGEN SUMMARIES ---")
for key, ps in data["pathogens"].items():
    n_app_gaps = len(ps.get("gap_applications", []))
    n_gene_gaps = len(ps.get("unstudied_genes", []))
    print(f"  {ps['name']:30s} | papers={ps['total_papers']:>5} | app_gaps={n_app_gaps} | gene_gaps={n_gene_gaps}")
    if ps.get("unstudied_genes"):
        print(f"    unstudied: {', '.join(ps['unstudied_genes'])}")
    if ps.get("gap_applications"):
        print(f"    app_gaps:  {', '.join(ps['gap_applications'])}")

print("\n--- CRITICAL QUERIES TO AUDIT ---")
# Show the actual PubMed queries used for landscape
for key, land in data["landscape"].items():
    if land["total_papers"] == 0:
        print(f"  ZERO-PAPER CLAIM: {key}")
        print(f"    Query: {land['query']}")

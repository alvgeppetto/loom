#!/usr/bin/env python3
"""
Artifact-backed claim verification for the merged CRISPR paper.

Checks that every verifiable number in the manuscript matches
the source artifact. Run as a pre-build gate.

Usage:
    python3 scripts/verify-paper-claims.py publications/papers/pangenomic-crispr-targets-merged.md

Exit codes:
    0 = all checks pass
    1 = one or more FAIL
"""
import json
import re
import sys
from pathlib import Path

FAIL = 0


def check(label: str, condition: bool, detail: str = ""):
    global FAIL
    if condition:
        print(f"  PASS  {label}")
    else:
        print(f"  FAIL  {label}")
        if detail:
            print(f"        {detail}")
        FAIL += 1


def main():
    if len(sys.argv) < 2:
        print("Usage: verify-paper-claims.py <paper.md>")
        sys.exit(2)

    paper = Path(sys.argv[1])
    if not paper.exists():
        print(f"ERROR: {paper} not found")
        sys.exit(2)

    text = paper.read_text()
    root = Path(".")

    print(f"\n=== Artifact-Backed Claim Verification ===")
    print(f"Paper: {paper}\n")

    # ── 1. Denominator consistency ──────────────────────────────────────
    print("[1] Denominator consistency")
    canon = len(re.findall(r"9,193,298|9193298", text))
    legacy = len(re.findall(r"9,419,528|9419528", text))
    check("Canon denominator present", canon > 0, f"Found {canon} occurrences")
    # Allow exactly 2 provenance mentions in §2.5 (dual-FASTA disclosure)
    check("No legacy denominator (except provenance)", legacy <= 2, f"Found {legacy} occurrences (max 2 for provenance)")
    check("No stale dedup phrasing",
          "stricter header-level genome deduplication" not in text)

    # ── 2. Same-corpus revalidation artifact ────────────────────────────
    print("\n[2] Same-corpus artifact")
    sc_path = root / "data/crispr_guides/novel_targets_same_corpus.json"
    if sc_path.exists():
        with open(sc_path) as f:
            sc = json.load(f)
        sc_denoms = set(t["total_genomes"] for t in sc)
        check("Same-corpus JSON denominator = 9,193,298",
              sc_denoms == {9193298},
              f"Found denominators: {sc_denoms}")
        check("Same-corpus has 52 entries", len(sc) == 52, f"Got {len(sc)}")
    else:
        check("Same-corpus artifact exists", False, f"{sc_path} missing")

    # ── 3. Novelty distance ─────────────────────────────────────────────
    print("\n[3] Novelty distance")
    nr_path = root / "data/crispr_guides/novelty_report.json"
    if nr_path.exists():
        with open(nr_path) as f:
            nr = json.load(f)
        dists = [t["closest_distance"] for t in nr["targets"]]
        min_dist = min(dists)
        # Extract what the paper claims: "≥N mismatches"
        m = re.search(r"closest match.*?≥(\d+)\s*mismatch", text)
        if m:
            claimed = int(m.group(1))
            check(f"Novelty distance claim (≥{claimed}) matches artifact (min={min_dist})",
                  claimed <= min_dist,
                  f"Paper says ≥{claimed}, artifact min = {min_dist}")
        else:
            check("Novelty distance claim found in text", False, "Could not parse ≥N pattern")

        check(f"52 novel targets in novelty report",
              nr["novel_count"] == 52,
              f"Got {nr['novel_count']}")
    else:
        check("Novelty report exists", False, f"{nr_path} missing")

    # ── 4. Gene gap count ───────────────────────────────────────────────
    print("\n[4] Gene gap count")
    ps_path = root / "web/data/pubmed-scan-v4.json"
    if ps_path.exists():
        with open(ps_path) as f:
            ps = json.load(f)
        gaps = ps.get("confirmed_gene_gaps", [])

        # Deduplicate: J2R=thymidine kinase, A56R=hemagglutinin
        unique = set()
        for g in gaps:
            gene = g["gene"]
            pathogen = g["pathogen"]
            if gene in ("J2R", "thymidine kinase"):
                unique.add(("mpox", "J2R"))
            elif gene in ("A56R", "hemagglutinin"):
                unique.add(("mpox", "A56R"))
            else:
                unique.add((pathogen, gene))

        # Extract paper claim: "N unique gene targets" or "N genuine gene-level gaps"
        m_gaps = re.search(r"\*\*(\d+)\s+genuine gene-level\s+gaps\*\*", text)
        m_unique = re.search(r"\((\d+)\s+unique\s+gene\s+targets?\)", text)

        if m_gaps:
            claimed_raw = int(m_gaps.group(1))
            # Raw gaps = 8 query experiments, 6 after dedup in the body table
            check(f"Gene-level gaps raw count in text ({claimed_raw}) ≤ artifact ({len(gaps)})",
                  claimed_raw <= len(gaps),
                  f"Paper says {claimed_raw}, artifact has {len(gaps)} entries")

        if m_unique:
            claimed_unique = int(m_unique.group(1))
            check(f"Unique gene targets ({claimed_unique}) matches deduplicated ({len(unique)})",
                  claimed_unique == len(unique),
                  f"Paper says {claimed_unique}, artifact deduplicates to {len(unique)}")
    else:
        check("PubMed scan v4 artifact exists", False, f"{ps_path} missing")

    # ── 5. Cross-reactivity (NGG + Cas12a) ──────────────────────────────
    print("\n[5] Cross-reactivity")
    cr_path = root / "web/data/cross-reactivity.json"
    cr12_path = root / "web/data/cross-reactivity-cas12a.json"
    if cr_path.exists():
        with open(cr_path) as f:
            cr = json.load(f)
        summary = cr["summary"]

        # NGG paper claims: "1,097 of 1,100" and "99.7%"
        check("NGG CR total in paper",
              "1,100" in text, "Expected '1,100' in paper")
        check("NGG CR specific in paper",
              "1,097" in text, "Expected '1,097' in paper")
        check("NGG CR rate 99.7% in paper",
              "99.7%" in text, "Expected '99.7%' in paper")

        # Verify against NGG artifact
        check(f"NGG artifact total = {summary['total_guides_checked']}",
              summary["total_guides_checked"] == 1100)
        check(f"NGG artifact specific = {summary['specific_guides']}",
              summary["specific_guides"] == 1097)
        check(f"NGG artifact rate = {summary['specificity_rate']}%",
              summary["specificity_rate"] == 99.7)
    else:
        check("NGG cross-reactivity artifact exists", False, f"{cr_path} missing")

    if cr12_path.exists():
        with open(cr12_path) as f:
            cr12 = json.load(f)
        s12 = cr12["summary"]

        # Cas12a paper claims: "970 of 972" and "99.8%"
        check("Cas12a CR total in paper",
              "972" in text, "Expected '972' in paper")
        check("Cas12a CR specific in paper",
              "970" in text, "Expected '970' in paper")
        check("Cas12a CR rate 99.8% in paper",
              "99.8%" in text, "Expected '99.8%' in paper")

        # Verify against Cas12a artifact
        check(f"Cas12a artifact total = {s12['total_guides_checked']}",
              s12["total_guides_checked"] == 972)
        check(f"Cas12a artifact specific = {s12['specific_guides']}",
              s12["specific_guides"] == 970)
        check(f"Cas12a artifact rate = {s12['specificity_rate']}%",
              s12["specificity_rate"] == 99.8)

        # Combined: paper claims "2,067 of 2,072" and "99.8%"
        check("Combined CR total in paper",
              "2,072" in text, "Expected '2,072' in paper")
        check("Combined CR specific in paper",
              "2,067" in text, "Expected '2,067' in paper")
        check("Combined CR rate 99.8% (combined) in paper",
              "99.8%" in text, "Expected '99.8%' in paper")
    else:
        check("Cas12a cross-reactivity artifact exists", False, f"{cr12_path} missing")

    # ── 6. Published guide conservation spot-checks ─────────────────────
    print("\n[6] Published guide conservation")
    pg_path = root / "data/crispr_guides/published_guides_same_corpus.json"
    if pg_path.exists():
        with open(pg_path) as f:
            pg = json.load(f)
        # Spot-check DETECTR_N and SHERLOCK_S
        pg_map = {t["guide_id"]: t["conservation_pct"] for t in pg}
        if "DETECTR_N" in pg_map:
            check(f"DETECTR_N conservation = {pg_map['DETECTR_N']:.2f}%",
                  "94.21" in text or f"{pg_map['DETECTR_N']:.2f}" in text,
                  f"Artifact: {pg_map['DETECTR_N']:.2f}%, expected in paper")
        if "SHERLOCK_S" in pg_map:
            check(f"SHERLOCK_S conservation = {pg_map['SHERLOCK_S']:.2f}%",
                  "63.28" in text or f"{pg_map['SHERLOCK_S']:.2f}" in text,
                  f"Artifact: {pg_map['SHERLOCK_S']:.2f}%, expected in paper")
    else:
        check("Published guides same-corpus artifact exists", False, f"{pg_path} missing")

    # ── 7. No placeholder/TODO markers ──────────────────────────────────
    print("\n[7] No placeholders")
    for marker in ["TODO", "FIXME", "TBD", "XXX", "PLACEHOLDER"]:
        count = len(re.findall(rf"\b{marker}\b", text))
        check(f"No '{marker}' markers", count == 0, f"Found {count}")

    # ── 8. Section numbering integrity ──────────────────────────────────
    print("\n[8] Section numbering")
    headings = re.findall(r"^###?\s+(\d+\.\d+)\s", text, re.MULTILINE)
    # Check no gaps in numbering within each major section
    by_major = {}
    for h in headings:
        major, minor = h.split(".")
        by_major.setdefault(major, []).append(int(minor))
    for major, minors in by_major.items():
        expected = list(range(1, max(minors) + 1))
        check(f"Section {major}.x numbering contiguous",
              sorted(set(minors)) == expected,
              f"Found {sorted(set(minors))}, expected {expected}")

    # ── 9. All Table S references resolve ───────────────────────────────
    print("\n[9] Table references")
    # Match both "Table S9" and "Table S9a" / "Tables S9a–b" style
    refs = set(re.findall(r"Table[s]?\s+S(\d+)", text, re.IGNORECASE))
    defs = set(re.findall(r"^###\s+S(\d+)\.", text, re.MULTILINE))
    for ref in sorted(refs, key=int):
        check(f"Table S{ref} has a definition section",
              ref in defs,
              f"S{ref} referenced but no '### S{ref}.' heading found")

    # ── 10. nsp10 not grouped with "0 PubMed" ──────────────────────────
    print("\n[10] nsp10 PubMed accuracy")
    # nsp10 has 1 PubMed paper (artifact-verified); it must NOT be directly
    # attributed "0 PubMed" (e.g., "nsp7 and nsp10 ... (0 PubMed ...")
    # Allow nsp10 on the same line if it's clearly in a different group.
    nsp10_zero = bool(re.search(
        r"nsp10[^.]*?\(0 PubMed|"             # nsp10 followed by (0 PubMed in same sentence
        r"\(0 PubMed[^)]*nsp10",              # (0 PubMed ... nsp10) inside parens
        text
    ))
    check("nsp10 NOT claimed as 0 PubMed (artifact: 1 paper)",
          not nsp10_zero,
          "nsp10 has 1 PubMed paper (narrow+broad) per v4 artifact")

    # ── 11. S4 Hamming distance matches artifact ───────────────────────
    print("\n[11] S4 Hamming distance")
    if nr_path.exists():
        s4_match = re.search(r"minimum\s+Hamming\s+distance\s+observed\s+was\s+(\d+)", text)
        if s4_match:
            claimed_hamming = int(s4_match.group(1))
            check(f"S4 minimum Hamming distance ({claimed_hamming}) = artifact min ({min_dist})",
                  claimed_hamming == min_dist,
                  f"Paper says {claimed_hamming}, artifact min = {min_dist}")
        else:
            check("S4 Hamming distance claim parseable", False)

    # ── 12. nsp14 position range matches Table S2 ──────────────────────
    print("\n[12] nsp14 position range")
    # Extract nsp14 positions from Table S2 in the paper itself
    s2_nsp14_positions = [
        int(m.group(1))
        for m in re.finditer(r"\|\s*\d+\s*\|\s*(\d+)\s*\|.*?nsp14", text)
    ]
    if s2_nsp14_positions:
        s2_min = min(s2_nsp14_positions)
        s2_max = max(s2_nsp14_positions)
        # Check prose range: "candidate windows (positions X–Y)" near nsp14 context
        range_match = re.search(
            r"(?:candidate\s+)?windows\s*\(positions\s+([\d,]+)[–-]([\d,]+)\)",
            text,
        )
        if range_match:
            lo = int(range_match.group(1).replace(",", ""))
            hi = int(range_match.group(2).replace(",", ""))
            check(f"nsp14 range low ({lo:,}) matches S2 min ({s2_min:,})",
                  lo == s2_min,
                  f"Prose says {lo:,}, table has {s2_min:,}")
            check(f"nsp14 range high ({hi:,}) matches S2 max ({s2_max:,})",
                  hi == s2_max,
                  f"Prose says {hi:,}, table has {s2_max:,}")
        else:
            check("nsp14 position prose parseable", False)

    # ── 13. Cross-reactivity computational verification ───────────────
    print("\n[13] Cross-reactivity (computational)")
    cr_path = root / "web/data/cross-reactivity.json"
    n52_path = root / "data/crispr_guides/novel_targets_52.json"
    if cr_path.exists():
        with open(cr_path) as f:
            cr_data = json.load(f)
        cr_pathogens = cr_data["pathogens"]
        cr_summary = cr_data["summary"]

        # All CR guides must be 23-mers ending in GG (NGG PAM)
        all_cr_seqs = [g["seq"] for guides in cr_pathogens.values() for g in guides]
        all_ngg = all(len(s) == 23 and s.endswith("GG") for s in all_cr_seqs)
        check("CR guides are all NGG 23-mers (artifact)",
              all_ngg,
              f"{sum(1 for s in all_cr_seqs if not (len(s)==23 and s.endswith('GG')))} non-NGG guides found")

        # Verify SARS-CoV-2 is NOT in the CR screen
        check("SARS-CoV-2 not in CR pathogens (artifact)",
              "sars-cov-2" not in cr_pathogens,
              "CR screen should cover 11 non-SARS-CoV-2 pathogens only")

        # Paper claims 11 pathogens
        m_path = re.search(r"(\d+)\s+pathogens\)", text)
        if m_path:
            claimed_path = int(m_path.group(1))
            check(f"CR pathogen count ({claimed_path}) matches artifact ({len(cr_pathogens)})",
                  claimed_path == len(cr_pathogens))

        # Paper must report Cas12a TTTV screen results (not just acknowledge gap)
        check("§3.4 reports Cas12a TTTV screen results",
              bool(re.search(r"TTTV-PAM.*screen|Cas12a.*TTTV.*screen", text, re.DOTALL)),
              "§3.4 should report Cas12a cross-reactivity results")

    if n52_path.exists():
        with open(n52_path) as f:
            n52_data = json.load(f)
        regions = n52_data["regions"]
        n_cas9 = sum(1 for r in regions if "Cas9" in r["pam_types"])
        n_cas12a = sum(1 for r in regions if "Cas12a" in r["pam_types"])

        # Paper claims "24 carry SpCas9 NGG" and "28 carry Cas12a TTTV"
        check(f"Cas9 (NGG) count = {n_cas9} matches paper '24'",
              "24" in text and n_cas9 == 24,
              f"Artifact: {n_cas9}")
        check(f"Cas12a (TTTV) count = {n_cas12a} matches paper '28'",
              "28" in text and n_cas12a == 28,
              f"Artifact: {n_cas12a}")

    # ── 14. SARS-CoV-2 gap result explained in abstract ────────────────
    print("\n[14] SARS-CoV-2 gap result framing")
    # The abstract must convey that SARS-CoV-2 has no outright gene-level gaps
    # but nsp13/nsp14 have only non-diagnostic CRISPR work.
    abstract_end = text.find("## 1.")  # abstract ends at §1
    if abstract_end == -1:
        abstract_end = 2000  # fallback
    abstract_text = text[:abstract_end]
    check("Abstract explains SARS-CoV-2 gap result",
          bool(re.search(r"no confirmed gaps|no diagnostic guide designs|unexploited opportunit", abstract_text)),
          "Abstract must explain SARS-CoV-2 gap result (no confirmed gaps, nsp13/nsp14 unexploited)")
    # Must NOT say "confirms" in patent context (overstates exact-match search)
    check("No FALSE-like contradiction (diagnostic opportunity after FALSE)",
          "unexploited diagnostic opportunit" in text,
          "Paper should explain the FALSE + diagnostic opportunity nuance")

    # ── 15. Patent language not overstated ──────────────────────────────
    print("\n[15] Patent language")
    # Main text should NOT overclaim patent novelty
    check("No overclaimed patent novelty",
          "no patented CRISPR guide" not in text.lower(),
          "Main text should not claim 'no patented guide overlap' without comprehensive search")
    # Patent section should acknowledge limitations
    check("Patent limitation acknowledged",
          bool(re.search(r"freedom-to-operate|patent counsel|not constitute", text)),
          "Patent section should acknowledge it is not a freedom-to-operate analysis")

    # ── 16. Europe PMC ≤30 threshold justified ─────────────────────────
    print("\n[16] Threshold justification")
    # The ≤30 threshold separating CONFIRMED from UNCERTAIN must have a
    # justification sentence (e.g., referencing the RSV SH reclassification)
    threshold_ctx = text[text.find("EuropePMC ≤ 30"):text.find("EuropePMC ≤ 30") + 600] if "EuropePMC ≤ 30" in text else ""
    check("≤30 threshold has justification",
          bool(re.search(r"≤\s*30.*?(conserv|RSV|minimize|false.*classif|validat)", threshold_ctx, re.DOTALL)),
          "CONFIRMED threshold (≤30 Europe PMC) should include justification (e.g., RSV SH case)")

    # ── 17. Table 2 tie-breaking documented ────────────────────────────
    print("\n[17] Table 2 tie-breaking")
    # Ranks 2–5 are all 99.1%; the tie-breaking rule must be stated
    check("Table 2 tie-breaking rule stated",
          bool(re.search(r"tied.*rank|ordering.*tied|tied.*99\.1|ascending|arbitrary.*order", text, re.IGNORECASE)),
          "Table 2 has tied ranks (99.1%); must state the tie-breaking rule")

    # ── 18. Abstract word count ─────────────────────────────────────────
    print("\n[18] Abstract word count")
    abstract_m = re.search(r"## Abstract\s*\n(.*?)(?=\n## )", text, re.DOTALL)
    if abstract_m:
        abstract_words = len(abstract_m.group(1).split())
        check("Abstract ≤ 300 words",
              abstract_words <= 300,
              f"Abstract has {abstract_words} words (target ≤ 300)")
    else:
        check("Abstract section found", False, "Could not locate ## Abstract heading")

    # ── 19. GC content analysis present ─────────────────────────────────
    print("\n[19] GC content analysis")
    gc_artifact = root / "data" / "crispr_guides" / "gc_content_52.json"
    check("GC content artifact exists", gc_artifact.exists())
    check("GC content referenced in paper",
          bool(re.search(r"gc_content_52\.json|GC.*percent|GC%", text, re.IGNORECASE)),
          "Paper should reference GC content analysis")
    check("Table 2 has GC% column",
          bool(re.search(r"GC%.*MFE|GC%", text)),
          "Table 2 should include a GC% column")

    # ── 20. RNAfold analysis present ────────────────────────────────────
    print("\n[20] RNAfold analysis")
    rnafold_artifact = root / "data" / "crispr_guides" / "rnafold_top8.json"
    check("RNAfold artifact exists", rnafold_artifact.exists())
    check("RNAfold/MFE referenced in paper",
          bool(re.search(r"rnafold_top8\.json|RNAfold|MFE.*kcal", text, re.IGNORECASE)),
          "Paper should reference RNA secondary structure analysis")

    # ── 21. Low-conservation tier marked in S1 ──────────────────────────
    print("\n[21] Low-conservation tier separation")
    check("Exploratory tier label in S1",
          bool(re.search(r"[Ee]xploratory tier|< ?92%.*conservation", text)),
          "S1 should visually separate ranks 42–52 as exploratory tier")

    # ── 22. Guide length asymmetry explained ────────────────────────────
    print("\n[22] Guide length asymmetry")
    check("20-mer vs 23-mer explanation present",
          bool(re.search(r"20-mer.*23-mer|23-mer.*20-mer|Cas12a.*23|spacer.*length", text, re.IGNORECASE)),
          "Paper should explain why SARS-CoV-2 uses 20-mers and multi-pathogen uses 23-mers")

    # ── 23. Seed-region analysis referenced ─────────────────────────────
    print("\n[23] Seed-region analysis")
    seed_artifact = root / "data" / "crispr_guides" / "offtarget_seed_analysis_top8.json"
    check("Seed-region artifact exists", seed_artifact.exists())
    check("Seed-region analysis in off-target section",
          bool(re.search(r"seed.region.*analysis|seed.*mismatch|PAM-proximal seed", text, re.IGNORECASE)),
          "Off-target section should include seed-region analysis for top 8")

    # ── 24. Unifying narrative in introduction ──────────────────────────
    print("\n[24] Unifying narrative")
    intro_end = text.find("## 2.")  # intro ends at §2
    if intro_end == -1:
        intro_end = 5000
    intro_text = text[:intro_end]
    check("Two-contribution bridge in intro",
          bool(re.search(r"linked|complementary|generaliz", intro_text)),
          "Introduction should explain why SARS-CoV-2 deep-dive + multi-pathogen belong together")

    # ── 25. Convergence claim softened ──────────────────────────────────
    print("\n[25] Convergence claim")
    check("No unsupported 'continues to converge' claim",
          "continues to converge" not in text,
          "§4.2 should soften or date-qualify the convergence claim")

    # ── 26. Hassan preprint backed by peer-reviewed ref ─────────────────
    print("\n[26] Hassan preprint mitigation")
    check("Peer-reviewed ExoN reference present",
          bool(re.search(r"Ogando|exoribonuclease.*fidelity|replication fidelity", text)),
          "Hassan preprint should be backed by a peer-reviewed reference (e.g., Ogando et al.)")

    # ── 27. ORF10 precedent citation ────────────────────────────────────
    print("\n[27] ORF10 precedent")
    check("ORF10 diagnostic precedent cited",
          bool(re.search(r"Broughton.*non-coding|irrespective.*coding|SHERLOCK.*DETECTR.*genomic loci", text)),
          "ORF10 section should cite established practice of CRISPR detecting RNA at non-coding loci")

    # ── 28. Reference numbering sanity ──────────────────────────────────
    print("\n[28] Reference numbering")
    check("No [9a] style sub-references",
          "[9a]" not in text,
          "References must use sequential integers, not sub-letters like [9a]")
    # Check sequential numbering in reference list (after ## References heading)
    ref_section_start = text.find("## References")
    if ref_section_start != -1:
        ref_section = text[ref_section_start:ref_section_start + 3000]
        ref_nums = [int(m) for m in re.findall(r"^\[(\d+)\]", ref_section, re.MULTILINE)]
        if ref_nums:
            expected = list(range(ref_nums[0], ref_nums[0] + len(ref_nums)))
            check("Reference numbers are sequential",
                  ref_nums == expected,
                  f"Expected {expected}, got {ref_nums}")

    # ── 29. Seed-region uses PAM-proximal numbering ─────────────────────
    print("\n[29] Seed-region numbering convention")
    check("Seed positions use PAM-proximal numbering",
          bool(re.search(r"PAM-proximal position \d+", text)),
          "Seed-region analysis should specify PAM-proximal position numbering")
    check("No raw 0-indexed positions in seed text",
          "mismatch at position 19" not in text and "mismatch at position 8" not in text,
          "Seed analysis should not use raw 0-indexed positions from the artifact")

    # ── 30. GC low-activity citation present ────────────────────────────
    print("\n[30] GC activity citation")
    check("GC low-activity claim has citation",
          bool(re.search(r"below 40%.*\[\d+\]|GC.*\[\d+\].*activity|cleavage.*below.*GC.*\[\d+\]", text)),
          "Claim about Cas activity below 40% GC needs a citation")

    # ── 31. Position 19218 low GC flagged ───────────────────────────────
    print("\n[31] Low-GC candidate flagged")
    check("Position 19218 low GC noted",
          bool(re.search(r"19218.*30%|19218.*lowest.*GC|19218.*empirical", text, re.DOTALL)),
          "Position 19218 (30% GC) should be flagged as requiring empirical activity confirmation")

    # ── 32. Top-8 2mm hit counts reported ───────────────────────────────
    print("\n[32] Top-8 2mm hit reporting")
    check("2mm hit range for top 8 stated",
          bool(re.search(r"all 8.*2-mismatch|2-mismatch.*all 8|10.*32.*loci", text, re.IGNORECASE)),
          "Paper should report that all 8 top candidates have 2mm human hits (range 10-32)")

    # ── 33. Unique 23-mer exhaustiveness data ─────────────────────────
    print("\n[33] 10K cutoff exhaustiveness")
    u23_path = root / "data" / "crispr_guides" / "unique_23mer_counts.json"
    check("unique_23mer_counts.json artifact exists",
          u23_path.exists(),
          "Missing unique_23mer_counts.json artifact")
    if u23_path.exists():
        u23 = json.loads(u23_path.read_text())
        check("Zika unique 23-mers in artifact (>100K)",
              u23.get("zika", {}).get("unique_23mers", 0) > 100_000,
              "Zika should have >100K unique 23-mers")
        check("Ebola unique 23-mers in artifact (>100K)",
              u23.get("ebola", {}).get("unique_23mers", 0) > 100_000,
              "Ebola should have >100K unique 23-mers")
    check("Unique 23-mer counts referenced in paper",
          "unique_23mer_counts" in text,
          "Paper should reference unique_23mer_counts.json artifact")
    check("Zika exhaustiveness data in paper",
          bool(re.search(r"Zika.*208", text)),
          "Paper should cite Zika's ~208K unique 23-mers")
    check("Limitation 6 updated with conservation context",
          bool(re.search(r"conservation threshold|conservation.depth threshold", text)),
          "Limitation 6 should describe cutoff as conservation threshold, not arbitrary")

    # ── 34. ORF naming consistency ──────────────────────────────────────
    print("\n[34] ORF naming consistency")

    # ── 35. Extended in silico characterization ─────────────────────────
    # Must reference the full-52 characterization artifact
    print("\n[35] Extended in silico characterization")
    insilico_path = root / "data" / "crispr_guides" / "insilico_characterization_52.json"
    check("insilico_characterization_52.json artifact exists",
          insilico_path.exists(),
          "Missing data/crispr_guides/insilico_characterization_52.json")
    if insilico_path.exists():
        isc = json.loads(insilico_path.read_text())
        check("In silico artifact has 52 targets",
              len(isc.get("targets", [])) == 52,
              f"Expected 52, got {len(isc.get('targets', []))}")
        check("Homopolymer data present in artifact",
              all("max_homopolymer_run" in t for t in isc["targets"]),
              "Missing homopolymer data in artifact")
        check("Tm data present in artifact",
              all("tm_celsius" in t for t in isc["targets"]),
              "Missing Tm data in artifact")
    check("Paper references insilico_characterization_52.json",
          "insilico_characterization_52" in text,
          "Paper should reference insilico_characterization_52.json")
    check("RNAfold extended to all 52 in paper",
          bool(re.search(r"all 52 candidates.*(?:ViennaRNA|RNAfold|MFE)", text)),
          "Paper should mention RNAfold MFE computed for all 52")
    check("Homopolymer screening mentioned in paper",
          "homopolymer" in text.lower(),
          "Paper should discuss homopolymer runs")
    check("Melting temperature mentioned in paper",
          bool(re.search(r"[Mm]elting temperature|T~m~", text)),
          "Paper should report melting temperature")

    # ── 36. Extended in silico: accessibility + poly-T + patent seq ─────
    print("\n[36] Extended in silico (accessibility, poly-T, patent seq)")
    ext_path = root / "data" / "crispr_guides" / "extended_insilico_52.json"
    check("extended_insilico_52.json artifact exists",
          ext_path.exists(),
          "Missing data/crispr_guides/extended_insilico_52.json")
    if ext_path.exists():
        ext = json.loads(ext_path.read_text())
        check("Extended artifact has 52 targets",
              len(ext.get("targets", [])) == 52,
              f"Expected 52, got {len(ext.get('targets', []))}")
        check("Accessibility data in artifact",
              all("target_site_accessibility" in t for t in ext["targets"]),
              "Missing target_site_accessibility in artifact")
        check("Poly-T data in artifact",
              all("pol3_termination_risk" in t for t in ext["targets"]),
              "Missing pol3_termination_risk in artifact")
    check("Target-site accessibility mentioned in paper",
          bool(re.search(r"[Tt]arget.site.*accessibility|unpaired probabil", text)),
          "Paper should discuss target-site RNA accessibility")
    check("Poly-T termination risk mentioned in paper",
          bool(re.search(r"[Pp]ol.III.*terminat|poly.T.*terminat", text)),
          "Paper should discuss Pol III poly-T termination risk")
    check("Patent sequence-level check mentioned in paper",
          bool(re.search(r"sequence.level.*patent|nucleotide.*patent.*search|patent.*nucleotide", text)),
          "Paper should mention sequence-level patent check")

    # 34a. Lowercase 'orf' outside code blocks and published assay names
    # Strip fenced code blocks and inline code before checking
    text_no_code = re.sub(r'```.*?```', '', text, flags=re.DOTALL)
    text_no_code = re.sub(r'`[^`]+`', '', text_no_code)
    # Allow 'Orf' only inside published assay names like SHERLOCK_Orf1ab
    text_no_assay = re.sub(r'SHERLOCK_Orf1ab', '', text_no_code)
    lowercase_orf = re.findall(r'\b[Oo]rf\d', text_no_assay)
    check("ORF always uppercase (except published assay names)",
          len(lowercase_orf) == 0,
          f"Found lowercase/mixed-case ORF: {lowercase_orf[:5]}")

    # 34b. nsp13–16 are ORF1b genes; must never be labelled as ORF1a
    orf1b_nsps = ['nsp12', 'nsp13', 'nsp14', 'nsp15', 'nsp16']
    bad_orf1a_labels = []
    for m in re.finditer(r'ORF1a\s*\(([^)]+)\)', text):
        parenthetical = m.group(1)
        for nsp in orf1b_nsps:
            if nsp in parenthetical:
                bad_orf1a_labels.append(f"ORF1a({parenthetical}) contains {nsp}")
    check("No ORF1b genes (nsp12–16) labelled as ORF1a",
          len(bad_orf1a_labels) == 0,
          f"ORF1b genes mislabelled as ORF1a: {bad_orf1a_labels[:3]}")

    # 34c. ORF1a and nsp13/14 listed together must mention ORF1b
    # Catches the error: "ORF1a, nsp13, nsp14" without ORF1b, which implies
    # nsp13/14 are part of ORF1a (they're in ORF1b). We check:
    #  (i) Comma-adjacent patterns: "ORF1a, nsp13" or "ORF1a and nsp13"
    #  (ii) Abstract and intro sections must mention ORF1b if nsp13/14 present
    prose_lines = [line for line in text.split('\n')
                   if not line.strip().startswith('|')]
    prose_text = '\n'.join(prose_lines)

    # (i) Direct comma/and adjacency: "ORF1a, nsp13" or "ORF1a, nsp14"
    comma_adjacent = re.findall(
        r'ORF1a\b(?:,?\s+(?:and\s+)?)?nsp1[34]\b', prose_text)
    missing_orf1b = []
    for match in comma_adjacent:
        # Check ±80 chars around the match for ORF1b
        idx = prose_text.find(match)
        window = prose_text[max(0, idx - 80):idx + len(match) + 80]
        if 'ORF1b' not in window and 'ORF1a/b' not in window:
            missing_orf1b.append(match.strip()[:80])

    # (ii) Abstract and intro: if nsp13 or nsp14 mentioned, ORF1b must be too
    intro_end = text.find("## 2.")
    if intro_end == -1:
        intro_end = 3000
    front_matter = text[:intro_end]
    has_nsp13_14_front = bool(re.search(r'nsp1[34]\b', front_matter))
    has_orf1b_front = bool(re.search(r'ORF1b', front_matter))
    if has_nsp13_14_front and not has_orf1b_front:
        missing_orf1b.append("Abstract/intro mentions nsp13/14 without ORF1b")
    check("nsp13/nsp14 near ORF1a always accompanied by ORF1b",
          len(missing_orf1b) == 0,
          f"ORF1a + nsp13/14 without ORF1b mention: {missing_orf1b[:3]}")

    # 34d. Gene column in markdown tables: ORF must be uppercase
    table_gene_cols = re.findall(r'\|\s*[Oo]rf\d[^\|]*\|', text)
    bad_table_orf = [c for c in table_gene_cols
                     if re.search(r'(?<!\w)[Oo]rf\d', c) and 'ORF' not in c
                     and 'SHERLOCK_Orf' not in c]
    check("Gene columns in tables use uppercase ORF",
          len(bad_table_orf) == 0,
          f"Table gene columns with lowercase ORF: {bad_table_orf[:3]}")

    # 34e. "outperform" without conservation qualifier is overconfident
    outperform_matches = re.findall(r'outperform.*?guide', text, re.IGNORECASE)
    check("No unqualified 'outperform' claims about guides",
          len(outperform_matches) == 0,
          f"Found unqualified 'outperform': {outperform_matches[:3]}")

    # 34f. Any "match or exceed" / "exceed" about guides should mention
    #       "conservation" within ~120 chars
    exceed_positions = [m.start() for m in
                        re.finditer(r'(?:match or )?exceed.*?(?:published|diagnostic).*?guide',
                                    text, re.IGNORECASE)]
    unqualified_exceed = []
    for pos in exceed_positions:
        window = text[max(0, pos - 30):pos + 150]
        if 'conservation' not in window.lower():
            snippet = text[pos:pos + 80].replace('\n', ' ').strip()
            unqualified_exceed.append(snippet)
    check("'exceed...guides' claims mention conservation metrics",
          len(unqualified_exceed) == 0,
          f"Exceed claim without 'conservation' qualifier: {unqualified_exceed[:3]}")

    # ── Summary ─────────────────────────────────────────────────────────
    print(f"\n{'=' * 50}")
    if FAIL == 0:
        print(f"CLAIM VERIFICATION PASSED. All checks passed.")
    else:
        print(f"CLAIM VERIFICATION FAILED. {FAIL} check(s) failed.")
    print(f"{'=' * 50}\n")
    sys.exit(1 if FAIL else 0)


if __name__ == "__main__":
    main()

"""Tests for PubMed query builder logic in scripts/pubmed_scan.py.

These are CRITICAL correctness tests. Bugs in query construction silently
produce false gap claims in medical literature, which can tarnish the
credibility of the entire project. See docs/lessons/pubmed-query-bugs.md.

Bug history:
  - 2026-03-04: "Mpox OR Monkeypox" was treated as a literal phrase instead
    of boolean OR, causing 0 results when 69 papers existed.
  - 2026-03-04: Gene queries used only primary names, missing synonyms like
    "enoyl-ACP reductase" for inhA, producing false gene gaps.
  - 2026-03-04: Duplicate gene names (NP + nucleoprotein for same gene)
    inflated gap counts.
  - 2026-03-04: MERS "S protein" listed separately from "spike" — same gene.
"""

import sys
import os
import importlib.util

# Import pubmed_scan.py as a module
_scan_path = os.path.join(os.path.dirname(__file__), "..", "scripts", "pubmed_scan.py")
_spec = importlib.util.spec_from_file_location("pubmed_scan", os.path.abspath(_scan_path))
pubmed_scan = importlib.util.module_from_spec(_spec)
# Prevent the module from executing main() on import
sys.modules["pubmed_scan"] = pubmed_scan
_spec.loader.exec_module(pubmed_scan)

# ---------------------------------------------------------------------------
# _name_clause: handles OR-separated pathogen names
# ---------------------------------------------------------------------------

class TestNameClause:
    """Verify _name_clause correctly splits OR-separated names into boolean PubMed queries."""

    def test_simple_name(self):
        result = pubmed_scan._name_clause("Ebola")
        assert result == '("Ebola"[Title/Abstract])'

    def test_or_separated_names(self):
        result = pubmed_scan._name_clause("Mpox OR Monkeypox")
        assert "OR" in result
        assert '"Mpox"[Title/Abstract]' in result
        assert '"Monkeypox"[Title/Abstract]' in result
        # Must NOT contain the literal phrase "Mpox OR Monkeypox"
        assert '"Mpox OR Monkeypox"' not in result

    def test_or_separated_produces_parenthesized_boolean(self):
        result = pubmed_scan._name_clause("Mpox OR Monkeypox")
        expected = '("Mpox"[Title/Abstract] OR "Monkeypox"[Title/Abstract])'
        assert result == expected

    def test_multi_word_name_no_or(self):
        result = pubmed_scan._name_clause("Mycobacterium tuberculosis")
        assert result == '("Mycobacterium tuberculosis"[Title/Abstract])'

    def test_vibrio_cholerae(self):
        result = pubmed_scan._name_clause("Vibrio cholerae")
        assert result == '("Vibrio cholerae"[Title/Abstract])'

    def test_three_way_or(self):
        """If we ever have 3 synonyms separated by OR."""
        result = pubmed_scan._name_clause("A OR B OR C")
        assert '"A"[Title/Abstract]' in result
        assert '"B"[Title/Abstract]' in result
        assert '"C"[Title/Abstract]' in result
        assert result.count(" OR ") == 2


# ---------------------------------------------------------------------------
# GENE_SYNONYMS: ensure critical genes have expanded queries
# ---------------------------------------------------------------------------

class TestGeneSynonyms:
    """Verify that genes known to have synonyms are properly mapped."""

    def test_tb_inha_has_synonyms(self):
        syns = pubmed_scan.GENE_SYNONYMS.get("inhA", [])
        assert len(syns) >= 1
        # Must include the enzyme name that PubMed actually indexes
        names_lower = [s.lower() for s in syns]
        assert any("enoyl" in n for n in names_lower), \
            f"inhA synonyms missing enoyl-ACP reductase: {syns}"

    def test_tb_gyrb_has_synonyms(self):
        syns = pubmed_scan.GENE_SYNONYMS.get("gyrB", [])
        assert len(syns) >= 1

    def test_tb_pnca_has_synonyms(self):
        syns = pubmed_scan.GENE_SYNONYMS.get("pncA", [])
        assert len(syns) >= 1
        names_lower = [s.lower() for s in syns]
        assert any("pyrazinam" in n for n in names_lower)

    def test_tb_rrs_has_synonyms(self):
        syns = pubmed_scan.GENE_SYNONYMS.get("rrs", [])
        assert len(syns) >= 1
        names_lower = [s.lower() for s in syns]
        assert any("16s" in n for n in names_lower)

    def test_tb_cfp10_has_synonyms(self):
        syns = pubmed_scan.GENE_SYNONYMS.get("CFP-10", [])
        assert len(syns) >= 1

    def test_cholera_ctxa_has_synonyms(self):
        syns = pubmed_scan.GENE_SYNONYMS.get("ctxA", [])
        assert len(syns) >= 1
        names_lower = [s.lower() for s in syns]
        assert any("cholera toxin" in n for n in names_lower)

    def test_ebola_np_has_synonyms(self):
        syns = pubmed_scan.GENE_SYNONYMS.get("NP", [])
        assert len(syns) >= 1
        assert "nucleoprotein" in syns

    def test_mers_orf1a_has_synonyms(self):
        syns = pubmed_scan.GENE_SYNONYMS.get("ORF1a", [])
        assert len(syns) >= 1


# ---------------------------------------------------------------------------
# PATHOGEN_GENES: structural integrity checks
# ---------------------------------------------------------------------------

class TestPathogenGenes:
    """Verify PATHOGEN_GENES dict is structurally valid and contains no
    known duplicate/synonym traps that inflate gap counts."""

    def test_all_pathogens_have_search_name(self):
        for key, cfg in pubmed_scan.PATHOGEN_GENES.items():
            assert "search_name" in cfg, f"{key} missing search_name"
            assert len(cfg["search_name"]) > 0

    def test_all_pathogens_have_genes(self):
        for key, cfg in pubmed_scan.PATHOGEN_GENES.items():
            assert "genes" in cfg, f"{key} missing genes list"
            assert len(cfg["genes"]) >= 1

    def test_all_pathogens_have_applications(self):
        for key, cfg in pubmed_scan.PATHOGEN_GENES.items():
            assert "applications" in cfg, f"{key} missing applications list"
            assert len(cfg["applications"]) >= 1

    def test_no_duplicate_genes_per_pathogen(self):
        """Each pathogen's gene list must not have exact duplicates."""
        for key, cfg in pubmed_scan.PATHOGEN_GENES.items():
            genes = cfg["genes"]
            assert len(genes) == len(set(genes)), \
                f"{key} has duplicate genes: {genes}"

    def test_no_duplicate_applications_per_pathogen(self):
        for key, cfg in pubmed_scan.PATHOGEN_GENES.items():
            apps = cfg["applications"]
            assert len(apps) == len(set(apps)), \
                f"{key} has duplicate applications: {apps}"

    def test_mers_no_spike_s_protein_duplicate(self):
        """MERS must not list both 'spike' and 'S protein' as separate genes.
        They are the same gene and would create a false double-gap."""
        mers = pubmed_scan.PATHOGEN_GENES.get("mers", {})
        genes = [g.lower() for g in mers.get("genes", [])]
        has_spike = "spike" in genes
        has_s_protein = "s protein" in genes
        assert not (has_spike and has_s_protein), \
            "MERS lists both 'spike' and 'S protein' — same gene, inflates gap count"

    def test_ebola_no_np_nucleoprotein_duplicate(self):
        """Ebola must not list both 'NP' and 'nucleoprotein' as separate genes."""
        ebola = pubmed_scan.PATHOGEN_GENES.get("ebola", {})
        genes = [g.lower() for g in ebola.get("genes", [])]
        has_np = "np" in genes
        has_nucleoprotein = "nucleoprotein" in genes
        assert not (has_np and has_nucleoprotein), \
            "Ebola lists both 'NP' and 'nucleoprotein' — same gene, inflates gap count"

    def test_rsv_no_n_protein_nucleoprotein_duplicate(self):
        """RSV must not list both 'N protein' and 'nucleoprotein' as separate genes."""
        rsv = pubmed_scan.PATHOGEN_GENES.get("rsv", {})
        genes = [g.lower() for g in rsv.get("genes", [])]
        has_n = "n protein" in genes
        has_nucleoprotein = "nucleoprotein" in genes
        assert not (has_n and has_nucleoprotein), \
            "RSV lists both 'N protein' and 'nucleoprotein' — same gene, inflates gap count"

    def test_human_grch38_not_in_pathogens(self):
        """Human GRCh38 is not a pathogen and should not generate gap claims."""
        assert "human" not in pubmed_scan.PATHOGEN_GENES, \
            "Human GRCh38 should not be in PATHOGEN_GENES — it inflates gap counts"

    def test_or_names_have_no_embedded_quotes(self):
        """Search names with OR must not contain embedded quotes that would
        break PubMed query syntax."""
        for key, cfg in pubmed_scan.PATHOGEN_GENES.items():
            name = cfg["search_name"]
            assert '"' not in name, \
                f"{key} search_name contains quotes: {name}"


# ---------------------------------------------------------------------------
# Gene query construction: expanded queries with synonyms
# ---------------------------------------------------------------------------

class TestGeneQueryExpansion:
    """Verify that gene_scan builds expanded queries using GENE_SYNONYMS."""

    def test_gene_with_synonyms_produces_or_query(self):
        """inhA query should include 'inhA' plus all its synonyms."""
        all_names = ["inhA"] + pubmed_scan.GENE_SYNONYMS.get("inhA", [])
        gene_clause = " OR ".join(
            f'"{g}"[Title/Abstract]' for g in all_names
        )
        # Must have multiple OR clauses
        assert gene_clause.count(" OR ") >= 1
        assert '"inhA"[Title/Abstract]' in gene_clause
        assert '"enoyl-ACP reductase"[Title/Abstract]' in gene_clause

    def test_gene_without_synonyms_uses_primary_only(self):
        """A gene not in GENE_SYNONYMS should still produce a valid query.
        Use a gene unlikely to appear in any ontology."""
        # prM is a dengue polyprotein cleavage site, not an NCBI gene symbol
        all_names = ["prM"] + pubmed_scan.GENE_SYNONYMS.get("prM", [])
        gene_clause = " OR ".join(
            f'"{g}"[Title/Abstract]' for g in all_names
        )
        assert gene_clause == '"prM"[Title/Abstract]'


# ---------------------------------------------------------------------------
# Regression: the exact bugs that caused false claims
# ---------------------------------------------------------------------------

class TestRegressions:
    """Regressions for specific bugs that produced false medical claims."""

    def test_mpox_not_literal_phrase(self):
        """BUG 2026-03-04: 'Mpox OR Monkeypox' was wrapped in quotes and
        treated as a literal phrase, returning 0 instead of 69 papers.

        The _name_clause function must split on ' OR ' and produce
        proper PubMed boolean, NOT a quoted phrase search."""
        clause = pubmed_scan._name_clause("Mpox OR Monkeypox")
        # The exact query that caused the bug:
        bad_query = '("Mpox OR Monkeypox"[Title/Abstract])'
        assert clause != bad_query, \
            "REGRESSION: Mpox query still produces literal phrase search"

    def test_inha_query_includes_enzyme_name(self):
        """BUG 2026-03-04: inhA query returned 0 because PubMed papers use
        'enoyl-ACP reductase' not 'inhA'. Synonym expansion is required."""
        syns = pubmed_scan.GENE_SYNONYMS.get("inhA", [])
        all_names = ["inhA"] + syns
        assert any("enoyl" in n.lower() for n in all_names), \
            "REGRESSION: inhA missing enoyl-ACP reductase synonym"


# ---------------------------------------------------------------------------
# Ontology-enhanced synonym generation
# ---------------------------------------------------------------------------

class TestOntologySynonyms:
    """Verify that _load_ontology_synonyms correctly maps NCBI gene
    annotations to our scanner gene names."""

    def test_ontology_synonyms_loaded(self):
        """Ontology synonym loader should produce non-empty results."""
        ont = pubmed_scan._ontology_synonyms
        assert len(ont) > 0, "No ontology synonyms were loaded"

    def test_sars_cov2_spike_gets_surface_glycoprotein(self):
        """SARS-CoV-2 S protein → 'surface glycoprotein' from NCBI."""
        syns = pubmed_scan.GENE_SYNONYMS.get("S protein", [])
        assert "surface glycoprotein" in syns, \
            f"S protein missing 'surface glycoprotein', got: {syns}"

    def test_ebola_vp35_gets_polymerase_complex(self):
        """Ebola VP35 → 'polymerase complex protein' from NCBI."""
        syns = pubmed_scan.GENE_SYNONYMS.get("VP35", [])
        assert "polymerase complex protein" in syns, \
            f"VP35 missing 'polymerase complex protein', got: {syns}"

    def test_rsv_f_protein_gets_fusion_glycoprotein(self):
        """RSV F protein → 'fusion glycoprotein' from NCBI gene name."""
        syns = pubmed_scan.GENE_SYNONYMS.get("F protein", [])
        assert "fusion glycoprotein" in syns, \
            f"F protein missing 'fusion glycoprotein', got: {syns}"

    def test_rsv_g_protein_gets_attachment_glycoprotein(self):
        """RSV G protein → 'attachment glycoprotein' from NCBI gene name."""
        syns = pubmed_scan.GENE_SYNONYMS.get("G protein", [])
        assert "attachment glycoprotein" in syns, \
            f"G protein missing 'attachment glycoprotein', got: {syns}"

    def test_tb_katg_gets_catalase_peroxidase(self):
        """TB katG → 'catalase-peroxidase' from NCBI gene name."""
        syns = pubmed_scan.GENE_SYNONYMS.get("katG", [])
        assert "catalase-peroxidase" in syns, \
            f"katG missing 'catalase-peroxidase', got: {syns}"

    def test_tb_rpob_gets_rna_polymerase(self):
        """TB rpoB → official NCBI enzyme name."""
        syns = pubmed_scan.GENE_SYNONYMS.get("rpoB", [])
        assert any("polymerase" in s.lower() for s in syns), \
            f"rpoB missing RNA polymerase synonym, got: {syns}"

    def test_manual_synonyms_preserved(self):
        """Manual synonyms must not be overwritten by ontology data."""
        # CFP-10 has manual synonyms (no ontology match expected)
        syns = pubmed_scan.GENE_SYNONYMS.get("CFP-10", [])
        assert "Rv3874" in syns, "Manual CFP-10 synonym 'Rv3874' was lost"
        assert "culture filtrate protein 10" in syns

    def test_no_self_synonym(self):
        """A gene should not appear as its own synonym."""
        for gene, syns in pubmed_scan.GENE_SYNONYMS.items():
            assert gene not in syns, \
                f"Gene '{gene}' listed as its own synonym"

    def test_no_duplicate_synonyms(self):
        """No synonym list should contain duplicates."""
        for gene, syns in pubmed_scan.GENE_SYNONYMS.items():
            lower_syns = [s.lower() for s in syns]
            assert len(lower_syns) == len(set(lower_syns)), \
                f"Gene '{gene}' has duplicate synonyms: {syns}"

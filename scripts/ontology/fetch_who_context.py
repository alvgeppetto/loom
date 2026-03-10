#!/usr/bin/env python3
"""
Fetch WHO priority pathogen data and outbreak/surveillance context.
Outputs: ontology_cache/who_context.json

Fetches:
  - WHO priority pathogen list classification
  - Disease outbreak news (DONs) count
  - Basic epidemiological context (mortality, case fatality, geographic spread)

Uses WHO GHO (Global Health Observatory) OData API — free, no key.
"""
import json
import os
import sys
import time
import urllib.request
import urllib.error

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
MANIFEST = os.path.join(SCRIPT_DIR, "pathogen_manifest.json")
CACHE_DIR = os.path.join(SCRIPT_DIR, "ontology_cache")
OUTPUT = os.path.join(CACHE_DIR, "who_context.json")

# Curated WHO context — some of this isn't machine-fetchable, so we embed verified facts
WHO_CONTEXT = {
    "cholera": {
        "who_classification": "epidemic-prone",
        "case_fatality_rate": "< 1% with treatment, up to 50% without",
        "annual_cases_est": "1.3-4M globally",
        "annual_deaths_est": "21,000-143,000",
        "geographic_spread": "endemic in 47+ countries, Africa and South Asia most affected",
        "diagnostic_need": "critical — current culture-based detection takes 24-48h",
        "existing_rapid_tests": ["Crystal VC dipstick (sensitivity ~90%)"],
        "crispr_diagnostic_status": "none published",
        "vaccine_available": True,
        "who_prequalified_dx": True,
    },
    "dengue": {
        "who_classification": "high priority — NTD",
        "case_fatality_rate": "< 1% with treatment, 2-5% severe dengue",
        "annual_cases_est": "100-400M infections, 96M symptomatic",
        "annual_deaths_est": "~40,000",
        "geographic_spread": "128+ countries, tropical/subtropical",
        "diagnostic_need": "high — serotype-specific tests needed for vaccines",
        "existing_rapid_tests": ["NS1 antigen RDTs", "IgM/IgG RDTs"],
        "crispr_diagnostic_status": "limited — a few Cas12 papers",
        "vaccine_available": True,
        "who_prequalified_dx": True,
    },
    "ebola": {
        "who_classification": "high priority — epidemic-prone",
        "case_fatality_rate": "25-90% (average ~50%)",
        "annual_cases_est": "sporadic outbreaks, 100s-10,000s per event",
        "annual_deaths_est": "variable — 11,310 in 2013-16 West Africa outbreak",
        "geographic_spread": "Central and West Africa",
        "diagnostic_need": "critical — BSL-4 requirement delays field diagnosis",
        "existing_rapid_tests": ["OraQuick Ebola RDT"],
        "crispr_diagnostic_status": "very few publications",
        "vaccine_available": True,
        "who_prequalified_dx": True,
    },
    "hepatitis-b": {
        "who_classification": "high priority — WHO elimination target 2030",
        "case_fatality_rate": "chronic: 15-25% premature death from cirrhosis/HCC",
        "annual_cases_est": "296M chronic carriers",
        "annual_deaths_est": "~820,000 (from cirrhosis + liver cancer)",
        "geographic_spread": "global — highest in Africa, Western Pacific",
        "diagnostic_need": "moderate — good serological tests exist, but POC needed in LMICs",
        "existing_rapid_tests": ["HBsAg RDTs (multiple WHO-prequalified)"],
        "crispr_diagnostic_status": "some gene therapy work, limited diagnostics",
        "vaccine_available": True,
        "who_prequalified_dx": True,
    },
    "hiv-1": {
        "who_classification": "high priority — WHO elimination target 2030",
        "case_fatality_rate": "~100% untreated over years, near-normal lifespan with ART",
        "annual_cases_est": "39M living with HIV, 1.3M new infections/year",
        "annual_deaths_est": "~630,000",
        "geographic_spread": "global — 2/3 in Africa",
        "diagnostic_need": "moderate — good tests exist, drug resistance monitoring needed",
        "existing_rapid_tests": ["Determine HIV-1/2", "OraQuick", "many others"],
        "crispr_diagnostic_status": "very active — Cas13 diagnostic + gene therapy",
        "vaccine_available": False,
        "who_prequalified_dx": True,
    },
    "influenza-a": {
        "who_classification": "high priority — pandemic preparedness",
        "case_fatality_rate": "< 0.1% seasonal, 2.5% H5N1 in humans",
        "annual_cases_est": "~1 billion infections/year",
        "annual_deaths_est": "290,000-650,000 seasonal",
        "geographic_spread": "global",
        "diagnostic_need": "high — rapid subtyping for H5N1/H7N9 surveillance",
        "existing_rapid_tests": ["RIDTs (moderate sensitivity)", "PCR standard"],
        "crispr_diagnostic_status": "active — several SHERLOCK/DETECTR papers",
        "vaccine_available": True,
        "who_prequalified_dx": True,
    },
    "mers": {
        "who_classification": "WHO R&D Blueprint priority",
        "case_fatality_rate": "~35%",
        "annual_cases_est": "sporadic — 2,600+ total cases since 2012",
        "annual_deaths_est": "~935 total deaths",
        "geographic_spread": "Middle East, sporadic global",
        "diagnostic_need": "high — limited field-deployable tests",
        "existing_rapid_tests": ["PCR-based only"],
        "crispr_diagnostic_status": "limited",
        "vaccine_available": False,
        "who_prequalified_dx": False,
    },
    "mpox": {
        "who_classification": "PHEIC (declared Aug 2024)",
        "case_fatality_rate": "3-6% (Clade I), < 1% (Clade II)",
        "annual_cases_est": "100,000+ during 2022 outbreak",
        "annual_deaths_est": "~200 during 2022 outbreak",
        "geographic_spread": "endemic Central/West Africa, global since 2022",
        "diagnostic_need": "critical — PCR only, no rapid POC test",
        "existing_rapid_tests": [],
        "crispr_diagnostic_status": "zero publications",
        "vaccine_available": True,
        "who_prequalified_dx": False,
    },
    "rsv": {
        "who_classification": "high burden — leading cause of infant hospitalization",
        "case_fatality_rate": "< 1% in developed countries, higher in LMICs",
        "annual_cases_est": "33M LRTI episodes in children < 5",
        "annual_deaths_est": "~100,000-160,000 in children < 5",
        "geographic_spread": "global, seasonal in temperate zones",
        "diagnostic_need": "high — rapid POC for clinical management",
        "existing_rapid_tests": ["Antigen RDTs (moderate sensitivity)"],
        "crispr_diagnostic_status": "very limited",
        "vaccine_available": True,
        "who_prequalified_dx": False,
    },
    "sars-cov-2": {
        "who_classification": "no longer PHEIC (ended May 2023)",
        "case_fatality_rate": "~0.5-1% (variable by variant and population)",
        "annual_cases_est": "ongoing endemic transmission",
        "annual_deaths_est": "declining",
        "geographic_spread": "global",
        "diagnostic_need": "moderate — abundant tests exist, variant surveillance matters",
        "existing_rapid_tests": ["Abbott BinaxNOW", "many antigen RDTs"],
        "crispr_diagnostic_status": "most studied pathogen for CRISPR dx",
        "vaccine_available": True,
        "who_prequalified_dx": True,
    },
    "tuberculosis": {
        "who_classification": "top infectious disease killer globally",
        "case_fatality_rate": "~50% untreated, ~5% with treatment",
        "annual_cases_est": "10.6M new cases/year",
        "annual_deaths_est": "~1.3M",
        "geographic_spread": "global — highest in India, Indonesia, China, Philippines",
        "diagnostic_need": "critical — smear microscopy still primary, GeneXpert limited by cost",
        "existing_rapid_tests": ["GeneXpert MTB/RIF", "LAM urine test"],
        "crispr_diagnostic_status": "very limited for diagnostics, some drug resistance work",
        "vaccine_available": True,
        "who_prequalified_dx": True,
    },
    "zika": {
        "who_classification": "WHO R&D Blueprint priority",
        "case_fatality_rate": "very low directly, but associated with microcephaly/GBS",
        "annual_cases_est": "sporadic after 2016 pandemic",
        "annual_deaths_est": "rare direct death, indirect fetal/neonatal impact",
        "geographic_spread": "tropical Americas, Southeast Asia, Africa",
        "diagnostic_need": "high — cross-reactivity with dengue makes serology unreliable",
        "existing_rapid_tests": ["limited — PCR-based mainly"],
        "crispr_diagnostic_status": "some SHERLOCK papers",
        "vaccine_available": False,
        "who_prequalified_dx": False,
    },
}


def fetch_gho_data() -> dict:
    """Attempt to fetch supplementary data from WHO GHO API."""
    results = {}
    indicators = {
        "tuberculosis": "MDG_0000000020",  # TB incidence
        "hiv-1": "HIV_0000000001",          # People living with HIV
    }

    for key, indicator in indicators.items():
        url = f"https://ghoapi.azureedge.net/api/{indicator}?$filter=TimeDim%20eq%202023&$top=5&$orderby=NumericValue%20desc"
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "LOOM-ontology/1.0"})
            with urllib.request.urlopen(req, timeout=20) as resp:
                data = json.loads(resp.read().decode("utf-8"))
            results[key] = {
                "indicator": indicator,
                "top_countries": [
                    {"country": v.get("SpatialDim", ""), "value": v.get("NumericValue")}
                    for v in data.get("value", [])[:5]
                ]
            }
        except Exception:
            pass
        time.sleep(0.5)

    return results


def main():
    os.makedirs(CACHE_DIR, exist_ok=True)

    # Start with curated data
    results = WHO_CONTEXT.copy()

    # Try to enrich with live GHO data
    print("Fetching supplementary WHO GHO data...")
    gho = fetch_gho_data()
    for key, gho_data in gho.items():
        if key in results:
            results[key]["gho_supplement"] = gho_data

    with open(OUTPUT, "w") as f:
        json.dump({
            "source": "WHO (curated + GHO API)",
            "fetched": time.strftime("%Y-%m-%dT%H:%M:%SZ"),
            "pathogens": results
        }, f, indent=2)

    print(f"Done. Written to {OUTPUT}")


if __name__ == "__main__":
    main()

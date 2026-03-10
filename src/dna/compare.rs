//! DNA benchmark comparison: plain-index vs transformer baseline.
//!
//! Loads two [`DnaBenchmarkReport`] JSON files produced by the `dna-benchmark`
//! command, computes metric deltas, and emits a comparison report in JSON and
//! Markdown formats.
//!
//! The "baseline" run is typically the deterministic plain-index (exact
//! substring search).  The "candidate" run is typically the transformer-scored
//! retrieval.  The comparison flags whether the candidate meets the production
//! KPI thresholds defined in the DNA pivot plan.
//!
//! # Usage
//!
//! ```no_run
//! use loom::dna::compare::{load_report, compare_reports};
//! use std::path::Path;
//!
//! let baseline = load_report(Path::new("plain_index_report.json")).unwrap();
//! let candidate = load_report(Path::new("transformer_report.json")).unwrap();
//! let cmp = compare_reports(&baseline, &candidate);
//! println!("{}", cmp.render_markdown());
//! ```

use std::fs;
use std::path::Path;

use anyhow::Context;
use serde::Serialize;

use crate::dna::benchmark::DnaBenchmarkReport;

// ---------------------------------------------------------------------------
// KPI thresholds (from DNA pivot plan Week-2 gates)
// ---------------------------------------------------------------------------

/// Minimum acceptable hit\@k for the candidate run.
pub const KPI_MIN_HIT_AT_K: f64 = 0.80;
/// Maximum acceptable regression in hit\@k vs baseline (−5 pp).
pub const KPI_MAX_HIT_REGRESSION: f64 = -0.05;
/// Minimum acceptable recall\@k.
pub const KPI_MIN_RECALL_AT_K: f64 = 0.70;
/// Maximum p95 query latency in milliseconds.
pub const KPI_MAX_P95_LATENCY_MS: f64 = 500.0;

// ---------------------------------------------------------------------------
// Report types
// ---------------------------------------------------------------------------

/// Delta for a single metric (candidate − baseline).
#[derive(Debug, Clone, Serialize)]
pub struct MetricDelta {
    /// Baseline value.
    pub baseline: f64,
    /// Candidate value.
    pub candidate: f64,
    /// Absolute difference (candidate − baseline).
    pub delta: f64,
    /// Relative change as a fraction of the baseline (NaN when baseline = 0).
    pub relative: f64,
}

impl MetricDelta {
    fn new(baseline: f64, candidate: f64) -> Self {
        let delta = candidate - baseline;
        let relative = if baseline == 0.0 { f64::NAN } else { delta / baseline };
        Self { baseline, candidate, delta, relative }
    }
}

/// KPI gate result for a single threshold check.
#[derive(Debug, Clone, Serialize)]
pub struct KpiCheck {
    /// Human-readable label.
    pub label: String,
    /// Whether the check passed.
    pub passed: bool,
    /// Observed value.
    pub observed: f64,
    /// Threshold that must be met.
    pub threshold: f64,
    /// Direction: `">="` (must be at least) or `"<="` (must be at most).
    pub direction: String,
}

impl KpiCheck {
    fn at_least(label: &str, observed: f64, threshold: f64) -> Self {
        Self {
            label: label.to_string(),
            passed: observed >= threshold,
            observed,
            threshold,
            direction: ">=".to_string(),
        }
    }

    fn at_most(label: &str, observed: f64, threshold: f64) -> Self {
        Self {
            label: label.to_string(),
            passed: observed <= threshold,
            observed,
            threshold,
            direction: "<=".to_string(),
        }
    }
}

/// Full comparison between a baseline and a candidate benchmark run.
#[derive(Debug, Clone, Serialize)]
pub struct ComparisonReport {
    /// Label for the baseline run (e.g. `"plain-index"`).
    pub baseline_label: String,
    /// Label for the candidate run (e.g. `"transformer"`).
    pub candidate_label: String,
    /// Deltas for each metric.
    pub hit_at_k: MetricDelta,
    pub precision_at_k: MetricDelta,
    pub recall_at_k: MetricDelta,
    pub ndcg_at_k: MetricDelta,
    /// Latency delta (p95 ms).
    pub p95_latency_ms: MetricDelta,
    /// KPI gate results.
    pub kpi_checks: Vec<KpiCheck>,
    /// Overall verdict: `"PASS"` if all KPI checks pass, else `"FAIL"`.
    pub verdict: String,
    /// Narrative recommendation.
    pub recommendation: String,
    /// ISO 8601 timestamp when this report was generated.
    pub generated_at: String,
}

impl ComparisonReport {
    /// Render a human-readable Markdown report.
    pub fn render_markdown(&self) -> String {
        let mut md = String::new();

        md.push_str("# DNA Benchmark Comparison Report\n\n");
        md.push_str(&format!(
            "**Baseline:** `{}`  \n**Candidate:** `{}`  \n**Generated:** {}\n\n",
            self.baseline_label, self.candidate_label, self.generated_at
        ));

        // Metrics table
        md.push_str("## Retrieval Metrics\n\n");
        md.push_str("| Metric | Baseline | Candidate | Delta | Relative |\n");
        md.push_str("|--------|----------|-----------|-------|----------|\n");
        for (name, d) in [
            ("hit@k", &self.hit_at_k),
            ("precision@k", &self.precision_at_k),
            ("recall@k", &self.recall_at_k),
            ("ndcg@k", &self.ndcg_at_k),
        ] {
            md.push_str(&format!(
                "| {} | {:.4} | {:.4} | {:+.4} | {:.1}% |\n",
                name,
                d.baseline,
                d.candidate,
                d.delta,
                if d.relative.is_nan() { 0.0 } else { d.relative * 100.0 }
            ));
        }
        md.push('\n');

        // Latency table
        md.push_str("## Latency (p95 ms)\n\n");
        md.push_str("| | Baseline | Candidate | Delta |\n");
        md.push_str("|--|----------|-----------|-------|\n");
        md.push_str(&format!(
            "| p95 | {:.3} | {:.3} | {:+.3} |\n\n",
            self.p95_latency_ms.baseline,
            self.p95_latency_ms.candidate,
            self.p95_latency_ms.delta
        ));

        // KPI gate
        md.push_str("## KPI Gate\n\n");
        md.push_str("| Check | Threshold | Observed | Result |\n");
        md.push_str("|-------|-----------|----------|--------|\n");
        for kpi in &self.kpi_checks {
            let result = if kpi.passed { "✅ PASS" } else { "❌ FAIL" };
            md.push_str(&format!(
                "| {} | {} {:.4} | {:.4} | {} |\n",
                kpi.label, kpi.direction, kpi.threshold, kpi.observed, result
            ));
        }
        md.push('\n');

        // Verdict
        md.push_str(&format!("## Verdict: **{}**\n\n", self.verdict));
        md.push_str(&format!("{}\n", self.recommendation));

        md
    }
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Load a [`DnaBenchmarkReport`] from a JSON file.
pub fn load_report(path: &Path) -> anyhow::Result<DnaBenchmarkReport> {
    let text = fs::read_to_string(path)
        .with_context(|| format!("Cannot read benchmark report {}", path.display()))?;
    let report: DnaBenchmarkReport = serde_json::from_str(&text)
        .with_context(|| format!("JSON parse error in {}", path.display()))?;
    Ok(report)
}

/// Compare a baseline and candidate [`DnaBenchmarkReport`].
pub fn compare_reports(
    baseline: &DnaBenchmarkReport,
    candidate: &DnaBenchmarkReport,
) -> ComparisonReport {
    compare_reports_labeled(baseline, "plain-index", candidate, "transformer")
}

/// Compare with explicit run labels.
pub fn compare_reports_labeled(
    baseline: &DnaBenchmarkReport,
    baseline_label: &str,
    candidate: &DnaBenchmarkReport,
    candidate_label: &str,
) -> ComparisonReport {
    let hit_delta = MetricDelta::new(baseline.hit_at_k, candidate.hit_at_k);
    let prec_delta = MetricDelta::new(baseline.precision_at_k, candidate.precision_at_k);
    let rec_delta = MetricDelta::new(baseline.recall_at_k, candidate.recall_at_k);
    let ndcg_delta = MetricDelta::new(baseline.ndcg_at_k, candidate.ndcg_at_k);
    let lat_delta = MetricDelta::new(baseline.latency.p95_ms, candidate.latency.p95_ms);

    let kpi_checks = vec![
        KpiCheck::at_least(
            &format!("candidate hit@{} >= {:.2}", candidate.k, KPI_MIN_HIT_AT_K),
            candidate.hit_at_k,
            KPI_MIN_HIT_AT_K,
        ),
        KpiCheck::at_least(
            &format!("hit@{} regression <= {:.0} pp", candidate.k, KPI_MAX_HIT_REGRESSION.abs() * 100.0),
            hit_delta.delta,
            KPI_MAX_HIT_REGRESSION,
        ),
        KpiCheck::at_least(
            &format!("candidate recall@{} >= {:.2}", candidate.k, KPI_MIN_RECALL_AT_K),
            candidate.recall_at_k,
            KPI_MIN_RECALL_AT_K,
        ),
        KpiCheck::at_most(
            &format!("candidate p95 latency <= {:.0} ms", KPI_MAX_P95_LATENCY_MS),
            candidate.latency.p95_ms,
            KPI_MAX_P95_LATENCY_MS,
        ),
    ];

    let all_pass = kpi_checks.iter().all(|c| c.passed);
    let verdict = if all_pass { "PASS".to_string() } else { "FAIL".to_string() };

    let recommendation = build_recommendation(all_pass, &hit_delta, &rec_delta, &lat_delta, candidate_label);

    let generated_at = {
        use std::time::SystemTime;
        let secs = SystemTime::now()
            .duration_since(SystemTime::UNIX_EPOCH)
            .map(|d| d.as_secs())
            .unwrap_or(0);
        crate::dna::benchmark::format_utc_iso(secs)
    };

    ComparisonReport {
        baseline_label: baseline_label.to_string(),
        candidate_label: candidate_label.to_string(),
        hit_at_k: hit_delta,
        precision_at_k: prec_delta,
        recall_at_k: rec_delta,
        ndcg_at_k: ndcg_delta,
        p95_latency_ms: lat_delta,
        kpi_checks,
        verdict,
        recommendation,
        generated_at,
    }
}

fn build_recommendation(
    all_pass: bool,
    hit: &MetricDelta,
    rec: &MetricDelta,
    lat: &MetricDelta,
    candidate_label: &str,
) -> String {
    if all_pass {
        format!(
            "The {} run meets all KPI thresholds. \
             Hit@k delta: {:+.4}, recall@k delta: {:+.4}, p95 latency delta: {:+.3} ms. \
             Recommend proceeding to the next milestone.",
            candidate_label, hit.delta, rec.delta, lat.delta
        )
    } else {
        let issues: Vec<String> = [
            if hit.candidate < KPI_MIN_HIT_AT_K {
                Some(format!("hit@k {:.4} below {:.2} threshold", hit.candidate, KPI_MIN_HIT_AT_K))
            } else {
                None
            },
            if hit.delta < KPI_MAX_HIT_REGRESSION {
                Some(format!(
                    "hit@k regression {:+.4} exceeds allowed {:+.2}",
                    hit.delta, KPI_MAX_HIT_REGRESSION
                ))
            } else {
                None
            },
            if rec.candidate < KPI_MIN_RECALL_AT_K {
                Some(format!("recall@k {:.4} below {:.2} threshold", rec.candidate, KPI_MIN_RECALL_AT_K))
            } else {
                None
            },
            if lat.candidate > KPI_MAX_P95_LATENCY_MS {
                Some(format!(
                    "p95 latency {:.3} ms above {:.0} ms limit",
                    lat.candidate, KPI_MAX_P95_LATENCY_MS
                ))
            } else {
                None
            },
        ]
        .into_iter()
        .flatten()
        .collect();

        format!(
            "The {} run failed {} KPI check(s): {}. Do not proceed until resolved.",
            candidate_label,
            issues.len(),
            issues.join("; ")
        )
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dna::latency::LatencyStats;

    fn make_report(hit: f64, prec: f64, rec: f64, ndcg: f64, p95: f64) -> DnaBenchmarkReport {
        DnaBenchmarkReport {
            corpus_size: 100,
            n_queries: 50,
            k: 10,
            query_len: 30,
            hit_at_k: hit,
            precision_at_k: prec,
            recall_at_k: rec,
            ndcg_at_k: ndcg,
            latency: LatencyStats {
                n: 50,
                cold_ms: p95,
                warm_median_ms: p95 * 0.5,
                p95_ms: p95,
                p99_ms: p95 * 1.1,
                mean_ms: p95 * 0.6,
                min_ms: 0.1,
                max_ms: p95 * 1.2,
            },
            generated_at: "2026-03-01T00:00:00Z".to_string(),
            seed: 42,
        }
    }

    #[test]
    fn compare_identical_reports_zero_delta() {
        let b = make_report(0.9, 0.8, 0.85, 0.75, 50.0);
        let c = make_report(0.9, 0.8, 0.85, 0.75, 50.0);
        let cmp = compare_reports(&b, &c);
        assert!((cmp.hit_at_k.delta).abs() < 1e-9);
        assert!((cmp.precision_at_k.delta).abs() < 1e-9);
        assert!((cmp.recall_at_k.delta).abs() < 1e-9);
        assert!((cmp.ndcg_at_k.delta).abs() < 1e-9);
    }

    #[test]
    fn all_kpi_pass_gives_pass_verdict() {
        let baseline = make_report(0.95, 0.85, 0.90, 0.80, 40.0);
        let candidate = make_report(0.91, 0.82, 0.88, 0.78, 45.0);
        let cmp = compare_reports(&baseline, &candidate);
        assert_eq!(cmp.verdict, "PASS");
        assert!(cmp.kpi_checks.iter().all(|c| c.passed));
    }

    #[test]
    fn low_hit_at_k_gives_fail_verdict() {
        let baseline = make_report(0.95, 0.85, 0.90, 0.80, 40.0);
        let candidate = make_report(0.50, 0.40, 0.45, 0.35, 45.0); // below 0.80 threshold
        let cmp = compare_reports(&baseline, &candidate);
        assert_eq!(cmp.verdict, "FAIL");
        assert!(cmp.kpi_checks.iter().any(|c| !c.passed));
    }

    #[test]
    fn large_regression_gives_fail_verdict() {
        let baseline = make_report(0.95, 0.85, 0.90, 0.80, 40.0);
        let candidate = make_report(0.88, 0.80, 0.85, 0.75, 45.0); // −0.07 regression > −0.05 limit
        let cmp = compare_reports(&baseline, &candidate);
        assert_eq!(cmp.verdict, "FAIL");
    }

    #[test]
    fn high_latency_gives_fail_verdict() {
        let baseline = make_report(0.92, 0.82, 0.88, 0.78, 40.0);
        let candidate = make_report(0.92, 0.82, 0.88, 0.78, 600.0); // p95 > 500 ms
        let cmp = compare_reports(&baseline, &candidate);
        assert_eq!(cmp.verdict, "FAIL");
    }

    #[test]
    fn render_markdown_contains_key_sections() {
        let b = make_report(0.95, 0.85, 0.90, 0.80, 40.0);
        let c = make_report(0.91, 0.82, 0.88, 0.78, 45.0);
        let cmp = compare_reports(&b, &c);
        let md = cmp.render_markdown();
        assert!(md.contains("# DNA Benchmark Comparison Report"));
        assert!(md.contains("Retrieval Metrics"));
        assert!(md.contains("KPI Gate"));
        assert!(md.contains("Verdict"));
    }

    #[test]
    fn metric_delta_relative_zero_baseline() {
        let d = MetricDelta::new(0.0, 0.5);
        assert!(d.relative.is_nan());
        assert!((d.delta - 0.5).abs() < 1e-9);
    }
}

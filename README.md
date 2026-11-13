# decidability-report

A small command-line tool that asks, for each cluster in a single-cell or spatial
transcriptomics dataset, **whether its cell-type label can be uniquely determined**
given the marker evidence you actually have — and writes out the alternatives that
could not be excluded.

Most annotation pipelines report a label and a confidence. This tool reports
something different: an *equivalence class* (the current label plus the alternative
identities that the measured genes cannot rule out), a three-level grade derived from
that class, and a per-alternative verdict with the reason. It is meant to be run on
your own data before the labels are used for downstream analysis.

Companion code for a manuscript on the decidability of cell-type annotation
(under review). The rules are fixed in the source and are not tuned per dataset.

## Install

```bash
pip install -r requirements.txt
```

Python >= 3.9. `scanpy` / `anndata` are only needed to read `.h5ad`; the optional
`requests` is only used with `--use-llm`.

## Quick start

```bash
python decidability_report.py --h5ad data.h5ad --label-col cell_type --outdir report/
```

With an existing clustering and a restricted gene panel:

```bash
python decidability_report.py --h5ad data.h5ad --label-col cell_type \
    --cluster-col leiden --genes panel_genes.txt --outdir report/
```

## Inputs

| flag | required | meaning |
|---|---|---|
| `--h5ad` | yes | AnnData object with log-normalised expression in `X` (or counts; detection rates are computed on non-zero values) |
| `--label-col` | yes | column in `adata.obs` holding the cluster labels to be tested |
| `--cluster-col` | no | an existing clustering column; if absent, Leiden (resolution 1) is computed |
| `--genes` | no | text file, one gene per line; restricts the evidence to a marker panel. Genes absent from the object are treated as "not measured" |
| `--min-cells` | no | minimum cells per cluster (default 30); smaller clusters are reported but marked |
| `--use-llm` | no | enumerate alternatives with an LLM instead of the built-in rule library; needs `DEEPSEEK_API_KEY` |
| `--tissue` | no | tissue hint passed to the LLM enumeration |
| `--outdir` | no | output directory (default `report/`) |

## Outputs

`decidability.csv` — one row per cluster:
`cluster, n_cells, consensus_label, n_alternatives, n_retained, equivalence_class_size, grade`

`alternatives.csv` — one row per (cluster, alternative):
`alternative, must_present, must_absent, axis, verdict, n_present_meas, n_absent_meas, measurable_frac`

`report.md` — human-readable summary with the rule thresholds and the limitations.

Grades: equivalence class of 1 → `decidable`; 2–3 → `partial`; 4 or more →
`undecidable`. A cluster whose discriminating genes are all outside the measured
panel is `unevaluable` and is kept out of the grade counts.

## The rules (fixed, not tuned)

- **must-present gene** in a cluster: detection ≥ 25% and rank ≤ 5 within the cluster
  → supports the alternative; 10–25% → inconclusive; < 10% → contradicts, the
  alternative is excluded.
- **must-absent gene**: detection ≥ 50% in the cluster → contradicts, excluded.
- all discriminating genes outside the measured panel → unmeasurable, stays in the
  equivalence class.

## Two alternative sets, and why the grade depends on them

The tool offers two ways to build the set of alternatives:

1. **Offline rule library** (default, no API key): a curated table of confusable pairs
   for common lineages (12 lineage families, 2–4 alternatives each, 4 per cluster on
   average). Fast, deterministic, reproducible.
2. **`--use-llm`**: the model enumerates alternatives per cluster (typically 14–24 per
   cluster, including state-versus-type and granularity variants). Slower, and the
   enumeration itself varies between runs.

The grade is **relative to the alternative set**, not an absolute property of the
cluster. On the same heart MERFISH dataset used in the manuscript, the offline library
yields 0 decidable / 17 partial / 3 undecidable clusters, while the LLM enumeration
yields 0 / 8 / 12. Report the set you used. This is the tool's main message rather
than a bug: what can be decided depends on what hypotheses you were willing to test.

## Limits

- Evidence is expression only. Protein, spatial context and functional experiments
  would change the equivalence class.
- The rule thresholds are a convention; they are written here and in the source so a
  reader can recompute with other thresholds.
- LLM enumeration is not deterministic; use the offline library when reproducibility
  across machines matters.

## License

MIT (see `LICENSE`). Citation details will be added when the companion manuscript is
accepted.

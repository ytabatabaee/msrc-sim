# Spatial reservoir reference prototype audit

See `docs/msrc_spatial_prototype_validation.md` for settings and interpretation.
The main validation script generates primary raw data and figures; run from
the repository root with the project .venv. The primary seed is 20261121.

- `marginals.csv`, `validation_summary.json`: 18 cases × 1,000 replicates.
- `refined_marginals_10000.csv`, `refined_summary_10000.json`: expanded four-taxon
  p=0.2, k=1 case, including the initial 1,000 draws.
- `selected_case_draws.csv`, `root_time_distributions.png`: four-taxon p=0.5,
  k=1; these repeat that primary case's same draws, not new evidence.
- `spatial_summary.csv`: chromosome-bootstrap correlation and Wilson topology
  intervals for positions 25 and 75, using the primary run.
- `rate_audit.csv`: 12,000 breakpoint/background marks per each of 18 settings.
- `benchmark.csv`, `profile.txt`: diagnostic short-chromosome execution costs.
- `runtime.csv`: total audit elapsed time per case across all three methods,
  including statistics extraction; use benchmark.csv for simulator-only costs.

The primary comparison family is Holm-corrected together. Practical intervals
are marginal 95% intervals; they are not simultaneous equivalence guarantees.

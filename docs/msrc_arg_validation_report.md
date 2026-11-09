# Two-background MSRC-ARG validation report

## Scope

This report covers the first structured-ARG implementation for two genetic
backgrounds. It is deliberately a simplified mechanistic model: physical
recombination is sampled from a joint position/background rate measure,
cross-background events are suppressed by `kappa_cross`, and an anchor marker
determines the background inherited on each side of a crossover. It is not a
complete inversion-heterokaryotype or gene-conversion model.

No MSRC-specific code was added to the neutral ARG's default mode.

## Reproducibility

The primary experiment used 1,000 independent chromosomes for each
`kappa_cross` value (`0, 0.01, 0.1, 0.5, 1.0`), seed `20261109`, four taxa,
sequence length 1,000 bp, constant recombination rate `1e-5` per bp per
generation, fixed terminal backgrounds `(1, 1, 0, 0)`, and a frozen frequency
history `p(t)=0.5` until a structured-origin age of 100 generations.

Run it with:

```bash
./.venv/bin/python validation/msrc_arg/validate_msrc_arg.py \
  --replicates 1000 --seed 20261109 \
  --output validation/msrc_arg/results/msrc_arg_validation_1000.json \
  --figure validation/msrc_arg/results/msrc_suppression_events.png
```

The compact result is in
`validation/msrc_arg/results/msrc_arg_validation_1000.json`; the generated
figure is `validation/msrc_arg/results/msrc_suppression_events.png`.

## Event-rate and inheritance checks

The implementation uses the same piecewise integration for hazard calculation
and breakpoint sampling. Inversion boundaries are included in the integration
partition, so a crossover interval cannot silently use the wrong suppression
factor. Unit tests cover a single active interval, disjoint active intervals,
gaps, stored-but-resolved ancestry, zero-suppression regions, nonempty parent
lineages, and anchor-oriented background inheritance.

For the synthetic map, the direct transition-rate targets are the expected
`p_g * integral(r0 * eta)` values. For example, with `p=0.37`, state 0, and a
100-bp focal-to-anchor distance outside the inversion, the target increases
from `5.55e-4` at κ=0 to `1.48e-3` at κ=1. The full target table is stored in
the JSON result.

## Empirical suppression experiment

| κ | mean physical recombinations | mean inside | mean outside | mean marginal trees | MRCA correlation |
|---:|---:|---:|---:|---:|---:|
| 0 | 2.768 | 1.046 | 1.722 | 2.618 | 0.456 |
| 0.01 | 2.793 | 1.040 | 1.753 | 2.636 | 0.449 |
| 0.1 | 2.805 | 1.129 | 1.676 | 2.689 | 0.360 |
| 0.5 | 3.010 | 1.347 | 1.663 | 2.896 | 0.339 |
| 1 | 3.276 | 1.673 | 1.603 | 3.205 | 0.331 |

The results show linked local genealogies in every setting, and increasing
within-inversion physical recombination as cross-background suppression is
relaxed. The correlation statistic is calculated across independent
chromosome replicates, not by treating windows from one chromosome as
independent observations.

The focal topology at 100 bp remained predominantly the species-tree topology
(0.970, 0.968, 0.954, 0.942, and 0.892 for the same κ sequence). These are
descriptive validation results, not evidence for a biological introgression
mechanism.

## Single-locus comparison

The validation script also runs the existing `simulate_genealogy()` with the
matched effective rate

`rho_eff = integral(r0 * kappa_cross)`

for the focal-to-anchor path. For κ values 0, 0.01, 0.1, 0.5, and 1, the
reference quartet frequencies were respectively:

```text
[0.961, 0.022, 0.017]
[0.961, 0.022, 0.017]
[0.956, 0.024, 0.020]
[0.938, 0.031, 0.031]
[0.910, 0.042, 0.048]
```

The direct ARG marginal frequencies are not expected to be numerically
identical in this first experiment because the ARG includes a finite physical
chromosome, explicit ancestry-changing recombination, species boundaries, and
the marker-only transfer mechanism, whereas the legacy single-locus simulator
uses an effective switching generator. The comparison is therefore a
consistency diagnostic, not an equivalence claim.

## Interpretation and limitations

The neutral limit and ordinary tree-sequence invariants remain covered by the
existing neutral ARG test suite. The new tests verify finite-copy
background-dependent coalescence, deterministic transition boundaries,
anchor-oriented inheritance, and the no-empty-parent invariant.

In the 1,000-replicate frozen-history run, the mean number of logged
`background_transfer` events was zero for every κ. This is not evidence that
the mechanism is absent: the direct gap/marker-transfer regression test
constructs the required ancestry geometry and observes the event. It means
that this particular full-chromosome setup did not generate such an
ancestry-invisible event at appreciable frequency. A dedicated high-rate,
marker-separated experiment is still needed to estimate its frequency and
validate its contribution to single-locus switching.

The current model also does not implement multiple rearrangements, gene
conversion, a full haplotype-frequency process, or inversion-specific viability
selection. Background inheritance is explicitly a documented effective model.
Consequently, readiness for biological MSRC inference is **CONDITIONAL PASS**:
the modular implementation and core event-rate tests are in place, but the
single-locus limit and background-transfer frequency require a larger targeted
validation before scientific use.

## Files

- Model specification: `docs/msrc_arg_model.md`
- Runnable configuration: `examples/arg_msrc_two_backgrounds.yaml`
- Unit tests: `tests/test_msrc_arg.py`
- Validation program: `validation/msrc_arg/validate_msrc_arg.py`
- Machine-readable results: `validation/msrc_arg/results/msrc_arg_validation_1000.json`
- Figure: `validation/msrc_arg/results/msrc_suppression_events.png`

# Neutral ARG validation report

## Scope and corrected bugs

This report records the post-fix validation of the repository's independent
neutral ARG engine. `msprime` is used only as a reference implementation.

The audit found and corrected two concrete bugs:

1. The recombination hazard was based on `active` ancestry, while breakpoint
   sampling used all stored segments. Both now use the same map-weighted
   eligible region: the half-open span from the leftmost to rightmost
   unresolved active interval. Gaps inside that span are eligible; stored
   resolved material outside it is not. Zero-rate map pieces contribute zero
   hazard and cannot be sampled.
2. Recombination edge recording did not split a stored segment at the sampled
   breakpoint, so an edge could incorrectly cover both parental intervals and
   produce unrooted marginal trees. Stored segments are now split before left
   and right parental edges are recorded.

Empty active parents are rejected, inactive lineages are retired after local
completion, and the event log remains separate from the marginal tree-sequence
view. The production tree sequence retains provenance nodes; validation
statistics simplify that tree sequence before counting marginal genealogies.

## Reproducible commands

```bash
./.venv/bin/python -m pytest -q tests/test_arg.py
MPLCONFIGDIR=/tmp/msrcsim-mpl ./.venv/bin/python \
  validation/arg/validate_neutral.py --replicates 1000 \
  --output validation/arg/results/neutral_validation_1000.json \
  --figures validation/arg/results/figures/neutral_1000
MPLCONFIGDIR=/tmp/msrcsim-mpl ./.venv/bin/python \
  validation/arg/benchmark.py --replicates 3 \
  --output validation/arg/results/benchmark.csv
./.venv/bin/python validation/arg/heterogeneous_check.py --replicates 3000 \
  --output validation/arg/results/heterogeneous_3000.json
```

The archived primary run used seed `20261108`, 1,000 independent chromosomes,
`N_e=100`, 100-kb sequence length, and `r=1e-7` per bp per generation.
Software versions were Python 3.14.5, msprime 1.4.5, tskit 1.0.4, and
NumPy 2.5.3. Machine-readable results are in
`validation/arg/results/neutral_validation_1000.json` and
`validation/arg/results/benchmark.csv`; the supplementary heterogeneous-map
check is in `validation/arg/results/heterogeneous_3000.json`.

## Statistical results

For two samples, the theoretical mean and variance are 200 and 40,000
generations. The simulator gave mean 200.932 (SE 6.259), msprime gave
200.867 (SE 6.155), and the independent implementation difference was
0.066 (SE 8.779; 95% CI [-17.140, 17.272]). Both estimates are consistent
with the exponential reference.

At separations 2, 5, 10, 25, and 50 kb, msrcsim versus msprime MRCA-time
correlations were respectively:

| Distance | msrcsim | msprime | Difference (95% CI) |
|---:|---:|---:|---:|
| 2 kb | 0.959 | 0.936 | 0.023 [-0.016, 0.064] |
| 5 kb | 0.913 | 0.850 | 0.063 [-0.006, 0.129] |
| 10 kb | 0.811 | 0.786 | 0.025 [-0.040, 0.101] |
| 25 kb | 0.628 | 0.586 | 0.042 [-0.048, 0.137] |
| 50 kb | 0.421 | 0.419 | 0.001 [-0.091, 0.105] |

Identical marginal-tree probabilities differed by -0.009, -0.008, -0.028,
-0.022, and -0.011 at those distances; all independent-sample 95% intervals
included zero. Correlations and identity probabilities decrease with distance,
showing linkage-induced dependence rather than independent-block behavior.

For the symmetric quartet tree used here, each of the AB and CD ancestral
branches has length `t`; therefore the theoretical concordant probability is
`1 - (2/3) exp(-2t)` and each discordant probability is `(1/3) exp(-2t)`.
At high, moderate, and low ILS, the expected concordant frequencies were
0.397, 0.755, and 0.998. The empirical frequencies and implementation
differences are in the JSON result, with all three-topology differences having
95% intervals covering zero.

For 8 and 16 samples, total branch-length differences were -22.9 (95% CI
[-66.8, 21.0]) and -14.7 ([-60.0, 30.6]) generations. Heterogeneous-map
marginal-tree-count difference was -0.095 ([-0.233, 0.043]). A supplementary
3,000-replicate heterogeneous-map check gave MRCA means 203.48 versus 197.95
and difference 5.53 (95% CI [-4.59, 15.65]), resolving the wider 1,000-run
fluctuation at one fixed position.

## Figures and performance

Figures are in `validation/arg/results/figures/neutral_1000/`:

- `two_locus_correlation.png`
- `quartet_topology_frequencies.png`

The benchmark separates simulation and tree-sequence construction for msrcsim
and records peak Python-traced memory, event counts, marginal-tree counts,
peak lineages, and peak segments. The current implementation is slower than
msprime for larger samples and recombination rates because it enumerates
same-population pairs and performs Python interval bookkeeping. This is an
expected correctness-first limitation, not a reason to replace the engine.

Across the three benchmark configurations and three timing replicates, mean
total runtime was 0.0067/0.0064 s for 4 taxa, 0.0084/0.0047 s for 8 taxa,
and 0.0210/0.0050 s for 16 taxa (msrcsim/msprime). These are small smoke
benchmarks, not production-scale throughput claims. msrcsim peak traced memory
averaged approximately 0.031, 0.013, and 0.023 MB for 4, 8, and 16 taxa;
Python allocation tracing is not a process-wide RSS measurement.

## Readiness assessment

**CONDITIONAL PASS.** The neutral engine now has consistent recombination
hazard/sampling, no empty active parental lineages, rooted marginal trees, and
quantitative neutral distributions consistent with theory and msprime within
the reported uncertainty. The remaining limitations are performance and the
fact that this validation uses a focused grid rather than every possible
species-tree/demography configuration. Proceeding to the two-background MSRC
extension is reasonable only if the existing regression suite remains green
and the heterogeneous-map supplementary result is retained with the archive.
No MSRC background or inversion-specific process was implemented here.

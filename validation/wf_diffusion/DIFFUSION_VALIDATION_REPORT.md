# Wright–Fisher diffusion validation report

Status: **EXPERIMENTAL — FURTHER VALIDATION REQUIRED**

The validation was run in `msrc-sim` at commit
`507aea5c097ea2c2f26ae07173a0985e200d07d8` with seed `20261011`. The
discrete Wright–Fisher and Moran implementations were not modified.

## Executive result

The corrected diffusion has the intended neutral scaling, bounded frequencies,
absorbing endpoints, and correct origin timing. At `Ne=500`, the selected
four-taxon genealogy comparisons were close to discrete WF within the reported
history-level Monte Carlo errors. The 30-taxon benchmark is computationally
practical: the frequency history took about 0.001 seconds and 500 complete
genealogies took 11–13 seconds on the validation machine.

The feature is not promoted to validated status. The beta transition is a
moment-matched approximation rather than an exact Wright–Fisher transition,
and the convergence grid still shows sampling-scale variation in finite-time
absorption and terminal-frequency statistics. The current evidence supports
large-scale exploratory use with convergence checks, not an unattended Mai
study.

## Method and acceptance rules

The process is

```text
dp_u = sqrt(p_u(1-p_u)) dW_u,       u=t/(2Ne).
```

Interior steps use a beta distribution with mean `p` and variance
`p(1-p)du`. Values in the boundary layer are projected to `{0, epsilon}` or
`{1-epsilon, 1}` with a mean-preserving probability. Origins are hard time
boundaries; the history stores absence immediately older than an interior
origin and the initial frequency exactly at the origin.

Before examining results, the screening acceptance rule was set to four Monte
Carlo standard errors for analytical means and heterozygosity, with an absolute
floor of 0.005 for common frequencies. Fixation probabilities used four
binomial standard errors. Genealogy topology differences were assessed against
the combined history-level standard error, using 0.05 as a practical effect
size. These are screening criteria, not claims of publication-level power.

## Analytical validation

The analytical stage used `du=0.002`, `tolerance=1e-5`, coalescent times
`u={0.01,0.1,0.5,1,3}`, and initial frequencies
`{0.001,0.01,0.1,0.5,0.9,0.99,0.999}`. Common initial frequencies used 1,000
trajectories; frequencies 0.01 and 0.99 used 5,000; frequencies 0.001 and
0.999 used 20,000. Eventual absorption trajectories were followed for up to
50 coalescent units.

The largest recorded mean-frequency error was 0.0194 at `p0=0.1,u=0.5`, with
MC SE 0.00536. The largest heterozygosity error was 0.0183 at the same setting,
with MC SE 0.00505. Both are below four MC SE only marginally and motivate
repeat runs for study-specific settings. At `p0=0.5,u=1`, observed
heterozygosity was 0.1781 versus theory 0.1839 (MC SE 0.00639).

Eventual fixation estimates were close to theory. Selected rows were:

| p0 | estimated fixation | expected | MC SE | histories |
|---:|---:|---:|---:|---:|
| 0.001 | 0.00125 | 0.001 | 0.00025 | 20,000 |
| 0.01 | 0.0104 | 0.01 | 0.00143 | 5,000 |
| 0.1 | 0.101 | 0.1 | 0.00953 | 1,000 |
| 0.5 | 0.501 | 0.5 | 0.0158 | 1,000 |
| 0.9 | 0.903 | 0.9 | 0.00936 | 1,000 |
| 0.999 | 0.999 | 0.999 | 0.00022 | 20,000 |

The complete machine-readable results are in
`results/analytical_moments.tsv`. The analytical tests and boundary tests pass
in CI, but the finite-time moment outlier means the diffusion should not be
treated as numerically certified from one seed.

## Numerical convergence

The convergence grid used 1,000 trajectories at `p0=0.5,u=1` over
`max_step_cu={0.05,0.02,0.01,0.005,0.002}` and
`tolerance={1e-2,1e-3,1e-4,1e-5}`. Mean terminal frequencies ranged roughly
from 0.478 to 0.528 across settings, with per-cell SE about 0.0126. Expected
heterozygosity was approximately 0.177–0.194, with SE about 0.0063, around the
theoretical value 0.18394.

No monotone plateau was established for finite-time absorption: estimates were
approximately 0.42–0.48 across the grid. This reflects both Monte Carlo noise
and sensitivity of absorption timing to the boundary-layer approximation.
The results support `max_step_cu=0.01,tolerance=1e-3` as a reasonable starting
setting, but do not prove it is optimal.

## WF and Moran comparisons

The matched process stage used 100 independent rearrangement histories per
process and five conditional gene trees per history on a balanced four-taxon
tree. Gene trees sharing a history were summarized within history before
calculating standard errors.

| Ne | process | mean root frequency | topology 0 | topology 1 | topology 2 |
|---:|---|---:|---:|---:|---:|
| 20 | WF | 0.2185 ± 0.0247 | 0.504 ± 0.023 | 0.250 ± 0.019 | 0.246 ± 0.018 |
| 20 | Moran | 0.2647 ± 0.0243 | 0.476 ± 0.023 | 0.296 ± 0.020 | 0.228 ± 0.018 |
| 20 | diffusion | 0.2707 ± 0.0302 | 0.502 ± 0.018 | 0.236 ± 0.018 | 0.262 ± 0.019 |
| 100 | WF | 0.2861 ± 0.0128 | 0.376 ± 0.023 | 0.320 ± 0.019 | 0.304 ± 0.021 |
| 100 | Moran | 0.2554 ± 0.0124 | 0.346 ± 0.021 | 0.348 ± 0.020 | 0.306 ± 0.020 |
| 100 | diffusion | 0.2539 ± 0.0128 | 0.378 ± 0.022 | 0.330 ± 0.023 | 0.292 ± 0.020 |
| 500 | WF | 0.2444 ± 0.0056 | 0.306 ± 0.020 | 0.362 ± 0.023 | 0.332 ± 0.021 |
| 500 | Moran | 0.2467 ± 0.0057 | 0.340 ± 0.018 | 0.338 ± 0.018 | 0.322 ± 0.020 |
| 500 | diffusion | 0.2441 ± 0.0051 | 0.330 ± 0.022 | 0.330 ± 0.020 | 0.340 ± 0.021 |

At `Ne=500`, diffusion and WF terminal means differed by 0.0003 and quartet
differences were within roughly one combined history-level SE. At `Ne=20`, the
finite-population processes differ more, as expected. These comparisons used a
moderate suppression fraction of 0.1; they are not a complete suppression grid.

## Genealogy and history discretization

The structured genealogy engine consumed the realized piecewise-constant
histories through its existing `frequency_at` and `next_frequency_boundary`
interface. The validation did not introduce a second genealogy simulator or
force terminal patterns. Full four-taxon trees completed, and the existing
30-taxon integrity checks passed.

The current convergence tables show no large systematic topology shift at
`Ne=500`, but the validation harness does not yet provide a high-power
history-discretization experiment with common random numbers and conditional
rare-persistence strata. That remains a limitation, especially near 0 and 1
where short frequency excursions affect switching and coalescence rates.

## Runtime benchmark

The original-scale example used 30 taxa and `Ne=10,000`. Measurements were:

| max_step_cu | history records | history s | 20 loci s | 500 loci s |
|---:|---:|---:|---:|---:|
| 0.05 | 104 | 0.0010 | 0.487 | 11.274 |
| 0.01 | 104 | 0.0008 | 0.414 | 12.285 |
| 0.002 | 106 | 0.0009 | 0.623 | 12.566 |

The genealogy engine is the bottleneck by several orders of magnitude. The
frequency simulation no longer scales with every raw generation. Peak resident
memory was approximately 50–55 MB in these runs. The benchmark is a simulator
runtime measurement, not a Mai D1 reproduction.

## Recommendation

Use `wf_diffusion` only as **experimental** for original-scale Mai-like
simulations. Start with `tolerance=0.001,max_step_cu=0.01`, run the convergence
grid on the exact species tree and suppression settings, and report
history-level uncertainty. Do not use it as the sole basis for biological
conclusions until a larger genealogy-level comparison, including rare and
persistent histories and multiple suppression strengths, is completed.

Reproduce the validation with:

```bash
PYTHONPATH=. python validation/wf_diffusion/scripts/run_validation.py \
  --stage all --replicates 100 --output validation/wf_diffusion
```

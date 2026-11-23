# Four-taxon chromosome biological validation

Conditional positive control: fixed p=0.5 on every species-tree branch,
terminal states ABCD=1010, Ne=1000, root age 200 generations, origin 1200,
marker 5 Mb, suppressed interval [3,7) Mb, r=1e-8 per bp per generation.
Model A is ordinary neutral uniform ARG, B neutral low-map ARG, C MSRC-ARG
with kappa=.01 and D the same structured model with kappa=1.
No ARG event generator or SMC code is changed.

Run from the simulator repository with its installed environment:

```bash
python validation/chromosome_biological_validation/theory_grid.py
python validation/chromosome_biological_validation/run.py
MPLCONFIGDIR=/tmp/msrc-mpl python validation/chromosome_biological_validation/analyze.py
```

The theory adapter loads the original full balanced prescribed-frequency
operator from Desktop/msrc-paper, computing terminal and sister population
transitions as well as structured root ancestry. Fixed terminal states on this
prescribed history are a conditional positive control, not a posterior sample
from stochastic Moran histories. Binomial daughter sampling is not applied to
an externally prescribed constant history.

The pilot chromosomes scale marker and region coordinates with chromosome
length while retaining the specified per-base rate and demography. They test
runtime/validity, not the full 10-Mb switching profile. Each run is isolated,
limited to 120 seconds and 1.5 GiB, and failures are retained. Initial B pilots
had a tuple/list input validation error before any simulated events; corrected
same-seed retries use distinct case names. No biological history is rejected.

Results include per-replicate interval tables and native tree sequences.
Uncertainty uses one observation per chromosome at each selected coordinate;
q intervals are exact Clopper-Pearson. Failure to reject a theory probability
is not an equivalence result. Window bands are observed chromosome ranges.
Branch-length results are simulation predictions in generations, with no
mutation-rate or empirical parameter fitting claim.

```bash
python validation/chromosome_biological_validation/sensitivity.py
MPLCONFIGDIR=/tmp/msrc-mpl python validation/chromosome_biological_validation/details.py
MSRC_THEORY_SOURCE=/Users/ytabatabaee/Desktop/msrc-paper python -m pytest -q -o addopts='' -p no:cacheprovider
```

`point_estimates.csv` retains all four individual 95% interval exclusions;
none is removed or relabeled as a passing equivalence test. Timing reflects
concurrent diagnostic runs, not a controlled performance comparison. Source
and seed provenance are in results/provenance.json. The numerical change in
msrc-study is commit 87362a2; existing unrelated working changes were preserved.
Large-M attempts with root duration 50 and budgets 30/120 seconds are separate
from the centered-marker duration 1000 attempt with a 30-second budget. Each
attempt's input and failed stage are explicit in its JSON.

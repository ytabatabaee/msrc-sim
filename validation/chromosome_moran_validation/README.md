# Stochastic Moran chromosome validation

The worker calls the production continuous-time Moran history simulator and
MSRC-ARG kernel. M=2000, Ne=1000, a centered initial count K=1000, nonzero
terminal/sister branches 100 generations, root age 200, structured root 1000,
origin 1200, r=1e-8, marker 5 Mb, region [3,7)Mb and kappa=.01 are the primary design.
No per-window history or genealogy is resampled. History and ARG RNGs are
SeedSequence([20261123, history_seed, 0/1]).

Run from the repository root using the installed environment:

```bash
python validation/chromosome_moran_validation/history_pool.py
python validation/chromosome_moran_validation/small_M.py
python validation/chromosome_moran_validation/assess.py
python validation/chromosome_moran_validation/conditional_extension.py
MPLCONFIGDIR=/tmp/msrc-mpl python validation/chromosome_moran_validation/analyze.py
MSRC_THEORY_SOURCE=/Users/ytabatabaee/Desktop/msrc-paper python -m pytest -q -o addopts='' -p no:cacheprovider
```

`config.json` records the target 100 chromosomes/condition, initial two-pilot
stage, original fixed 64-prior-draw audit, time/memory limits and window widths.
The initial fixed 64 trials estimate acceptance probability. Their only 1010
is seed 30. The separately recorded sequential rejection extension stops at the
next 1010 or 64 extra trials; it found seed 64 immediately. The extension is not
pooled into a fixed-n binomial probability estimate. Both accepted histories
are selected by sampled terminal states before any genealogy. No genealogy
failure is replaced and no survival filtering is used.

An individual example uses the same job format archived in each summary:

```bash
python validation/chromosome_moran_validation/worker.py '{"length":10000000,"model":"Moran_A","duration":1000,"kappa":0.01,"pattern":"sampled","seed":0,"history_seed":0,"seconds":600}' /tmp/msrc-moran-A0
python validation/chromosome_moran_validation/worker.py '{"length":10000000,"model":"Moran_B","duration":1000,"kappa":0.01,"pattern":"sampled","seed":30,"history_seed":30,"condition":"1010","seconds":600}' /tmp/msrc-moran-B30
```

Use fresh output directories: worker locks prevent overlapping duplicate jobs.
The seed0 180-second attempt, seed1 180-second attempt, seed0 600-second replay,
seed 30/64 conditional pilots and the explicitly separate T200 sensitivity are
retained. Seed0's instrumented 30-second profile is diagnostic, not additional
replication. Some jobs ran concurrently; timing is a resource diagnostic, not
a controlled comparative speed benchmark. Failed jobs have no complete tree
counts; partial generation/carriers are separate fields. Replays of the same
seed are deduplicated in statistical analyses.

`windows.py` intersects every tree span with each window. Native tree times
are generations. Per-window `.csv.gz` files retain each chromosome separately
at 40/200/400/1000 windows; neither trees nor windows are independent replicate
units. Summary/bootstrap curves are empirical and remain unvalidated against
large-M exact theory. For n=1 no interval is drawn; n=2 bootstrap intervals do
not establish adequate precision. The earlier fixed-p and neutral chromosomes
are loaded as unchanged references (20 per condition).

The full large-M theory attempt is separately profiled; its exact target is
not replaced with a constant-frequency prediction. `small_M.json` instead
uses the original full stochastic balanced recursion with nonzero populations,
4000 independent prior draws and all accepted 1010 draws (100), with no losses
or genealogy errors discarded. `small_M_comparisons.csv` retains unresolved
criteria. History traces are representative downsampled queries; they are not
the complete histories. Full histories can be regenerated from exact seeds.

A progressive cohort driver is available (the100-replicate setting was not
executed because the pilots already exceed practical resource budgets):

```bash
python validation/chromosome_moran_validation/reproduce.py --replicates 2 --seconds 600 --output /tmp/msrc-moran-pilots
```

For larger cohorts it continues rejection draws from seed65, with explicit
trial/history-time limits, saves all rejections and the accepted cohort plan,
and never replaces a failed genealogy. An adaptive rejection extension is not
used as a fixed-n binomial probability estimate. The archived original64 draws
and seed64 extension establish the first two conditional draws.

All three completed chromosomes (two primary examples, one distinct T200
sensitivity) passed native/coordinate/Newick export audit: 1448 intervals,
expected ABCD labels, ultrametric Newick and generation metadata. Large complete
event-by-event frequency-history pickles are not archived; exact histories are
regenerated from seeds. Downsampled per-population traces are clearly labeled.
The planned 100-chromosome cohorts and broad sensitivity grid were not executed.
No published file is described as a 100-replicate dataset.

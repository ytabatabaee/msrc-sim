# Marker-aware MSRC-SMC′ validation status

## Readiness: FAIL — no simulator kernel implemented

This status is deliberate. The first minimal-state candidate was assessed
against the requested conditional-transition definition and found
under-specified. No local genealogy correlation, marginal, or runtime result
is claimed. The existing neutral SMC′ validation remains in
[`msrc_smc_validation.md`](msrc_smc_validation.md).

## Analytical transition counterexample

For background frequency `p_h=1/2` and effective size `N`, one compatible
pair coalesces at `1/(2N p_h)=1/N`. If a detached focal lineage has `m_h`
compatible marker-associated guide lineages for `2N` generations, then

```
P(no coalescence | m_h) = exp(-m_h * (1/N) * 2N) = exp(-2*m_h).
```

Thus `m_h=1` gives 0.1353 and `m_h=2` gives 0.0183. A focal tree plus a
single linked/separated marker bit can be identical in both hidden states,
while the probability of reaching each re-coalescence outcome differs by
more than seven-fold. The proposed state cannot determine a unique transition
kernel.

Adding the marker tree and its current branch association identifies
`m_h` at the present transition, but a scan to a third coordinate requires a
conditional distribution for that new coordinate's marker/focal phase and
ancestry. The repository equations specify the one-locus switch rates
`q_(h->g)` but not this multi-position conditional distribution. The
full-segment particle option would require an explicit proposal and weighting
scheme; none is available to validate here.

## Audit outcomes

- The proposed marginal target is the legacy single-locus structured
  coalescent: background-specific pair coalescence and
  `q01=rho_eff*p`, `q10=rho_eff*(1-p)`, conditional on one shared frequency
  trajectory and a required explicit origin age.
- The current MSRC-ARG's marker-mediated event measure yields the intended
  one-focal switching intensity when the marker and active focal ancestry are
  retained under the specified geometry. This rate identity does not establish
  the joint law of neighboring genealogies.
- A binary marker-linkage flag fails the re-coalescence target-count check
  above. A marker tree plus current branch mapping still leaves the next-site
  conditional ancestry unspecified.
- Existing fixation/loss floor and missing-origin defaults require explicit
  harmonization before marginal comparisons. The old experiments are not
  changed.

## Tests and simulations

No new runtime tests or simulations were run because there is no proposed
transition kernel to test. No statistical validation table or figure is
reported. The analytic calculation above is the worked conditional
counterexample; it is not represented as empirical validation.

## Next implementable research step

Define and calibrate a conditional hidden-state law, for example a particle
proposal over marker-linked segment ancestries conditioned on the retained
local tree, background histories, and linkage summary. Specify the proposal,
weights, resampling, and particle error diagnostics first. Then a small
two-taxon/fixed-frequency prototype can measure both one-site marginal error
and two-site dependence against the legacy MSRC and full MSRC-ARG, with the
ARG treated as a comparison model rather than unquestioned truth. Until that
step is complete, Phase B should not start.

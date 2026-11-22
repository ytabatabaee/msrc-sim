# Exact Moran quartet theory: compatibility audit

Audited against `msrc-paper` commit
`ba260f8d4303b5d3bc884fb76a9c5a5032b7e477` and simulator base
`956e20add44859903b99f00ac1d3b648aca67da5`. The five theory source hashes
are recorded with the numerical results. No theory or production simulator
was changed. This audit concerns the hierarchical MSRC model, not an ancestral
dual of the same microscopic Moran replacement events.

| Parameter or convention | Existing exact theorem / single-locus simulator | Reservoir ARG | Structured SMC′ |
|---|---|---|---|
| Population size | Branch M=2Ne chromosome copies | Same Ne; integer Ne interface | Same |
| Moran clock | Neutral q+=q-=K(M-K)/(2M); attempted replacement rate M/2 per generation-equivalent unit | Same supplied production history | Same supplied history |
| Selection | Native genic fitness 1+s, p_sel=p(1+s)/(1+sp); s>-1 | Same generator, not a second selection clock | Same |
| Frequency state | K(t)/M, absorbing endpoints, discrete jumps; forward histories read backward by genealogy | Same recorded history, shared across chromosome | Same |
| Initial law | pi_old at older root boundary, not a stationary or terminal law | Sample K from pi_old once per chromosome before history generation | Same |
| Speciation | Conditionally independent binomial daughter counts given common parental count; default unbiased split | Same history generator; background labels persist backward at node | Same |
| Branch sizes/durations | Full balanced A1/A2/A4 or pectinate A1/A2/A3/A4 recursion; separate branch M and duration | Explicit species populations and divergence times | Same embedding of dated clades |
| Terminal samples | Bernoulli per tip conditional on its branch-end frequency; conditional results use posterior history given Z | Explicit states or Bernoulli through unified API; validation supplies matched states | Explicit states supplied in validation |
| Backward switching | 0→1: rho p, 1→0: rho(1-p) | rho(x)=integral_marker_to_x r0(y) k_cross(y) dy under reservoir exchangeability | rho(x*) fixed throughout scan; physical cut rate still varies with x |
| Pair coalescence | 1/max(K,1), 1/max(M-K,1), no cross-background merger during structure | 1/(M p_h) for occupied positive p_h | Same, and requires occupied p_h>=1/M for initialization compatibility |
| Fixation/loss | Floor defines even globally incompatible states | Occupied absent state errors, never silently floored | Same; invalid sub-one-copy occupied frequency also errors |
| Origin | Erase structure at explicit origin; no unique-founder merger | Same erasure, explicit positive origin required | Same; augmented root stem retained to origin |
| Older phase | Ordinary same-population MSC pair rate 1/M; exchangeable unresolved quartet mass allocated 1/3 | Same | Same |

## Global boundary difference versus reachable-state agreement

The finite-copy floor and strict boundary policy define different generators
on the entire formal state space. Equality cannot be inferred by substituting
mean p(t), treating an arbitrary prescribed path as Moran, or filtering failed
simulations. A boundary bridge is still needed for externally prescribed paths
that strand occupied backgrounds, sub-one-copy frequencies, mutation or other
resurrection mechanisms, or different founder conventions.

There is, however, a narrower compatibility argument for the **existing
absorbing, no-mutation Moran model with binomial inheritance**:

1. If K=0 at a forward branch time, it stays zero on the younger remainder.
   Every daughter receives count zero with probability one. Recursing down
   that subtree gives only background-0 terminal samples. Background-0
   lineages cannot switch to 1 while p=0; all physical parental marks also
   have background 0. Thus a backward lineage cannot enter the missing
   background at that boundary. The K=M case is identical by complementation.
2. On every occupied positive background in a genuine finite-M history, its
   count is an integer at least one. The finite-copy floor is inactive:
   `1/max(K,1)=1/K` and similarly for background 0.
3. Binomial daughter inheritance cannot create a background missing in its
   parent. Backward species transitions consequently preserve this support.
4. The origin itself erases background labels in both constructions; it is
   not an additional structured interval at count zero. The production
   nextafter origin annotation and 1e-10 age tolerance are retained. Freezing
   a history for the conditional operator excludes origin annotation slivers
   within that tolerance, rather than creating a spurious absent epoch.

These statements give a reachable-state generator match for the reservoir
marginal under its exchangeable parental-reservoir assumption. They do not
make the two policies interchangeable for arbitrary inputs. Every attempted
simulation in the empirical audit is retained or explicitly marked failed;
finite-grid success is not a validation of all demographic histories.

The root-origin comparisons initialize K at the older root boundary, use
ordinary MSC older than that boundary, and retain the literal speciation
inheritance events. Even zero-duration branches would retain binomial draws;
they must not be replaced by the collapsed single-population theorem.

## Conditioning distinctions

- The unconditional finite-M theory averages frequency histories **and** tip
  backgrounds, with their actual probabilities and shared parental counts.
- Fixed count/frequency experiments use the prescribed-frequency multispecies
  operator. They have no stochastic daughter-count inheritance; M sets
  coalescence rates and can differ across populations.
- One realized Moran history is evaluated with the prescribed-history operator,
  conditional on its complete path. Its mean frequency is not substituted into
  the unconditional theorem.
- For a stochastic history conditional on terminal Z, the history law is the
  posterior proportional to its prior probability times P(Z|history). Fixing
  tip bits on independently drawn prior histories gives the wrong law.
  The dedicated audit uses a fixed number of prior/terminal draws and selects
  Z before any genealogy simulation. No trajectory is selected for survival
  or computational success, and failed selected runs are never replaced.

## Legacy ARG and SMC′ qualifications

The legacy ARG shares the theorem's floor but loses marker-to-focal exchange
when marker ancestry disappears. Its formal floor alone does not make its
spatial marginal the theorem's target. This workflow uses reservoir_v1 as the
ARG comparator and the independent original single-locus simulator as a control.

SMC′ initializes from the conditional single-locus model and targets that
path density at rho(x*) through its existing Metropolis correction. Its
frozen-reference target is compared separately from rho(x)-specific theory.
A match at x* does not establish whole-chromosome marginal correctness,
spatial joint correctness, or accuracy of its inherited-state parental closure.

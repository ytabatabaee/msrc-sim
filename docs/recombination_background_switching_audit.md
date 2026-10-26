# Recombination and genetic-background switching audit

## Executive summary

The simulator is a mathematically defined two-stage MSRC model. A forward population process generates a trajectory for the frequency of one binary arrangement/background, and a backward structured coalescent generates focal-locus genealogies conditional on that trajectory.

The backward process uses the effective intensity

\[
\rho = r_{\mathrm{base}}\kappa,
\]

where baseline_rate is the configured baseline and effective_cross_arrangement_fraction, \(\kappa\), is the retained cross-background exchange fraction. For a lineage in background 0 the switch rate is \(\rho p(t)\); for background 1 it is \(\rho[1-p(t)]\).

The implementation does not simulate two-locus haplotypes or recombinant offspring. It simulates arrangement-frequency evolution and then applies an effective ancestral switching model. Therefore it is:

- **B** as a conditional effective structured-coalescent simulator;
- **C** for a literal biological recombination interpretation, because the forward process does not define a focal locus, a background-defining locus, linkage phase, or recombinant haplotypes.

No interior-rate implementation bug was found. The main limitation is a modeling and terminology assumption. Boundary behavior for lineages assigned to an absent background is an explicit finite-copy convention, not a derivation from a two-locus forward model.

## Code inventory

| Component | Current implementation | Expected interpretation | Status | Recommendation |
|---|---|---|---|---|
| msrcsim/wright_fisher.py:61-102 | Binomial frequency updates with \(p_{\rm sel}=p(1+s)/(1+sp)\) | Discrete-generation drift and genic selection of arrangement frequency | Verified | Describe as frequency evolution |
| msrcsim/moran.py:30-50,157-205 | CTMC on copy count \(k\), \(M=2N_e\) | Continuous-time Moran frequency process | Verified | State generation-equivalent time scale |
| msrcsim/wf_diffusion.py:72-217 | Bounded beta transitions and boundary-layer snapping | Neutral diffusion approximation | Verified as approximation | Label as experimental approximation |
| msrcsim/structured_coalescent.py:22-112,176-347 | Backward switching and same-background coalescence conditional on frequency history | Effective two-background structured coalescent | Verified for stated model | Call the parameter an effective switching rate |
| msrcsim/experiments.py:110-116,123-196 | Multiplies baseline rate by suppression fraction and samples terminal labels | Coupling of forward history and focal genealogies | Verified | Document as a two-stage construction |
| msrcsim/speciation_split.py:203-275 | Binomial daughter initialization, with optional frequency distortion | Arrangement-frequency inheritance at speciation | Verified | Treat split offsets as extra assumptions |
| msrcsim/population_process.py:22-88 | Dispatches WF, Moran, and diffusion | Population-process selection | Verified | Preserve model-specific metadata |

## Forward-time model

The population process tracks one binary state, A1, interpreted as the focal arrangement/background. Frequency records contain A1/A0 copy counts, total copies, and status. They do not contain phased haplotypes, a second locus, linkage phase, crossover locations, or recombinant chromosomes.

Wright--Fisher uses

\[
K_{t+1}\sim\operatorname{Binomial}(M,p_{\rm sel}),\qquad
p_{\rm sel}=\frac{p_t(1+s)}{1+s p_t},
\]

with \(M=2N_e\). Moran uses

\[
q_+(k)=\frac{M}{2}p_{\rm sel}(1-p),\qquad
q_-(k)=\frac{M}{2}(1-p_{\rm sel})p,
\]

with \(p=k/M\), and absorbing count states \(k=0,M\). Under neutrality,
\(\operatorname{Var}(dp)/dt=p(1-p)/M=p(1-p)/(2N_e)\).

The diffusion mode is neutral only. It uses bounded beta transitions with branch-local coalescent-time increments and snaps a numerical boundary layer to 0 or 1.

Offspring are not formed by choosing parental chromosomes and recombining segments. Forward events are frequency drift, selection, origin, and daughter initialization. The optional split policy changes daughter target frequencies and samples daughter counts; it is not a crossover model.

## Backward switching equations

experiments.py sets

\[
\rho=r_{\rm base}\kappa.
\]

structured_coalescent.py implements, with state order \((0,1)\),

\[
\lambda_{0\rightarrow1}(t)=\rho p(t),\qquad
\lambda_{1\rightarrow0}(t)=\rho[1-p(t)],
\]

and therefore

\[
Q(t)=
\begin{pmatrix}
-\rho p(t)&\rho p(t)\\
\rho[1-p(t)]&-\rho[1-p(t)]
\end{pmatrix}.
\]

These are the conventional random-exchange rates for a focal locus whose background-defining state 1 has frequency \(p(t)\). At \(p=1/2\) the rates are equal. At fixed p the switching component has stationary distribution \((1-p,p)\).

A switch changes only the background label of a backward focal-locus lineage. It does not create a third state, alter the forward arrangement count, or change an individual's physical arrangement. The parameter is not migration between species or populations.

## Coalescence and time scaling

For a structured branch with effective population size \(N\),

\[
c_0(t)=\frac{1}{2N\max(1-p(t),1/(2N))},\qquad
c_1(t)=\frac{1}{2N\max(p(t),1/(2N))}.
\]

For interior frequencies this is

\[
c_0=\frac{1}{M(1-p)},\qquad c_1=\frac{1}{Mp},\qquad M=2N.
\]

Different-background lineages cannot coalesce in the structured phase. They must switch into a common background or cross the rearrangement-origin boundary. At and beyond that boundary, background labels are discarded and the ordinary MSC rate \(1/(2N)\) is used for every same-population pair.

At species boundaries, lineages move to the parent population and retain their background label. The forward history is continued through daughter branches by the specified initialization kernel.

The time conventions are: WF generations; Moran generation-equivalent continuous time with attempted replacement rate \(M/2\); diffusion updates in branch-local coalescent increments \(dt/(2N_e)\) while recorded ages remain on the configured time scale.

## Boundary and origin behavior

At \(p=0\), the code gives

\[
\lambda_{0\to1}=0,\qquad \lambda_{1\to0}=\rho.
\]

At \(p=1\),

\[
\lambda_{0\to1}=\rho,\qquad \lambda_{1\to0}=0.
\]

Thus a lineage in an absent background can switch one way into the present background, and its coalescence rate is limited by the finite-copy floor. A generated terminal sample has probability zero of being in the absent state, so this is normally reachable only through an explicitly supplied or conditioned state. This is a boundary convention, not a uniquely derived two-locus process.

When \(\rho=0\), no switches occur. When \(\kappa=0\), the same no-switching limit is obtained. When \(\kappa=1\), the configured baseline rate is used. The spatial simulator can use distinct effective fractions inside and outside a rearrangement.

The structured process stops at the recorded rearrangement origin. Older time is ordinary unstructured MSC, so arrangement labels do not persist before the rearrangement existed.

## Forward-backward consistency

The construction is coherent as an effective conditional model:

1. Generate a frequency history.
2. Sample terminal arrangement labels from terminal frequencies.
3. Generate focal genealogies conditional on the history and labels.
4. Allow background association to switch along each backward lineage.

It is compatible with a conventional structured coalescent if a distinct background-defining locus exists, recombination can be summarized by a memoryless effective rate, the frequency trajectory is exogenous to the focal genealogy, and inverse background copy count is the coalescence approximation.

The forward code does not generate recombinant haplotypes, so the backward rates cannot be uniquely derived from it alone. Many two-locus mechanisms could produce the same frequency history and effective switch rate. Overall classification is **C** for a literal biological recombination model, while the conditional mathematical process itself is coherent as **B**. It is not **D** based on the audited implementation.

## Implications for quartet theory

The quartet calculations depend on frequency/count trajectories, \(\rho\), inverse-frequency coalescence, the origin boundary, and species-boundary inheritance. The simulator's interior rates match the finite-\(M\) structured-coalescent conventions used by the theory when \(M=2N_e\) and the same frequency definition is used.

Exact identity still requires matching the boundary floor, count-history versus fixed-frequency assumptions, daughter inheritance kernels, time units, and origin-boundary behavior. The quartet implementation is in the separate msrc-study repository; this audit did not inspect it, so a final direct generator/operator comparison remains necessary.

## Edge-case validation

The focused existing suite was run with the SciPy environment:

    20 passed

It covered Moran neutral/selected rates, absorbing count boundaries, frequency-history boundaries, reproducible Moran histories, optimized/reference backward genealogy invariants, analytic switching identities, structured-origin behavior, and ordinary MSC behavior after the origin.

Static inspection additionally verified no switches at \(\rho=0\), equal switch rates at \(p=1/2\), frequency-boundary updates, WF/Moran loss and fixation, and diffusion boundary treatment. The suite does not prove that the absent-background finite-copy floor is biologically correct; that remains an explicit modeling convention.

## Terminology recommendations

- **Genetic background:** binary ancestral association state carried by the focal locus; not the entire chromosome ancestry.
- **Arrangement frequency:** population frequency of A1.
- **Background switching:** backward-time change in focal-locus association.
- **Effective background-switching rate:** \(\rho=r_{\rm base}\kappa\).
- **Recombination suppression:** reduction in effective exchange between focal and background-defining loci.
- Avoid **arrangement switching** for backward events; it implies physical rearrangement change.
- Reserve **migration** for movement between populations/species-tree branches.
- **Haploid model** is accurate for sampled gene copies and chromosome-copy count \(M=2N_e\), but does not imply clonal whole-genome reproduction.
- Use **recombination rate** only with biological units or explicitly qualify baseline_rate as effective.

A haploid two-background model is adequate as a first-order effective model for inversion-associated ancestry when the focal locus is neutral and exchange is memoryless. It omits phase, crossover locations, gene conversion, multiple breakpoints, and multilocus selection.

## Recommendations

1. Qualify baseline_rate and \(\kappa\) as an effective switching model in Methods.
2. Replace backward “arrangement switching” with “background switching”.
3. Compare the quartet-theory generator and boundary operators directly with structured_coalescent.py.
4. Decide whether externally supplied samples in absent backgrounds should be rejected or retained under the current one-way-switch convention.
5. Add an explicit two-locus forward model only if the manuscript needs claims about recombinant haplotype formation.

No production behavior change is required by this audit.

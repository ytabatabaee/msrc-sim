# Grouped structured-coalescent implementation

`msrcsim.structured_coalescent.simulate_genealogy_reference` preserves the
pre-optimization lineage-level event enumeration. `simulate_genealogy` uses
the same stochastic model with grouped hazards:

- switching hazards are aggregated by population and background;
- coalescence hazards use `choose(n, 2)` within each same-background group;
- ordinary ancestry groups all backgrounds within a population;
- species boundaries and the explicit structured-origin boundary are handled
  before event generation;
- frequency history is queried once per occupied population per event step.

The event category is selected from the exact aggregate hazards, then a
lineage or unordered pair is selected uniformly in that category. The random
number stream therefore changes, but the transition kernel does not.

The benchmark is reproducible with:

```bash
python scripts/benchmark_structured_coalescent.py \
  --config path/to/mai30_config.yaml --replicates 10 --rho 0.02 \
  --output benchmark.json
```

On Mai replicate 007, the measured mean was 0.276 seconds per genealogy for
the reference and 0.128 seconds for the grouped implementation (10
replicates, 30 taxa, rho=0.02). Focused tests compare topology distributions,
complete coalescence, nonnegative times, and origin-boundary behavior.

The implementation does not replace Wright--Fisher histories, approximate
waiting times, or alter the ordinary ancestral fallback.

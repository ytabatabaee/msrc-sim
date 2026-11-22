"""Primary versioned MSRC-ARG API; all versions use the shared ARG scheduler."""
import math
import numpy as np
from .engine import ARGConfig, simulate_arg
from .backgrounds import create_background_model, BACKGROUND_MODEL_VERSIONS


def _seed_metadata(value):
    """JSON-safe description without narrowing numpy.default_rng's seed API."""
    if isinstance(value, np.random.SeedSequence):
        return {"kind": "SeedSequence", **_seed_metadata(value.state)}
    if isinstance(value, np.ndarray): return _seed_metadata(value.tolist())
    if isinstance(value, np.generic): return value.item()
    if isinstance(value, dict): return {k: _seed_metadata(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)): return [_seed_metadata(v) for v in value]
    if value is None or isinstance(value, (int, float, str, bool)): return value
    return {"kind": type(value).__name__}


def simulate_msrc_arg(tree, rec_map, history, marker, suppressed_region,
                      kappa=1.0, background_model="reservoir_v1", seed=1,
                      sampled_backgrounds=None, *, record_events=True, rng=None):
    """Simulate one chromosome conditional on one shared frequency history.

    Versions are 'reservoir_v1' (new default) and 'legacy'. Terminal states
    are supplied or sampled once by Bernoulli in tree.taxa order. Providing
    rng continues that generator; seed is then unused. Historical low-level
    model-object threshold sampling remains unchanged in simulate_arg.
    Neither version silently retries histories that strand an absent state.
    """
    if background_model not in BACKGROUND_MODEL_VERSIONS:
        raise ValueError(f"unknown MSRC background version: {background_model!r}; choose {BACKGROUND_MODEL_VERSIONS}")
    # Keep the original reservoir wrapper guards version-specific. Legacy
    # retains the map/marker acceptance of its historical low-level adapter.
    if background_model == "reservoir_v1":
        if not 0 <= marker <= rec_map.length:
            raise ValueError("marker must lie inside or at an endpoint of the chromosome")
        a, b = suppressed_region
        if not 0 <= a < b <= rec_map.length:
            raise ValueError("suppressed region must lie on the chromosome")
        if any(branch.effective_population_size <= 0 for branch in tree.branches.values()):
            raise ValueError("effective population sizes must be positive")
        points = (0.0,) + tuple(rec_map.breakpoints) + (rec_map.length,)
        if (not all(math.isfinite(x) for x in points)
                or any(b <= a for a, b in zip(points, points[1:]))
                or not all(math.isfinite(r) and r >= 0 for r in rec_map.rates)):
            raise ValueError("map requires finite ordered coordinates and nonnegative finite rates")
    supplied_rng = rng is not None
    rng = np.random.default_rng(seed) if rng is None else rng
    if sampled_backgrounds is None:
        sampled_backgrounds = {}
        for taxon in tree.taxa:
            p = float(history.terminal_frequency(taxon))
            if background_model == "reservoir_v1" and (not np.isfinite(p) or not 0 <= p <= 1):
                raise ValueError("invalid terminal background frequency")
            sampled_backgrounds[taxon] = int(rng.random() < p)
    if set(sampled_backgrounds) != set(tree.taxa):
        raise ValueError("one background assignment is required for every taxon")
    if any(h not in (0, 1) for h in sampled_backgrounds.values()):
        raise ValueError("backgrounds must be 0 or 1")
    model = create_background_model(background_model, history, marker,
                                    suppressed_region, kappa, sampled_backgrounds)
    result = simulate_arg(tree, ARGConfig(rec_map.length, rec_map, model, record_events,
                                          history, sampled_backgrounds), rng)
    source = ("caller_generator" if supplied_rng else "seed_generator" if isinstance(seed, (np.random.Generator, np.random.BitGenerator))
              else "system_entropy" if seed is None else "seed")
    result.metadata.update(seed=None if supplied_rng else _seed_metadata(seed),
                           rng_source=source,
                           marker=marker, suppressed_region=list(suppressed_region),
                           cross_background_suppression=kappa)
    return result

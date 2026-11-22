"""Compatibility imports for the consolidated MSRC-ARG implementation.

New callers should import simulate_msrc_arg from msrcsim.arg. These names
preserve the original reservoir model, defaults, sampling order and RNG stream.
"""
from ..arg.backgrounds.reservoir import ReservoirBackground
from ..arg.msrc import simulate_msrc_arg


def simulate_reservoir_arg(tree, rec_map, history, marker, suppressed_region,
                           kappa=1.0, seed=1, sampled_backgrounds=None):
    return simulate_msrc_arg(tree, rec_map, history, marker, suppressed_region,
                             kappa=kappa, background_model="reservoir_v1",
                             seed=seed, sampled_backgrounds=sampled_backgrounds)

from .base import BackgroundModel
from .neutral import NeutralBackground
from .msrc_two_state import MSRCTwoBackground
from .reservoir import ReservoirBackground

BACKGROUND_MODEL_VERSIONS = ("legacy", "reservoir_v1")


def create_background_model(version, history, marker, suppressed_region,
                            kappa=1.0, sampled_backgrounds=None):
    """Construct a versioned adapter without drawing random numbers."""
    models = {"legacy": MSRCTwoBackground, "reservoir_v1": ReservoirBackground}
    if version not in models:
        raise ValueError(f"unknown MSRC background version: {version!r}")
    return models[version](history, marker, suppressed_region, kappa, sampled_backgrounds)


__all__ = ["BackgroundModel", "NeutralBackground", "MSRCTwoBackground",
           "ReservoirBackground", "create_background_model", "BACKGROUND_MODEL_VERSIONS"]

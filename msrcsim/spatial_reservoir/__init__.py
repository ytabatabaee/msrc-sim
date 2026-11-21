"""Experimental trajectory-conditioned reservoir ARG; not an SMC kernel."""
from .model import ReservoirBackground, simulate_reservoir_arg

__all__ = ["ReservoirBackground", "simulate_reservoir_arg"]

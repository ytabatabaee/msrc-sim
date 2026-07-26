import numpy as np
from msrcsim.model_fitting import (
    compare_models, fit_all_msc, fit_all_networks, msc_probabilities,
    network_probabilities, off_arm_statistics,
)


def test_exact_msc_fit():
    counts = (800, 100, 100)
    best, _ = fit_all_msc(counts)
    assert best.topology == 0
    assert np.allclose((best.q1, best.q2, best.q3), (0.8, 0.1, 0.1))
    assert best.distance < 1e-12


def test_off_arm_signal():
    stats = off_arm_statistics((169, 758, 73))
    assert stats["nearest_msc_topology"] == 1
    assert abs(abs(stats["off_arm_difference"]) - 0.096) < 1e-12
    assert stats["distance_to_nearest_msc_arm"] > 0.06
    assert stats["off_arm_p_value"] < 1e-6


def test_network_represents_interior_point():
    best, _ = fit_all_networks((169, 758, 73))
    assert best.representation_error < 1e-5
    assert best.nondegenerate


def test_msc_vector_stays_on_arm():
    result = compare_models((100, 800, 100))
    assert result["nearest_msc_topology"] == 1
    assert result["distance_to_nearest_msc_arm"] < 1e-12

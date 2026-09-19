from __future__ import annotations

import numpy as np

from msrcsim.prevalence_grid import (
    Cell,
    classify_absorption,
    corrected_moran_cells,
    deterministic_seed,
    discordant_2to2,
    exact_absorption,
    exact_moran_initial_drift,
    exact_wf_initial_drift,
    heterozygosity,
    native_generations_from_tau,
    process_specific_scaled_selection,
    raw_s_from_gamma,
    raw_s_from_gamma_target,
    scaled_selection,
    simulate_cell,
    strict_2to2,
    wf_moran_comparison_rows,
)


def test_absorption_components_sum_to_one_and_schema_for_processes():
    fields = None
    for process in ("wright_fisher", "moran"):
        row = simulate_cell(Cell(process, 8, 0.5, 0.1, 0.0, 50, "test"))
        assert np.isclose(row["P_loss"] + row["P_persist"] + row["P_fix"], 1.0)
        assert row["n_loss"] + row["n_persist"] + row["n_fix"] == 50
        if fields is None:
            fields = set(row)
        else:
            assert set(row) == fields


def test_persistence_and_boundaries():
    assert classify_absorption(0, 20) == "loss"
    assert classify_absorption(20, 20) == "fix"
    assert classify_absorption(7, 20) == "persist"


def test_heterozygosity():
    assert heterozygosity(0.0) == 0.0
    assert heterozygosity(1.0) == 0.0
    assert np.isclose(heterozygosity(0.5), 0.5)


def test_strict_and_discordant_2to2_classification():
    assert strict_2to2("0011")
    assert strict_2to2("0101")
    assert not strict_2to2("0111")
    assert not discordant_2to2("0011")
    assert discordant_2to2("0101")
    assert discordant_2to2("1001")


def test_conditional_unconditional_bookkeeping():
    row = simulate_cell(Cell("wright_fisher", 10, 0.5, 0.2, 0.0, 200, "bookkeeping"))
    assert row["P_discordant_2to2"] <= row["P_strict_2to2"] + 1e-12
    if row["n_persist"]:
        assert 0.0 <= row["P_strict_2to2_given_persist"] <= 1.0
        assert row["P_discordant_2to2_given_persist"] <= row["P_strict_2to2_given_persist"] + 1e-12


def test_deterministic_seed_and_cell_reproducibility():
    a = deterministic_seed(1, "wright_fisher", 20, 0.5, 0.6, 0)
    b = deterministic_seed(1, "wright_fisher", 20, 0.5, 0.6, 0)
    c = deterministic_seed(1, "moran", 20, 0.5, 0.6, 0)
    assert a == b
    assert a != c
    cell = Cell("wright_fisher", 10, 0.5, 0.2, 0.0, 100, "seed")
    assert simulate_cell(cell) == simulate_cell(cell)


def test_raw_scaled_conversion():
    assert raw_s_from_gamma(5, 100) == 0.025
    assert scaled_selection(0.025, 100) == 5.0
    assert native_generations_from_tau(0.6, 20) == 24


def test_process_specific_gamma_target_conversions():
    assert raw_s_from_gamma_target(5, 100, "wright_fisher") == 0.025
    assert raw_s_from_gamma_target(5, 100, "moran") == 0.05
    assert raw_s_from_gamma_target(0, 100, "wright_fisher") == 0.0
    assert raw_s_from_gamma_target(0, 100, "moran") == 0.0
    assert raw_s_from_gamma_target(-2, 20, "wright_fisher") < 0.0
    assert raw_s_from_gamma_target(-2, 20, "moran") < 0.0
    assert process_specific_scaled_selection(0.025, 100, "wright_fisher") == 5.0
    assert process_specific_scaled_selection(0.05, 100, "moran") == 5.0


def test_matched_initial_weak_selection_drift():
    Ne = 1000
    gamma = 2.0
    p0 = 0.37
    s_wf = raw_s_from_gamma_target(gamma, Ne, "wright_fisher")
    s_moran = raw_s_from_gamma_target(gamma, Ne, "moran")
    assert np.isclose(s_moran, 2.0 * s_wf)
    assert abs(exact_wf_initial_drift(p0, s_wf) - exact_moran_initial_drift(p0, s_moran)) < 2e-7


def test_corrected_output_schema_and_old_provenance_comparison():
    wf = simulate_cell(Cell("wright_fisher", 8, 0.5, 0.1, 2.0, 50, "primary_scaled"))
    moran = simulate_cell(Cell("moran", 8, 0.5, 0.1, 2.0, 50, "moran_corrected_matched"))
    for row in (wf, moran):
        assert "gamma_target" in row
        assert "raw_s" in row
        assert "process_specific_scaled_selection" in row
    assert wf["process_specific_scaled_selection"] == moran["process_specific_scaled_selection"] == 2.0
    old_moran = dict(moran)
    old_moran["grid"] = "moran_matched"
    old_moran["P_persist"] = 0.0
    comp = wf_moran_comparison_rows([wf], [moran], [old_moran])
    assert comp
    assert any(r["metric"] == "P_persist" and r["old_absolute_difference"] >= 0.0 for r in comp)


def test_corrected_moran_grid_includes_all_target_gammas():
    gammas = {cell.gamma for cell in corrected_moran_cells()}
    assert gammas == {-5.0, -2.0, 0.0, 2.0, 5.0}
    assert len(corrected_moran_cells()) == 2 * 2 * 3 * 5


def test_exact_vs_monte_carlo_small_ne():
    cell = Cell("wright_fisher", 6, 0.5, 0.2, 0.0, 1200, "exact")
    row = simulate_cell(cell)
    exact = exact_absorption(cell.process, cell.Ne, cell.p0, cell.tau, cell.gamma)
    assert abs(row["P_persist"] - exact["P_persist_exact"]) < 0.07
    moran = Cell("moran", 6, 0.5, 0.2, 0.0, 600, "exact")
    mrow = simulate_cell(moran)
    mexact = exact_absorption(moran.process, moran.Ne, moran.p0, moran.tau, moran.gamma)
    assert abs(mrow["P_persist"] - mexact["P_persist_exact"]) < 0.09


def test_exact_vs_monte_carlo_selected_moran_cell():
    cell = Cell("moran", 6, 0.5, 0.2, 2.0, 1600, "selected_exact")
    row = simulate_cell(cell)
    exact = exact_absorption(cell.process, cell.Ne, cell.p0, cell.tau, cell.gamma)
    assert np.isclose(row["P_loss"] + row["P_persist"] + row["P_fix"], 1.0)
    assert abs(row["P_loss"] - exact["P_loss_exact"]) < 0.08
    assert abs(row["P_persist"] - exact["P_persist_exact"]) < 0.08
    assert abs(row["P_fix"] - exact["P_fix_exact"]) < 0.08


def test_focal_parameter_cell_runs():
    row = simulate_cell(Cell("wright_fisher", 20, 0.5, 0.6, 0.0, 200, "focal_test"))
    assert row["Ne"] == 20
    assert row["p0"] == 0.5
    assert row["branch_duration_scaled"] == 0.6
    assert row["scaled_selection"] == 0.0


def test_root_and_legacy_persistence_labels_are_distinct():
    row = simulate_cell(Cell("wright_fisher", 20, 0.5, 0.6, 0.0, 50, "label_test"))
    assert "P_persist" in row
    assert "P_legacy_daughter_persist" not in row

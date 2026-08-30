import csv
import subprocess
import sys

import pytest


FIELDNAMES = [
    "rearrangement_id",
    "branch_id",
    "parent_branch_id",
    "forward_generation",
    "absolute_age",
    "copy_count_A1",
    "copy_count_A0",
    "population_chromosomes",
    "frequency_A1",
    "frequency_A0",
    "status",
    "is_origin",
    "is_branch_start",
    "is_branch_end",
    "selection_coefficient",
    "effective_population_size",
]


def _write_run_dir(path):
    path.mkdir()
    rows = [
        ["inv_1", "ROOT", "", 0, 20, 4, 16, 20, 0.2, 0.8, "newly_originated", True, True, False, 0.0, 10],
        ["inv_1", "ROOT", "", 1, 10, 10, 10, 20, 0.5, 0.5, "segregating", False, False, True, 0.0, 10],
        ["inv_1", "1", "ROOT", 0, 10, 10, 10, 20, 0.5, 0.5, "segregating", False, True, False, 0.0, 10],
        ["inv_1", "1", "ROOT", 1, 0, 20, 0, 20, 1.0, 0.0, "fixed", False, False, True, 0.0, 10],
        ["inv_1", "2", "ROOT", 0, 10, 10, 10, 20, 0.5, 0.5, "segregating", False, True, False, 0.0, 10],
        ["inv_1", "2", "ROOT", 1, 0, 0, 20, 20, 0.0, 1.0, "lost", False, False, True, 0.0, 10],
    ]
    with (path / "frequency_history.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(FIELDNAMES)
        writer.writerows(rows)
    with (path / "sampled_arrangements.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["taxon", "arrangement", "terminal_frequency_A1"])
        writer.writerow(["1", "1", "1.0"])
        writer.writerow(["2", "0", "0.0"])
    (path / "config.resolved.yaml").write_text("mode: mechanistic\n")


def _require_matplotlib():
    try:
        import matplotlib

        matplotlib.use("Agg")
    except Exception as exc:
        pytest.skip(f"Matplotlib is not importable in this environment: {exc}")


def test_cli_run_dir_autodiscovery(tmp_path):
    _require_matplotlib()
    run_dir = tmp_path / "run"
    _write_run_dir(run_dir)
    result = subprocess.run(
        [sys.executable, "-m", "msrcsim.plot_history_cli", "--run-dir", str(run_dir)],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "Wrote" in result.stdout
    assert (run_dir / "wright_fisher_history.png").exists()
    assert (run_dir / "wright_fisher_history.figure.json").exists()


def test_cli_format_smoke_tests(tmp_path):
    _require_matplotlib()
    run_dir = tmp_path / "run"
    _write_run_dir(run_dir)
    for fmt in ("png", "pdf", "svg"):
        out = tmp_path / f"history.{fmt}"
        subprocess.run(
            [
                sys.executable,
                "-m",
                "msrcsim.plot_history_cli",
                "--frequency-history",
                str(run_dir / "frequency_history.csv"),
                "--sampled-arrangements",
                str(run_dir / "sampled_arrangements.csv"),
                "--output",
                str(out),
                "--format",
                fmt,
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        assert out.exists() and out.stat().st_size > 0
        assert (tmp_path / "history.figure.json").exists()


def test_mechanistic_cli_can_make_history_plot(tmp_path):
    _require_matplotlib()
    out = tmp_path / "mechanistic"
    config = tmp_path / "config.yaml"
    config.write_text(
        """
mode: mechanistic
seed: 9
num_loci: 2
species_tree:
  newick: "((1:8,2:8)A:4,(3:8,4:8)B:4)ROOT;"
  root_extension: 12
  default_effective_population_size: 20
rearrangement:
  id: inv_1
  type: inversion
  origin_branch: ROOT
  origin_time_from_branch_start: 3
  initial_copy_count: 20
  selection:
    coefficient: 0.0
recombination:
  baseline_rate: 0.01
  effective_cross_arrangement_fraction: 0.05
output:
  directory: {out}
  record_resolved_config: true
  record_frequency_history: true
  record_sampled_arrangements: true
  record_gene_trees: false
  record_coalescence_times: true
  record_backward_events: false
  make_history_plot: true
  history_plot:
    filename: auto_history.png
    glyphs_per_row: 8
    max_rows_per_branch: 10
""".format(out=out)
    )
    subprocess.run([sys.executable, "-m", "msrcsim.cli", "--config", str(config)], check=True, capture_output=True, text=True)
    assert (out / "frequency_history.csv").exists()
    assert (out / "auto_history.png").exists()
    assert (out / "auto_history.figure.json").exists()

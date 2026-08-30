from __future__ import annotations

import csv
import json


def write_csv(path, rows, fieldnames=None):
    if fieldnames is None:
        fieldnames = list(rows[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def msrc_dir(tmp_path, *, windows=True):
    out = tmp_path / "msrc"
    out.mkdir(parents=True)
    loci = [
        {"locus_id": i, "position": i * 100, "inside_rearrangement": 2 <= i <= 5, "topology_index": 1 if 2 <= i <= 5 else 0, "topology": "13|24" if 2 <= i <= 5 else "12|34"}
        for i in range(10)
    ]
    write_csv(out / "spatial_loci.csv", loci)
    if windows:
        write_csv(out / "spatial_windows.csv", [
            {"window_id": 0, "start_position": 0, "end_position": 300, "center_position": 150, "n_loci": 4, "q1": 0.5, "q2": 0.5, "q3": 0.0, "dominant_topology": 0, "distance_to_nearest_msc_arm": 0.1, "off_arm_difference": 0.5},
            {"window_id": 1, "start_position": 400, "end_position": 700, "center_position": 550, "n_loci": 4, "q1": 0.5, "q2": 0.5, "q3": 0.0, "dominant_topology": 0, "distance_to_nearest_msc_arm": 0.1, "off_arm_difference": 0.5},
        ])
    (out / "spatial_summary.json").write_text(json.dumps({
        "genome_length": 1000,
        "overall_q1": 0.6,
        "overall_q2": 0.4,
        "overall_q3": 0.0,
        "rearrangement_interval": {"start": 200, "end": 500, "interval_id": "inv1"},
        "topology_names": ["12|34", "13|24", "14|23"],
    }))
    return out


def hyb_dir(tmp_path, *, windows=True):
    out = tmp_path / "hyb"
    out.mkdir(parents=True)
    loci = [
        {"locus_id": i, "position_bp": i * 100, "ancestry_state": 1 if i in {1, 2, 6, 7} else 0, "ancestry_label": "introgressed" if i in {1, 2, 6, 7} else "major", "tract_id": i // 2, "topology_index": 1 if i in {1, 2, 6, 7} else 0, "topology_label": "13|24" if i in {1, 2, 6, 7} else "12|34", "model": "pulse_hybridization"}
        for i in range(12)
    ]
    write_csv(out / "hybridization_loci.csv", loci)
    write_csv(out / "hybridization_tracts.csv", [
        {"tract_id": 0, "start_bp": 0, "end_bp": 100, "length_bp": 100, "ancestry_state": 0, "ancestry_label": "major"},
        {"tract_id": 1, "start_bp": 100, "end_bp": 300, "length_bp": 200, "ancestry_state": 1, "ancestry_label": "introgressed"},
        {"tract_id": 2, "start_bp": 300, "end_bp": 600, "length_bp": 300, "ancestry_state": 0, "ancestry_label": "major"},
        {"tract_id": 3, "start_bp": 600, "end_bp": 800, "length_bp": 200, "ancestry_state": 1, "ancestry_label": "introgressed"},
        {"tract_id": 4, "start_bp": 800, "end_bp": 1200, "length_bp": 400, "ancestry_state": 0, "ancestry_label": "major"},
    ])
    if windows:
        write_csv(out / "hybridization_windows.csv", [
            {"window_id": 0, "start_bp": 0, "end_bp": 300, "center_bp": 150, "num_loci": 4, "introgressed_fraction": 0.5, "q1": 0.5, "q2": 0.5, "q3": 0.0, "dominant_topology": 0, "distance_to_nearest_msc_arm": 0.1, "off_arm_difference": 0.5},
            {"window_id": 1, "start_bp": 400, "end_bp": 900, "center_bp": 650, "num_loci": 6, "introgressed_fraction": 0.3333333333, "q1": 0.6666666667, "q2": 0.3333333333, "q3": 0.0, "dominant_topology": 0, "distance_to_nearest_msc_arm": 0.05, "off_arm_difference": 0.3333333333},
        ])
    (out / "hybridization_summary.json").write_text(json.dumps({
        "chromosome_length_bp": 1200,
        "observed_marginal_q": [2 / 3, 1 / 3, 0.0],
        "expected_marginal_q": [0.65, 0.35, 0.0],
    }))
    return out

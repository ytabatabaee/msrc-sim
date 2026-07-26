from msrcsim.experiments import simulate_replicates


def _config(seed=7):
    return {
        "seed": seed,
        "species_tree": {
            "newick": "((1:5,2:5)A:2,(3:5,4:5)B:2)ROOT;",
            "root_extension": 10,
            "default_effective_population_size": 20,
        },
        "rearrangement": {
            "id": "inv", "type": "inversion", "origin_branch": "ROOT",
            "origin_time_from_branch_start": 0, "initial_copy_count": 20,
            "selection": {"coefficient": 0.0},
        },
        "recombination": {
            "baseline_rate": 0.02,
            "effective_cross_arrangement_fraction": 0.1,
        },
        "conditioning": {"mode": "none"},
        "experiment": {"replicates": 3, "loci_per_replicate": 20, "asymmetry_threshold": 0.1},
    }


def test_replicates_are_reproducible():
    a, sa = simulate_replicates(_config())
    b, sb = simulate_replicates(_config())
    assert a == b
    assert sa == sb
    assert sa["accepted_replicates"] == 3

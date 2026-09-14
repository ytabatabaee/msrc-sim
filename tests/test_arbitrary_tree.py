from __future__ import annotations

import json

import numpy as np

from msrcsim.arbitrary import run_arbitrary_simulation
from msrcsim.quartet import classify_quartet_newick, summarize_quartet, validate_newick_taxa
from msrcsim.rearrangement import Rearrangement
from msrcsim.species_tree import SpeciesTree
from msrcsim.structured_coalescent import simulate_genealogy
from msrcsim.wright_fisher import simulate_frequency_history


def test_species_tree_accepts_more_than_four_taxa():
    tree = SpeciesTree("((A:5,B:5)AB:10,((C:5,D:5)CD:5,(E:5,F:5)EF:5)Y:5)ROOT;", 50, 20)
    assert len(tree.taxa) == 6
    assert tree.parent_branch("AB") == "ROOT"
    assert tree.parent_branch("F") == "EF"


def test_arbitrary_genealogy_newick_contains_each_taxon_once():
    tree = SpeciesTree("((A:5,B:5)AB:10,((C:5,D:5)CD:5,(E:5,F:5)EF:5)Y:5)ROOT;", 50, 20)
    rearr = Rearrangement("inv", "inversion", "ROOT", 5, 20, 0.0)
    rng = np.random.default_rng(7)
    history = simulate_frequency_history(tree, rearr, rng)
    sampled = {taxon: int(rng.random() < history.terminal_frequency(taxon)) for taxon in tree.taxa}
    result = simulate_genealogy(0, tree, history, sampled, 0.01, 0.05, rng)
    validate_newick_taxa(result.newick, tree.taxa)
    assert result.topology_index == -1
    assert len(result.coalescences) == len(tree.taxa) - 1


def test_quartet_extraction_matches_legacy_ordering():
    trees = [
        "((A:1,B:1):1,(C:1,D:1):1);",
        "((A:1,C:1):1,(B:1,D:1):1);",
        "((A:1,D:1):1,(B:1,C:1):1);",
    ]
    assert [classify_quartet_newick(t, ["A", "B", "C", "D"]) for t in trees] == [0, 1, 2]
    summary = summarize_quartet(trees, ["A", "B", "C", "D"])
    assert (summary["n1"], summary["n2"], summary["n3"]) == (1, 1, 1)


def test_arbitrary_runner_is_deterministic(tmp_path):
    config = {
        "mode": "arbitrary_tree",
        "seed": 123,
        "num_loci": 4,
        "species_tree": {
            "newick": "((A:5,B:5)AB:10,((C:5,D:5)CD:5,(E:5,F:5)EF:5)Y:5)ROOT;",
            "root_extension": 20,
            "default_effective_population_size": 40,
        },
        "rearrangement": {
            "id": "inv",
            "type": "inversion",
            "origin_branch": "ROOT",
            "origin_time_from_branch_start": 5,
            "initial_copy_count": 20,
        },
        "recombination": {"baseline_rate": 0.01, "effective_cross_arrangement_fraction": 0.1},
        "quartets": [{"taxa": ["A", "B", "C", "D"]}],
        "output": {"directory": str(tmp_path / "run1")},
    }
    out1 = run_arbitrary_simulation(config)
    config["output"] = {"directory": str(tmp_path / "run2")}
    out2 = run_arbitrary_simulation(config)
    assert (out1 / "local_gene_trees.nwk").read_text() == (out2 / "local_gene_trees.nwk").read_text()
    summary = json.loads((out1 / "simulation_summary.json").read_text())
    assert summary["n_taxa"] == 6
    assert summary["gene_tree_validity"]["all_expected_taxa_once"]

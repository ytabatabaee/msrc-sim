import numpy as np
from msrcsim.conditioning import evaluate_conditioning
from msrcsim.species_tree import SpeciesTree
from msrcsim.rearrangement import Rearrangement
from msrcsim.wright_fisher import simulate_frequency_history


def test_terminal_pattern_conditioning_accepts_requested_pattern():
    tree = SpeciesTree("((1:5,2:5)A:2,(3:5,4:5)B:2)ROOT;", 20, 10)
    rr = Rearrangement("inv", "inversion", "ROOT", 0, 20, 0.0)
    history = simulate_frequency_history(tree, rr, np.random.default_rng(1))
    sampled = {"1":1, "2":0, "3":1, "4":0}
    result = evaluate_conditioning(
        {"mode":"terminal_pattern", "accepted_patterns":["1010"]},
        history, sampled, tree.taxa,
    )
    assert result.accepted

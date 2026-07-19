import numpy as np
from msrcsim.species_tree import BalancedQuartetTree
from msrcsim.rearrangement import Rearrangement
from msrcsim.wright_fisher import simulate_frequency_history,sample_terminal_arrangements
from msrcsim.structured_coalescent import simulate_genealogy

def test_mechanistic_reproducibility():
    tree=BalancedQuartetTree.create(10,20,40,50)
    r=Rearrangement(initial_copies=20,recombination_rate=.1,suppression_factor=.1)
    def run(seed):
        rng=np.random.default_rng(seed); h=simulate_frequency_history(tree,r,rng); z=sample_terminal_arrangements(h,tree,rng)
        gs=[simulate_genealogy(tree,h,z,r.effective_recombination_rate,rng) for _ in range(20)]
        return z,[g.newick for g in gs],[g.topology_index for g in gs]
    assert run(7)==run(7)

def test_no_cross_population_coalescence_before_root_possible():
    tree=BalancedQuartetTree.create(5,10,20,10)
    r=Rearrangement(initial_copies=0,recombination_rate=0)
    rng=np.random.default_rng(3); h=simulate_frequency_history(tree,r,rng)
    z={str(i):0 for i in range(1,5)}
    g=simulate_genealogy(tree,h,z,0,rng)
    # A discordant first cherry can only occur after root age.
    if g.topology_index in (1,2): assert g.coalescence_times[0] >= tree.root_age

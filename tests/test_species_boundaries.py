from msrcsim.species_tree import BalancedQuartetTree

def test_tree_boundaries():
    t=BalancedQuartetTree.create(10,20,40,100)
    assert t.parent_of("t1")=="left"
    assert t.parent_of("left")=="root"
    assert t.boundary_age("t1")==10
    assert t.boundary_age("left")==20

from msrcsim.species_tree import SpeciesTree

def test_balanced_and_unbalanced():
    b=SpeciesTree('((1:100,2:100)A:50,(3:100,4:100)B:50)ROOT;',1000,500)
    assert b.root.age==150
    u=SpeciesTree('(((1:100,2:100)A:50,3:150)B:50,4:200)ROOT;',1000,500)
    assert u.root.age==200 and u.parent_branch('1')=='A' and u.parent_branch('A')=='B'

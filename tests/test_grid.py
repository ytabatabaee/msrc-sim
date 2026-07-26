from msrcsim.experiments import parameter_combinations


def test_parameter_combinations_cartesian_product():
    rows = list(parameter_combinations({"a":[1,2], "b":[3,4,5]}))
    assert len(rows) == 6
    assert {tuple(sorted(x.items())) for x in rows} == {
        (("a",1),("b",3)), (("a",1),("b",4)), (("a",1),("b",5)),
        (("a",2),("b",3)), (("a",2),("b",4)), (("a",2),("b",5)),
    }

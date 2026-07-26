from msrcsim.wright_fisher import FrequencyHistory, FrequencyRecord


def rec(age, p):
    k = int(round(200*p))
    return FrequencyRecord('x','B',None,int(age),float(age),k,200-k,200,p,1-p,'segregating',False,False,False,0.0,100)


def test_piecewise_frequency_and_next_boundary():
    h = FrequencyHistory([], {'B':[rec(2,0.8),rec(1,0.4),rec(0,0.1)]})
    assert h.frequency_at('B',0.5) == 0.1
    assert h.frequency_at('B',1.5) == 0.4
    assert h.frequency_at('B',3.0) == 0.8
    assert h.next_frequency_boundary('B',0.5) == 1.0
    assert h.next_frequency_boundary('B',1.0) == 2.0

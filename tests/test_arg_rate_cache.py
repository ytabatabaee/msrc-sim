import pytest
from msrcsim.arg.rates import RateCache
from msrcsim.arg.backgrounds.reservoir import ReservoirBackground
from msrcsim.arg.backgrounds.neutral import NeutralBackground
from msrcsim.arg.recombination_map import PiecewiseConstantRecombinationMap as Map
from msrcsim.arg.lineages import ARGLineage
from msrcsim.arg.segments import Segment


class History:
    structured_origin_age=10.
    def frequency_at(self,pop,t):return .2 if t<5 else .8
    def next_frequency_boundary(self,pop,t):return 5. if t<5 else 10. if t<10 else float('inf')
    def terminal_frequency(self,taxon):return .2


def test_memoized_rates_are_bit_identical_across_geometry_time_and_eviction():
    bg=ReservoirBackground(History(),0,(25,75),.01);rec=Map(100,(25,75),(.01,.0001,.005));cache=RateCache(bg,rec,capacity=2)
    for t in (0,1,6,8,11,12,2,6):
        cache.epoch()
        for active in (((0,100),),((20,30),(70,80)),((20,80),),((80,90),)):
            for h in (0,1):
                a=ARGLineage(0,'ROOT',(Segment(0,100,0),),0,h,active)
                b=ARGLineage(1,'ROOT',(Segment(0,100,1),),0,1-h,active)
                direct=bg.recombination_hazard('ROOT',h,a.active,a.segments,t,rec)
                assert cache.recombination(a,t)==direct
                assert cache.recombination(a,t)==direct
                assert cache.coalescence(a,b,20,t)==bg.coalescence_rate('ROOT',h,1-h,20,t)
                assert cache.coalescence(a,a,20,t)==bg.coalescence_rate('ROOT',h,h,20,t)
                assert cache.boundary('ROOT',t)==bg.next_boundary('ROOT',t)
                assert len(cache.hazards)<=2


def test_absent_background_errors_are_not_hidden_by_cached_geometry():
    class Loss(History):
        def frequency_at(self,pop,t):return .5 if t<5 else 0.
    bg=ReservoirBackground(Loss(),0,(25,75),.01);cache=RateCache(bg,Map.constant(100,.001))
    a=ARGLineage(0,'ROOT',(Segment(0,100,0),),0,1)
    assert cache.recombination(a,0)>0
    cache.epoch()
    with pytest.raises(ValueError,match='absent'):cache.recombination(a,6)
    with pytest.raises(ValueError,match='absent'):cache.coalescence(a,a,20,6)


def test_custom_adapters_and_maps_do_not_use_builtin_memoization():
    class CustomNeutral(NeutralBackground):pass
    class CustomMap(Map):pass
    assert not RateCache.supported(CustomNeutral(),Map.constant(100,.001))
    assert not RateCache.supported(NeutralBackground(),CustomMap.constant(100,.001))

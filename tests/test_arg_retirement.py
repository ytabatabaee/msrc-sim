"""Coverage sweep versus the original elementary-interval retirement rule."""
import numpy as np
from msrcsim.arg.coalescence import retire_completed_material
from msrcsim.arg.lineages import ARGLineage
from msrcsim.arg.segments import Segment


def reference(lineages):
    boundaries=sorted({x for l in lineages for interval in l.active for x in interval})
    out=[]
    for lineage in lineages:
        keep=[]
        for left,right in zip(boundaries,boundaries[1:]):
            if any(a<=left and right<=b for a,b in lineage.active):
                if sum(any(a<=left and right<=b for a,b in other.active) for other in lineages)>1:
                    keep.append((left,right))
        replacement=lineage.with_active(keep)
        if replacement.active:out.append(replacement)
    return out


def test_sweep_matches_original_on_fragmented_material():
    rng=np.random.default_rng(20261122)
    for _ in range(200):
        points=sorted(set([0.,100.,*rng.uniform(0,100,15)]));lineages=[]
        for i in range(int(rng.integers(2,15))):
            active=[(a,b) for a,b in zip(points,points[1:]) if rng.random()<.3]
            if active:lineages.append(ARGLineage(i,'ROOT',(Segment(0,100,i),),0,0,active))
        assert retire_completed_material(lineages)==reference(lineages)


def test_sweep_preserves_half_open_and_near_adjacent_boundaries():
    lineages=[ARGLineage(i,'ROOT',(Segment(0,10,i),),0,0,active) for i,active in enumerate([
        ((0,1),(2,3)),((1,2),(3,4)),((.5,1.5),(2-5e-13,3.5)),((5,6),)])]
    assert retire_completed_material(lineages)==reference(lineages)
    assert not retire_completed_material([lineages[-1]])

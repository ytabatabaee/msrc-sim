from __future__ import annotations

from .lineages import ARGLineage
from .segments import Segment, interval_union, normalize


def union_lineage(a: ARGLineage, b: ARGLineage, parent_node: int | None, time: float, new_id: int) -> ARGLineage:
    """Merge carriers without changing already resolved/disjoint local roots.

    A common-ancestor event changes the local genealogy only where both
    lineages carry material. Disjoint material is carried forward with its
    existing node, avoiding artificial unary ancestry above a completed root.
    """
    bounds = sorted({x for s in a.segments + b.segments for x in (s.left, s.right)})
    segments = []
    for left, right in zip(bounds, bounds[1:]):
        in_a = [s for s in a.segments if s.left <= left and right <= s.right]
        in_b = [s for s in b.segments if s.left <= left and right <= s.right]
        if in_a and in_b and parent_node is not None:
            segments.append(Segment(left, right, parent_node))
        elif in_a:
            segments.append(Segment(left, right, in_a[0].node))
        elif in_b:
            segments.append(Segment(left, right, in_b[0].node))
    active = interval_union(tuple(a.active) + tuple(b.active))
    return ARGLineage(new_id, a.population, normalize(segments, max(s.right for s in segments)), time, a.background, active)


def retire_completed_material(lineages):
    """Retire exactly intervals carried by a single lineage, using one sweep.

    A coverage count plus XOR of carrier indices identifies the sole carrier
    when the count is one. Only those carriers require updates. Canonical
    adjacency tolerance is applied per carrier by with_active, as before.
    """
    changes = {}
    toggles = {}
    for i, lineage in enumerate(lineages):
        for left, right in lineage.active:
            changes[left] = changes.get(left, 0) + 1
            changes[right] = changes.get(right, 0) - 1
            toggles[left] = toggles.get(left, 0) ^ i
            toggles[right] = toggles.get(right, 0) ^ i
    points = sorted(changes)
    if len(points) < 2:
        return list(lineages)
    removals = {}
    count = sole = 0
    for left, right in zip(points, points[1:]):
        count += changes[left]
        sole ^= toggles[left]
        if count == 1:
            removals.setdefault(sole, []).append((left, right))
    replacements = []
    for i, lineage in enumerate(lineages):
        if i not in removals:
            replacements.append(lineage)
            continue
        keep = []
        remove = removals[i]
        j = 0
        for left, right in lineage.active:
            cursor = left
            while j < len(remove) and remove[j][1] <= left:
                j += 1
            while j < len(remove) and remove[j][0] < right:
                a, b = remove[j]
                if a > cursor:
                    keep.append((cursor, min(a, right)))
                cursor = max(cursor, b)
                j += 1
            if cursor < right:
                keep.append((cursor, right))
        replacement = lineage.with_active(keep)
        if replacement.active:
            replacements.append(replacement)
    return replacements

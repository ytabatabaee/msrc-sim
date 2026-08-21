from __future__ import annotations
from typing import Any
import numpy as np
from .wright_fisher import FrequencyHistory


def summarize_frequency_history(history: FrequencyHistory) -> list[dict[str, Any]]:
    rows=[]
    for branch_id, recs in sorted(history.by_branch.items()):
        recs=sorted(recs, key=lambda r:r.absolute_age, reverse=True)
        ps=np.array([r.frequency_A1 for r in recs], dtype=float)
        # Each record represents approximately one generation; the final boundary has zero width.
        widths=np.array([max(0.0, recs[i].absolute_age-recs[i+1].absolute_age) for i in range(len(recs)-1)]+[0.0])
        duration=float(widths.sum())
        integral=float(np.dot(ps,widths))
        segregating=float(sum(w for r,w in zip(recs,widths) if 0.0<r.frequency_A1<1.0))
        fixed1=float(sum(w for r,w in zip(recs,widths) if r.frequency_A1>=1.0))
        fixed0=float(sum(w for r,w in zip(recs,widths) if r.status in {"lost","not_present"} or r.frequency_A1<=0.0))
        rows.append({
            "branch_id":branch_id,
            "parent_branch_id":recs[0].parent_branch_id,
            "start_age":recs[0].absolute_age,
            "end_age":recs[-1].absolute_age,
            "duration":duration,
            "frequency_at_branch_start":float(recs[0].frequency_A1),
            "frequency_at_branch_end":float(recs[-1].frequency_A1),
            "mean_frequency":float(integral/duration) if duration else float(ps[-1]),
            "variance_frequency":float(np.average((ps-(integral/duration if duration else ps[-1]))**2,weights=np.where(widths>0,widths,0))) if duration else 0.0,
            "time_segregating":segregating,
            "time_fixed_A1":fixed1,
            "time_fixed_A0_or_absent":fixed0,
            "integrated_A1_frequency":integral,
            "integrated_A0_frequency":duration-integral,
        })
    return rows

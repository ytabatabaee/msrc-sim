from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from .analytic import TOPOLOGY_PAIRS, compute_Hm, all_configurations


@dataclass(frozen=True)
class ConditionalResult:
    configuration: tuple[int, int, int, int]
    counts: np.ndarray
    empirical_probabilities: np.ndarray
    exact_probabilities: np.ndarray

    @property
    def absolute_error(self):
        return np.abs(self.empirical_probabilities - self.exact_probabilities)


def _parse_configuration(value):
    if isinstance(value, tuple) and len(value) == 4 and all(x in (0, 1) for x in value):
        return value
    s = str(value)
    if len(s) != 4 or any(c not in '01' for c in s):
        raise ValueError('configuration must be a four-character binary string')
    return tuple(int(c) for c in s)


def _events(z, m01, m10, lambda0, lambda1):
    events = []
    for i in range(4):
        rate = m01 if z[i] == 0 else m10
        if rate > 0:
            events.append(('switch', (i,), rate))
    for topology, pairs in TOPOLOGY_PAIRS.items():
        for i, j in pairs:
            if z[i] == z[j]:
                rate = lambda0 if z[i] == 0 else lambda1
                if rate > 0:
                    events.append(('coalescence', (i, j, topology), rate))
    return events


def simulate_one(configuration, duration, m01, m10, lambda0, lambda1, rng):
    z = _parse_configuration(configuration)
    t = 0.0
    while t < duration:
        events = _events(z, m01, m10, lambda0, lambda1)
        total = sum(e[2] for e in events)
        if total <= 0:
            break
        wait = rng.exponential(1.0 / total)
        if t + wait >= duration:
            break
        t += wait
        u = rng.random() * total
        cumulative = 0.0
        for kind, payload, rate in events:
            cumulative += rate
            if u <= cumulative:
                if kind == 'switch':
                    i = payload[0]
                    zz = list(z)
                    zz[i] = 1 - zz[i]
                    z = tuple(zz)
                else:
                    return int(payload[2])
                break
    return int(rng.integers(0, 3))


def simulate_conditional(config):
    si = config['structured_interval']
    mig = si['migration']
    coal = si['coalescence']
    duration = float(si['duration'])
    m01 = float(mig['m01'])
    m10 = float(mig['m10'])
    lambda0 = float(coal['lambda0'])
    lambda1 = float(coal['lambda1'])
    z = _parse_configuration(si['configuration'])

    rng = np.random.default_rng(int(config['seed']))
    counts = np.zeros(3, dtype=int)
    for _ in range(int(config['num_loci'])):
        counts[simulate_one(z, duration, m01, m10, lambda0, lambda1, rng)] += 1

    H, states = compute_Hm(duration, m01, m10, lambda0, lambda1)
    exact = H[states.index(z)]
    empirical = counts / int(config['num_loci'])
    return ConditionalResult(z, counts, empirical, exact)

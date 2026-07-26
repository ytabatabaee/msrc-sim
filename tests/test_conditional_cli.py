import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import yaml

from msrcsim.conditional import simulate_conditional


def test_conditional_simulation_matches_exact():
    config = {
        'mode': 'conditional',
        'seed': 7,
        'num_loci': 50000,
        'structured_interval': {
            'duration': 1.0,
            'configuration': '1010',
            'migration': {'m01': 0.05, 'm10': 0.05},
            'coalescence': {'lambda0': 1.0, 'lambda1': 1.0},
        },
        'output': {'directory': 'unused'},
    }
    result = simulate_conditional(config)
    assert np.max(result.absolute_error) < 0.015


def test_cli_accepts_conditional_mode(tmp_path):
    out = tmp_path / 'out'
    config = {
        'mode': 'conditional',
        'seed': 3,
        'num_loci': 1000,
        'structured_interval': {
            'duration': 0.5,
            'configuration': '1010',
            'migration': {'m01': 0.1, 'm10': 0.1},
            'coalescence': {'lambda0': 1.0, 'lambda1': 1.0},
        },
        'output': {'directory': str(out), 'record_resolved_config': True},
    }
    path = tmp_path / 'config.yaml'
    path.write_text(yaml.safe_dump(config))
    subprocess.run([sys.executable, '-m', 'msrcsim.cli', '--config', str(path)], check=True)
    assert (out / 'summary.json').exists()
    data = json.loads((out / 'summary.json').read_text())
    assert data['mode'] == 'conditional'
    assert data['configuration'] == '1010'

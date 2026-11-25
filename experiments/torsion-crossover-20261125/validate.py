"""Uncapped relation-set equality replay, separate from timed campaign."""
import importlib.util
import json
from pathlib import Path

here = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('campaign', here/'run.py')
campaign = importlib.util.module_from_spec(spec)
spec.loader.exec_module(campaign)
rows = []
for n in (8, 16, 32, 64, 128, 256, 512, 1024):
    xs, fb = campaign.fixture(n)
    panel = campaign.targets(n, fb)
    direct = campaign.direct_scan(fb, panel)
    quotient, status, prefix, orbits, _ = campaign.quotient_table(xs, panel, None)
    assert status == 'complete' and prefix == 288 and quotient == direct
    rows.append({'n': n, 'status': 'matched', 'targets': 288,
                 'factor_base_points': len(fb), 'quotient_orbits': orbits,
                 'relation_pairs': sum(map(len, direct))})
    print(n, 'matched', flush=True)
(here/'validation.json').write_text(json.dumps({'rows': rows}, indent=2)+'\n')

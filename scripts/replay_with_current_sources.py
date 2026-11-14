"""Replay a --check-only run whose recorded source hashes predate later edits of the sources.

Some drivers (scripts/run_ic_chain_depth8.py and the depth-9 run made with it) refuse to replay when
a source file listed under `source_sha256` in the record's verification.json has changed since the
run. For experiments/ic-chain-depth8 that is the case: three sources were edited after the run (see
its NOTES.md). This wrapper leaves the record and the driver untouched. It copies the record
directory to a temporary directory, replaces only the `source_sha256` entries by the hashes of the
current files (printing each one that changed), keeps every `data_sha256` entry, and runs the driver
with `--check-only --out <temporary directory>`. The driver then checks the data hashes and replays
every system, comparing each result with the saved one.

    python scripts/replay_with_current_sources.py experiments/ic-chain-depth8 -- --workers 2
"""
from pathlib import Path
import argparse, hashlib, json, shutil, subprocess, sys, tempfile

ROOT = Path(__file__).resolve().parents[1]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('record', type=Path, help='record directory, e.g. experiments/ic-chain-depth8')
    ap.add_argument('--driver', default='scripts/run_ic_chain_depth8.py')
    ap.add_argument('driver_args', nargs=argparse.REMAINDER, help='arguments after -- go to the driver')
    a = ap.parse_args()
    extra = a.driver_args[1:] if a.driver_args[:1] == ['--'] else a.driver_args
    src = ROOT / a.record
    report = json.loads((src / 'verification.json').read_text())
    with tempfile.TemporaryDirectory() as tmp:
        dst = Path(tmp) / src.name
        shutil.copytree(src, dst)
        for name, recorded in report['source_sha256'].items():
            now = hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
            if now != recorded:
                print(json.dumps(dict(source_changed_since_run=name, recorded=recorded, current=now)), flush=True)
                report['source_sha256'][name] = now
        (dst / 'verification.json').write_text(json.dumps(report, indent=2) + '\n')
        cmd = [sys.executable, str(ROOT / a.driver), '--check-only', '--out', str(dst)] + extra
        sys.exit(subprocess.run(cmd, cwd=ROOT).returncode)


if __name__ == '__main__':
    main()

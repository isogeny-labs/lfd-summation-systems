"""Rebuild every project module without project artifacts, then run the saved verifier.

Mathlib's pinned dependency cache may be populated first. Run only in a fresh
checkout: this script refuses an existing project build directory and deletes nothing.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import signal
import subprocess
import sys
import time

from verify_lean_fall_degrees import lean_code


def module_order(formal):
    files = sorted(formal.glob('*.lean')) + sorted((formal / 'ICFallDegrees').glob('*.lean'))
    modules = {'.'.join(p.relative_to(formal).with_suffix('').parts): p for p in files}
    order, active, done = [], set(), set()

    def visit(name):
        if name in done:
            return
        if name in active:
            raise ValueError(f'cyclic project imports: {name}')
        active.add(name)
        for line in lean_code(modules[name].read_text()).splitlines():
            if line.strip().startswith('import '):
                for dep in line.strip()[7:].split():
                    if dep.startswith('ICFallDegrees'):
                        if dep not in modules:
                            raise ValueError(f'missing project import: {dep}')
                        visit(dep)
        active.remove(name)
        done.add(name)
        order.append(name)

    for name in sorted(modules):
        visit(name)
    return order


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    root = Path(__file__).resolve().parents[1]
    formal = root / 'formal'
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    receipt = dict(status='running', project_build_directory_absent=not (formal / '.lake/build').exists(),
                   platform=platform.platform(), dependency_cache='pinned Mathlib only', stages=[])

    def save():
        (out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')

    def run(name, argv, cwd=formal):
        print(f'Running {name}', flush=True)
        started = time.monotonic()
        stage = dict(name=name, argv=list(map(str, argv)), log=f'{name}.log')
        receipt['stages'].append(stage)
        save()
        with (out / stage['log']).open('w') as log:
            proc = subprocess.Popen(argv, cwd=cwd, env=env, stdout=log,
                                    stderr=subprocess.STDOUT, start_new_session=True)
            try:
                stage['exit_code'] = proc.wait(timeout=7200)
            except BaseException:
                os.killpg(proc.pid, signal.SIGTERM)
                try:
                    proc.wait(timeout=30)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid, signal.SIGKILL)
                    proc.wait()
                raise
            finally:
                stage['elapsed_seconds'] = round(time.monotonic() - started, 3)
                save()
        if stage['exit_code']:
            print((out / stage['log']).read_text()[-12000:], flush=True)
            raise RuntimeError(f'{name} exited {stage["exit_code"]}')

    env = dict(os.environ, MATHLIB_NO_CACHE_ON_UPDATE='1', LEAN_NUM_THREADS='2')
    try:
        if not receipt['project_build_directory_absent']:
            raise RuntimeError('Project build artifacts exist; use a fresh checkout, not a cleanup.')
        baseline = json.loads((root / 'experiments/ic-lean/verification.json').read_text())
        for name, digest in baseline['source_sha256'].items():
            if hashlib.sha256((root / name).read_bytes()).hexdigest() != digest:
                raise RuntimeError(f'Source differs from reviewed audit: {name}')
        receipt['source_sha256'] = baseline['source_sha256']
        receipt['commit'] = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root, text=True).strip()
        receipt['github'] = {k: os.environ[k] for k in
                             ('GITHUB_RUN_ID', 'GITHUB_RUN_ATTEMPT', 'GITHUB_SHA') if k in os.environ}
        lake = subprocess.check_output(['elan', 'which', 'lake'], cwd=formal, text=True).strip()
        receipt['lean_version'] = subprocess.check_output([lake, 'env', 'lean', '--version'],
                                                         cwd=formal, env=env, text=True).strip()
        if '4.19.0' not in receipt['lean_version']:
            raise RuntimeError('Wrong Lean version')
        manifest = json.loads((formal / 'lake-manifest.json').read_text())
        actual = {}
        for package in manifest['packages']:
            actual[package['name']] = subprocess.check_output(
                ['git', 'rev-parse', 'HEAD'], cwd=formal / '.lake/packages' / package['name'], text=True).strip()
            if actual[package['name']] != package['rev']:
                raise RuntimeError(f'Dependency revision mismatch: {package["name"]}')
        receipt['dependency_revisions'] = actual
        receipt['module_order'] = module_order(formal)
        save()
        # Separate Lake calls bound project compilation to one module at a time.
        for i, name in enumerate(receipt['module_order']):
            run(f'{i:03d}-{name}', [lake, '--no-cache', 'build', name])
        run('verifier', [sys.executable, str(root / 'scripts/verify_lean_fall_degrees.py'),
                         '--lake', lake], root)
        for name in ('verification.json', 'build.log', 'axioms.log', 'sparse-demo.log'):
            shutil.copy2(root / 'experiments/ic-lean' / name, out / name)
        fresh = json.loads((out / 'verification.json').read_text())
        if fresh['audited_declarations'] != baseline['audited_declarations']:
            raise RuntimeError('Fresh declaration audit differs from reviewed audit')
        if fresh['source_sha256'] != baseline['source_sha256']:
            raise RuntimeError('Source hashes changed during rebuild')
        receipt.update(status='passed', audited_declarations=len(fresh['audited_declarations']),
                       source_hashes_match=True, declaration_audits_match=True)
    except BaseException as exc:
        receipt.update(status='failed', error=str(exc))
        raise
    finally:
        save()
    print(json.dumps({k: receipt[k] for k in
                      ('status', 'commit', 'audited_declarations', 'source_hashes_match', 'declaration_audits_match')}, indent=2))


if __name__ == '__main__':
    main()

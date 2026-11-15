"""Build every Lean module of formal/ and audit the axioms of every audited declaration.

    python3 scripts/lean_check_all.py --out <directory>

Unlike scripts/verify_lean_fall_degrees.py, which builds the default target (the modules
imported by formal/ICFallDegrees.lean) and audits formal/ICFallDegrees/Audit.lean, this script

  1. scans every .lean file under formal/ for placeholders and escape hatches;
  2. builds every module formal/ICFallDegrees/*.lean, whether or not the root file imports it;
  3. elaborates every audit file (formal/ICFallDegrees/Audit.lean, formal/N8*Audit.lean) and
     every module that prints axioms itself, and collects each `#print axioms` answer;
  4. fails unless every answer lists only propext, Classical.choice and Quot.sound.

It needs the pinned toolchain (formal/lean-toolchain) and the Mathlib build cache
(`lake exe cache get` in formal/; on a machine where that binary does not start, use
`lake env lean --run .lake/packages/mathlib/Cache/Main.lean get`). It writes build.log,
axioms.log and summary.json to --out and changes nothing else in the repository.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent))
from verify_lean_fall_degrees import lean_code  # noqa: E402

ALLOWED = {'propext', 'Classical.choice', 'Quot.sound'}
PROHIBITED = re.compile(r'\b(?:sorry|admit|axiom|unsafe|native_decide|implemented_by|extern)\b')
# Built alone before the rest: with parallel jobs it exhausted 15 GB in the cold build of
# 2026-11-12 (experiments/ic-lean/BUILD-COLD-20261112.md).
HEAVY_FIRST = ['ICFallDegrees.GeneratedPrimeUpper']
HEADLINE = [
    'ICFallDegrees.GrowthFixedE6PrimeSequence.exists_actual_growing_prime_theta_sequence',
]


def run(command, cwd, log, timeout):
    start = time.time()
    proc = subprocess.run(command, cwd=cwd, text=True, capture_output=True, timeout=timeout)
    text = proc.stdout + proc.stderr
    log.write(f'$ {" ".join(command)}\n{text}\n[exit {proc.returncode}, {time.time() - start:.0f} s]\n\n')
    log.flush()
    return proc.returncode, text


def parse_axioms(text):
    found = {name: sorted(set(re.findall(r'[A-Za-z_][A-Za-z0-9_.]*', body)))
             for name, body in re.findall(r"'([^']+)' depends on axioms: \[([^\]]*)\]", text, re.S)}
    found.update({name: [] for name in re.findall(r"'([^']+)' does not depend on any axioms", text)})
    return found


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', type=Path, required=True, help='directory for build.log, axioms.log, summary.json')
    ap.add_argument('--retries', type=int, default=3, help='build attempts (a killed job resumes incrementally)')
    ap.add_argument('--skip-build', action='store_true', help='audit an already built tree')
    args = ap.parse_args()
    root = Path(__file__).resolve().parents[1]
    formal = root / 'formal'
    args.out.mkdir(parents=True, exist_ok=True)

    library = sorted((formal / 'ICFallDegrees').glob('*.lean'))
    top = sorted(formal.glob('*.lean'))
    for path in library + top:
        match = PROHIBITED.search(lean_code(path.read_text()))
        if match:
            raise SystemExit(f'Unapproved source construct {match.group()} in {path}')
    modules = ['ICFallDegrees'] + [f'ICFallDegrees.{p.stem}' for p in library]

    with (args.out / 'build.log').open('w') as log:
        if not args.skip_build:
            for heavy in HEAVY_FIRST:
                code, _ = run(['lake', 'build', heavy], formal, log, 4 * 3600)
                if code:
                    raise SystemExit(f'lake build {heavy} failed with exit code {code}; see build.log')
            for attempt in range(1, args.retries + 1):
                code, text = run(['lake', 'build'] + modules, formal, log, 6 * 3600)
                if code == 0:
                    break
                if 'error:' in text:  # a Lean error, not a killed job: retrying cannot help
                    break
            if code:
                print(text[-6000:])
                raise SystemExit(f'lake build failed with exit code {code}; see build.log')

    # Every `#print axioms` in the project: the audit files, and modules that print their own.
    printing = [p for p in library + top
                if re.search(r'#print\s+axioms', lean_code(p.read_text()))]
    expected, answers, files = {}, {}, {}
    with (args.out / 'axioms.log').open('w') as log:
        for path in printing:
            names = re.findall(r'#print\s+axioms\s+([A-Za-z_][A-Za-z0-9_.]*)', lean_code(path.read_text()))
            code, text = run(['lake', 'env', 'lean', str(path.relative_to(formal))], formal, log, 3600)
            if code:
                raise SystemExit(f'{path.relative_to(root)} does not elaborate (exit {code}); see axioms.log')
            if 'sorryAx' in text or 'Lean.ofReduceBool' in text:
                raise SystemExit(f'{path.relative_to(root)}: sorryAx or native evaluation in the audit output')
            found = parse_axioms(text)
            replies = len(re.findall(r"depends on axioms: \[|does not depend on any axioms", text))
            if replies != len(names):
                raise SystemExit(f'{path.relative_to(root)}: {len(names)} #print axioms commands, '
                                 f'{replies} answers')
            files[str(path.relative_to(root))] = replies
            expected[str(path.relative_to(root))] = names
            answers.update(found)

    offending = {name: axioms for name, axioms in answers.items() if not set(axioms) <= ALLOWED}
    if offending:
        raise SystemExit(f'Declarations with other axioms: {json.dumps(offending, indent=1)}')
    # `#print axioms` prints the fully qualified name; a command inside a namespace gives a short one.
    missing = [h for h in HEADLINE if h not in answers]
    if missing:
        raise SystemExit(f'Headline declarations not audited: {missing}')

    version = subprocess.check_output(['lake', 'env', 'lean', '--version'], cwd=formal, text=True).strip()
    manifest = json.loads((formal / 'lake-manifest.json').read_text())
    mathlib = next(p for p in manifest['packages'] if p['name'] == 'mathlib')['rev']
    hashed = library + top + [formal / 'lean-toolchain', formal / 'lakefile.toml', formal / 'lake-manifest.json']
    summary = dict(
        lean_version=version, mathlib_revision=mathlib,
        build='skipped (--skip-build)' if args.skip_build else 'passed',
        modules_built=len(modules), library_files=len(library), top_level_files=len(top),
        files_with_print_axioms=files, audited_declarations=len(answers),
        allowed_foundational_axioms=sorted(ALLOWED), other_axioms=[],
        headline={h: answers[h] for h in HEADLINE},
        axioms_by_declaration=dict(sorted(answers.items())),
        source_sha256={str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest() for p in hashed},
    )
    (args.out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps({k: summary[k] for k in ('lean_version', 'mathlib_revision', 'build', 'modules_built',
                                              'audited_declarations', 'other_axioms', 'headline')}, indent=2))


if __name__ == '__main__':
    main()

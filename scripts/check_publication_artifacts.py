"""Read-only consistency checks of saved publication artifacts; no proof replay."""
from pathlib import Path
from collections import Counter
import hashlib
import json
import re
import subprocess

from verify_lean_fall_degrees import lean_code


def main():
    root = Path(__file__).resolve().parents[1]
    failures = []
    checks = {}

    def check(name, condition):
        checks[name] = bool(condition)
        if not condition:
            failures.append(name)

    reports = {}
    for name in ('ic-lean', 'ic-curve-paper'):
        report = json.loads((root / 'experiments' / name / 'verification.json').read_text())
        reports[name] = report
        for group in ('source_sha256', 'data_sha256'):
            for path, expected in report.get(group, {}).items():
                file = root / path
                check(f'{name}:{group}:{path}', file.is_file() and
                      hashlib.sha256(file.read_bytes()).hexdigest() == expected)

    audit = (root / 'experiments/ic-lean/axioms.log').read_text()
    dependencies = {name: sorted(set(re.findall(r'[A-Za-z_][A-Za-z0-9_.]*', body)))
                    for name, body in re.findall(
                        r"'([^']+)' depends on axioms: \[([^\]]*)\]", audit, re.S)}
    dependencies.update({name: [] for name in re.findall(
        r"'([^']+)' does not depend on any axioms", audit)})
    expected = set(re.findall(r'#print\s+axioms\s+([A-Za-z_][A-Za-z0-9_.]*)', lean_code(
        (root / 'formal/ICFallDegrees/Audit.lean').read_text())))
    check('audit declarations match source', set(dependencies) == expected)
    check('audit log matches verification report',
          dependencies == reports['ic-lean']['audited_declarations'])
    allowed = {'propext', 'Classical.choice', 'Quot.sound'}
    check('saved audit has only foundational axioms',
          all(set(axioms) <= allowed for axioms in dependencies.values()))
    check('paper build has no warnings or overfull boxes',
          reports['ic-curve-paper']['latex_warnings'] == [] and
          reports['ic-curve-paper']['overfull_boxes'] == 0)

    sweep = root / 'experiments/ic-curve-sweep'
    summary = json.loads((sweep / 'summary.json').read_text())
    cases = [json.loads(line) for line in (sweep / 'cases.jsonl').read_text().splitlines()]
    rows = [json.loads(line) for line in (sweep / 'results.jsonl').read_text().splitlines()]
    systems = {c['system_hash'] for c in cases}
    check('sweep case count', len(cases) == summary['cases'])
    check('sweep has exactly one result per system', len(rows) == len(systems) ==
          len({r['system_hash'] for r in rows}) == summary['unique_systems'])
    check('sweep result identities', {r['system_hash'] for r in rows} == systems)
    check('sweep LFD histogram', dict(Counter(str(r['exact_lfd']) for r in rows)) ==
          summary['lfd_histogram'])

    paper = root / 'paper'
    tex = (paper / 'curve-lfd-bounds.tex').read_text()
    abstract = re.search(r'\\begin\{abstract\}(.*?)\\end\{abstract\}', tex, re.S).group(1)
    abstract = abstract.replace(r'\SWEEPSYSTEMS{}', str(summary['unique_systems']))
    solving = json.loads((root / 'experiments/ic-solving-degree/summary.json').read_text())
    check('solving-degree macros match summary', '\\newcommand{\\SDSYSTEMS}{%s}' %
          f"{solving['exact_lfd_systems']:,}".replace(',', '{,}') in
          (paper / 'solving-degree-data.tex').read_text())
    abstract = abstract.replace(r'\SDSYSTEMS{}', str(solving['exact_lfd_systems']))
    abstract = abstract.replace(r'\"o', 'ö')
    abstract = ' '.join(abstract.replace('~', ' ').split()) + '\n'
    check('submission abstract matches manuscript', abstract ==
          (paper / 'eprint-curve-lfd-bounds/eprint-abstract.txt').read_text())
    sheet = (paper / 'eprint-curve-lfd-bounds/EPRINT.md').read_text()
    recorded = re.findall(r'^\| `([^`]+)` \| `([0-9a-f]{64})` \|$', sheet, re.M)
    check('submission sheet records the manuscript digests', {path for path, _ in recorded} == {
        'paper/curve-lfd-bounds.pdf', 'paper/curve-lfd-bounds.tex',
        'paper/eprint-curve-lfd-bounds/eprint-abstract.txt'})
    for path, expected in recorded:
        file = root / path
        check(f'submission sheet digest:{path}', file.is_file() and
              hashlib.sha256(file.read_bytes()).hexdigest() == expected)
    check('submission sheet abstract length',
          re.search(r'\((\d+) chars', sheet).group(1) == str(len(abstract) - 1))
    citations = {key for group in re.findall(r'\\cite(?:\[[^\]]*\])?\{([^}]+)\}', tex)
                 for key in group.split(',')}
    check('citation keys resolve', citations <= set(re.findall(r'\\bibitem\{([^}]+)\}', tex)))
    pdf = paper / 'curve-lfd-bounds.pdf'
    info = subprocess.check_output(['pdfinfo', str(pdf)], text=True)
    title = 'Last fall degrees of elliptic curve summation systems: machine-checked bounds and exact values'
    check('PDF title metadata', re.search(r'^Title:\s+(.+)$', info, re.M).group(1) == title)
    check('PDF uses A4 or US letter pages', bool(re.search(r'^Page size:.*\((A4|letter)\)', info, re.M)))
    check('submission sheet page count', re.search(r'\((\d+) pages\)', sheet).group(1) ==
          re.search(r'^Pages:\s+(\d+)$', info, re.M).group(1))
    first_page = subprocess.check_output(['pdftotext', '-f', '1', '-l', '1', str(pdf), '-'], text=True)
    check('author and contact on first page', 'Adam Buran' in first_page and
          'a.buran28@gmail.com' in first_page)
    print(json.dumps(dict(status='failed' if failures else 'passed',
                         scope='artifact consistency only; not a fresh Lean build or closure replay',
                         checks=len(checks), audited_declarations=len(dependencies),
                         sweep_systems=len(systems), failures=failures), indent=2))
    raise SystemExit(bool(failures))


if __name__ == '__main__':
    main()

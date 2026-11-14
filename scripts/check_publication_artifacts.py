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
    clean = root / 'experiments/ic-lean/clean-20261114'
    if clean.is_dir():
        from clean_lean_audit import module_order
        receipt = json.loads((clean / 'receipt.json').read_text())
        fresh = json.loads((clean / 'verification.json').read_text())
        check('clean build passed from empty project directory',
              receipt['status'] == 'passed' and receipt['project_build_directory_absent'] is True)
        order = module_order(root / 'formal')
        check('clean build compiled every project module', receipt['module_order'] == order and
              [s['name'] for s in receipt['stages']] ==
              [f'{i:03d}-{name}' for i, name in enumerate(order)] + ['verifier'] and
              all(s['exit_code'] == 0 for s in receipt['stages']))
        check('clean build source hashes match reviewed sources',
              fresh['source_sha256'] == receipt['source_sha256'] == reports['ic-lean']['source_sha256'])
        check('fresh declaration audit matches reviewed audit',
              fresh['audited_declarations'] == dependencies)
        check('fresh audit log matches original audit log',
              (clean / 'axioms.log').read_bytes() == (root / 'experiments/ic-lean/axioms.log').read_bytes())
        for name, digest in json.loads((clean / 'artifact-manifest.json').read_text()).items():
            check(f'clean build artifact:{name}', hashlib.sha256((clean / name).read_bytes()).hexdigest() == digest)
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
    abstract = ' '.join(abstract.replace('~', ' ').split()) + '\n'
    check('submission abstract matches manuscript', abstract ==
          (paper / 'eprint-curve-lfd-bounds/eprint-abstract.txt').read_text())
    citations = {key for group in re.findall(r'\\cite(?:\[[^\]]*\])?\{([^}]+)\}', tex)
                 for key in group.split(',')}
    check('citation keys resolve', citations <= set(re.findall(r'\\bibitem\{([^}]+)\}', tex)))
    pdf = paper / 'curve-lfd-bounds.pdf'
    info = subprocess.check_output(['pdfinfo', str(pdf)], text=True)
    title = 'Last fall degrees of elliptic curve summation systems: machine-checked bounds and exact values'
    check('PDF title metadata', re.search(r'^Title:\s+(.+)$', info, re.M).group(1) == title)
    check('PDF uses A4 or US letter pages', bool(re.search(r'^Page size:.*\((A4|letter)\)', info, re.M)))
    first_page = subprocess.check_output(['pdftotext', '-f', '1', '-l', '1', str(pdf), '-'], text=True)
    check('author and contact on first page', 'Adam Buran' in first_page and
          'a.buran28@gmail.com' in first_page)
    print(json.dumps(dict(status='failed' if failures else 'passed',
                         scope='read-only artifact consistency; no new proof or closure replay',
                         checks=len(checks), audited_declarations=len(dependencies),
                         sweep_systems=len(systems), failures=failures), indent=2))
    raise SystemExit(bool(failures))


if __name__ == '__main__':
    main()

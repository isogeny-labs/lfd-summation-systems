"""Build the pinned Lean project and audit source/proof dependencies."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess


def lean_code(text):
    """Remove nested comments and strings before scanning Lean source tokens."""
    result=[]
    i=0
    depth=0
    quoted=False
    while i<len(text):
        if depth:
            if text.startswith('/-',i):depth+=1;i+=2
            elif text.startswith('-/',i):depth-=1;i+=2
            else:i+=1
        elif quoted:
            if text[i]=='\\':i+=2
            elif text[i]=='"':quoted=False;i+=1
            else:i+=1
        elif text.startswith('/-',i):depth=1;i+=2;result.append(' ')
        elif text.startswith('--',i):
            end=text.find('\n',i);i=len(text) if end<0 else end;result.append(' ')
        elif text[i]=='"':quoted=True;i+=1;result.append(' ')
        else:result.append(text[i]);i+=1
    return ''.join(result)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--lake',type=Path,help='Explicit Lake executable; otherwise use PATH.')
    args=ap.parse_args()
    root=Path(__file__).resolve().parents[1]
    formal=root/'formal'
    out=root/'experiments/ic-lean'
    out.mkdir(parents=True,exist_ok=True)
    lake=args.lake or (Path(shutil.which('lake')) if shutil.which('lake') else None)
    if lake is None:
        raise SystemExit('Lake not found. Install the pinned lean-toolchain, or pass --lake.')
    lake=lake.resolve()
    env=os.environ.copy()
    env['PATH']=str(lake.parent)+os.pathsep+env.get('PATH','')
    env['MATHLIB_NO_CACHE_ON_UPDATE']='1'
    sources=sorted(formal.glob('*.lean'))+sorted((formal/'ICFallDegrees').glob('*.lean'))
    prohibited=re.compile(r'\b(?:sorry|admit|axiom|unsafe|native_decide|implemented_by|extern)\b')
    for path in sources:
        match=prohibited.search(lean_code(path.read_text()))
        if match:raise SystemExit(f'Unapproved source construct {match.group()} in {path}')
    if not (formal/'.lake/packages/mathlib/lakefile.lean').exists():
        subprocess.run([str(lake),'update'],cwd=formal,env=env,check=True)
    build=subprocess.run([str(lake),'--no-cache','build','ICFallDegrees'],cwd=formal,
                         env=env,text=True,capture_output=True,timeout=7200)
    (out/'build.log').write_text(build.stdout+build.stderr)
    if build.returncode:
        print((build.stdout+build.stderr)[-15000:])
        raise SystemExit(build.returncode)
    audit=subprocess.run([str(lake),'env','lean','ICFallDegrees/Audit.lean'],cwd=formal,
                         env=env,text=True,capture_output=True,check=True,timeout=1800)
    audit_text=audit.stdout+audit.stderr
    (out/'axioms.log').write_text(audit_text)
    allowed={'propext','Classical.choice','Quot.sound'}
    declarations=re.findall(r"'([^']+)' depends on axioms: \[([^\]]*)\]",audit_text,re.S)
    no_axioms=re.findall(r"'([^']+)' does not depend on any axioms",audit_text)
    expected=set(re.findall(r'#print\s+axioms\s+([A-Za-z_][A-Za-z0-9_.]*)',
                            lean_code((formal/'ICFallDegrees/Audit.lean').read_text())))
    assert len(expected)>=424
    dependencies={name:sorted(set(re.findall(r'[A-Za-z_][A-Za-z0-9_.]*',body)))
                  for name,body in declarations}
    dependencies.update({name:[] for name in no_axioms})
    assert set(dependencies)==expected,(expected-set(dependencies),set(dependencies)-expected)
    for name,axioms in dependencies.items():
        assert set(axioms)<=allowed,(name,axioms)
    assert 'sorryAx' not in audit_text and 'Lean.ofReduceBool' not in audit_text
    demo=subprocess.run([str(lake),'env','lean','ICFallDegrees/SparseDemo.lean'],cwd=formal,
                        env=env,text=True,capture_output=True,check=True,timeout=1800)
    demo_text=demo.stdout+demo.stderr
    (out/'sparse-demo.log').write_text(demo_text)
    assert '[some 0, some 0, none, none]' in demo_text,demo_text
    version=subprocess.check_output([str(lake),'env','lean','--version'],cwd=formal,env=env,text=True).strip()
    assert '4.19.0' in version,version
    manifest=json.loads((formal/'lake-manifest.json').read_text())
    mathlib=next(package for package in manifest['packages'] if package['name']=='mathlib')
    expected=mathlib['rev']
    actual=subprocess.check_output(['git','-C',str(formal/'.lake/packages/mathlib'),'rev-parse','HEAD'],text=True).strip()
    assert actual==expected
    hashed=sources+[formal/'lean-toolchain',formal/'lakefile.toml',formal/'lake-manifest.json',Path(__file__)]
    report=dict(lean_version=version,mathlib_revision=actual,build='passed',
                audited_declarations=dependencies,allowed_foundational_axioms=sorted(allowed),
                project_axioms=[],placeholders=[],native_decide_used=False,
                executable_sparse_checker='passed: F5 proof trees and F257 shared DAG / containing spaces / affine strict-fall witness',
                checked_prime_gap=dict(prime=257,curve_a=105,curve_b=196,target=[42,56],k=3,variables=6,
                    encoding='six quadratic power-chain equations plus normalized S3',
                    affine_first_strict_fall=3,homogeneous_first_fall=3,
                    homogeneous_convention='ordinary polynomial syzygies modulo the full Koszul submodule',
                    last_fall_degree=4,lower_space_rows=58,upper_dag_nodes=177),
                checked_binary_dimension_two=[
                    dict(name='ordinary-binary',field_size=8,modulus_bits=11,a=0,b=2,target=[5,1],
                         factor_base_basis=[1,2],factor_base_dimension=2,variables=4,m=2,lfd=3,
                         frobenius_invariant=False,lower_rows=10,upper_dag_nodes=25),
                    dict(name='koblitz-noninvariant',field_size=16,modulus_bits=19,a=0,b=1,target=[7,0],
                         factor_base_basis=[1,2],factor_base_dimension=2,variables=4,m=2,lfd=3,
                         frobenius_invariant=False,lower_rows=10,upper_dag_nodes=36)],
                homogeneous_ffd_family=dict(m=2,k=3,variables=6,
                    condition='3*t^2+a=0',ffd=3,coefficient_ring='every nontrivial commutative ring',
                    theorem='ICFallDegrees.PrimeM2Homogeneous.derivative_zero_ffd_three'),
                general_prime_m2=dict(depth='every k>=2',prime='every p>3',summands=2,
                    ordinary_ffd='3 if 3*t^2+a=0, otherwise 4',
                    truncated_ffd='same exact 3/4 value in Fp[z]/(z_i^p), modulo Koszul and Frobenius rows',
                    ordinary_theorem='ICFallDegrees.PrimeFieldConclusions.ordinary_ffd',
                    truncated_theorem='ICFallDegrees.PrimeFieldConclusions.truncated_ffd'),
                growing_prime_sharp_family=dict(curve='y^2=x^3+a*x',target=[0,0],summands=2,
                    depth='every k>=2',hypotheses=['prime p>3','2^k divides p-1','a!=0','a^(2^k)!=1'],
                    normalized_lfd=4,ideal='unit',
                    theorem='ICFallDegrees.PrimeFieldConclusions.sharp_lfd'),
                exact_transport=dict(
                    scalar_extension='Full closure reflection/preservation and exact LFD for every finite coefficient-field extension with an actual basis',
                    moore_matrix='Invertibility proved from an actual finite binary-field basis, without an assumed rank or inverse certificate',
                    inverse_coordinates='Arbitrary-dimensional reconstruction/Frobenius identities, square/inverse rows at2, raw normalization at4',
                    affine='Exact LFD transport with affine images in both directions; same-budget input replacement',
                    scope='Generic transport and forward coordinate links; full inverse normalization/lower/filtered-input specialization remains open'),
                broad_binary_bounds=dict(
                    raw_split='LFD <= max(r+1,4), every target/extension/basis/S3 incidence graph',
                    normalized_split='LFD <= max(r+1,3), every target/extension/basis/S3 incidence graph',
                    coefficient_degree='Basis coordinate extraction does not increase total degree',
                    local_criterion='Bounded algebra section plus local low-degree reductions implies a full-ideal LFD bound',
                    inverse_family_status='All-n parity equality and scaled-subfield extension written; carry rigidity/local rewrites checked in Lean; complete bridge not formalized',
                    binary_fiber='Actual S3 quadratic fiber identity and at-most-two-roots theorem checked in Lean',
                    cyclic_carries='Constructed from modular equality; odd uniqueness/even alternating rigidity checked without a supplied carry witness',
                    even_bridge='All-length alternating product bridge at max(3,d+1), with charged ordinary products',
                    split_lower='Every base-F2 normalized split chain with m>=4 and either target bit has LFD>=4; upper4 remains written'),
                source_sha256={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in hashed},
                limitations=['The graded-lex GB-to-normal-form bridge, finite-quotient theorem, grid/Boolean corollaries, and growing-dimension Boolean sharpness are formalized.',
                             'The concrete prime PKM chains are kernel-checked to be reduced-ideal trivial (ChainQuotient.reduced_chain_ideal_trivial), giving LFD <= max(mk+1, delta) for every m and k.',
                             'The F5 example proves presentation dependence; the F257 normalized curve system proves ordinary Koszul FFD 3, affine first strict fall 3, and exact LFD 4.',
                             'The F2689 FFD4/LFD5 example and the satisfiable Frobenius-invariant F8 case remain Python-checked; the two inconsistent dimension-two binary descents are now kernel-checked.',
                             'The new characteristic-two sharp family has two summands, a dimension-one coordinate base, and an order-two target.',
                             'Binary homogeneous/filtered FFD is not inferred. Prime m=2, all k>=2 now have exact ordinary and truncated FFD classification.',
                             'The semantic checker remains noncomputable; a separate executable sparse checker now has proved soundness and kernel-replayed prime examples.',
                             'Shared-DAG and containing-space checks now replay F257 and two actual four-variable binary descents; the larger F2689 and satisfiable invariant binary imports remain pending.',
                             'The growing-depth prime ECC family has exact LFD 4 but uses an order-two target and inconsistent systems; it does not prove the general mk+1 bound sharp.',
                             'No sharp bound for arbitrary targets/m>=3, growing-dimensional binary/Koblitz bases, or publication novelty is established.'])
    (out/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(build=report['build'],lean=version,mathlib=actual,
                         audited_declarations=len(dependencies),project_axioms=[],placeholders=[]),indent=2))


if __name__=='__main__':
    main()

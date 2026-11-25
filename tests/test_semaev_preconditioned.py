"""The fast shared input builder and preconditioned solver against independent controls."""
import random

from csd import ic_dreg as I
from csd import semaev_matched as M
from csd import semaev_preconditioned as P


def test_fast_builder_equals_audited_inputs():
    for n in (2,3,4,5,7):
        E=I.BinaryCurve(I.binary_field(n),1,3)
        basis=[1,2]
        assert P.lifted_points(E,range(1<<n))==M.rational_points(E)
        for R in P.lifted_points(E,range(1<<n))[:2]:
            for form in ('original','split'):
                assert P.build_system(E,basis,R[0],form)==M.build_system(E,basis,R[0],form)


def test_preconditioned_and_baseline_solve_match_point_enumeration():
    rng=random.Random(20261124)
    for n in (2,3,4,5):
        E=I.BinaryCurve(I.binary_field(n),rng.randrange(2),rng.randrange(1,1<<n))
        basis=[1,2]
        oracle=M.decomposition_map(E,I.span(basis))
        points=P.lifted_points(E,range(1<<n))
        targets=[points[0],points[-1]]
        for R in targets:
            for form in ('original','split'):
                for pre in (False,True):
                    got,metrics=P.solve(E,basis,R,form,pre)
                    assert got==oracle.get(R,set())
                    assert metrics['residual_variables']<=metrics['input']['variables']


def test_saved_preconditioning_record_and_warm_confirmation():
    from pathlib import Path
    import subprocess
    import sys
    root=Path(__file__).resolve().parents[1]
    subprocess.run([sys.executable,str(root/'scripts/check_semaev_preconditioned.py')],
                   cwd=root,check=True,capture_output=True,text=True)
    subprocess.run([sys.executable,str(root/'scripts/run_semaev_preconditioned_warm.py'),'check'],
                   cwd=root,check=True,capture_output=True,text=True)

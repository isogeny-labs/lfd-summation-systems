"""Equality gate of csd.ic_bounds_split_gpu against csd.ic_bounds_split_m (small cells, CPU and
fake-gpu), its negative controls, the symmetry checks, the batch echelon and checkpoints."""
import json
import os
import pathlib
import sys
import tempfile
import unittest

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

from csd import ic_bounds_big as B              # noqa: E402
from csd import ic_bounds_split_gpu as G        # noqa: E402
from csd import ic_bounds_split_m as S          # noqa: E402
import ic_split_gpu_gate as GATE                # noqa: E402

FIELDS = GATE.FIELDS


class EqualityGate(unittest.TestCase):
    def check(self, m, k, D, arms):
        p, gens, prof = GATE.instance(m, k)
        ref = S.closure(gens, m, k, p, D)
        self.assertEqual({f: ref[f] for f in FIELDS}, {f: prof[D][f] for f in FIELDS}, 'reference vs its record')
        for arm in arms:
            rec = GATE.run_candidate(gens, m, k, p, D, arm, ref_gens=gens, batch=256)
            self.assertEqual(GATE.compare(ref, rec), [], f'{(m, k, D, arm)}')
            if arm.endswith(':split'):
                self.assertTrue(rec['symmetry'].startswith('split:'), rec['symmetry'])
                self.assertEqual(len(rec['part_rank']), 2)

    def test_m4_k2_all_degrees_all_arms(self):
        for D in (3, 4, 5, 6):          # certified value 6: closed without 1 through 5, 1 enters at 6
            self.check(4, 2, D, ['cpu:split', 'cpu:plain', 'fake-gpu:split', 'fake-gpu:plain'])

    def test_m5_k2_without_one(self):
        for D in (4, 5):
            self.check(5, 2, D, ['cpu:split', 'fake-gpu:split', 'fake-gpu:plain'])

    def test_m4_k3_without_one(self):
        self.check(4, 3, 5, ['cpu:split', 'fake-gpu:split'])

    def test_exact_value_ladder_m4_k2(self):
        p, gens, _ = GATE.instance(4, 2)
        value = None
        for D in range(3, 8):
            rec, _ = G.closure(gens, 4, 2, p, D, device='fake-gpu')
            if rec['contains_one']:
                value = D
                break
            self.assertTrue(rec['verified_closed'])
        self.assertEqual(value, 6)


class NegativeControls(unittest.TestCase):
    def setUp(self):
        self.m, self.k, self.D = 4, 2, 5
        self.p, self.gens, _ = GATE.instance(self.m, self.k)
        self.ref = S.closure(self.gens, self.m, self.k, self.p, self.D)

    def run_arm(self, gens, arm='cpu:split', **kw):
        rec = GATE.run_candidate(gens, self.m, self.k, self.p, self.D, arm, ref_gens=self.gens, **kw)
        return rec, GATE.compare(self.ref, rec)

    def test_no_fault_passes(self):
        self.assertEqual(self.run_arm(self.gens)[1], [])

    def test_dropped_operator_is_detected_by_dimension_and_by_verification(self):
        rec, diff = self.run_arm(self.gens, fault='drop-operator')
        self.assertIn('dim_V_D', diff)
        self.assertEqual(rec['status'], 'verification-failed')

    def test_wrong_gather_is_detected(self):
        for arm in ('cpu:split', 'cpu:plain'):
            self.assertTrue(self.run_arm(self.gens, arm, fault='wrong-gather')[1])

    def test_perturbed_coefficient_is_seen_only_by_the_space_comparison(self):
        bad = GATE.perturb(self.gens, self.p)
        for arm in ('cpu:split', 'cpu:plain'):
            rec, diff = self.run_arm(bad, arm)
            self.assertEqual(diff, ['space'])        # same dimensions, another space

    def test_asymmetric_input_refuses_the_split(self):
        g0 = dict(self.gens[0])
        key = next(key for key in sorted(g0) if G._swap_generator({key: 1}, self.k) != {key: 1})
        g0[key] = (g0[key] + 1) % self.p or 1
        rec, _ = G.closure([g0] + self.gens[1:], self.m, self.k, self.p, 4)
        self.assertTrue(rec['symmetry'].startswith('refused:'), rec['symmetry'])
        self.assertEqual(len(rec['part_rank']), 1)
        ref = S.closure([g0] + self.gens[1:], self.m, self.k, self.p, 4)
        self.assertEqual({f: ref[f] for f in FIELDS}, {f: rec[f] for f in FIELDS})


class DeviceDiscipline(unittest.TestCase):
    def test_fake_device_rejects_host_arrays(self):
        B.use_device('fake-gpu')
        try:
            d = B.DV.to_dev(np.arange(6, dtype=np.float64).reshape(2, 3))
            with self.assertRaises(TypeError):
                d[:, np.array([0, 1])]
            with self.assertRaises(TypeError):
                d + np.ones((2, 3))
        finally:
            B.use_device('cpu')

    def test_batch_echelon_against_host_elimination(self):
        rng = np.random.default_rng(7)
        p = 7681
        for dev in ('cpu', 'fake-gpu'):
            B.use_device(dev)
            for trial in range(5):
                A = rng.integers(0, p, size=(40, 12)) @ rng.integers(0, p, size=(12, 90)) % p   # rank 12
                A[:, :30] = 0
                A[::3, :60] = 0                       # rows of different leading degree
                fdeg = np.repeat([3, 2, 1], 30)
                st = dict(rounds=0, host_cells=0, madds_rref=0)
                L, lead = G.rref_batch(B.DV.to_dev(A.astype(np.float64)), fdeg, p, B.max_inner(p), 8, st)
                Lh = B.DV.to_host(L).astype(np.int64)
                _, r, _ = G.host_rref_transform(A, p)
                self.assertEqual(len(lead), r)
                self.assertTrue(np.array_equal(Lh[:, lead], np.eye(r, dtype=np.int64)))
                for i in range(r):                    # pivot in the highest degree block of its row
                    self.assertEqual(fdeg[lead[i]], fdeg[np.flatnonzero(Lh[i])].max())
                both = np.concatenate([A, Lh])
                self.assertEqual(G.host_rref_transform(both, p)[1], r)   # same row space
        B.use_device('cpu')

    def test_checkpoint_resume_gives_the_same_closure(self):
        p, gens, prof = GATE.instance(4, 2)
        with tempfile.TemporaryDirectory() as tmp:
            ck = os.path.join(tmp, 'c.npz')
            eng = G.SplitClosure(gens, 4, 2, p, 5, device='fake-gpu', batch=64, checkpoint=ck, budget_s=1e-9)
            eng.add_inputs()
            self.assertTrue(eng.saturate().startswith('stopped'))
            self.assertTrue(os.path.exists(ck))
            rec, _ = G.closure(gens, 4, 2, p, 5, device='fake-gpu', batch=64, checkpoint=ck, resume=True)
            self.assertEqual({f: rec[f] for f in FIELDS}, {f: prof[5][f] for f in FIELDS})
            self.assertTrue(rec['verified_closed'])
            other = GATE.perturb(gens, p)
            with self.assertRaises(ValueError):
                G.closure(other, 4, 2, p, 5, device='fake-gpu', checkpoint=ck, resume=True)


if __name__ == '__main__':
    unittest.main()

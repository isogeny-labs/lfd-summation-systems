#!/usr/bin/env python3
"""Re-run the PRIOR producer's ell=131 construction (../construct_ppav.py, unmodified on disk) in memory
and dump its raw objects for an INDEPENDENT audit in PARI (audit_prior_ppav131.gp).
Nothing in the parent directory is written. Output: prior_ppav131_dump.txt (gp-readable)."""
import sys, pathlib, time
here = pathlib.Path(__file__).resolve().parent
parent = here.parent
sys.path.insert(0, str(parent))
src = (parent / "construct_ppav.py").read_text()
# expose internals: return the omega (form on O_K basis), principal (form on I basis), v, pivot, degree
src = src.replace("    return record\n", "    return record, omega, principal, v, pivot, degree\n")
ns = {"__name__": "prior_construct", "__file__": str(parent / "construct_ppav.py")}
t = time.time()
exec(compile(src, "construct_ppav.py(in-memory)", "exec"), ns)
record, omega, principal, v, pivot, degree = ns["construct"](None)
print("construct() ran in %.1fs" % (time.time() - t))
def gpmat(M):
    return "[" + ";".join(",".join(str(int(M[i, j])) for j in range(M.ncols())) for i in range(M.nrows())) + "]"
with open(here / "prior_ppav131_dump.txt", "w") as f:
    f.write("OMEGA = %s;\n" % gpmat(omega))
    f.write("PRINC = %s;\n" % gpmat(principal))
    f.write("VNUM = %s~;\n" % ("[" + ",".join(map(str, v)) + "]"))
    f.write("PIVOT = %d;\nDEG = %d;\n" % (pivot, degree))
    f.write("SHA = \"%s\";\n" % record["polarization"]["principal_matrix_sha256_decimal_rows"])
print("pivot", pivot, "degree", degree, "sha", record["polarization"]["principal_matrix_sha256_decimal_rows"])

#!/usr/bin/env python3
"""Heuristic expected number of Jacobians in the isogeny class of A_ell over F_2.

E[#Jac] ~ N_ppav(class) * #M_g(F_2)/#A_g(F_2),  #M_g ~ 2^(3g-3), #A_g ~ 2^(g(g+1)/2)  (Lang-Weil leading terms).
N_ppav(class) ~ sum_R |Pic R|/|Pic R+| * |S+(R)|, dominated by R0 = Z[pi,pibar]:
   |Pic R0|/|Pic R0+| ~ h^-(K) * [O_K:R0]/[O_K+:R0+] * (Euler factors) / ([O_K^*:R0^*]/[O_K+^*:R0+^*])
Inputs (computed in this folder): h^-(K) exact (hminus.out, Q=2 forced by integrality), indices (index_plus.out).
Cross-check: DiPippo-Howe count of q-Weil polynomials v_g q^(g(g+1)/4) => average ppav per class.
"""
import json, math
L2 = math.log2
out = {}
# --- inputs from gp outputs (copied; regenerate with hminus.gp / index_plus.gp)
data = {
  11: dict(hminus=5, log2_minus_index=22.618, unit_ratio=(11, 22), splus_OK=0),
  17: dict(hminus=45, log2_minus_index=64.657, unit_ratio=(17, 34), splus_OK=0),
  19: dict(ppav_exist=False, hminus=148, log2_minus_index=79.524, unit_ratio=(19, 38), splus_OK=2),
  131: dict(hminus=1128247071505617895582943540743356840380563482194685535482598575,  # hminus131.out
            log2_minus_index=4229.073, unit_ratio=(131, 262), splus_OK=None),
}
def DH_log2_vg(g):
    s = g - sum(L2(k) for k in range(1, g+1))
    s += sum((g-j)*L2(2*j/(2*j-1)) for j in range(1, g+1))
    return s
rows = {}
for ell, d in data.items():
    g = ell-1
    lo_units, hi_units = L2(d["unit_ratio"][0]), L2(d["unit_ratio"][1])
    central = L2(d["hminus"]) + d["log2_minus_index"] - (lo_units+hi_units)/2
    # band: units ratio (+-0.5), norm cokernel (-1..0), Euler factors (-1.5..0), sum over orders (0..+1), S+ (0..g-1)
    lo = central - 0.5 - 1 - 1.5
    hi = central + 0.5 + 1 + (g-1)
    mg, ag = 3*g-3, g*(g+1)//2
    avg = ag - (g*(g+1)/4 + DH_log2_vg(g))
    if d.get("ppav_exist", True) is False:
        rows[ell] = {"g": g, "note": "lead 1 parity rule: no ppav in the class (eps=-1, only split-type primes); "
                     "the class-number estimate %.1f bits must be gated to 0 -> E[#Jac]=0 (negative control of the gate)" % central}
        continue
    rows[ell] = {
      "g": g, "log2_hminus": round(L2(d["hminus"]), 3), "log2_minus_index": d["log2_minus_index"],
      "log2_Nppav_central(S+=1)": round(central, 2), "log2_Nppav_band": [round(lo, 1), round(hi, 1)],
      "log2_DH_average_ppav_per_class": round(avg, 1),
      "dim_M_g": mg, "dim_A_g": ag,
      "log2_E_Jac_central": round(central + mg - ag, 1),
      "log2_E_Jac_band": [round(lo + mg - ag, 1), round(hi + mg - ag, 1)],
    }
out["rows"] = rows
out["h_minus_K131_exact"] = str(data[131]["hminus"])
out["notes"] = "S+ (polarization units mod norms) observed 0,0,2 bits at ell=11,17,19 for O_K; band allows up to g-1 bits."
print(json.dumps(out, indent=1))

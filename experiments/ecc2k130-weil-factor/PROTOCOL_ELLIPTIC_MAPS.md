# Elliptic-map obstruction probe (2026-11-24)

Input frozen before the probe: explicit_ppav.json and the principal
matrix with SHA-256
4ab7dd31504657e5ad4f180a96dbd9e5838d4cfaeb9c7fe7fccde26b988b94a4.

Question: does the constructed principal polarization admit an elliptic
map of degree in {1,2,3,5,7,11}? This is a test of one polarized
representative, not of every member of the isogeny class.

The Hom(E,B) lattice is J=O+Z*(v mod 263)/263; it is O_K-stable,
unlike I at ell. With Omega the form on O and alpha's multiplication
matrix A, its symmetric Gram matrix is Omega*(2*A+1), restricted to J.
For x in J, half the norm is the elliptic degree.

Use exact-Gram FLINT LLL, delta=0.99, eta=0.51. Inspect all resulting
basis vectors. Independently check every reported vector's norm in the
original Gram matrix. A degree 1 contradicts genus 130 directly. For a
prime degree d<=11, x and zeta*x are independent maps, of the same
degree, and Castelnuovo--Severi gives g<=d^2+1<130. Their function fields
generate the curve field: otherwise the common degree is d and the
maps differ only through a genus-one common factor, contradicting the
M-linear independence of x and zeta*x. Degree 2 is also separable here:
a purely inseparable degree-two map between smooth curves would
preserve genus.

A found vector excludes only the specified principal polarization.
No found vector, or a timeout, is inconclusive. LLL is not a minimum
certificate. Resource cap: 120 seconds for the LLL worker, one attempt.
Save the Gram identity, smallest observed degrees, raw vector if a
criterion is met, runtime and status in elliptic_map_probe.json.

## Adaptive alternative after the exact-Gram timeout

The exact-Gram worker reached the 120-second cap without a result.
Keep that outcome. Run one materially different arithmetic mode,
FLINT approximate-Gram LLL, with a 60-second cap. This changes the
internal reduction arithmetic, not the acceptance criterion: the final
transform identity and any small-vector norm must still be checked
exactly against the input integer Gram matrix. Save this separately
as elliptic_map_probe_approx.json. No further LLL retry is planned in
this bounded probe if this also times out.

import ICFallDegrees.GrowthPrimeIntegralUnit
import Mathlib

/-!
# N8: integral Bezout identities specialize to unit ideals

The actual E6 chain input polynomials have integer coefficients.
If a nonzero integer constant d is a combination of the integer
inputs, then the reduced inputs generate the unit ideal for every
prime not dividing d. The abstract proof is stated for an arbitrary
commutative-ring homomorphism, and therefore works for maps of
multivariate polynomial coefficient rings as soon as an explicit
integral Bezout certificate is supplied.

Unlike a statement about rational unit ideals, this theorem explicitly
records the denominator/exceptional-prime obstruction. In particular
it CANNOT justify an all-primes E6 unit claim (p = 17,k = 2 fails).

No claim of a rational/integral certificate for every E6 depth, nor
of prime-field LFD growth, is made here.
-/

noncomputable section
namespace ICFallDegrees.GrowthIntegralBezoutTransfer

variable {R S : Type*} [CommRing R] [CommRing S]
variable {ι : Type*} [Fintype ι]

/-- A finite Bezout identity with a unit image transports the unit ideal.
This is the exact algebraic interface needed after clearing the
denominators of a rational E6 unit certificate. -/
theorem unit_image_ideal_of_bezout
    (φ : R →+* S) (f q : ι → R) (d : R)
    (h : ∑ i, q i * f i = d) (hd : IsUnit (φ d)) :
    Ideal.span (Set.range (fun i : ι => φ (f i))) = ⊤ := by
  classical
  let I : Ideal S := Ideal.span (Set.range (fun i : ι => φ (f i)))
  have hfi (i : ι) : φ (f i) ∈ I := by
    apply Ideal.subset_span
    exact ⟨i, rfl⟩
  have hdmem : φ d ∈ I := by
    have hsum : (∑ i, φ (q i) * φ (f i)) ∈ I := by
      apply Submodule.sum_mem
      intro i hi
      exact I.mul_mem_left (φ (q i)) (hfi i)
    have hmap : φ d = ∑ i, φ (q i) * φ (f i) := by
      rw [← h]
      simp
    rw [hmap]
    exact hsum
  change I = ⊤
  apply le_antisymm le_top
  intro x _hx
  obtain ⟨v, hv⟩ := (isUnit_iff_dvd_one).mp hd
  have hone : (1 : S) ∈ I := by
    have hmul : φ d * v ∈ I := by
      simpa only [mul_comm] using I.mul_mem_left v hdmem
    rw [hv]
    exact hmul
  simpa only [mul_one] using I.mul_mem_left x hone

/-- The same transport statement for a nonzero image in a field. -/
theorem field_image_ideal_of_bezout
    {K : Type*} [Field K]
    (φ : R →+* K) (f q : ι → R) (d : R)
    (h : ∑ i, q i * f i = d) (hd : φ d ≠ 0) :
    Ideal.span (Set.range (fun i : ι => φ (f i))) = ⊤ := by
  exact unit_image_ideal_of_bezout φ f q d h (isUnit_iff_ne_zero.mpr hd)

#print axioms unit_image_ideal_of_bezout
#print axioms field_image_ideal_of_bezout

end ICFallDegrees.GrowthIntegralBezoutTransfer

# Draft bounded critical-pair criterion

Status: proof draft for review. It is not a Lean theorem or an accepted
large-instance certificate. The current software reports membership only;
its provisional queue completion must not be promoted to exactness.

## Setting

Work over the recorded prime field in the ordinary polynomial ring on
chain variables x_(b,j) and ordinary auxiliary variables. Give every
variable degree one, with graded lexicographic order extending the normal
monomial `(degree, encoded_code)` order. The fixed chain relations are
x_(b,j)^2 - x_(b,j+1), and x_(b,k-1)^2 - 1. Their leading monomials are
pairwise coprime squares; their reductions do not increase total degree.
The chain-normal monomials have squarefree chain masks and unrestricted
auxiliary exponents. All degree bounds used here are below 63, so the
six-bit auxiliary encoding does not overflow.

Let G be a finite monic list of chain-normal polynomials of degree at most
D. Let T_D be the image, in the chain quotient, of the span of all ordinary
monomial multiples h*g satisfying deg(h)+deg(g) <= D, for g in G. Include
the corresponding degree-bounded multiples of the fixed chain relations
before passing to the quotient.

## Conditions to check

1. Every original generator of degree at most D reduces to zero.
2. For every pair g_i,g_j, check its S-polynomial when the degree of the
   leading-monomial lcm is at most D.
3. For each chain variable occurring in the leading monomial of g, check
   the chain-normal product x*g when deg(g)+1 <= D.
4. All reductions in these checks must use products of degree no greater
   than the input pair's degree. A nonzero remainder or a resource cap is
   a failed/incomplete condition check.

Reduction by g uses a normal monomial multiplier disjoint from its leading
chain mask, with the required differences in auxiliary exponents. Its
product has the expected leading monomial: a carry in a lower term lowers
total degree, while carry-free terms preserve the graded lexicographic
comparison. Normalize through the fixed chain relations after each step.

## Proposed argument

First establish the degree-bounded version of the usual critical-pair
argument. On ordinary monomials of degree at most D, rewriting by the
monic G and fixed chain relations terminates because the monomial order
strictly decreases. Any overlap within this degree range has lcm degree
at most D. The G/G overlaps are condition 2. A G/chain overlap that shares
a chain variable is condition 3 after chain normalization. If the leading
monomials are coprime, the usual product reduction resolves the overlap
within its lcm degree. Chain/chain overlaps are of this latter type.

The degree bound must be retained in the overlap argument: multiplying a
critical reduction at lcm L by a context monomial q introduces no term of
degree above deg(q)+deg(L), which is the degree of the overlapping word
and is at most D. Thus the critical-pair proof can be restricted to this
finite degree range without introducing a forbidden product. This is the
technical lemma that needs a careful independent audit/formalization.

With that lemma, reduction on this finite degree range is confluent. Its
normal-form map is linear. Its kernel is exactly the bounded span above:
each reduction changes a polynomial by an allowed multiple; conversely an
allowed multiple of a defining polynomial reduces to zero. Confluence
then applies to sums of these multiples, including cancellations that
lower degree.

If f in T_D has normal total degree below D, division gives a zero normal
form and a representation by multiples whose degrees do not exceed
deg(f). Multiplication of this representation by one variable stays within
D. Chain normalization also stays within D. Therefore T_D is closed under
the degree-admissible variable multiplications in the definition of V_D.
Condition 1 puts every required input generator in T_D, so minimality gives
V_D subset T_D.

For this upper containment, the polynomials in G need not themselves be
proved to belong to V_D. A proper, checked upper space would suffice to
exclude 1 from V_D. A CAS output that fails any condition is not such a
certificate.

For equality, additionally replay a degree-bounded derivation of each g in
G. Then g belongs to V_D, and all its allowed monomial multiples belong
to V_D. This gives T_D subset V_D and hence equality. A directly replayed
constant polynomial proves an upper bound even without this criterion.

If equality is established, leading monomials determine the exact finite
span dimension by counting their upward closure among normal monomials
through D. This counting step must be checked separately; matching a
published dimension alone does not prove input containment.

## Checks and remaining gates

The Python candidate has been compared with direct finite-space closures
on 200 generated controls and all 625 coefficient choices for one
polynomial over the four-element cyclic grid at each of D=1,2,3. Its row
membership proofs were independently replayed. Missing-carry and
missing-pair controls fail as intended. These finite checks support
falsification; they do not prove the technical lemma at arbitrary depth.

Before claiming a large lower bound, finish the bounded-overlap proof,
implement an independent checker for all finite conditions, and bind its
input and polynomial-list hashes to the original instance. Before using
it for the depth14 lift, also establish membership/equality and the
actual archived input binding. None of those gates is waived by this draft.

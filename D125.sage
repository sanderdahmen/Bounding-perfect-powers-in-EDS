### Setup
from modular_method.diophantine_equations.conditions import apply_to_conditional_value
from modular_method.number_fields.galois_group import cyclotomic_galois_isomorphism
load('frey_curves.sage')

### Setting up $E_D(\Q)$
D = 125
E125 = EllipticCurve([0, 0, 0, D, 0])
P, = E125.gens(); T = E125.torsion_points()[0]

### Determining the points on $E_D(\Q)$ with B not divisible by a prime >3
S_integral_points = [E125(P.Coordinates().sage()) for P in magma(E125).SIntegralPoints([2, 3])]
S_integral_points += [-P for P in S_integral_points]
compare = [
    P,
    T,
]
compare += [-P for P in compare]
assert set(S_integral_points) == set(compare)

### The case a = 1
# Corresponding to points m*P with m an integer
a = 1
E1zw = Frey_curve_of_divisibility_sequence(a, D, precision=1)

### Conductor computation for a = 1
e2 = apply_to_conditional_value(lambda E: E.conductor_exponent(2, verbose=True), E1zw)
# Taking into account that 2 \mid B and B is at least a square
z, w = E1zw[0][0].parameters()
con_extra = CongruenceCondition(w^2 - a*z^4, 2^8)
e2 = ConditionalValue([
    (e, con) for e, con in e2
    if not (con & con_extra).pAdic_tree(pAdics=pAdicBase(QQ, 2)).is_empty()
]); e2
# 0 if ('z', 'w') is 1 of 256 possibilities mod 256
# 1 if ('z', 'w') is 1 of 256 possibilities mod 256
assert all(
    mod(ww^2 - a*zz^4, 2^9) == 2^8
    for zz, ww in e2[0][1].pAdic_tree().give_as_congruence_condition()[0]
)
assert all(
    mod(ww^2 - a*zz^4, 2^9) == 0
    for zz, ww in e2[1][1].pAdic_tree().give_as_congruence_condition()[0]
)
assert apply_to_conditional_value(lambda E: E.conductor_exponent(5), E1zw) == 1

### Newform computations for a = 1
nfs1 = apply_to_conditional_value(
    lambda E: E.newform_candidates(bad_primes=[2,5], verbose=2),
    E1zw,
)
nfs1 = ConditionalValue([
    (nfs, con) for nfs, con in nfs1
    if not (con & con_extra).pAdic_tree(pAdics=pAdicBase(QQ, 2)).is_empty()
]); nfs1
# [] if ('z', 'w') is 1 of 128 possibilities mod 128

### The case a = 125
# Corresponding to points m*P + T with m an integer
a = 125
E125zw = Frey_curve_of_divisibility_sequence(a, D, precision=1)

### Q-curve data computations for a = 125
E125zw.splitting_character()
# Dirichlet character modulo 20 of conductor 20 mapping 11 |--> -1, 17 |--> zeta4
K = E125zw.decomposition_field()
L.<zeta40> = CyclotomicField(40)
assert K.is_isomorphic(L.subfield(zeta40 + zeta40^(-1))[0])
Gval = [n for n in range(20) if gcd(n, 40) == 1]
G = [cyclotomic_galois_isomorphism(n, N=40) for n in Gval]
matrix([[E125zw.c(s, t) for t in G] for s in G])
# [ 1  1  1  1  1  1  1  1]
# [ 1 -2 -2  1  1 -2 -2  1]
# [ 1  2  2  1  1  2  2  1]
# [ 1  1  1  1  1  1  1  1]
# [ 1 -1 -1  1  1 -1 -1  1]
# [ 1 -2 -2  1  1 -2 -2  1]
# [ 1  2  2  1  1  2  2  1]
# [ 1 -1 -1  1  1 -1 -1  1]
matrix([[E125zw.c_splitting_map(s, t) for t in G] for s in G])
# Warning: The restriction of scalars of this Q-curve over the decomposition field does not decompose into abelian varieties of GL_2-type. Use the method decomposable_twist to find a twist that does.
# [ 1  1  1  1  1  1  1  1]
# [ 1 -2  2  1  1  2 -2  1]
# [ 1  2  2 -1 -1  2  2  1]
# [ 1  1 -1 -1 -1 -1  1  1]
# [ 1  1 -1 -1 -1 -1  1  1]
# [ 1  2  2 -1 -1  2  2  1]
# [ 1 -2  2  1  1  2 -2  1]
# [ 1  1  1  1  1  1  1  1]
alpha_val = {
    1: 1,
    -1: 1,
    19: 1,
    -19: 1,
    3: zeta40^17 + zeta40^(-17),
    -3: zeta40^17 + zeta40^(-17),
    17: zeta40^17 + zeta40^(-17),
    -17: zeta40^17 + zeta40^(-17),
    7: (zeta40 + zeta40^(-1))^(-1),
    -7: (zeta40 + zeta40^(-1))^(-1),
    13: (zeta40 + zeta40^(-1))^(-1),
    -13: (zeta40 + zeta40^(-1))^(-1),
    9: (zeta40^3 + zeta40^(-3))*(zeta40^9 + zeta40^(-9)),
    -9: (zeta40^3 + zeta40^(-3))*(zeta40^9 + zeta40^(-9)),
    11: (zeta40^3 + zeta40^(-3))*(zeta40^9 + zeta40^(-9)),
    -11: (zeta40^3 + zeta40^(-3))*(zeta40^9 + zeta40^(-9)),
}
alpha = {cyclotomic_galois_isomorphism(key, N=40): value for key, value in alpha_val.items()}
assert all(alpha[s] * s(alpha[t]) * alpha[s*t]^(-1) == E125zw.c(s, t) / E125zw.c_splitting_map(s, t)
           for t in L.galois_group() for s in L.galois_group())
gamma = product(zeta40^k + zeta40^(-k) for k in [1, 2, 3])
assert all(s(gamma) == alpha[s]^2 * gamma for s in L.galois_group())
E125zwg = E125zw.twist(gamma)
assert E125zwg.does_decompose()
K0 = E125zwg.definition_field()
assert K0.is_isomorphic(L.subfield(zeta40^2 + zeta40^(-2))[0])
assert K0 == E125zwg.decomposition_field()

### Newform levels and characters for a = 125
z, w = E125zwg.parameters()
z, w = z.change_ring(QQ), w.change_ring(QQ)
E125zwg._condition = E125zwg._condition & ~CongruenceCondition(w, 5)
N125g = E125zwg.conductor(additive_primes=K0.primes_above(2*5)); N125g
# (64)*Rad_P( ((2141250000000*zeta400^3 - 5032500000000*zeta400^2 - 31000000000000*zeta400 + 72875000000000)) * (z^2 + (1/50*zeta400^2 - 1/5)*w) * (z^2 + (-1/50*zeta400^2 + 1/5)*w)^2 )
E125zwg.splitting_image_field('conjugacy')
# (Number Field in zeta80 with defining polynomial x^2 + 2*x + 2 with zeta80 = -1 - 1*I,
#  Number Field in zeta80 with defining polynomial x^2 + 2*x + 2 with zeta80 = -1 - 1*I)
E125zwg.newform_levels(bad_primes=K0.primes_above(2*5))
# [(1280, 6400), (6400, 1280)]
E125zwg.splitting_character('conjugacy')
# (Dirichlet character modulo 20 of conductor 20 mapping 11 |--> -1, 17 |--> zeta4,
#  Dirichlet character modulo 20 of conductor 20 mapping 11 |--> -1, 17 |--> zeta4)

### Newform computation for a = 125
nfs125 = E125zwg.newform_candidates(bad_primes=K0.primes_above(2*5), algorithm='magma')
assert len(nfs125) == 144
primes = [p for p in prime_range(50) if p != 2 and p != 5]
nfs125 = eliminate_by_traces(E125zwg, nfs125, condition=CoprimeCondition([z, w]),
                             primes=primes, verbose=True)
assert sum(1 for nf in nfs125 if nf[-1] == 0) == 24
assert lcm(nf[-1] for nf in nfs125 if nf[-1] != 0).prime_factors() == [2, 3, 5, 13, 17]
assert all(nf[0].coefficient_field().is_isomorphic(QuadraticField(-1))
           for nf in nfs125 if nf[-1] == 0)

### Considering odd multiples of P1 = P + T
P1 = P + T
assert P1.xy()[0].denominator().prime_factors() == [11]
nfs125P = eliminate_by_trace(E125zwg, nfs125, 11,
                             condition=(CoprimeCondition([z, w]) &
                                        CongruenceCondition(w^2 - a*z^4, 11)),
                             verbose=True)
assert lcm(nf[-1] for nf in nfs125P).prime_factors() == [2, 3, 11]

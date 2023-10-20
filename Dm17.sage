### Setup
from modular_method.diophantine_equations.conditions import apply_to_conditional_value
load("frey_curves.sage")

### Setting up $E_D(\Q)$
D = -17
Em17 = EllipticCurve([0, 0, 0, D, 0])
P, Q = Em17.gens(); T = Em17.torsion_points()[0]

### Determining the points of $E_D(\Q)$ with B not divisible by a prime >3
S_integral_points = [Em17(P.Coordinates().sage()) for P in magma(Em17).SIntegralPoints([2, 3])]
S_integral_points += [-P for P in S_integral_points]
compare = [
    P,
    Q,
    T,
    2*P,
    2*Q,
    P + T,
    Q + T,
    2*Q + T,
    P + Q,
    P - Q,
    2*P - 2*Q,
    P - Q + T,
    P - 2*Q + T,
]
compare += [-P for P in compare]
assert set(S_integral_points) == set(compare)
perfect_powers = [
    P for P in compare
    if sqrt(P.xy()[0].denominator()).is_perfect_power() and
    sqrt(P.xy()[0].denominator()) != 1
]
assert set(perfect_powers) == set([
    2*P,
    -2*P,
    2*Q,
    -2*Q,
    2*Q + T,
    -2*Q + T,
])
assert sqrt((2*P).xy()[0].denominator()) == 2^2
assert sqrt((-2*P).xy()[0].denominator()) == 2^2
assert sqrt((2*Q).xy()[0].denominator()) == 2^2
assert sqrt((-2*Q).xy()[0].denominator()) == 2^2
assert sqrt((2*Q + T).xy()[0].denominator()) == 3^2
assert sqrt((-2*Q + T).xy()[0].denominator()) == 3^2

### The case a = 1
# Corresponding points m*P + n*Q with m + n even
a = 1
E1zw = Frey_curve_of_divisibility_sequence(a, D, precision=1)

### Conductor computation for a = 1
e2 = apply_to_conditional_value(lambda E: E.conductor_exponent(2), E1zw)
# Limit it to the only relevant cases
z, w = E1zw[0][0].parameters()
con_extra = (CongruenceCondition(w^2 - a*z^4 - 1, 2) |
             CongruenceCondition(w^2 - a*z^4, 2^8))
e2 = ConditionalValue([
    (e, con) for e, con in e2
    if not (con & con_extra).pAdic_tree(pAdics=pAdicBase(QQ, 2)).is_empty()
]); e2
# 8 if ('z', 'w') == (1, 0), (3, 0) mod 4
# 0 if ('z', 'w') is 1 of 256 possibilities mod 256
# 1 if ('z', 'w') is 1 of 256 possibilities mod 256
assert all(mod(ww^2 - a*zz^4, 2) == 1
           for zz, ww in e2[0][1].pAdic_tree().give_as_congruence_condition()[0])
assert all(mod(ww^2 - a*zz^4, 2^9) == 2^8
           for zz, ww in e2[1][1].pAdic_tree().give_as_congruence_condition()[0])
assert all(mod(ww^2 - a*zz^4, 2^9) == 0
           for zz, ww in e2[2][1].pAdic_tree().give_as_congruence_condition()[0])
# To save computation time uncomment the next line to set the conductor at 17
# apply_to_conditional_value(lambda E: E.conductor_exponent.set_cache(1, 17), E1zw)
assert apply_to_conditional_value(lambda E: E.conductor_exponent(17), E1zw) == 1

### Newform computation for a = 1
Enfs1 = apply_to_conditional_value(
    lambda E: apply_to_conditional_value(
        lambda nfs: (E, nfs),
        E.newform_candidates(
            bad_primes=[2,17],
            algorithm='magma',
        ),
    ),
    E1zw,
)
Enfs1 = ConditionalValue([
    (Enfs, con) for Enfs, con in Enfs1
    if not (con & con_extra).pAdic_tree(pAdics=pAdicBase(QQ, 2)).is_empty()
])
apply_to_conditional_value(lambda Enfs: len(Enfs[1]), Enfs1)
# 33 if ('z', 'w') == (1, 0), (3, 0) mod 4
# 1  if ('z', 'w') is 1 of 128 possibilities mod 128
Enfs1 = apply_to_conditional_value(
    lambda Enfs: (Enfs[0], eliminate_by_traces(
        Enfs[0],
        Enfs[1],
        condition=CoprimeCondition([z, w]),
        primes=[p for p in prime_range(50) if not p.divides(2*17)],
        verbose=True,
    )),
    Enfs1,
)
apply_to_conditional_value(lambda Enfs: sum(1 for nf in Enfs[1] if nf[-1] == 0), Enfs1)
# 8 if ('z', 'w') == (1, 0), (3, 0) mod 4
# 1 if ('z', 'w') is 1 of 128 possibilities mod 128
apply_to_conditional_value(
    lambda Enfs: lcm(nf[-1] for nf in Enfs[1] if nf[-1] != 0).prime_factors(),
    Enfs1,
)
# [2, 3, 5, 7] if ('z', 'w') == (1, 0), (3, 0) mod 4
# []           if ('z', 'w') is 1 of 128 possibilities mod 128

### Finding the (pseudo)-solutions corresponding to not eliminated newforms for a = 1
bad_nfs = [nf[0] for Enfs, _ in Enfs1 for nf in Enfs[1] if nf[-1] == 0]
assert all(nf.coefficient_field() == QQ for nf in bad_nfs)
bad_j = set(nf._f.EllipticCurve().sage().j_invariant() for nf in bad_nfs)
E1zwj = E1zw[0][0].j_invariant()
polys = set(poly for Efj in bad_j for poly, _ in (E1zwj - Efj).numerator().factor()
            if poly.degree(w) == 1)
# Possible values of w / z^2
wdivzsq = set(-poly(1, 0) / poly(0, 1) for poly in polys)
zw = [(sqrt(val.denominator()), val.numerator()) for val in wdivzsq if val.denominator().is_square()]
assert all(any(poly(z_, w_) == 0 for poly in polys) for z_, w_ in zw)
def tmp(val):
    try:
        return val.nth_root(4)
    except ValueError:
        return 'pseudo'
[(z_, w_, tmp((w_^2 - z_^4) / (-17))) for z_, w_ in zw]
# [(3, -8, 1),
#  (15, 353, 'pseudo'),
#  (3, 8, 1),
#  (12, -145, 'pseudo'),
#  (12, 145, 'pseudo'),
#  (23, 495, 'pseudo')]

### Considering multiples of P1 = 2*P + 2*Q
P1 = 2*P + 2*Q; P1.xy()
# (3568321/451584, 5750178337/303464448)
P1.xy()[0].denominator().prime_factors()
# [2, 3, 7]
C2 = CongruenceCondition(w^2 - a*z^4, 2)
C3 = CongruenceCondition(w^2 - a*z^4, 3)
C7 = CongruenceCondition(w^2 - a*z^4, 7)
Enfs1P = ConditionalValue([
    (Enfs, C & C2) for Enfs, C in Enfs1
    if not (C & C2).pAdic_tree(pAdics=pAdicBase(QQ, 2)).is_empty()
])
Enfs1P = apply_to_conditional_value(
    lambda Enfs: (Enfs[0], eliminate_by_traces(
        Enfs[0],
        Enfs[1],
        condition=CoprimeCondition([z, w]) & C3 & C7,
        primes=[3, 7],
        verbose=True,
    )),
    Enfs1P,
)
apply_to_conditional_value(lambda Enfs: [nf[-1].prime_factors() for nf in Enfs[1]], Enfs1P)
# [[2, 3]]

### The case a = -17
# Corresponding points m*P + n*Q + T with m + n even
a = -17
Em17zw = Frey_curve_of_divisibility_sequence(a, D, precision=1)

### Q-curve data computations for a = -17
K = Em17zw.decomposition_field()
assert K.is_isomorphic(QQ[sqrt(2), sqrt(-17)])
assert Em17zw.splitting_character().conductor() == 1
G = K.galois_group()
sqrt2, sqrtm17 = sqrt(K(2)), sqrt(K(-17))
s2 = next(s for s in G if s != G(1) and s(sqrt2) == sqrt2)
s17 = next(s for s in G if s != G(1) and s(sqrtm17) == sqrtm17)
Gls = [G(1), s2, s17, s2*s17]
matrix([[Em17zw.c(s, t) for t in Gls] for s in Gls])
# [ 1  1  1  1]
# [ 1  2  1  2]
# [ 1 -1  1 -1]
# [ 1 -2  1 -2]
matrix([[Em17zw.c_splitting_map(s, t) for t in Gls] for s in Gls])
# Warning: The restriction of scalars of this Q-curve over the decomposition field does not decompose into abelian varieties of GL_2-type. Use the method decomposable_twist to find a twist that does.
# [1 1 1 1]
# [1 2 1 2]
# [1 1 1 1]
# [1 2 1 2]
alpha = {
    G(1): 1,
    s2: -1,
    s17: (1 - 3*sqrt2) / sqrtm17,
    s2*s17: (1 - 3*sqrt2) / sqrtm17,
}
assert all(alpha[s] * s(alpha[t]) * alpha[s*t]^(-1) == Em17zw.c(s, t) / Em17zw.c_splitting_map(s, t)
           for t in G for s in G)
gamma = 1 + 3*sqrt2
assert all(s(gamma) == alpha[s]^2 * gamma for s in G)
Em17zwg = Em17zw.twist(gamma)
assert Em17zwg.does_decompose()
assert K == Em17zwg.definition_field()
assert K == Em17zwg.decomposition_field()

### Newform levels for a = -17
Set([P.smallest_integer() for P in Em17zwg.primes_of_possible_additive_reduction()])
# Warning: Assuming that (-8160/19*lu^3 - 89760/19*lu - 51680)*z^2 + (48*lu^3 - 288*lu^2 + 2352*lu - 4320)*w and (-8716055040/19*lu^3 - 95876605440/19*lu - 27576944128)*z^6 + (-42688768*lu^3 + 512709120*lu^2 - 2091749632*lu + 7690636800)*z^4*w + (-512709120/19*lu^3 - 5639800320/19*lu - 1622173184)*z^2*w^2 + (-2511104*lu^3 + 30159360*lu^2 - 123044096*lu + 452390400)*w^3 are coprime outside ('(2, 1/76*lu^3 - 1/4*lu^2 + 11/76*lu - 11/4)', '(17, 1/76*lu^3 - 1/4*lu^2 + 11/76*lu - 3/4)', '(17, 1/76*lu^3 - 1/4*lu^2 + 11/76*lu - 27/4)').
# {17, 2}
Nm17 = Em17zwg.conductor(additive_primes=K.primes_above(2*17)); Nm17
# (2, 1/76*lu^3 - 1/4*lu^2 + 11/76*lu - 11/4)^n0*(17, 1/38*lu^3 + 11/38*lu + 6)*(17, 1/38*lu^3 + 11/38*lu - 6)*Rad_P( ((-8716055040/19*lu^3 - 95876605440/19*lu - 27576944128)) * (z^2 + (-1/646*lu^3 - 49/646*lu)*w) * (z^2 + (1/646*lu^3 + 49/646*lu)*w)^2 )
#  where 
# n0 = 16 if ('z', 'w') == (1, 0), (3, 0) mod 4
#      6  if ('z', 'w') is 1 of 8 possibilities mod 8
assert (Nm17.left().right() * Nm17.left().left().right() == K.ideal(17) and
        Nm17.left().left().left().left() == K.prime_above(2))
z, w = Em17zw.parameters()
z = z.change_ring(QQ); w = w.change_ring(QQ)
assert 17*product(p for p, _ in Em17zwg.discriminant().factor()) == w^2 - a*z^4
Nm17.left().left().left().right()[1]
# (6,
#  The condition that ('z', 'w') == (0, 1), (0, 7), (2, 1), (2, 7), (4, 1), (4, 7), (6, 1), (6, 7) mod 8)
levels = Em17zwg.newform_levels(bad_primes=K.primes_above(2*17)); levels
# [(73984, 73984)]               if ('z', 'w') == (1, 0), (3, 0) mod 4
# [(9248, 18496), (18496, 9248)] if ('z', 'w') is 1 of 8 possibilities mod 8
assert (73984 == 2^8 * 17^2 and
        18496 == 2^6 * 17^2 and
        9248 == 2^5 * 17^2)

### Newform computations for a = -17
nfsm17 = Em17zwg.newform_candidates(algorithm='magma', conjugates=False)
# Comment the line above and uncomment the line below to load newforms from a file instead
# nfsm17 = Em17zw.newform_candidates(algorithm='file', path='Dm17am17.nfs', conjugates=False)
primes = [3, 5, 7, 11, 13, 19, 29, 31, 37, 41, 43, 47, 59, 67, 73, 97, 113]
nfsm17 = eliminate_by_traces(Em17zwg, nfsm17, condition=CoprimeCondition([z, w]),
                             primes=primes, verbose=True, use_minpoly=True)
apply_to_conditional_value(lambda nfs: sum(1 for nf in nfs if nf[-1] == 0), nfsm17)
# 6 if ('z', 'w') == (1, 0), (3, 0) mod 4
# 4 if ('z', 'w') is 1 of 8 possibilities mod 8
apply_to_conditional_value(lambda nfs: lcm(nf[-1] for nf in nfs if nf[-1] != 0).prime_factors(), nfsm17)
# [2, 3, 5, 7, 13, 17, 23, 31] if ('z', 'w') == (1, 0), (3, 0) mod 4
# [2, 7, 13, 17]               if ('z', 'w') is 1 of 8 possibilities mod 8

### Considering odd multiples of P1 = P + Q + T
P1 = P + Q + T
assert P1.xy()[0].denominator().prime_factors() == [7]
nfsm17P = eliminate_by_trace(Em17zwg, nfsm17, 7,
                             condition=(CoprimeCondition([z, w]) &
                                        CongruenceCondition(w^2 - a*z^4, 7)),
                             verbose=True, use_minpoly=True )
apply_to_conditional_value(lambda nfs: product(nf[-1] for nf in nfs).prime_factors(), nfsm17P)
# [2, 5, 7, 13, 17] if ('z', 'w') == (1, 0), (3, 0) mod 4
# [2, 7, 17]        if ('z', 'w') is 1 of 8 possibilities mod 8

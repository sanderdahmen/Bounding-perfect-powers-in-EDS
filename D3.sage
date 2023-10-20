### Setup
from modular_method.diophantine_equations.conditions import apply_to_conditional_value
load('frey_curves.sage')

### Setting up $E_D(\Q)$
D = 3
E3 = EllipticCurve([0, 0, 0, D, 0])
P, = E3.gens(); T = E3.torsion_points()[0]
trace_primes = [p for p in prime_range(50) if not p.divides(2*D)]

### Determining the points on $E_D(\Q)$ with B not divisible by a prime >3
S_integral_points = [E3(P.Coordinates().sage()) for P in magma(E3).SIntegralPoints([2, 3])]
S_integral_points += [-P for P in S_integral_points]
compare = [
    P,
    2*P,
    3*P,
    T,
    P + T,
    2*P + T,
]
compare += [-P for P in compare]
assert set(S_integral_points) == set(compare)

### The case a = 1
# Corresponding to points m*P with m an integer
a = 1
E1zw = Frey_curve_of_divisibility_sequence(a, D, precision=1)

### Conductor computation for a = 1
e2 = apply_to_conditional_value(lambda E: E.conductor_exponent(2, verbose=True), E1zw)
# Taking into account that B is at least a square
z, w = E1zw[0][0].parameters()
con_extra = (CongruenceCondition(w^2 - a*z^4, 2^8) |
             ~CongruenceCondition(w^2 - a*z^4, 2))
e2 = ConditionalValue([
    (e, con) for e, con in e2
    if not (con & con_extra).pAdic_tree(pAdics=pAdicBase(QQ, 2)).is_empty()
]); e2
# 8 if ('z', 'w') == (1, 2), (3, 2) mod 4
# 0 if ('z', 'w') is 1 of 256 possibilities mod 256
# 1 if ('z', 'w') is 1 of 256 possibilities mod 256
assert all(
    mod(ww^2 - a*zz^4, 2) == 1
    for zz, ww in e2[0][1].pAdic_tree().give_as_congruence_condition()[0]
)
assert all(
    mod(ww^2 - a*zz^4, 2^9) == 2^8
    for zz, ww in e2[1][1].pAdic_tree().give_as_congruence_condition()[0]
)
assert all(
    mod(ww^2 - a*zz^4, 2^9) == 0
    for zz, ww in e2[2][1].pAdic_tree().give_as_congruence_condition()[0]
)
assert apply_to_conditional_value(lambda E: E.conductor_exponent(3), E1zw) == 1

### Newform computations for a = 1
Enfs1 = apply_to_conditional_value(
    lambda E: apply_to_conditional_value(
        lambda nfs: (E, nfs),
        E.newform_candidates(
            bad_primes=[2,3],
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
# 12 if ('z', 'w') == (1, 2), (3, 2) mod 4
# 0  if ('z', 'w') is 1 of 640 possibilities mod 128
Enfs1 = apply_to_conditional_value(
    lambda Enfs: (Enfs[0], eliminate_by_traces(
        Enfs[0],
        Enfs[1],
        primes=trace_primes,
        condition=CoprimeCondition([z, w]),
        verbose=True,
    )),
    Enfs1,
)
apply_to_conditional_value(
    lambda Enfs: sum(1 for nf in Enfs[1] if nf[-1] == 0),
    Enfs1,
)
# 8 if ('z', 'w') == (1, 2), (3, 2) mod 4
# 0 if ('z', 'w') is 1 of 640 possibilities mod 128
apply_to_conditional_value(
    lambda Enfs: lcm(nf[-1] for nf in Enfs[1] if nf[-1] != 0).prime_factors(),
    Enfs1,
)
# [2, 3, 7] if ('z', 'w') == (1, 2), (3, 2) mod 4
# []        if ('z', 'w') is 1 of 640 possibilities mod 128

### Considering multiples of P1 = 2*P
P1 = 2*P; P1.xy()
# (1/4, -7/8)
assert P1.xy()[0].denominator().prime_factors() == [2]
Enfs1P = ConditionalValue([
    (Enfs, con) for Enfs, con in Enfs1
    if not (con & CongruenceCondition(w^2 - a*z^4, 2)).pAdic_tree(pAdics=pAdicBase(QQ, 2)).is_empty()
])
assert apply_to_conditional_value(
    lambda Enfs: len(Enfs[1]),
    Enfs1P,
) == 0

### The case a = 3
# Corresponding to points m*P + T with m an integer
a = 3
E3zw = Frey_curve_of_divisibility_sequence(a, D, precision=1)

### Q-curve data computations for a = 3
E3zwg = E3zw.decomposable_twist()
K = E3zwg.definition_field()
assert K == E3zwg.decomposition_field()

### Newform computation for a = 3
nfs3 = E3zwg.newform_candidates(bad_primes=K.primes_above(2*D), algorithm='magma')
apply_to_conditional_value(len, nfs3)
# 32 if ('z', 'w') == (1, 2), (3, 2) mod 4 and ('z', 'w') is 1 of 6 possibilities mod 3
# 28 if ('z', 'w') == (1, 2), (3, 2) mod 4 and ('z', 'w') == (1, 0), (2, 0) mod 3
# 4  if ('z', 'w') is 1 of 4 possibilities mod 8 and ('z', 'w') is 1 of 6 possibilities mod 3 or ('z', 'w') is 1 of 4 possibilities mod 8 and ('z', 'w') == (1, 0), (2, 0) mod 3 or ('z', 'w') is 1 of 4 possibilities mod 8 and ('z', 'w') is 1 of 6 possibilities mod 3
# 0  if ('z', 'w') is 1 of 4 possibilities mod 8 and ('z', 'w') == (1, 0), (2, 0) mod 3
z, w = E3zwg.parameters()
z, w = z.change_ring(QQ), w.change_ring(QQ)
nfs3 = eliminate_by_traces(
    E3zwg,
    nfs3,
    condition=CoprimeCondition([z, w]),
    primes=trace_primes,
    verbose=True,
)
apply_to_conditional_value(
    lambda nfs: sum(1 for nf in nfs if nf[-1] == 0),
    nfs3,
)
# 16 if ('z', 'w') == (1, 2), (3, 2) mod 4 and ('z', 'w') is 1 of 6 possibilities mod 3
# 12 if ('z', 'w') == (1, 2), (3, 2) mod 4 and ('z', 'w') == (1, 0), (2, 0) mod 3
# 4  if ('z', 'w') is 1 of 4 possibilities mod 8 and ('z', 'w') is 1 of 6 possibilities mod 3 or ('z', 'w') is 1 of 4 possibilities mod 8 and ('z', 'w') is 1 of 6 possibilities mod 3
# 0  if ('z', 'w') is 1 of 4 possibilities mod 8 and ('z', 'w') == (1, 0), (2, 0) mod 3 or ('z', 'w') is 1 of 4 possibilities mod 8 and ('z', 'w') == (1, 0), (2, 0) mod 3
apply_to_conditional_value(
    lambda nfs: lcm(nf[-1] for nf in nfs if nf[-1] != 0).prime_factors(),
    nfs3,
)
# [2]               if ('z', 'w') == (1, 2), (3, 2) mod 4 and ('z', 'w') is 1 of 6 possibilities mod 3
# [2, 3, 7, 11, 17] if ('z', 'w') == (1, 2), (3, 2) mod 4 and ('z', 'w') == (1, 0), (2, 0) mod 3
# []                if ('z', 'w') is 1 of 4 possibilities mod 8 and ('z', 'w') is 1 of 6 possibilities mod 3 or ('z', 'w') is 1 of 4 possibilities mod 8 and ('z', 'w') is 1 of 6 possibilities mod 3 or ('z', 'w') is 1 of 4 possibilities mod 8 and ('z', 'w') == (1, 0), (2, 0) mod 3
# [2, 3, 5]         if ('z', 'w') is 1 of 4 possibilities mod 8 and ('z', 'w') == (1, 0), (2, 0) mod 3

### Considering odd multiples of P1 = 3*P + T
P1 = 3*P + T; P1.xy()
# (27/121, 1098/1331)
assert P1.xy()[0].denominator().prime_factors() == [11]
nfs3P = eliminate_by_trace(E3zwg, nfs3, 11,
                             condition=(CoprimeCondition([z, w]) &
                                        CongruenceCondition(w^2 - a*z^4, 11)),
                             verbose=True)
apply_to_conditional_value(
    lambda nfs: lcm(nf[-1] for nf in nfs).prime_factors(),
    nfs3P
)
# [2, 3, 11]           if ('z', 'w') == (1, 2), (3, 2) mod 4 and ('z', 'w') is 1 of 6 possibilities mod 3 or ('z', 'w') is 1 of 4 possibilities mod 8 and ('z', 'w') is 1 of 6 possibilities mod 3 or ('z', 'w') is 1 of 4 possibilities mod 8 and ('z', 'w') is 1 of 6 possibilities mod 3
# [2, 3, 5, 7, 11, 17] if ('z', 'w') == (1, 2), (3, 2) mod 4 and ('z', 'w') == (1, 0), (2, 0) mod 3
# [2, 3]               if ('z', 'w') is 1 of 4 possibilities mod 8 and ('z', 'w') == (1, 0), (2, 0) mod 3
# []                   if ('z', 'w') is 1 of 4 possibilities mod 8 and ('z', 'w') == (1, 0), (2, 0) mod 3

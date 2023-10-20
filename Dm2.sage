### Setup
from modular_method.diophantine_equations.conditions import apply_to_conditional_value
load('frey_curves.sage')

### Setting up $E_D(\Q)$
D = -2
Em2 = EllipticCurve([0, 0, 0, D, 0])
P, = Em2.gens(); T = Em2.torsion_points()[0]
trace_primes = [p for p in prime_range(50) if not p.divides(2*D)]

### Determining the points on $E_D(\Q)$ with B not divisible by a prime >3
S_integral_points = [Em2(P.Coordinates().sage()) for P in magma(Em2).SIntegralPoints([2, 3])]
S_integral_points += [-P for P in S_integral_points]
compare = [
    P,
    2*P,
    T,
    P + T,
    2*P + T,
    3*P + T,
]
compare += [-P for P in compare]
assert set(S_integral_points) == set(compare)

### The case a = 1
# Corresponding to points m*P with m even
a = 1
E1zw = Frey_curve_of_divisibility_sequence(a, D, precision=1)

### Conductor computation for a = 1
e2 = apply_to_conditional_value(lambda E: E.conductor_exponent(2, verbose=True), E1zw)
# Taking into account that B is at least a square
z, w = E1zw[0][0].parameters()
con_extra = (CongruenceCondition(w^2 - a*z^4, 2^9) |
             ~CongruenceCondition(w^2 - a*z^4, 2^2))
e2 = ConditionalValue([
    (e, con) for e, con in e2
    if not (con & con_extra).pAdic_tree(pAdics=pAdicBase(QQ, 2)).is_empty()
]); e2
# 1 if ('z', 'w') is 1 of 256 possibilities mod 256
assert all(
    mod(ww^2 - a*zz^4, 2^9) == 0
    for zz, ww in e2[0][1].pAdic_tree().give_as_congruence_condition()[0]
)

### Newform computations for a = 1
Enfs1 = apply_to_conditional_value(
    lambda E: apply_to_conditional_value(
        lambda nfs: (E, nfs),
        E.newform_candidates(
            bad_primes=(2*D).prime_factors(),
            algorithm='magma',
        ),
    ),
    E1zw,
)
Enfs1 = ConditionalValue([
    (Enfs, con) for Enfs, con in Enfs1
    if not (con & con_extra).pAdic_tree(pAdics=pAdicBase(QQ, 2)).is_empty()
])
assert apply_to_conditional_value(lambda Enfs: len(Enfs[1]), Enfs1) == 0

### The case a = -1
# Corresponding to points m*P with m odd
a = -1
Em1zw = Frey_curve_of_divisibility_sequence(a, D, precision=1)

### Q-curve data computations for a = -1
Em1zwg = Em1zw.decomposable_twist()
K = Em1zwg.definition_field()
assert K == Em1zwg.decomposition_field()

### Newform computation for a = -1
nfsm1 = Em1zwg.newform_candidates(bad_primes=K.primes_above(2*D), algorithm='magma')
assert apply_to_conditional_value(len, nfsm1) == 16
z, w = Em1zwg.parameters()
z, w = z.change_ring(QQ), w.change_ring(QQ)
nfsm1 = eliminate_by_traces(
    Em1zwg,
    nfsm1,
    condition=CoprimeCondition([z, w]),
    primes=trace_primes,
    verbose=True,
)
assert apply_to_conditional_value(
    lambda nfs: sum(1 for nf in nfs if nf[-1] == 0),
    nfsm1,
) == 4
assert apply_to_conditional_value(
    lambda nfs: lcm(nf[-1] for nf in nfs if nf[-1] != 0).prime_factors(),
    nfsm1,
) == [2, 3, 7]

### Considering odd multiples of P1 = 3*P
P1 = 3*P; P1.xy()
# (-1/169, 239/2197)
assert P1.xy()[0].denominator().prime_factors() == [13]
nfsm1P = eliminate_by_trace(
    Em1zwg,
    nfsm1,
    13,
    condition=(CoprimeCondition([z, w]) &
               CongruenceCondition(w^2 - a*z^4, 13)),
    verbose=True,
)
assert apply_to_conditional_value(
    lambda nfs: lcm(nf[-1] for nf in nfs).prime_factors(),
    nfsm1P
) == [2, 5, 7, 13]
# l-th power must be multiple of 13*P1
assert [
    p for p in prime_range(100)
    if p.divides((13*P1).xy()[0].denominator())
] == [13, 37, 41]
nfsm1P = eliminate_by_trace(
    Em1zwg,
    nfsm1P,
    37,
    condition=(CoprimeCondition([z, w]) &
               CongruenceCondition(w^2 - a*z^4, 37)),
    verbose=True,
)
nfsm1P = eliminate_by_trace(
    Em1zwg,
    nfsm1P,
    41,
    condition=(CoprimeCondition([z, w]) &
               CongruenceCondition(w^2 - a*z^4, 41)),
    verbose=True,
)
assert apply_to_conditional_value(
    lambda nfs: lcm(nf[-1] for nf in nfs).prime_factors(),
    nfsm1P
) == [2]

### The case a = 2
# Corresponding to points m*P + T with m odd
a = 2
E2zw = Frey_curve_of_divisibility_sequence(a, D, precision=1)

### Q-curve data computations for a = 2
E2zwg = E2zw.decomposable_twist()
K = E2zwg.definition_field()
assert K == E2zwg.decomposition_field()

### Newform computation for a = 2
nfs2 = E2zwg.newform_candidates(bad_primes=K.primes_above(2*D), algorithm='magma')
assert apply_to_conditional_value(len, nfs2) == 28
z, w = E2zwg.parameters()
z, w = z.change_ring(QQ), w.change_ring(QQ)
nfs2 = eliminate_by_traces(
    E2zwg,
    nfs2,
    condition=CoprimeCondition([z, w]),
    primes=trace_primes,
    verbose=True,
)
assert apply_to_conditional_value(
    lambda nfs: sum(1 for nf in nfs if nf[-1] == 0),
    nfs2,
) == 12
assert apply_to_conditional_value(
    lambda nfs: lcm(nf[-1] for nf in nfs if nf[-1] != 0).prime_factors(),
    nfs2,
) == [2]

### Considering odd multiples of P1 = 5*P + T
P1 = 5*P + T; P1.xy()
# (4651250/1803649, -8388283850/2422300607)
assert P1.xy()[0].denominator().prime_factors() == [17, 79]
nfs2P = eliminate_by_trace(
    E2zwg,
    nfs2,
    17,
    condition=(CoprimeCondition([z, w]) &
               CongruenceCondition(w^2 - a*z^4, 17)),
    verbose=True,
)
nfs2P = eliminate_by_trace(
    E2zwg,
    nfs2,
    79,
    condition=(CoprimeCondition([z, w]) &
               CongruenceCondition(w^2 - a*z^4, 79)),
    verbose=True,
)
assert apply_to_conditional_value(
    lambda nfs: lcm(nf[-1] for nf in nfs).prime_factors(),
    nfs2P
) == [2, 5, 79]

### The case a = -2
# Corresponding to points m*P + T with m even
a = -2
Em2zw = Frey_curve_of_divisibility_sequence(a, D, precision=1)

### Q-curve data computations for a = -2
Em2zwg = Em2zw.decomposable_twist()
K = Em2zwg.definition_field()
assert K == Em2zwg.decomposition_field()

### Newform computation for a = -2
nfsm2 = Em2zwg.newform_candidates(bad_primes=K.primes_above(2*D), algorithm='magma')
assert apply_to_conditional_value(len, nfsm2) == 28
z, w = Em2zwg.parameters()
z, w = z.change_ring(QQ), w.change_ring(QQ)
nfsm2 = eliminate_by_traces(
    Em2zwg,
    nfsm2,
    condition=CoprimeCondition([z, w]),
    primes=trace_primes,
    verbose=True,
)
assert apply_to_conditional_value(
    lambda nfs: sum(1 for nf in nfs if nf[-1] == 0),
    nfsm2,
) == 4
assert apply_to_conditional_value(
    lambda nfs: lcm(nf[-1] for nf in nfs if nf[-1] != 0).prime_factors(),
    nfsm2,
) == [2, 3, 7]

### Considering odd multiples of P1 = 2*P + T
P1 = 2*P + T; P1.xy()
# (-8/9, -28/27)
assert P1.xy()[0].denominator().prime_factors() == [3]
nfsm2P = eliminate_by_trace(
    Em2zwg,
    nfsm2,
    3,
    condition=(CoprimeCondition([z, w]) &
               CongruenceCondition(w^2 - a*z^4, 3)),
    verbose=True,
)
assert apply_to_conditional_value(
    lambda nfs: lcm(nf[-1] for nf in nfs).prime_factors(),
    nfsm2P
) == [2, 3, 7]
# l-th power must be multiple of 3*P1
assert [
    p for p in prime_range(100)
    if p.divides((3*P1).xy()[0].denominator())
] == [3, 11]
nfsm2P = eliminate_by_trace(
    Em2zwg,
    nfsm2P,
    11,
    condition=(CoprimeCondition([z, w]) &
               CongruenceCondition(w^2 - a*z^4, 11)),
    verbose=True,
)
assert apply_to_conditional_value(
    lambda nfs: lcm(nf[-1] for nf in nfs).prime_factors(),
    nfsm2P
) == [2, 3]

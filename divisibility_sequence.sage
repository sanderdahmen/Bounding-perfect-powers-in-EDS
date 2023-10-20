### Doing the long conductor computations for a = 1 once:
R.<z, w> = QQ[]
DB4 = w^2 - z^4
C = OrderCondition(DB4, 3) | (OrderCondition(z, 0) & PowerCondition(DB4, 8))
E1 = FreyCurve([0, 4*z, 0, 2*(z^2 - w), 0], condition=C)
Em1 = FreyCurve([0, -4*z, 0, 2*(z^2 - w), 0], condition=C)
E2 = FreyCurve([0, -4*z, 0, 2*(z^2 + w), 0], condition=C)
Em2 = FreyCurve([0, 4*z, 0, 2*(z^2 + w), 0], condition=C)
### Conductor exponent computation for E1
%time N1 = E1.conductor_exponent(2); N1
# 8 if ('z', 'w') == (0, 1), (1, 0) mod 2
# 7 if ('z', 'w') is 1 of 8 possibilities mod 8
# 6 if ('z', 'w') is 1 of 1088 possibilities mod 128
# 5 if ('z', 'w') is 1 of 4 possibilities mod 8
# 0 if ('z', 'w') is 1 of 64 possibilities mod 256
# 1 if ('z', 'w') is 1 of 64 possibilities mod 256
# 4 if ('z', 'w') is 1 of 32 possibilities mod 128
# Total time: 6h 10min 48s
# Human readable format:
# 0 if z == 1 (mod 4) and w - z^2 == 128 (mod 256) 
# 1 if z == 1 (mod 4) and w - z^2 == 0   (mod 256)
# 4 if z == 3 (mod 4) and w - z^2 == 0   (mod 128)
# 5 if z == 0 (mod 2) and w - z^2 == 6   (mod 8)
# 6 if z == 0 (mod 2) and w - z^2 == 2   (mod 8)
#   or z == 1 (mod 2) and w + z^2 == 0   (mod 128)
# 7 if z == 1 (mod 2) and w - z^2 == 2   (mod 8) # Combined
#   or z == 1 (mod 2) and w - z^2 == 4   (mod 8) # w^2 - z^4 == 8 (mod 16)
# 8 if z == 0 (mod 1) and w - z^2 == 1   (mod 2)
### Conductor exponent computation for Em1
%time Nm1 = Em1.conductor_exponent(2); Nm1
# 8 if ('z', 'w') == (0, 1), (1, 0) mod 2
# 7 if ('z', 'w') is 1 of 8 possibilities mod 8
# 6 if ('z', 'w') is 1 of 1088 possibilities mod 128
# 5 if ('z', 'w') is 1 of 4 possibilities mod 8
# 0 if ('z', 'w') is 1 of 64 possibilities mod 256
# 1 if ('z', 'w') is 1 of 64 possibilities mod 256
# 4 if ('z', 'w') is 1 of 32 possibilities mod 128
# Total time: 6h 18min 25s
# Human readable format:
# 0 if z == 3 (mod 4) and w - z^2 == 128 (mod 256)
# 1 if z == 3 (mod 4) and w - z^2 == 0   (mod 256)
# 4 if z == 1 (mod 4) and w - z^2 == 0   (mod 128)
# 5 if z == 0 (mod 2) and w - z^2 == 6   (mod 8)
# 6 if z == 0 (mod 2) and w - z^2 == 2   (mod 8)
#   or z == 1 (mod 2) and w + z^2 == 0   (mod 128)
# 7 if z == 1 (mod 2) and w - z^2 == 2   (mod 8) # Combined
#   or z == 1 (mod 2) and w - z^2 == 4   (mod 8) # w^2 - z^4 == 8 (mod 16)
# 8 if z == 0 (mod 1) and w - z^2 == 1   (mod 2)
### Conductor exponent computation for E2
%time N2 = E2.conductor_exponent(2); N2
# 8 if ('z', 'w') == (0, 1), (1, 0) mod 2
# 7 if ('z', 'w') is 1 of 8 possibilities mod 8
# 6 if ('z', 'w') is 1 of 1088 possibilities mod 128
# 5 if ('z', 'w') is 1 of 4 possibilities mod 8
# 0 if ('z', 'w') is 1 of 64 possibilities mod 256
# 1 if ('z', 'w') is 1 of 64 possibilities mod 256
# 4 if ('z', 'w') is 1 of 32 possibilities mod 128
# Total time: 6h 5min 48s
# Human readable format:
# 0 if z == 3 (mod 4) and w + z^2 == 128 (mod 256)
# 1 if z == 3 (mod 4) and w + z^2 == 0   (mod 256)
# 4 if z == 1 (mod 4) and w + z^2 == 0   (mod 128)
# 5 if z == 0 (mod 2) and w + z^2 == 2   (mod 8)
# 6 if z == 0 (mod 2) and w + z^2 == 6   (mod 8)
#   or z == 1 (mod 2) and w - z^2 == 0   (mod 128)
# 7 if z == 1 (mod 2) and w + z^2 == 4   (mod 8) # Combined
#   or z == 1 (mod 2) and w + z^2 == 6   (mod 8) # w^2 - z^4 == 8 (mod 16)
# 8 if z == 0 (mod 1) and w + z^2 == 1   (mod 2)
### Conductor exponent computation for Em2
%time Nm2 = Em2.conductor_exponent(2); Nm2
# 8 if ('z', 'w') == (0, 1), (1, 0) mod 2
# 7 if ('z', 'w') is 1 of 8 possibilities mod 8
# 6 if ('z', 'w') is 1 of 1088 possibilities mod 128
# 5 if ('z', 'w') is 1 of 4 possibilities mod 8
# 0 if ('z', 'w') is 1 of 64 possibilities mod 256
# 1 if ('z', 'w') is 1 of 64 possibilities mod 256
# 4 if ('z', 'w') is 1 of 32 possibilities mod 128
# Total time: 6h 24min 24s
# Human readable format:
# 0 if z == 1 (mod 4) and w + z^2 == 128 (mod 256)
# 1 if z == 1 (mod 4) and w + z^2 == 0   (mod 256)
# 4 if z == 3 (mod 4) and w + z^2 == 0   (mod 128)
# 5 if z == 0 (mod 2) and w + z^2 == 2   (mod 8)
# 6 if z == 0 (mod 2) and w + z^2 == 6   (mod 8)
#   or z == 1 (mod 2) and w - z^2 == 0   (mod 128)
# 7 if z == 1 (mod 2) and w + z^2 == 4   (mod 8) # Combined
#   or z == 1 (mod 2) and w + z^2 == 6   (mod 8) # w^2 - z^4 == 8 (mod 16)
# 8 if z == 0 (mod 1) and w + z^2 == 1   (mod 2)
### Overall bests:
# 0 if z == 1 (mod 4) and w - z^2 == 128 (mod 256) with E1
#   or z == 1 (mod 4) and w + z^2 == 128 (mod 256) with Em2
#   or z == 3 (mod 4) and w - z^2 == 128 (mod 256) with Em1
#   or z == 3 (mod 4) and w + z^2 == 128 (mod 256) with E2
# 1 if z == 1 (mod 4) and w - z^2 == 0   (mod 256) with E1
#   or z == 1 (mod 4) and w + z^2 == 0   (mod 256) with Em2
#   or z == 3 (mod 4) and w - z^2 == 0   (mod 256) with Em1
#   or z == 3 (mod 4) and w + z^2 == 0   (mod 256) with E2
# 5 if z == 0 (mod 2) and w + z^2 == 2   (mod 8)   with E2/Em2 # Can not happen
#   or z == 0 (mod 2) and w - z^2 == 6   (mod 8)   with E1/Em1 # as z and w can't
# 6 if z == 0 (mod 2) and w - z^2 == 2   (mod 8)   with E1/Em1 # simultaneously be
#   or z == 0 (mod 2) and w + z^2 == 6   (mod 8)   with E2/Em2 # even
# 7 if z == 1 (mod 2) and w - z^2 == 2   (mod 8)   with E1/Em1 # All the same condition:
#   or z == 1 (mod 2) and w - z^2 == 4   (mod 8)   with E1/Em1 # w^2 - z^4 == 8 (mod 16)
#   or z == 1 (mod 2) and w + z^2 == 4   (mod 8)   with E2/Em2 # ""
#   or z == 1 (mod 2) and w + z^2 == 6   (mod 8)   with E2/Em2 # ""
# 8 if z == 0 (mod 1) and w + z^2 == 1   (mod 2)   with E1/E2/Em1/Em2

### Setup
from modular_method.number_fields.field_constructors import field_with_root
from modular_method.diophantine_equations.conditions import ConditionalValue
from modular_method.padics.pAdic_base import pAdicBase
from modular_method.diophantine_equations.conditions import TreeCondition

def Frey_curve_of_divisibility_sequence_1(R, C):
    z, w = R.gens()
    C0_1 = CongruenceCondition(z - 1, 4) & CongruenceCondition(w - z^2 - 128, 256)
    C0_m2 = CongruenceCondition(z - 1, 4) & CongruenceCondition(w + z^2 - 128, 256)    
    C0_m1 = CongruenceCondition(z - 3, 4) & CongruenceCondition(w - z^2 - 128, 256)
    C0_2 = CongruenceCondition(z - 3, 4) & CongruenceCondition(w + z^2 - 128, 256)
    C1_1 = CongruenceCondition(z - 1, 4) & CongruenceCondition(w - z^2, 256)
    C1_m2 = CongruenceCondition(z - 1, 4) & CongruenceCondition(w + z^2, 256)    
    C1_m1 = CongruenceCondition(z - 3, 4) & CongruenceCondition(w - z^2, 256)
    C1_2 = CongruenceCondition(z - 3, 4) & CongruenceCondition(w + z^2, 256)
    C7_1 = CongruenceCondition(z - 1, 2) & CongruenceCondition(w^2 - z^4 - 8, 16)
    C8_1 = CongruenceCondition(w + z^2 - 1, 2)
    pAdics = pAdicBase(QQ, 2)
    result = []
    C_1 = (C0_1 | C1_1 | C7_1 | C8_1) & C
    E1 = FreyCurve([0, 4*z, 0, 2*(z^2 - w), 0], condition=C_1)
    N2_1 = []
    T0_1 = TreeCondition((C0_1 & C).pAdic_tree(pAdics=pAdics))
    if (not T0_1.never()): N2_1.append((0, T0_1))
    T1_1 = TreeCondition((C1_1 & C).pAdic_tree(pAdics=pAdics))
    if (not T1_1.never()): N2_1.append((1, T1_1))
    T7_1 = TreeCondition((C7_1 & C).pAdic_tree(pAdics=pAdics))
    if (not T7_1.never()): N2_1.append((7, T7_1))
    T8_1 = TreeCondition((C8_1 & C).pAdic_tree(pAdics=pAdics))
    if (not T8_1.never()): N2_1.append((7, T8_1))
    if len(N2_1) > 0:
        E1.conductor_exponent.set_cache(ConditionalValue(N2_1), 2)
        result.append((E1, C_1))
    C_m1 = (C0_m1 | C1_m1) & C
    Em1 = FreyCurve([0, -4*z, 0, 2*(z^2 - w), 0], condition=C_m1)
    N2_m1 = []
    T0_m1 = TreeCondition((C0_m1 & C).pAdic_tree(pAdics=pAdics))
    if (not T0_m1.never()): N2_m1.append((0, T0_m1))
    T1_m1 = TreeCondition((C1_m1 & C).pAdic_tree(pAdics=pAdics))
    if (not T1_m1.never()): N2_m1.append((1, T1_m1))
    if len(N2_m1) > 0:
        Em1.conductor_exponent.set_cache(ConditionalValue(N2_m1), 2)
        result.append((Em1, C_m1))
    C_2 = (C0_2 | C1_2) & C
    E2 = FreyCurve([0, -4*z, 0, 2*(z^2 + w), 0], condition=C)
    N2_2 = []
    T0_2 = TreeCondition((C0_2 & C).pAdic_tree(pAdics=pAdics))
    if (not T0_2.never()): N2_2.append((0, T0_2))
    T1_2 = TreeCondition((C1_2 & C).pAdic_tree(pAdics=pAdics))
    if (not T1_2.never()): N2_2.append((1, T1_2))
    if len(N2_2) > 0:
        E2.conductor_exponent.set_cache(ConditionalValue(N2_2), 2)
        result.append((E2, C_2))
    C_m2 = (C0_m2 | C1_m2) & C
    Em2 = FreyCurve([0, 4*z, 0, 2*(z^2 + w), 0], condition=C)
    N2_m2 = []
    T0_m2 = TreeCondition((C0_m2 & C).pAdic_tree(pAdics=pAdics))
    if (not T0_m2.never()): N2_m2.append((0, T0_m2))
    T1_m2 = TreeCondition((C1_m2 & C).pAdic_tree(pAdics=pAdics))
    if (not T1_m2.never()): N2_m2.append((1, T1_m2))
    if len(N2_m2) > 0:
        Em2.conductor_exponent.set_cache(ConditionalValue(N2_m2), 2)
        result.append((Em2, C_m2))
    return ConditionalValue(result)

def Frey_curve_of_divisibility_sequence(a, D, precision=1):
    R.<z, w, B> = QQ[]
    d = ZZ(D / a)
    C = ExistsCondition(a*z^4 + d*B^4 - w^2, [z, w], precision=precision)
    C2 = TreeCondition(ExistsCondition(a*z^4 + d*B^4 - w^2, [z, w],
                                       precision=precision+3).pAdic_tree(pAdics=pAdicBase(QQ, 2)))
    z, w = QQ[z, w].gens()
    C = CoprimeCondition([z, w]) & C & C2
    if a == 1:
        return Frey_curve_of_divisibility_sequence_1(z.parent(), C)
    K = field_with_root(QQ, a)
    sqrta = sqrt(K(a))
    invariants = [0, 4*sqrta*z, 0, 2*(a * z^2 + sqrta * w), 0]
    if K == QQ:
        return FreyCurve(invariants, condition=C)
    return FreyQcurve(invariants, condition=C, guessed_degrees=[2])

def possible_squarefree_parts(D, use_curve=True):
    if not use_curve:
        return [d for d in D.divisors() if d.is_squarefree()]
    E = EllipticCurve([0, 0, 0, D, 0])
    ls = [1]
    gens = E.gens()
    if len(gens) == 0:
        return []
    gens.extend([E(P) for P in E.torsion_subgroup().gens()])
    for P in gens:
        a = (D if P.xy()[0] == 0 else P.xy()[0].squarefree_part())
        if a not in ls:
            ls.append(a)
    n = len(ls)
    for i in range(n):
        for j in range(i+1, n):
            c = (ls[i]*ls[j]).squarefree_part()
            if c not in ls:
                ls.append(c)
    return ls

def possible_a(D, use_curve=True):
    if not use_curve:
        return [d for d in D.divisors()]
    E = EllipticCurve([0, 0, 0, D, 0])
    ls = [1]
    gens = E.gens()
    if len(gens) == 0:
        return [D]
    result = []
    gens.extend([E(P) for P in E.torsion_subgroup().gens() if 2.divides(P.order())])
    for d in mrange([2 for P in gens]):
        P = sum(d[i] * gens[i] for i in range(len(gens)))
        result.append(1 if P == E(0)
                      else (D if P == E((0, 0))
                            else sign(P.xy()[0].numerator())*gcd(P.xy()[0].numerator(), D)))
    return list(Set(result))

### helper functions:
def common(ls):
    result = 1
    for item in ls:
        result = lcm(result, item[-1])
    return result.prime_factors()

def minrepr(x):
    if x == 0:
        return x
    else:
        return product(x.prime_factors())

def last(ls):
    result = []
    for item in ls:
        result.append(minrepr(item[-1]))
    return result

def cm_test(ls, cap=200):
    primes = prime_range(cap)
    return [sum(nf[0].coefficient(p) == 0 for p in primes)/len(primes) for nf in ls]

from modular_method.padics.pAdic_base import pAdicBase
from modular_method.diophantine_equations.conditions import TreeCondition
def premake_tree(condition, primes, field=QQ, **kwds):
    result = None
    for p in primes:
        C = TreeCondition(condition.pAdic_tree(pAdics=pAdicBase(field, p),
                                               **kwds))
        if result is None:
            result = C
        else:
            result = result & C
    return result

### Examples without additional integral points
def is_good_candidate(D):
    if D == 0 or any(e >= 4 for p, e in ZZ(D).factor()):
        return False
    E = EllipticCurve([0, 0, 0, D, 0])
    return (E.rank(proof=False) > 0 and
            len(E.integral_points()) <= 1 + 2*(ZZ(-D).is_square()))

def check(N):
    for D in range(-N, N):
        if is_good_candidate(D):
            print(D, "=", ZZ(D).factor())

def advanced_check():
    D = 0
    while True:
        D = (-D if D > 0 else 1 - D)
        try:
            if is_good_candidate(D):
                print("Trying for D =", D, "=", D.factor())
                for a in possible_a(D):
                    if a.is_square():
                        continue
                    print("Case for a =", a)
                    Ea = Frey_curve_of_divisibility_sequence(a, D)
                    Ea = Ea.decomposable_twist()
                    Ka = Ea.decomposition_field()
                    for p in D.prime_factors():
                        if p.divides(6):
                            continue
                        if 4.divides(Ka(D).valuation(Ka.prime_above(p))):
                            print(p, "could appear good in newform level")
                        else:
                            print(p, "will appear quadratic in newform level")
        except Exception as e:
            if isinstance(e, KeyboardInterrupt):
                return

# check(30) # -23, 13, 21, 29

### Example for D = 3
from modular_method.diophantine_equations.conditions import apply_to_conditional_value
from modular_method.diophantine_equations.conditions import conditional_product
from modular_method.padics.pAdic_base import pAdicBase
D = 3
print("Possible a", possible_a(D)) # [1, 3]
E = EllipticCurve([0, 0, 0, D, 0])
# Rank 1 with one 2-torsion point
P0 = E.gens()[0]; T = E.torsion_points()[0]
# Case a = 1; corresponding to points m*P0 for m any integer
a = 1; d = ZZ(D / a)
E1 = Frey_curve_of_divisibility_sequence(a, D)
print("Additive primes for E1:", apply_to_conditional_value(lambda E1_ : E1_.primes_of_possible_additive_reduction(), E1)) # [2]
z, w = E1[0][0].parameters()
nfs1 = apply_to_conditional_value(lambda E1_ : E1_.newform_candidates(bad_primes=[2, 3], algorithm='magma'), E1)
# nfs1 has two cases, one of which has no newforms at all, so we replace it
print(len(nfs1) == 2 and len(nfs1[0][0]) == 0) # True
nfs1, C1 = nfs1[1]
# Only one of the curves corresponds to C1, so we replace it with only that one
print([not (C_ & C1).pAdic_tree(pAdics=pAdicBase(QQ, 2)).is_empty() for E1_, C_ in E1])
E1, C = E1[0]
C = premake_tree(C, primes=prime_range(5, 20), precision_cap=1, verbose=-1)
nfs1 = eliminate_by_traces(E1, nfs1, primes=prime_range(5, 20), condition=C)
# Newforms remain even ones without CM. Should choose a point to work this out.
# Case P = 2*P0:
P = 2*P0
print ("Prime factors in denominator of P:",
       P.xy()[0].denominator().prime_factors()) # [2]
print ("Prime factors in denominator of 2*P:",
       (2*P).xy()[0].denominator().prime_factors()) # [2, 7]
nfs1P = eliminate_by_trace(E1, nfs1, 7, condition=(C & CongruenceCondition(w^2 - a*z^4, 7)))
print("Remaining primes for P:")
print(common(nfs1P)) # [2, 3, 5, 7]
# Case P = 3*P0:
P = 3*P0
print ("Prime factors in denominator of P:",
       P.xy()[0].denominator().prime_factors()) # [3]
print ("Prime factors in denominator of 3*P:",
       (3*P).xy()[0].denominator().prime_factors()) # [3, 2521, 21961]
print ("Prime factors in denominator of 3^2*P:",
       (3^2*P).xy()[0].denominator().prime_factors()) # [3, 107, ...]
nfs1P = eliminate_by_trace(E1, nfs1, 107, condition=CoprimeCondition([z, w]) & CongruenceCondition(w^2 - a*z^4, 107))
print("Remaining primes for P:")
print(common(nfs1P)) # [2, 3, 5, 7, 13, 107]

# Case for a = 3; corresponding to points m*P0 + T for m any integer
a = 3; d = ZZ(D / a)
E3 = Frey_curve_of_divisibility_sequence(a, D)
E3 = E3.decomposable_twist()
print("Additive primes for E3", [P.smallest_integer() for P in
                                 E3.primes_of_possible_additive_reduction()]) # [2, 3]
z, w = E3.parameters()
z = z.change_ring(QQ); w = w.change_ring(QQ)
nfs3 = E3.newform_candidates(algorithm='magma')
C = premake_tree(E3._condition, prime_range(5, 20), precision_cap=1, verbose=-1)
nfs3 = eliminate_by_traces(E3, nfs3, condition=C, primes=prime_range(5, 20))
# Can not eliminate all newforms in general
# Some are not eliminated and also definitely non-cm
# If we choose a point we can solve the problem
# Case P = 3*P0 + T
P = 3*P0 + T
print("Prime factors in denominator of P:",
       P.xy()[0].denominator().prime_factors()) # [11]
nfs3P = eliminate_by_trace(E3, nfs3, 11, condition=(C & CongruenceCondition(w^2 - a*z^4, 11)))
print("Primes of no elimination for P:",
      apply_to_conditional_value(common, nfs3P)) # [2, 3, 11]

### Example for D = -2
from modular_method.diophantine_equations.conditions import apply_to_conditional_value
from modular_method.diophantine_equations.conditions import conditional_product
D = -2
print("Squarefree parts", possible_squarefree_parts(D)) # [1, -1, 2, -2]
E = EllipticCurve([0, 0, 0, D, 0])
# Rank 1 with one 2-torsion point
P0 = E.gens()[0]; T = E.torsion_points()[0]
# Case a = 1; corresponding to points m*P0 for m even
a = 1; d = ZZ(D / a)
E1 = Frey_curve_of_divisibility_sequence(a, D)
print("Additive primes for E1:", apply_to_conditional_value(lambda E1_ : E1_.primes_of_possible_additive_reduction(), E1)) # [2]
z, w = E1[0][0].parameters()
nfs1 = apply_to_conditional_value(lambda E1_ : E1_.newform_candidates(algorithm='magma'), E1)
C = premake_tree(E1[0][1]._left, primes=prime_range(3, 10), precision_cap=1, verbose=-1)
E1_ = ConditionalValue([(val, C & con._right) for val, con in E1])
Enfs1 = conditional_product(E1_, nfs1)
# Note that the group of these points is generated by P = 2*P0
P = 2*P0
print("Prime factors in denominator of P:",
      P.xy()[0].denominator().prime_factors()) # [2]
# Prime powers thus appear in multiples of 2*P
print("Prime factors in denominator of 2*P:",
      (2*P).xy()[0].denominator().prime_factors()) # [2, 3, 7]
# We start eliminating at 3 and 7
nfs1 = apply_to_conditional_value((lambda Enfs_, con : eliminate_by_traces(Enfs_[0], Enfs_[1], primes=[3, 7], condition=(con & CongruenceCondition(w^2 - a*z^4, 3*7)))), Enfs1, use_condition=True)
print("Primes of no elimination for a=1:",
      apply_to_conditional_value(common, nfs1)) # [2, 3]

# Case for a = -1; corresponding to points m*P0 for m odd
a = -1; d = ZZ(D / a)
Em1 = Frey_curve_of_divisibility_sequence(a, D)
Em1 = Em1.decomposable_twist()
print("Bad primes for Em1", [P.smallest_integer() for P in
                             Em1.primes_of_possible_additive_reduction()]) # [2]
z, w = Em1.parameters()
z = z.change_ring(QQ); w = w.change_ring(QQ)
nfsm1 = Em1.newform_candidates(algorithm='magma')
C = premake_tree(Em1._condition, prime_range(3, 15), precision_cap=1, verbose=-1)
nfsm1 = eliminate_by_traces(Em1, nfsm1, condition=C, primes=prime_range(3, 15))
# Can not eliminate all newforms in general
# Some are not eliminated and also definitely non-cm
# If we choose a point we can solve the problem
# Case P = 3*P0
P = 3*P0
print("Prime factors in denominator of P:",
       P.xy()[0].denominator().prime_factors()) # [13]
nfsm1P = eliminate_by_trace(Em1, nfsm1, 13, condition=(C & CongruenceCondition(w^2 - a*z^4, 13)))
print("Primes of no elimination for P:",
      apply_to_conditional_value(common, nfsm1P)) # [2, 5, 13]
# We can do better as all powers appear in multiples of 13*P
# hence in multiples of 13*P0, so we use primes there
print("Prime factors in denominator of 13*P0:",
      (13*P0).xy()[0].denominator().prime_factors()) # [37, 41, ...]
nfsm1P = eliminate_by_trace(Em1, nfsm1P, 37, condition=(CoprimeCondition([z, w]) & CongruenceCondition(w^2 - a*z^4, 37)))
print("Primes of no elimination for P:",
      apply_to_conditional_value(common, nfsm1P)) # [2]

# Case for a = 2; corresponding to points m*P0 + T for m odd
a = 2; d = ZZ(D / a)
E2 = Frey_curve_of_divisibility_sequence(a, D)
E2 = E2.decomposable_twist()
print("Bad primes for E2", [P.smallest_integer() for P in
                            E2.primes_of_possible_additive_reduction()]) # [2]
z, w = E2.parameters()
z = z.change_ring(QQ); w = w.change_ring(QQ)
nfs2 = E2.newform_candidates(algorithm='magma')
C = premake_tree(E2._condition, prime_range(3, 15), precision_cap=1, verbose=-1)
nfs2 = eliminate_by_traces(E2, nfs2, condition=C, primes=prime_range(3, 15))
# Can not eliminate all newforms in general
# Some are not eliminated and also definitely non-cm
# If we choose a point we can solve the problem
# Case P = 5*P0 + T
P = 5*P0 + T
print("Prime factors in denominator of P:",
       P.xy()[0].denominator().prime_factors()) # [17, 79]
nfs2P = eliminate_by_trace(E2, nfs2, 17, condition=(CoprimeCondition([z, w]) & CongruenceCondition(w^2 - a*z^4, 17)))
nfs2P = eliminate_by_trace(E2, nfs2P, 79, condition=(CoprimeCondition([z, w]) & CongruenceCondition(w^2 - a*z^4, 79))) # long time
print("Primes of no elimination for P:",
      apply_to_conditional_value(common, nfs2P)) # [2, 5]

# Case for a = -2; corresponding to points m*P0 + T for m even
a = -2; d = ZZ(D / a)
Em2 = Frey_curve_of_divisibility_sequence(a, D)
Em2 = Em2.decomposable_twist()
print("Bad primes for Em2", [P.smallest_integer() for P in
                            Em2.primes_of_possible_additive_reduction()]) # [2]
z, w = Em2.parameters()
z = z.change_ring(QQ); w = w.change_ring(QQ)
nfsm2 = Em2.newform_candidates(algorithm='magma')
C = premake_tree(Em2._condition, prime_range(3, 15), precision_cap=1, verbose=-1)
nfsm2 = eliminate_by_traces(Em2, nfsm2, condition=C, primes=prime_range(3, 15))
# Can not eliminate all newforms, but those that remain probably have CM
# If we choose a point we can solve the problem
# Case P = 2*P0 + T
P = 2*P0 + T
print("Prime factors in denominator of P:",
       P.xy()[0].denominator().prime_factors()) # [3]
nfsm2P = eliminate_by_trace(Em2, nfsm2, 3, condition=(CoprimeCondition([z, w]) & CongruenceCondition(w^2 - a*z^4, 3)))
print("Primes of no elimination for P:",
      apply_to_conditional_value(common, nfsm2P)) # [2, 3, 7]
# Also primes in 3*P must divide the prime power
print("Prime factors in denominator of 3*P:",
       (3*P).xy()[0].denominator().prime_factors()) # [3, 11, 577]
nfsm2P = eliminate_by_trace(Em2, nfsm2P, 11, condition=(CoprimeCondition([z, w]) & CongruenceCondition(w^2 - a*z^4, 11)))
print("Primes of no elimination for P:",
      apply_to_conditional_value(common, nfsm2P)) # [2, 3]

### Example for D = -17
from modular_method.diophantine_equations.conditions import apply_to_conditional_value
from modular_method.diophantine_equations.conditions import conditional_product
D = -17
print("Squarefree parts", possible_squarefree_parts(D)) # [1, -1, 17, -17]
E = EllipticCurve([0, 0, 0, D, 0])
# Rank 2 with one 2-torsion point
P0, Q0 = E.gens(); T = E.torsion_points()[0]
# Integral points are P0, Q0, T, Q0 + T, P0 - Q0
# Case a = 1; corresponding to points m*P0 + n*Q0 with m + n even
# This case includes the integral point P0 - Q0
a = 1; d = ZZ(D / a)
E1 = Frey_curve_of_divisibility_sequence(a, D)
print("Additive primes for E1:", apply_to_conditional_value(lambda E1_ : E1_.primes_of_possible_additive_reduction(), E1)) # [2]
z, w = E1[0][0].parameters()
nfs1 = apply_to_conditional_value(lambda E1_ : E1_.newform_candidates(bad_primes=[2, 17], algorithm='magma'), E1)
matches = [(i, j) for i in range(len(E1)) for j in range(len(nfs1))
           if not (E1[i][1] & nfs1[j][1]).pAdic_tree(pAdics=pAdicBase(QQ, 2)).is_empty()]
Enfs1 = ConditionalValue([((E1[i][0], nfs1[j][0]), E1[i][1] & nfs1[j][1]) for i, j in matches])
nfs1 = apply_to_conditional_value((lambda Enfs_, con : eliminate_by_traces(Enfs_[0], Enfs_[1], primes=prime_range(3, 15), condition=con)), Enfs1, use_condition=True)
print("Primes of no elimination for a=1:",
      apply_to_conditional_value(common, nfs1)) # [2, 3]

# Case for a = -1; corresponding to points m*P0 + n*Q0 for m + n odd
# This case includes the integral points P0 and Q0
a = -1; d = ZZ(D / a)
Em1 = Frey_curve_of_divisibility_sequence(a, D)
Em1 = Em1.decomposable_twist()
print("Bad primes for Em1", [P.smallest_integer() for P in
                             Em1.primes_of_possible_additive_reduction()]) # [2]
z, w = Em1.parameters()
z = z.change_ring(QQ); w = w.change_ring(QQ)
nfsm1 = Em1.newform_candidates(algorithm='magma')
C = premake_tree(Em1._condition, prime_range(3, 15), precision_cap=1, verbose=-1)
nfsm1 = eliminate_by_traces(Em1, nfsm1, condition=C, primes=prime_range(3, 15))
# Can not eliminate all newforms in general
# Some are not eliminated and also definitely non-cm
# If we choose a point we can solve the problem
# Case P = 3*P0
P = 3*P0
print("Prime factors in denominator of P:",
       P.xy()[0].denominator().prime_factors()) # [13]
nfsm1P = eliminate_by_trace(Em1, nfsm1, 13, condition=(C & CongruenceCondition(w^2 - a*z^4, 13)))
print("Primes of no elimination for P:",
      apply_to_conditional_value(common, nfsm1P)) # [2, 5, 13]
# We can do better as all powers appear in multiples of 13*P
# hence in multiples of 13*P0, so we use primes there
print("Prime factors in denominator of 13*P0:",
      (13*P0).xy()[0].denominator().prime_factors()) # [37, 41, ...]
nfsm1P = eliminate_by_trace(Em1, nfsm1P, 37, condition=(CoprimeCondition([z, w]) & CongruenceCondition(w^2 - a*z^4, 37)))
print("Primes of no elimination for P:",
      apply_to_conditional_value(common, nfsm1P)) # [2]

### Case for a = 17; corresponding to points m*P0 + n*Q0 + T for m + n odd
# This case includes the integral point Q0 + T
a = 17; d = ZZ(D / a)
E17 = Frey_curve_of_divisibility_sequence(a, D)
E17 = E17.decomposable_twist()
print("Bad primes for E17", [P.smallest_integer() for P in
                             E17.primes_of_possible_additive_reduction()]) # [2, 17]
z, w = E17.parameters()
z = z.change_ring(QQ); w = w.change_ring(QQ)
# Doing the computation for the conductor at 2 in a smart way and
# storing it beforehand
K17 = E17.definition_field()
P17, Q17 = K17.primes_above(2)
Csq = CongruenceCondition(w - z^2, 4)
Cmsq = CongruenceCondition(w + z^2, 4)
C0 = (~Csq) & (~Cmsq)
P2 = ConditionalValue([(P17, C0), (P17, Csq), (Q17, Cmsq)])
N17_2 = apply_to_conditional_value((lambda prime, con: E17.conductor_exponent(prime, condition=(E17._condition & con))), P2, use_condition=True) # Time consuming: about 1.5 hours
E17.conductor_exponent.set_cache(N17_2, P17)
E17.conductor_exponent.set_cache(N17_2, Q17)
# # Twisting the curve by -1 gives a lower conductor for
# #    z = 1 mod 2 and w = 23*z^2 + 24*z      mod 32
# #    z = 1 mod 2 and w =  9*z^2 + 24*z      mod 32
# #    z = 1 mod 2 and w =  7*z^2 + 28*z + 30 mod 32
# #    z = 1 mod 2 and w = 29*z^2 + 28*z + 30 mod 32
# # in which case conductor exponent 4 changes into 0, 1, 2, or 3
# E17_ = E17.twist(QQ(-1))
# K17_ = E17_.definition_field()
# P17_, Q17_ = K17_.primes_above(2)
# P2_ = ConditionalValue([(P17_, C0), (P17_, Csq), (Q17_, Cmsq)])
# N17_2_ = apply_to_conditional_value((lambda prime, con: E17_.conductor_exponent(prime, condition=(E17_._condition & con))), P2_, use_condition=True) # Time consuming: about 1.5 hours
# E17_.conductor_exponent.set_cache(N17_2, P17_)
# E17_.conductor_exponent.set_cache(N17_2, Q17_)
# Now for the real computations
N17 = E17.conductor(additive_primes=K17.primes_above(2*17))
levels = E17.newform_levels(bad_primes=K17.primes_above(2*17))
# We do the computation for all levels except 2^8 * 17^2
# which corresponds to w ~= z^2 (mod 4) and w ~= -z^2 (mod 4)
lvls = ConditionalValue([(levels[i][0], levels[i][1]) for i in range(2, len(levels))])
nfs17 = apply_to_conditional_value((lambda lvl, con:
                                    E17._newform_candidates(lvl,
                                                            con &
                                                            E17._condition,
                                                            "magma",
                                                            None,
                                                            False)),
                                   lvls, use_condition=True)
nfs17 = eliminate_by_traces(E17, nfs17, primes=prime_range(3, 15), verbose=True)
### Can not eliminate all newforms in general
# Some are not eliminated and also definitely non-cm
# If we choose a point we can solve the problem
# Case P = 2*P0 - Q0 + T
P = 2*P0 - Q0 + T
print("Prime factors in denominator of P:",
       P.xy()[0].denominator().prime_factors()) # [19]
nfs17P = eliminate_by_trace(E17, nfs17, 19, condition=(CoprimeCondition([z, w]) & CongruenceCondition(w^2 - a*z^4, 19)))
print("Primes of no elimination for P:",
      apply_to_conditional_value(common, nfs17P)) # [2, 3, 5, 7, 11, 13, 17, 19]

# Case for a = -17; corresponding to points m*P0 + n*Q0 + T for m + n even
# This case includes the integral point T
a = -17; d = ZZ(D / a)
Em17 = Frey_curve_of_divisibility_sequence(a, D, precision=1)
# The default decomposable twist function chooses a twist that adds
# additional primes above 7 to the discriminant, which we do not
# want. We can avoid them by choosing our own twist
K = Em17.decomposition_field()
P17 = K.prime_above(17)
gamma = (P17^2).gens_reduced()[0]
Em17 = Em17.twist(gamma)
print(Em17.does_decompose()) # True
print("Bad primes for Em17", [P.smallest_integer() for P in
                            Em17.primes_of_possible_additive_reduction()]) # [2, 7, 17]
z, w = Em17.parameters()
z = z.change_ring(QQ); w = w.change_ring(QQ)
nfsm17 = Em17.newform_candidates(algorithm='magma')
C = premake_tree(Em17._condition, prime_range(3, 15), precision_cap=1, verbose=-1)
nfsm17 = eliminate_by_traces(Em17, nfsm17, condition=C, primes=prime_range(3, 15))
# Can not eliminate all newforms, but those that remain probably have CM
# If we choose a point we can solve the problem
# Case P = 2*P0 + T
P = 2*P0 + T
print("Prime factors in denominator of P:",
       P.xy()[0].denominator().prime_factors()) # [3]
nfsm17P = eliminate_by_trace(Em17, nfsm17, 3, condition=(CoprimeCondition([z, w]) & CongruenceCondition(w^2 - a*z^4, 3)))
print("Primes of no elimination for P:",
      apply_to_conditional_value(common, nfsm17P)) # [2, 3, 7]
# Also primes in 3*P must divide the prime power
print("Prime factors in denominator of 3*P:",
       (3*P).xy()[0].denominator().prime_factors()) # [3, 11, 577]
nfsm17P = eliminate_by_trace(Em17, nfsm17P, 11, condition=(CoprimeCondition([z, w]) & CongruenceCondition(w^2 - a*z^4, 11)))
print("Primes of no elimination for P:",
      apply_to_conditional_value(common, nfsm17P)) # [2, 3]

### Example for D = 13
from modular_method.diophantine_equations.conditions import apply_to_conditional_value
from modular_method.diophantine_equations.conditions import conditional_product
from modular_method.padics.pAdic_base import pAdicBase
D = 13
print("Possible a", possible_a(D)) # [1, 13]
E = EllipticCurve([0, 0, 0, D, 0])
# Rank 1 with one 2-torsion point
P0 = E.gens()[0]; T = E.torsion_points()[0]
# Case a = 1; corresponding to points m*P0 for m any integer
a = 1; d = ZZ(D / a)
E1 = Frey_curve_of_divisibility_sequence(a, D)
print("Additive primes for E1:", apply_to_conditional_value(lambda E1_ : E1_.primes_of_possible_additive_reduction(), E1)) # [2]
z, w = E1[0][0].parameters()
nfs1 = apply_to_conditional_value(lambda E1_ : E1_.newform_candidates(bad_primes=[2, 13], algorithm='magma'), E1)
# nfs1 has five cases, one of which has no newforms at all, the others
# correspond to one of the E1's each, so we replace it accordingly.
print(len(nfs1) == 5 and len(nfs1[0][0]) == 0) # True
nfs1 = ConditionalValue([nfs1[i] for i in range(1, len(nfs1))])
# The order of the newforms is the same as the elliptic curves
print([[not (C_ & C__).pAdic_tree(pAdics=pAdicBase(QQ, 2)).is_empty() for E1_, C_ in E1] for nfs1_, C__ in nfs1])
Enfs1 = ConditionalValue([((E1[i][0], nfs1[i][0]), E1[i][1] & nfs1[i][1]) for i in range(4)])
nfs1 = apply_to_conditional_value((lambda Enfs1_, C_ : eliminate_by_trace(Enfs1_[0], Enfs1_[1], prime=3, condition=C_)), Enfs1, use_condition=True)
print("Primes that can't be eliminated:", apply_to_conditional_value(common, nfs1)) # [3, 5, 7]
# Using a point we can do even better!
P = P0
print("Prime factors in denominator of P:",
      P.xy()[0].denominator().prime_factors()) # [2]
# Since the power is at least a cube
# we may look at points that are multiples of 4*P
print("Prime factors in denominator of 4*P:",
      (4*P).xy()[0].denominator().prime_factors()) # [2, 3, 7, 17, 127, 21559]
# We note that the additional elimination at 7 kills 5
Enfs1P = ConditionalValue([((E1[i][0], nfs1[i][0]), E1[i][1] & nfs1[i][1]) for i in range(4)])
nfs1P = apply_to_conditional_value((lambda Enfs1_ : eliminate_by_trace(Enfs1_[0], Enfs1_[1], prime=7, condition=(CoprimeCondition([z, w]) & CongruenceCondition(w^2 - z^4, 7)))), Enfs1P)
print("Primes that can't be eliminated:", apply_to_conditional_value(common, nfs1P)) # [3, 5, 7]
# Since 3 appeared, we can also look at primes in 6*P
print("Prime factors in denominator of 4*P:",
      (6*P).xy()[0].denominator().prime_factors()) # [2, 3, 7, 17, 127, 21559]
# 11 gets rid of some additional primes
Enfs1P = ConditionalValue([((E1[i][0], nfs1P[i][0]), E1[i][1] & nfs1P[i][1]) for i in range(4)])
nfs1P = apply_to_conditional_value((lambda Enfs1_ : eliminate_by_trace(Enfs1_[0], Enfs1_[1], prime=11, condition=(CoprimeCondition([z, w]) & CongruenceCondition(w^2 - z^4, 11)))), Enfs1P)
print("Primes for each newform:", apply_to_conditional_value(last, nfs1P)) # [3, 7]
# The remaining newforms are rational and the corresponding elliptic curves
# have 3-torsion and 7-torsion respectively.
print("Torsion order and remaining primes per newform:",
      apply_to_conditional_value((lambda nfs : [(nf[0]._f.EllipticCurve().sage().torsion_order(), nf[1]) for nf in nfs]),
                                 nfs1P))

# Case for a = 13; corresponding to points m*P0 + T for m any integer
a = 13; d = ZZ(D / a)
E13 = Frey_curve_of_divisibility_sequence(a, D)
E13 = E13.decomposable_twist()
print("Additive primes for E13", [P.smallest_integer() for P in
                                 E13.primes_of_possible_additive_reduction()]) # [2, 13]
z, w = E13.parameters()
z = z.change_ring(QQ); w = w.change_ring(QQ)
# Add the condition that neither z nor w may be divisible by 13
E13._condition = E13._condition & ~CongruenceCondition(z*w, 13)
print("Newform levels for E13:", E13.newform_levels()) # [(3328, 3328, 43264, 43264), (43264, 43264, 3328, 3328)]
print("Corresponding characters:", [eps^(-1) for eps in E13.splitting_character('conjugacy')])
nfs13 = E13.newform_candidates(algorithm='magma')
nfs13 = eliminate_by_traces(E13, nfs13, condition=CoprimeCondition([z, w]), primes=[3, 5, 7, 11])
# The first 8 newforms have CM by sqrt(-1) over which they are also defined
nfs13bad = nfs13[:8]
print(all((nf[0].coefficient_field().degree() == 2 and
           nf[0].coefficient_field().discriminant() == -4)
          for nf in nfs13bad))
eps = DirichletGroup(4).gens()[0]
bound = CuspForms((E13.splitting_character()^(-1)).extend(3328)).sturm_bound()
%time have_cm = all(nf[0].coefficient(p) == QQ(eps(p)) * nf[0].coefficient(p) for p in prime_range(bound) for nf in nfs13bad)
print(have_cm)
# We try some more elimination for the good ones using Kraus at 7, 5 and 3
# nfs13good = nfs13[8:]
# orig = last(nfs13good)
# %time nfs13good = kraus_method(E13, nfs13good, 7, w^2 - 13*z^4, primes=100, condition=CoprimeCondition([z, w]), verbose=2)
# %time nfs13good = kraus_method(E13, nfs13good, 5, w^2 - 13*z^4, primes=100, condition=CoprimeCondition([z, w]), verbose=2)
# %time nfs13good = kraus_method(E13, nfs13good, 3, w^2 - 13*z^4, primes=100, condition=CoprimeCondition([z, w]), verbose=2)
# print(common(nfs13good))
# We ran this, it eliminated some newforms, but by far not all for 2, 3, 5, and 7

# The bad cm curves arise from
#  z = 0, w = 1 (nfs13[4], nfs13[5], nfs13[6], nfs13[7])
# and
#  z = 0, w = -13 (nfs13[0], nfs13[1], nfs13[2], nfs13[3])

### Example for D = -23
from modular_method.diophantine_equations.conditions import apply_to_conditional_value
from modular_method.diophantine_equations.conditions import conditional_product
from modular_method.padics.pAdic_base import pAdicBase
D = -23
print("Possible a", possible_a(D)) # [1, -23]
E = EllipticCurve([0, 0, 0, D, 0])
# Rank 1 with one 2-torsion point
P0 = E.gens()[0]; T = E.torsion_points()[0]
# Case a = 1; corresponding to points m*P0 for m any integer
a = 1; d = ZZ(D / a)
E1 = Frey_curve_of_divisibility_sequence(a, D)
print("Additive primes for E1:", apply_to_conditional_value(lambda E1_ : E1_.primes_of_possible_additive_reduction(), E1)) # [2]
z, w = E1[0][0].parameters()
Enfs1 = apply_to_conditional_value(lambda E1_ : apply_to_conditional_value(lambda nfs : (E1_, nfs), E1_.newform_candidates(bad_primes=[2, 23], algorithm='magma')), E1)
# nfs1 has nine cases
Enfs1 = apply_to_conditional_value((lambda Enfs1_, C_ : (Enfs1_[0], eliminate_by_traces(Enfs1_[0], Enfs1_[1], primes=prime_range(3, 8), condition=Enfs1_[0]._condition))), Enfs1, use_condition=True)
# We start using a point to refine our search
P = P0
print("Prime factors in denominator of P:",
      P.xy()[0].denominator().prime_factors()) # [5]
Enfs1P = apply_to_conditional_value((lambda Enfs1_, C_ : (Enfs1_[0], eliminate_by_trace(Enfs1_[0], Enfs1_[1], prime=5, condition=(CoprimeCondition([z, w]) & CongruenceCondition(w^2 - a*z^4, 5))))), Enfs1, use_condition=True)
print("Prime factors in denominator of 5*P:",
      (5*P).xy()[0].denominator().prime_factors()) # [5, 17, 61, ...]
Enfs1P = apply_to_conditional_value((lambda Enfs1_, C_ : (Enfs1_[0], eliminate_by_trace(Enfs1_[0], Enfs1_[1], prime=17, condition=(CoprimeCondition([z, w]) & CongruenceCondition(w^2 - a*z^4, 17))))), Enfs1P, use_condition=True)
nfs1P = apply_to_conditional_value(lambda x : x[1], Enfs1P)
print("Primes that can't be eliminated:", apply_to_conditional_value(common, nfs1P)) # [2, 5, 11]

# Case for a = -23; corresponding to points m*P0 + T for m any integer
a = -23; d = ZZ(D / a)
Em23 = Frey_curve_of_divisibility_sequence(a, D)
Em23 = Em23.decomposable_twist()
print("Additive primes for Em23", [P.smallest_integer() for P in
                                   Em23.primes_of_possible_additive_reduction()]) # [2, 23]
z, w = Em23.parameters()
z = z.change_ring(QQ); w = w.change_ring(QQ)
# Since -23 = 1 mod 8 we need to do the computation of the conductor exponent
# at 2 in a smart way (or else run out of memory)
Km23 = Em23.definition_field()
Pm23, Qm23 = Km23.primes_above(2)
Csq = CongruenceCondition(w - z^2, 4)
Cmsq = CongruenceCondition(w + z^2, 4)
C0 = (~Csq) & (~Cmsq)
P2 = ConditionalValue([(Pm23, C0), (Pm23, Cmsq), (Qm23, Csq)])
Nm23_2 = apply_to_conditional_value((lambda prime, con: Em23.conductor_exponent(prime, condition=(Em23._condition & con))), P2, use_condition=True) # Time consuming: about 1.5 hours
Em23.conductor_exponent.set_cache(Nm23_2, Pm23)
Em23.conductor_exponent.set_cache(Nm23_2, Qm23)
print("Newform levels for Em23:", Em23.newform_levels()) # All contain a 23^2
print("Corresponding characters:", [eps^(-1) for eps in Em23.splitting_character('conjugacy')])
lvls = ConditionalValue([(val, con)
                         for val, con
                         in [(val, TreeCondition(con.pAdic_tree(pAdics=pAdicBase(QQ, 2))))
                             for val, con in Em23.newform_levels()]
                         if not con.pAdic_tree().is_empty()])
# # Takes too long
# nfsm23 = apply_to_conditional_value((lambda lvl, con:
#                                      Em23._newform_candidates(lvl, con & Em23._condition,
#                                                               "magma", None, False)),
#                                     lvls, use_condition=True)
# nfsm23 = eliminate_by_traces(Em23, nfsm23, condition=CoprimeCondition([z, w]), primes=[3, 5, 7, 11])

### Example for D = 21
from modular_method.diophantine_equations.conditions import apply_to_conditional_value
from modular_method.diophantine_equations.conditions import conditional_product
from modular_method.padics.pAdic_base import pAdicBase
D = 21
print("Possible a", possible_a(D)) # [1, -23]
E = EllipticCurve([0, 0, 0, D, 0])
# Rank 1 with one 2-torsion point
P0 = E.gens()[0]; T = E.torsion_points()[0]
# Case a = 1; corresponding to points m*P0 for m any integer
a = 1; d = ZZ(D / a)
E1 = Frey_curve_of_divisibility_sequence(a, D)
print("Additive primes for E1:", apply_to_conditional_value(lambda E1_ : E1_.primes_of_possible_additive_reduction(), E1)) # [2]
z, w = E1[0][0].parameters()
Enfs1 = apply_to_conditional_value(lambda E1_ : apply_to_conditional_value(lambda nfs : (E1_, nfs), E1_.newform_candidates(bad_primes=[2, 3, 7], algorithm='magma')), E1)
# Enfs1 has 8 cases
# We start using a point to refine our search
P = P0
print("Prime factors in denominator of P:",
      P.xy()[0].denominator().prime_factors()) # [2]
print("Prime factors in denominator of 2*P:",
      (2*P).xy()[0].denominator().prime_factors()) # [2, 5, 31]
Enfs1P = apply_to_conditional_value((lambda Enfs1_, C_ : (Enfs1_[0], eliminate_by_trace(Enfs1_[0], Enfs1_[1], prime=5, condition=(CoprimeCondition([z, w]) & CongruenceCondition(w^2 - a*z^4, 5))))), Enfs1, use_condition=True)
nfs1P = apply_to_conditional_value(lambda x : x[1], Enfs1P)
print("Primes that can't be eliminated:", apply_to_conditional_value(common, nfs1P)) # [2, 5]

# Case for a = 21; corresponding to points m*P0 + T for m any integer
a = 21; d = ZZ(D / a)
E21 = Frey_curve_of_divisibility_sequence(a, D)
E21 = E21.decomposable_twist()
print("Additive primes for E21", [P.smallest_integer() for P in
                                  E21.primes_of_possible_additive_reduction()]) # [2, 3, 7]
z, w = E21.parameters()
z = z.change_ring(QQ); w = w.change_ring(QQ)
E21._condition = E21._condition & (~CongruenceCondition(w, 3)) & (~CongruenceCondition(w, 7))
print("Newform levels for E21:", E21.newform_levels()) # 2^8 * 3^2 * 7^2
print("Corresponding characters:", [eps^(-1) for eps in E21.splitting_character('conjugacy')])
# # Takes too long
# nfsm23 = apply_to_conditional_value((lambda lvl, con:
#                                      Em23._newform_candidates(lvl, con & Em23._condition,
#                                                               "magma", None, False)),
#                                     lvls, use_condition=True)
# nfsm23 = eliminate_by_traces(Em23, nfsm23, condition=CoprimeCondition([z, w]), primes=[3, 5, 7, 11])

### Example for D = 125
from modular_method.diophantine_equations.conditions import apply_to_conditional_value
from modular_method.diophantine_equations.conditions import conditional_product
from modular_method.padics.pAdic_base import pAdicBase
D = 125
print("Possible a", possible_a(D)) # [1, 125]
E = EllipticCurve([0, 0, 0, D, 0])
# Rank 1 with one 2-torsion point
P0 = E.gens()[0]; T = E([0, 0])
# Case for a = 1; corresponding to points m*P0 for m any integer
a = 1; d = ZZ(D / a)
E1 = Frey_curve_of_divisibility_sequence(a, D)
print("Additive primes for E1:", apply_to_conditional_value(lambda E1_ : E1_.primes_of_possible_additive_reduction(), E1)) # [2]
z, w = E1[0][0].parameters()
Enfs1 = apply_to_conditional_value(lambda E1_ : apply_to_conditional_value(lambda nfs : (E1_, nfs), E1_.newform_candidates(bad_primes=[2, 5], algorithm='magma')), E1)
nfs1 = apply_to_conditional_value(lambda Enfs_ : Enfs_[1], Enfs1)
print(nfs1)
# No newforms at all

# Case for a = 125; corresponding to points m*P0 + T for m any integer
a = 125; d = ZZ(D / a)
E125 = Frey_curve_of_divisibility_sequence(a, D)
z, w = E125.parameters()
z = z.change_ring(QQ); w = w.change_ring(QQ)
E125._condition = E125._condition & ~CongruenceCondition(w, 5)
E125 = E125.decomposable_twist()
K125 = E125.definition_field()
nfs125 = E125.newform_candidates(bad_primes=K125.primes_above(10), algorithm='magma')
nfs125 = eliminate_by_traces(E125, nfs125, primes=[3, 7, 11], verbose=True)
# Choose a point to do better
P = P0 + T
print("Prime factors in denominator of P:",
      P.xy()[0].denominator().prime_factors()) # [11]
nfs125P = eliminate_by_trace(E125, nfs125, prime=11, condition=(CoprimeCondition([z, w]) & CongruenceCondition(w^2 - a*z^4, 11)))
print("Primes of no elimination for P:", common(nfs125P))

### Example for D = 31
from modular_method.diophantine_equations.conditions import apply_to_conditional_value
from modular_method.diophantine_equations.conditions import conditional_product
from modular_method.padics.pAdic_base import pAdicBase
D = 31
print("Possible a", possible_a(D)) # [1, 31]
E = EllipticCurve([0, 0, 0, D, 0])
# Rank 1 with one 2-torsion point
P0 = E.gens()[0]; T = E([0, 0])
# Case for a = 1; corresponding to points m*P0 for m any integer
a = 1; d = ZZ(D / a)
E1 = Frey_curve_of_divisibility_sequence(a, D)
print("Additive primes for E1:", apply_to_conditional_value(lambda E1_ : E1_.primes_of_possible_additive_reduction(), E1)) # [2]
z, w = E1[0][0].parameters()
C31 = CoprimeCondition([z, w]) & CongruenceCondition(w^2 - a*z^4, 31)
for E1_, con in E1:
    val = E1_.conductor_exponent(31, condition=C31)
    E1_.conductor_exponent.set_cache(val, 31)
Enfs1 = apply_to_conditional_value(lambda E1_ : apply_to_conditional_value(lambda nfs : (E1_, nfs), E1_.newform_candidates(bad_primes=[2, 31], algorithm='magma')), E1)
Enfs1 = apply_to_conditional_value((lambda Enfs1_: (Enfs1_[0], eliminate_by_traces(Enfs1_[0], Enfs1_[1], primes=[3, 5, 7, 11]))), Enfs1)
# We use the information of the point
P = P0
print("Prime factors in denominator of P:",
      P.xy()[0].denominator().prime_factors()) # [3]
Enfs1P = apply_to_conditional_value((lambda Enfs1_ : (Enfs1_[0], eliminate_by_trace(Enfs1_[0], Enfs1_[1], prime=3, condition=(CoprimeCondition([z, w]) & CongruenceCondition(w^2 - a*z^4, 3))))), Enfs1)
# And then denominators of 3*P0
P = 3*P0
print("Prime factors in denominator of P:",
      P.xy()[0].denominator().prime_factors()) # [3, 11, 71, 457]
Enfs1P = apply_to_conditional_value((lambda Enfs1_ : (Enfs1_[0], eliminate_by_trace(Enfs1_[0], Enfs1_[1], prime=11, condition=(CoprimeCondition([z, w]) & CongruenceCondition(w^2 - a*z^4, 11))))), Enfs1P)
Enfs1P = apply_to_conditional_value((lambda Enfs1_ : (Enfs1_[0], eliminate_by_trace(Enfs1_[0], Enfs1_[1], prime=71, condition=(CoprimeCondition([z, w]) & CongruenceCondition(w^2 - a*z^4, 71))))), Enfs1P)
print("Primes of no elimination for P:", apply_to_conditional_value(lambda Enfs: common(Enfs[1]), Enfs1P))

# Case for a = 31; corresponding to points m*P0 + T for m any integer
a = 31; d = ZZ(D / a)
E31 = Frey_curve_of_divisibility_sequence(a, D)
z, w = E31.parameters()
z = z.change_ring(QQ); w = w.change_ring(QQ)
E31 = E31.decomposable_twist()
K31 = E31.definition_field()
C31 = CoprimeCondition([z, w]) & ~CongruenceCondition(w, 31)
e31 = E31.conductor_exponent(K31.prime_above(31), condition=C31)
for P in K31.primes_above(31):
    E31.conductor_exponent.set_cache(e31, P)
nfs31 = E31.newform_candidates(bad_primes=K31.primes_above(2*31), algorithm='magma')
nfs31 = eliminate_by_traces(E31, nfs31, primes=[3, 7, 11], verbose=True)
# Choose a point to do better
# P = P0 + T
# print("Prime factors in denominator of P:",
#       P.xy()[0].denominator().prime_factors()) # [5]
# nfs31P = eliminate_by_trace(E125, nfs125, prime=11, condition=(CoprimeCondition([z, w]) & CongruenceCondition(w^2 - a*z^4, 11)))
# print("Primes of no elimination for P:", common(nfs125P))

### Other cases:
# D = -169, highest level 2^8*13^2, with character

### Necessary imports
from modular_method.diophantine_equations.conditions import CoprimeCondition
from modular_method.elliptic_curves.tates_algorithm import tates_algorithm_multiple
from modular_method.padics.pAdic_base import pAdicBase

### Computing conductors
R.<z, w> = QQ[]
C = CoprimeCondition([z, w])
pAdics = pAdicBase(QQ, 2)
N = {}
for gamma in [1, -1, 2, -2]:
    E1 = EllipticCurve([0, 4*z*gamma, 0, 2*(z^2 + w)*gamma^2, 0])
    E2 = EllipticCurve([0, -8*z*gamma, 0, 8*(z^2 - w)*gamma^2, 0])
    N[gamma] = tates_algorithm_multiple(
        (E1, E2),
        pAdics=(pAdics, pAdics),
        initial_values=C.pAdic_tree(pAdics=pAdics),
        only_calculate=['conductor'],
        verbose=True,
        precision_cap=25,
    )

### Functions for verification
def elements_modulo(N, con):
    """Iterator over the z, w modulo `N` that satisfy `con`"""
    ls, N_ = con.pAdic_tree().give_as_congruence_condition()
    N_ = N_.gens_reduced()[0]
    dif = ceil(N / N_)
    for z, w in ls:
        for zdif in range(dif):
            for wdif in range(dif):
                yield z + zdif*N_, w + wdif*N_

def satisfy_congruences(congruences, con):
    """Give `True` iff any of the z, w satisfying `con` also satisfy the
    `congruences`

    The argument `congruences` is a list of tuples. Each tuple
    contains a polynomial $f$ and an integer $N$, and represents the
    congruence $f(z, w) \equiv 0$ (mod $N$).

    """
    maxN = lcm(N for poly, N in congruences)
    return any(
        all(mod(poly(z, w), N) == 0 for poly, N in congruences)
        for z, w in elements_modulo(maxN, con)
    )

def possible_values(conditional_value, congruences):
    """Give the values of the `conditional_value` for which the
    `congruences` are satisfied

    """
    return [val for val, con in conditional_value if satisfy_congruences(congruences, con)]

def possible_N(congruences):
    """Give a dictionary of possible values satisfying the `congruences`
    for each twist

    """
    return {key: possible_values(value, congruences) for key, value in N.items()}
    
### Verifying the table
assert ({1: [8], -1: [8], 2: [8], -2: [8]} == possible_N([(w^2 - z^4 - 1, 2)]) and
        {1: [7], -1: [7], 2: [7], -2: [7]} == possible_N([(w^2 - z^4 - 8, 16)]) and
        {1: [4], -1: [3], 2: [6], -2: [6]} == possible_N([(w + z^2 - 8, 32), (z - 1, 4)]) and
        {1: [3], -1: [4], 2: [6], -2: [6]} == possible_N([(w + z^2 - 8, 32), (z - 3, 4)]) and
        {1: [2], -1: [4], 2: [6], -2: [6]} == possible_N([(w + z^2 - 24, 32), (z - 1, 4)]) and
        {1: [4], -1: [2], 2: [6], -2: [6]} == possible_N([(w + z^2 - 24, 32), (z - 3, 4)]) and
        {1: [6], -1: [6], 2: [4], -2: [2]} == possible_N([(w - z^2 - 8, 32), (z - 1, 4)]) and
        {1: [6], -1: [6], 2: [2], -2: [4]} == possible_N([(w - z^2 - 8, 32), (z - 3, 4)]) and
        {1: [6], -1: [6], 2: [3], -2: [4]} == possible_N([(w - z^2 - 24, 32), (z - 1, 4)]) and
        {1: [6], -1: [6], 2: [4], -2: [3]} == possible_N([(w - z^2 - 24, 32), (z - 3, 4)]) and
        {1: [5], -1: [5], 2: [6], -2: [6]} == possible_N([(w + z^2 - 2^4, 2^5)]) and
        {1: [6], -1: [6], 2: [5], -2: [5]} == possible_N([(w - z^2 - 2^4, 2^5)]) and
        {1: [3], -1: [4], 2: [6], -2: [6]} == possible_N([(w + z^2 - 2^5, 2^6), (z - 1, 4)]) and
        {1: [4], -1: [3], 2: [6], -2: [6]} == possible_N([(w + z^2 - 2^5, 2^6), (z - 3, 4)]) and
        {1: [6], -1: [6], 2: [4], -2: [3]} == possible_N([(w - z^2 - 2^5, 2^6), (z - 1, 4)]) and
        {1: [6], -1: [6], 2: [3], -2: [4]} == possible_N([(w - z^2 - 2^5, 2^6), (z - 3, 4)]) and
        {1: [3], -1: [4], 2: [6], -2: [6]} == possible_N([(w + z^2 - 2^6, 2^7), (z - 1, 4)]) and
        {1: [4], -1: [3], 2: [6], -2: [6]} == possible_N([(w + z^2 - 2^6, 2^7), (z - 3, 4)]) and
        {1: [6], -1: [6], 2: [4], -2: [3]} == possible_N([(w - z^2 - 2^6, 2^7), (z - 1, 4)]) and
        {1: [6], -1: [6], 2: [3], -2: [4]} == possible_N([(w - z^2 - 2^6, 2^7), (z - 3, 4)]) and
        {1: [0], -1: [4], 2: [6], -2: [6]} == possible_N([(w + z^2 - 2^7, 2^8), (z - 1, 4)]) and
        {1: [4], -1: [0], 2: [6], -2: [6]} == possible_N([(w + z^2 - 2^7, 2^8), (z - 3, 4)]) and
        {1: [6], -1: [6], 2: [4], -2: [0]} == possible_N([(w - z^2 - 2^7, 2^8), (z - 1, 4)]) and
        {1: [6], -1: [6], 2: [0], -2: [4]} == possible_N([(w - z^2 - 2^7, 2^8), (z - 3, 4)]) and
        {1: [1], -1: [4], 2: [6], -2: [6]} == possible_N([(w + z^2, 2^8), (z - 1, 4)]) and
        {1: [4], -1: [1], 2: [6], -2: [6]} == possible_N([(w + z^2, 2^8), (z - 3, 4)]) and
        {1: [6], -1: [6], 2: [4], -2: [1]} == possible_N([(w - z^2, 2^8), (z - 1, 4)]) and
        {1: [6], -1: [6], 2: [1], -2: [4]} == possible_N([(w - z^2, 2^8), (z - 3, 4)]))

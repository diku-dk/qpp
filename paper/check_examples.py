"""Reproduce the paper's small algebraic examples using exact arithmetic.

These are checks of the paper, not validation of a QPP exact backend.
The script uses only the Python standard library.
"""
from fractions import Fraction as F
from itertools import permutations
from math import gcd, isqrt, lcm


def transpose(a):
    return list(map(list, zip(*a)))


def matmul(a, b):
    return [[sum(x * y for x, y in zip(row, col))
             for col in transpose(b)] for row in a]


def determinant(a):
    n = len(a)
    total = 0
    for p in permutations(range(n)):
        inversions = sum(p[i] > p[j] for i in range(n) for j in range(i + 1, n))
        term = (-1) ** inversions
        for i in range(n):
            term *= a[i][p[i]]
        total += term
    return total


def inverse(a):
    n = len(a)
    rows = [[F(x) for x in row] + [F(i == j) for j in range(n)]
            for i, row in enumerate(a)]
    for j in range(n):
        p = next(i for i in range(j, n) if rows[i][j] != 0)
        rows[j], rows[p] = rows[p], rows[j]
        pivot = rows[j][j]
        rows[j] = [x / pivot for x in rows[j]]
        for i in range(n):
            if i != j:
                factor = rows[i][j]
                rows[i] = [x - factor * y for x, y in zip(rows[i], rows[j])]
    return [row[n:] for row in rows]


# Fixed field Q[x]/(x^8-x^4+1), used only to check the worked example.
# This small reference arithmetic is not the proposed Haskell scalar module.
class K24:
    def __init__(self, value=0):
        self.coeffs = (tuple(map(F, value)) if isinstance(value, (tuple, list))
                       else (F(value),) + (F(0),) * 7)
        assert len(self.coeffs) == 8

    @staticmethod
    def coerce(value):
        return value if isinstance(value, K24) else K24(value)

    def __add__(self, other):
        other = self.coerce(other)
        return K24([a + b for a, b in zip(self.coeffs, other.coeffs)])

    __radd__ = __add__

    def __neg__(self):
        return K24([-a for a in self.coeffs])

    def __sub__(self, other):
        return self + -self.coerce(other)

    def __rsub__(self, other):
        return self.coerce(other) + -self

    def __mul__(self, other):
        other = self.coerce(other)
        coeffs = [F(0)] * 15
        for i, a in enumerate(self.coeffs):
            for j, b in enumerate(other.coeffs):
                coeffs[i+j] += a*b
        for k in range(14, 7, -1):
            coeffs[k-4] += coeffs[k]
            coeffs[k-8] -= coeffs[k]
        return K24(coeffs[:8])

    __rmul__ = __mul__

    def __pow__(self, power):
        assert power >= 0
        result = K24(1)
        for _ in range(power):
            result = result * self
        return result

    def __eq__(self, other):
        return self.coeffs == self.coerce(other).coeffs

    def conjugate(self):
        return sum(a * zeta**((-j) % 24)
                   for j, a in enumerate(self.coeffs))


zeta = K24([0, 1, 0, 0, 0, 0, 0, 0])
assert zeta**8 - zeta**4 + 1 == 0
assert zeta**24 == 1 and zeta**12 == -1
h = F(1, 2) * (zeta**3 + zeta**21)
z = zeta**4
assert h == F(1, 2) * (zeta**3 - zeta**9)  # the K8 coordinates
assert h == K24([0, F(1, 2), 0, F(1, 2), 0, F(-1, 2), 0, 0])
assert h*h == F(1, 2) and h.conjugate() == h
assert h*(2*h) == 1  # the inverse example in Section 4
assert z + z.conjugate() == 1 and z*z.conjugate() == 1
bell = [[h*x for x in row]
        for row in [[1, 0, 1, 0], [0, 1, 0, 1],
                    [0, 1, 0, -1], [1, 0, -1, 0]]]
phase = [[K24(i == j) * [1, z, 1, z][i] for j in range(4)]
         for i in range(4)]
unitary = matmul(phase, bell)
assert unitary == [[h*x for x in row]
                   for row in [[1, 0, 1, 0], [0, z, 0, z],
                               [0, 1, 0, -1], [z, 0, -z, 0]]]
phased_bell = [row[0] for row in unitary]
assert phased_bell == [h, 0, 0, h*z]
assert sum(a.conjugate()*a for a in phased_bell) == 1
bell_dagger = [[a.conjugate() for a in row] for row in transpose(bell)]
observed = [row[0] for row in matmul(bell_dagger, [[a] for a in phased_bell])]
assert observed == [F(1, 2)*(1+z), 0, F(1, 2)*(1-z), 0]
assert [a.conjugate()*a for a in observed] == [F(3, 4), 0, F(1, 4), 0]
identity4 = [[K24(i == j) for j in range(4)] for i in range(4)]
unitary_dagger = [[a.conjugate() for a in row] for row in transpose(unitary)]
assert matmul(unitary_dagger, unitary) == identity4

# Complex tensor example: the row action uses transpose, not adjoint.
hh = [[h, h], [h, -h]]
tt = [[K24(1), K24(0)], [K24(0), zeta**3]]
values = [1, 2, 3, 5]
ht = [[hh[i][j]*tt[k][l] for j in range(2) for l in range(2)]
      for i in range(2) for k in range(2)]
assert [r[0] for r in matmul(ht, [[v] for v in values])] == [
    h*4, h*zeta**3*7, h*(-2), h*zeta**3*(-3)]

# Worked non-Pauli rotation, with independently written expected coordinates.
def rotation_coefficients(angle):
    angle = F(angle)
    order = 4 * angle.denominator
    assert 24 % order == 0
    w = zeta**((24 // order * angle.numerator) % 24)
    return F(1, 2)*(w + w.conjugate()), F(1, 2)*(w.conjugate() - w)


imaginary = zeta**6

# Linear combinations, inner products, and zero-qubit scalar evolution.
def inner(xs, ys):
    return sum(K24.coerce(x).conjugate() * y for x, y in zip(xs, ys))


v1, v2, v3 = [1, imaginary], [imaginary, 2], [z, h]
combination = [z*x + h*y for x, y in zip(v1, v2)]
assert inner(combination, v3) == z.conjugate()*inner(v1, v3) + h*inner(v2, v3)
scalar_state = 2 + imaginary
evolved_scalar = imaginary * scalar_state  # Phase (1/2), width zero
assert evolved_scalar == -1 + 2*imaginary
assert inner([evolved_scalar], [evolved_scalar]) == inner([scalar_state], [scalar_state]) == 5

rc, rd = rotation_coefficients(F(1, 3))
assert rc == K24([0, 0, 1, 0, 0, 0, F(-1, 2), 0])
assert rd == F(-1, 2)*imaginary
rh = [[rc + rd*h, rd*h], [rd*h, rc - rd*h]]
rotation_state = [rc + rd*h, rd*h]
assert rotation_state == [
    K24([0, F(1, 4), 1, F(-1, 4), 0, F(-1, 4), F(-1, 2), 0]),
    K24([0, F(1, 4), 0, F(-1, 4), 0, F(-1, 4), 0, 0])]
assert [r[0] for r in matmul(rh, [[1], [0]])] == rotation_state
assert [a.conjugate()*a for a in rotation_state] == [F(7, 8), F(1, 8)]
rh_dagger = [[a.conjugate() for a in row] for row in transpose(rh)]
assert matmul(rh_dagger, rh) == [[1, 0], [0, 1]]
# These expected values detect sign/half-angle errors that round trips can miss.
c1, d1 = rotation_coefficients(1)
assert [[c1+d1, K24(0)], [K24(0), c1-d1]] == [
    [-imaginary, 0], [0, imaginary]]  # R Z 1 = -i Z
c2, d2 = rotation_coefficients(2)
assert c2 + d2 == -1  # R (Id 0) 2

# Count coprime residues independently of the displayed totient product.
large_index = lcm(8, 4*101, 4*103)
large_degree = sum(gcd(k, large_index) == 1 for k in range(1, large_index+1))
assert large_index == 83224 and large_degree == 40800
assert 2**10 * large_degree == 41779200
assert gcd(101+103, 101*103) == 1  # the product of these primitive roots

# Squared singular values of M are roots of the characteristic polynomial.
m = [[1, 1, 0], [0, 1, 1], [0, 0, 2]]
g = matmul(transpose(m), m)
trace = sum(g[i][i] for i in range(3))
minor_sum = sum(g[i][i] * g[j][j] - g[i][j] * g[j][i]
                for i in range(3) for j in range(i + 1, 3))
assert (trace, minor_sum, determinant(g)) == (8, 15, 4)
assert sum(x*x for row in m for x in row) == 8
f = lambda x: x**3 - 8*x*x + 15*x - 4
assert all(f(x) != 0 for x in (-4, -2, -1, 1, 2, 4))
a, b, c = -8, 15, -4
disc = a*a*b*b - 4*b**3 - 4*a**3*c - 27*c*c + 18*a*b*c
assert disc == 916 and isqrt(disc)**2 != disc

# Rank-factorization formula with redundant columns and nonleading rows.
t = [[0, 0, 0, 0], [1, 2, 3, 4], [2, 1, 3, 5]]
cols, rows = [0, 1], [1, 2]
bmat = [[row[j] for j in cols] for row in t]
cmat = matmul(inverse([bmat[i] for i in rows]), [t[i] for i in rows])
assert matmul(bmat, cmat) == t

# beta^2 = 7 + 4 sqrt(2), so beta^4 - 14 beta^2 + 17 = 0.
# Represent Q(sqrt(2)) as rational pairs.
def qadd(x, y):
    return tuple(a+b for a, b in zip(x, y))


def qmul(x, y):
    a, b = x
    c, d = y
    return a*c + 2*b*d, a*d + b*c


radicand = (7, 4)
assert qadd(qadd(qmul(radicand, radicand), (-98, -56)), (17, 0)) == (0, 0)

# Prepare the SVD-obstruction state using H and controlled H, followed
# by a reversible computational-basis permutation. Amplitudes lie in Q(sqrt(2)).
def controlled_h(state, target, control=None):
    out = list(state)
    mask = 1 << (3 - target)
    for j in range(16):
        if j & mask or (control is not None and not j & (1 << (3 - control))):
            continue
        a, b = state[j], state[j | mask]
        out[j] = qmul(qadd(a, b), (0, F(1, 2)))
        out[j | mask] = qmul(qadd(a, tuple(-x for x in b)), (0, F(1, 2)))
    return out


state = [(0, 0)] * 16
state[0] = (1, 0)
state = controlled_h(state, 0)
state = controlled_h(state, 1, control=0)
state = controlled_h(state, 2, control=0)
specified = {0: 10, 8: 0, 10: 1, 12: 5, 14: 6}
remaining = iter(j for j in range(16) if j not in specified.values())
permutation = [specified[j] if j in specified else next(remaining) for j in range(16)]
assert sorted(permutation) == list(range(16))
permuted = [(0, 0)] * 16
for j, dest in enumerate(permutation):
    permuted[dest] = state[j]
padded = [row + [0] for row in m] + [[0] * 4]
expected = [qmul((x, 0), (0, F(1, 4))) for row in padded for x in row]
assert permuted == expected
norm = (0, 0)
for amp in permuted:
    norm = qadd(norm, qmul(amp, amp))
assert norm == (1, 0)

# Direct tensor action: rows receive B; columns then receive A.
aa = [[1, 2], [3, 4]]
bb = [[0, 1], [-1, 0]]
v = [1, 2, 3, 5]
tensor = [[aa[i][j]*bb[k][l] for j in range(2) for l in range(2)]
          for i in range(2) for k in range(2)]
view = [v[:2], v[2:]]
row_action = matmul(view, transpose(bb))
result = matmul(aa, row_action)
assert [x for row in result for x in row] == [x[0] for x in matmul(tensor, [[x] for x in v])]

print("Exact paper examples passed: cyclotomic coordinates, linear combinations and inner products, zero-qubit phases, field growth, phased Bell circuit and probabilities, non-Pauli rotation and expected-value cases, complex tensor action, SVD obstruction and state preparation, rank factorization, norm polynomial, tensor indexing.")

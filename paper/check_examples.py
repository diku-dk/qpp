"""Reproduce the paper's small algebraic examples using exact arithmetic.

These are checks of the paper, not validation of a QPP exact backend.
The script uses only the Python standard library.
"""
from fractions import Fraction as F
from itertools import permutations
from math import isqrt


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

print("Exact paper examples passed: SVD obstruction and state preparation, rank factorization, norm polynomial, tensor indexing.")

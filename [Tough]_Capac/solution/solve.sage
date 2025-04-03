#!/usr/bin/env sage

from Crypto.Util.number import inverse, long_to_bytes
import itertools
import ast


with open('output.txt', 'r') as f:
    output = f.read()
pubkey, enc = [ast.literal_eval(line.split(' = ')[-1]) for line in output.splitlines()]

nbit, r, s = 128, 4, 6


def add(P, Q, c, n):
    # Add two points P and Q on curve x^3 + c*y^3 + c^2*z^3 - 3*c*x*y*z = 1 in Zmod(n)
    (x1, y1, z1) = P
    (x2, y2, z2) = Q
    x3 = (x1 * x2 + c * (y2 * z1 + y1 * z2)) % n
    y3 = (x2 * y1 + x1 * y2 + c * z1 * z2) % n
    z3 = (y1 * y2 + x2 * z1 + x1 * z2) % n
    return (x3, y3, z3)


def mul(P, g, c, n):
    # Scalar multiplication on curve
    (x1, y1, z1) = P
    (x2, y2, z2) = (1, 0, 0)
    for b in bin(g)[2:]:
        (x2, y2, z2) = add((x2, y2, z2), (x2, y2, z2), c, n)
        if b == '1':
            (x2, y2, z2) = add((x2, y2, z2), (x1, y1, z1), c, n)
    return (x2, y2, z2)


n, e = pubkey
x, y, z = enc

p = 242195390295637766135570468161415180323
q = 171955359080932502551172606124599656543
assert p ^ r * q ^ s == n

E = [
    p ^ (2 * (r - 1)) * q ^ (2 * (s - 1)) * (p ^ 2 + p + 1) * (q ^ 2 + q + 1),
    p ^ (2 * (r - 1)) * q ^ (2 * (s - 1)) * (p - 1) ^ 2 * (q - 1) ^ 2,
    p ^ (2 * (r - 1)) * q ^ (2 * (s - 1)) * (p ^ 2 + p + 1) * (q - 1) ^ 2,
    p ^ (2 * (r - 1)) * q ^ (2 * (s - 1)) * (p - 1) ^ 2 * (q ^ 2 + q + 1),
]
D = [inverse(e, i) for i in E]

P.<t> = PolynomialRing(ZZ)
f = x ^ 3 + t * y ^ 3 + t ^ 2 * z ^ 3 - 3 * t * x * y * z - 1
cp = f.change_ring(Zp(p, prec=r)).roots(multiplicities=False)
cq = f.change_ring(Zp(q, prec=s)).roots(multiplicities=False)
c_list = [crt(list(map(ZZ, i)), [p ^ r, q ^ s]) for i in itertools.product(cp, cq)]
for c, d in itertools.product(c_list, D):
    x, y, z = mul(enc, d, c, n)
    if z != 0:
        continue
    flag = long_to_bytes(x) + long_to_bytes(y)
    print(flag)

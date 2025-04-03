#!/usr/bin/env sage

from Crypto.Util.number import long_to_bytes
import ast


with open('output.txt', 'r') as f:
    output = f.read()
e, p, PT = ast.literal_eval(output.split(' = ')[-1])


P.<x> = PolynomialRing(GF(p))
f = P.lagrange_polynomial(PT)
c = f[0]

for m in c.nth_root(e, all=True):
    flag = long_to_bytes(int(m))
    if flag.startswith(b'CCTF{'):
        print(flag)

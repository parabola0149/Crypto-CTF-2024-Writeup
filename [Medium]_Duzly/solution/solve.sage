#!/usr/bin/env sage

import string
import time
from sage.misc.banner import version_dict

assert version_dict()['major'] == 9


C = [1, 17761542461647558231, 13293668011354679701, 9204760597720472707, 8540722934676348527, 3568330912555059249]
pash = bytes.fromhex('29779e083d3a77a1a78a1ce145c0890c5996bd91d11cdf9f049de5895fae1de612f237cc6488eab27965d4cb2e3139affb83c79035265e08b50eed5eed8ad3741eca2a42814057e92ec5d857c6c17e4a95e7fec39e6dea49a9e369354ba9f72c138d5e0a5a813741080e87982382932ffe51e025ad4db71147a6cfde90821d143cb5e88f5742834b2a527cb6313e9281404d5ffbd5becd54b2f74aeb89f7079c23cf9c9e2b6da2896313aeef9282567fadb310b61f309efbfab5d62cbe650c81')


p = 2 ^ 64 - 59
E = [2 ^ 24 + 17, 2 ^ 24 + 3, 3, 2, 1, 0]

t0 = time.time()
P.<x> = PolynomialRing(GF(p), sparse=True)
pash = [int.from_bytes(pash[i:i + 8], 'big') for i in range(0, len(pash), 8)]
f0 = P(dict(zip(E, C)))

flag = b''
for idx, h in enumerate(pash):
    t1 = time.time()
    f1 = f0 - h
    roots = gp(f'polrootsmod({f1}, {p})')
    roots = [Integer(i.lift()) for i in roots]
    pt = [int(i).to_bytes(8, 'big') for i in roots]
    t = time.time() - t1
    print(f'{idx + 1:2}/{len(pash):2}: t = {t:.1f}, pt = {pt}')
    flag += max(pt, key=lambda b: sum(1 for i in b if i in string.printable.encode()))
print(flag)
t = time.time() - t0
print(t)

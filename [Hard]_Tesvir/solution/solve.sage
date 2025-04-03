#!/usr/bin/env sage

import ast


with open('output_updated.txt', 'r') as f:
    output = f.read()
r, pubkey, enc = [ast.literal_eval(line.split(' = ')[-1]) for line in output.splitlines() if line]


def solve_CVP(B, x):
    B, T = B.LLL(transformation=True)
    G, _ = B.gram_schmidt()
    v = x
    t = vector(ZZ, [0] * B.nrows())
    for i in reversed(range(G.nrows())):
        i -= G.nrows()
        c = round((v * G[i]) / (G[i] * G[i]))
        v -= c * B[i]
        t[i] = c
    return x - v, t * T


p, h = 223, 24
c = pubkey[:h] + pubkey[-h:]
n = len(c)
flag = b''
for s in enc:
    K = 2
    B = block_matrix(ZZ, [
        [matrix([c]).T, K * 1],
        [p ^ h - 1    , K * 0],
    ])
    x = vector(ZZ, [s] + [K // 2] * n)
    v, t = solve_CVP(B, x)
    
    m = sum(v * 2 ^ i for i, v in enumerate(reversed(t[:h])))
    flag += int(m).to_bytes(h // 8, 'big')

if not flag.endswith(b'}'):
    flag += b'}'
print(flag)

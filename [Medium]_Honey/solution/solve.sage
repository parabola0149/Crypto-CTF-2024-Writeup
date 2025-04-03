#!/usr/bin/env sage

from Crypto.Util.number import long_to_bytes
import ast


with open('params_enc.txt', 'r') as f:
    output = f.read()
p, Q, R, S, C = [ast.literal_eval(line.split(' = ')[-1]) for line in output.splitlines()]


def solve_CVP(B, x):
    B, T = B.LLL(transformation=True)
    # G, _ = B.gram_schmidt()
    nbits = max(i.nbits() for i in B.list())
    G = B.change_ring(RealField(2 * nbits))
    for i in range(G.nrows() - 1):
        G[i + 1:] -= ((G[i + 1:] * G[i:i + 1].T) / (G[i] * G[i])) * G[i:i + 1]
    v = x
    t = vector(ZZ, [0] * B.nrows())
    for i in reversed(range(G.nrows())):
        i -= G.nrows()
        c = ((v * G[i]) / (G[i] * G[i])).round()
        v -= c * B[i]
        t[i] = c
    return x - v, t * T


d = len(Q)
K1 = p // 2 ^ d
K2 = p
B = block_matrix(ZZ, [
    [K2 * matrix([Q])       , 1, K1 * 0, K1 * 0],
    [K2 * diagonal_matrix(R), 0, K1 * 1, K1 * 0],
    [K2 * diagonal_matrix(S), 0, K1 * 0, K1 * 1],
    [K2 * p                 , 0, K1 * 0, K1 * 0],
])
x = vector(ZZ, [K2 * i for i in C] + [p // 2] * (2 * d + 1))
v, t = solve_CVP(B, x)

m = t[0]
print(long_to_bytes(m))

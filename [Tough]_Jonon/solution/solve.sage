#!/usr/bin/env sage

import itertools
from string import printable as prn


with open('output.txt', 'r') as f:
    output = f.read()
pkey, C = output.removeprefix('pkey = (').split(')\nC = ')
pA, psA, pE = pkey.split(']], [')
pA = [matrix(ZZ, [[int(k) for k in j.strip('[]').split() if k] for j in i.splitlines()]) for i in pA.split(', ')]
psA = [matrix(ZZ, [[int(k) for k in j.strip('[]').split() if k] for j in i.splitlines()]) for i in psA.split(', ')]
pE = matrix(ZZ, [[int(k) for k in j.strip('[]').split() if k] for j in pE.splitlines()])
C = matrix(ZZ, [[int(k) for k in j.strip('[]').split() if k] for j in C.splitlines()])

k, l, _B = 5, 19, 63


def make_permutation(A, M):
    permutation_list = [([], M)]
    for _ in range(len(A)):
        permutation_list_2 = []
        for (L2, M2), (idx, A2) in itertools.product(permutation_list, enumerate(A)):
            if idx in L2:
                continue
            M3 = A2.solve_right(M2)
            if not all(i.is_integer() for i in M3.list()):
                continue
            permutation_list_2.append((L2 + [idx], M3))
        permutation_list = permutation_list_2
    [(L, _)] = permutation_list
    
    assert prod([A[i] for i in L]) == M
    return L


A, E, F = psA, pE, pA

M = (F[0] * E).list()
K = round(k ^ (l + 1) * (_B / 2) ^ (l + 2))
B = block_matrix(ZZ, [
    [K, matrix([M[1:]])],
    [0, -M[0]          ],
])
B, T = B.LLL(transformation=True)
idx = min(range(B.nrows()), key=lambda x: B[x] * B[x])
v, t = B[idx], T[idx]
if t[0] < 0:
    v, t = -v, -t

p1 = min(i // j for i, j in zip(M, t))
for i in itertools.count():
    p2 = p1 - i
    if not p2.is_prime():
        continue
    iE = E.change_ring(GF(p2)) ^ -1
    iA0 = A[0].change_ring(GF(p2)) ^ -1
    D2 = iA0 * iE * F[0] * E
    if all(i <= _B for i in D2.list()):
        p, D = p2, D2
        break
F = [i.change_ring(GF(p)) for i in F]
A = [i.change_ring(GF(p)) for i in A]
E = E.change_ring(GF(p))

C1 = E ^ -1 * C * E * D ^ -1
C2_list = [(D ^ -1 * i ^ -1 * C1).change_ring(ZZ) for i in A]
idx = min(range(l), key=lambda x: sum(C2_list[x].list()))
C2 = C2_list[idx]
D = D.change_ring(ZZ)
A = [i.change_ring(ZZ) for i in A]

C3 = A[idx] * D * C2 * D
A2 = [i * D for i in A]
_p = make_permutation(A2, C3)

flag = 'CCTF{' + ''.join([prn[10 + _p[i]] for i in range(l)]) + '}'
print(flag)

#!/usr/bin/env sage

from pwn import remote
import itertools


# host = '00.cr.yp.toc.tf'
host = 'localhost'
port = 31117


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


with remote(host, port) as conn:
    step = 12
    for level in range(step):
        k, l = 3, 12 + level
        conn.recvuntil(b'[Q]uit\n')
        
        conn.sendline(b'G')
        data = conn.recvuntil(b'[Q]uit\n')
        A = [matrix(ZZ, [[int(i3) for i3 in i2.strip('┃ []'.encode()).split(b' ') if i3] for i2 in data.splitlines()[k * i1:k * (i1 + 1)]]) for i1 in range(l)]
        
        conn.sendline(b'P')
        data = conn.recvuntil(b'[Q]uit\n')
        M = matrix(ZZ, [[int(i3) for i3 in i2.strip('┃ M = []'.encode()).split(b' ') if i3] for i2 in data.splitlines()[:k]])
        
        conn.sendline(b'S')
        L = make_permutation(A, M)
        print(L)
        conn.sendline(','.join(map(str, L)).encode())
    
    data = conn.recvall()
    print(data)

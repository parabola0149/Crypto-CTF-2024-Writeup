#!/usr/bin/env sage

from pwn import remote


# host = '00.cr.yp.toc.tf'
host = 'localhost'
port = 37771


def select_word(b, nl, nw):
    return b.splitlines()[nl].split(b' ')[nw]


def solve_CVP(B, x):
    # https://cims.nyu.edu/~regev/teaching/lattices_fall_2004/ln/cvp.pdf
    B, T = B.LLL(transformation=True)
    G, _ = B.gram_schmidt()
    v = x
    t = vector(ZZ, [0] * B.nrows())
    for i in reversed(range(G.nrows())):
        i -= G.nrows()
        c = ((v * G[i]) / (G[i] * G[i])).round(mode='even')
        v -= c * B[i]
        t[i] = c
    return x - v, t * T


def make_polynomial(n):
    d = ceil(2 * log(n)) - 1
    d *= 2
    P.<x> = PolynomialRing(ZZ)
    f1 = x ^ 2 + x + 2
    B1 = matrix(ZZ, d + 1, 2)
    for i in range(d + 1):
        f2 = x ^ i % f1
        B1[i] = [f2[0], f2[1]]
    C = d + 1
    B = block_matrix(ZZ, [[C * B1, 1]])
    x = vector(ZZ, [C * n, C * 0] + [0] * (d + 1))
    v, t = solve_CVP(B, x)
    
    f2 = P(list(t))
    return f2


with remote(host, port) as conn:
    for _ in range(11):
        data = conn.recvuntil(b'please send the polynomial f:\n')
        n = int(select_word(data, -1, 9).removesuffix(b','))
        f = make_polynomial(n)
        print(n, f)
        conn.sendline(f'{f}'.encode())
    data = conn.recvall()
    print(data)

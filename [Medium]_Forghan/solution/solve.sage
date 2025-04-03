#!/usr/bin/env sage

from pwn import remote
from Crypto.Util.number import long_to_bytes
import string


# host = '00.cr.yp.toc.tf'
host = 'localhost'
port = 13337


def select_word(b, nl, nw):
    return b.splitlines()[nl].split(b' ')[nw]


nbit = 256
p, q = [random_prime(2 ^ nbit, lbound = 2 ^ (nbit - 1)) for _ in range(2)]
print(f'{p = }')
print(f'{q = }')
n_factors = factor(1)
for i in [p - 1, p + 1, q - 1, q + 1]:
    i_factors = factor(i)
    print(i_factors)
    n_factors *= i_factors
phi = prod(p ^ (e - 1) * (p - 1) for p, e in n_factors)
n = n_factors.value()

with remote(host, port) as conn:
    conn.recvuntil(b'[Q]uit\n')
    conn.sendline(b's')
    conn.recvuntil(b'Send your desired prime numbers separated by comma: \n')
    conn.sendline(f'{p},{q}'.encode())
    
    mp, mq = None, None
    while True:
        conn.recvuntil(b'[Q]uit\n')
        conn.sendline(b'g')
        data = conn.recvuntil(b'Options: \n')
        cp = int(select_word(data, 0, -1))
        cq = int(select_word(data, 1, -1))
        
        conn.recvuntil(b'[Q]uit\n')
        conn.sendline(b'p')
        data = conn.recvuntil(b'Options: \n')
        yp = int(select_word(data, 2, -1))
        yq = int(select_word(data, 3, -1))
        
        if mp is None and gcd(yp, phi) == 1:
            mp = pow(cp, inverse_mod(yp, phi), n)
            print(f'{mp = }')
        if mq is None and gcd(yq, phi) == 1:
            mq = pow(cq, inverse_mod(yq, phi), n)
            print(f'{mq = }')
        
        if mp is not None and mq is not None:
            conn.recvuntil(b'[Q]uit\n')
            break

chr_list = string.printable.encode()
m_list = []
for m0 in [mp, mq]:
    g = prod(p ^ (e - 1) for p, e in n_factors if m0 % p == 0)
    l = []
    for i in range(g):
        m = (m0 + i * n // g) % n
        pt = long_to_bytes(int(m))
        if all(i in chr_list for i in pt):
            l.append(m)
    [m] = l
    m_list.append(m)
    print(f'{m = }')

flag = b''.join(long_to_bytes(int(m)) for m in m_list)
print(flag)

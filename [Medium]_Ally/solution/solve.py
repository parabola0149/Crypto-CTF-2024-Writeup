#!/usr/bin/env python3

from pwn import remote
from Crypto.Util.number import getPrime


# host = '00.cr.yp.toc.tf'
host = 'localhost'
port = 13777


def select_word(b, nl, nw):
    return b.splitlines()[nl].split(b' ')[nw]


with remote(host, port) as conn:
    for _ in range(20):
        data = conn.recvuntil(b'-bit prime:  \n')
        nbit = int(select_word(data, -1, 4).split(b'-')[0])
        
        while True:
            p = getPrime(nbit)
            if p % 4 == 1:
                break
        y = (p - 1) // 4
        x = 2 * y + 1
        print(p, x, y)
        
        conn.sendline(str(p).encode())
        data = conn.recvuntil(b' * (x - y)**3 = (x**2 + y) * (x + y**2)\n')
        conn.sendline(f'{x},{y}'.encode())
    
    data = conn.recvall()
    print(data)

#!/usr/bin/env python3

from pwn import remote


# host = '00.cr.yp.toc.tf'
host = 'localhost'
port = 37733


def select_word(b, nl, nw):
    return b.splitlines()[nl].split(b' ')[nw]


with remote(host, port) as conn:
    conn.recvuntil(b'[Q]uit\n')
    
    conn.sendline(b'G')
    data = conn.recvuntil(b'[Q]uit\n')
    q, m, d, u = [int(select_word(data, 0, i).strip(b'(,)')) for i in range(-4, -1 + 1)]
    
    conn.sendline(b'O')
    data = conn.recvuntil(b'[Q]uit\n')
    U = [int(i) for i in data.splitlines()[0].strip('┃ The oracle output U = []'.encode()).split(b' ') if i]
    print(U)
    
    conn.sendline(b'V')
    for _ in range(m + d):
        conn.recvuntil(b':\n')
        conn.sendline(','.join(map(str, U)).encode())
    
    data = conn.recvall()
    print(data)

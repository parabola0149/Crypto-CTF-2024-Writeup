#!/usr/bin/env python3

from pwn import remote


# host = '00.cr.yp.toc.tf'
host = 'localhost'
port = 13377


def select_word(b, nl, nw):
    return b.splitlines()[nl].split(b' ')[nw]


with remote(host, port) as conn:
    data = conn.recvuntil('┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛\n'.encode())
    print(data.decode(), end='')
    x1, y1, x2, y2 = [int(select_word(data, i, j).strip(b'f(),')) for i, j in [(4, 2), (4, 4), (5, 4), (5, 6)]]
    a = (y1 - y2) // (x1 - x2)
    b = y1 - a * x1
    
    while True:
        data = conn.recvline()
        print(data.decode(), end='')
        x = int(select_word(data, 0, -1).strip(b'f():'))
        
        y = a * x + b
        
        data = str(y)
        print(data)
        conn.sendline(data.encode())
        
        data = conn.recvline()
        print(data.decode(), end='')
        if not data.startswith('┃ Good job, try the next step '.encode()):
            break

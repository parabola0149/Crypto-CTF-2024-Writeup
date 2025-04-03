#!/usr/bin/env python3

from pwn import remote


# host = '00.cr.yp.toc.tf'
host = 'localhost'
port = 17113


def select_word(b, nl, nw):
    return b.splitlines()[nl].split(b' ')[nw]


with remote(host, port) as conn:
    while True:
        data = conn.recvuntil(b'Please send x, y separated by comma: \n')
        print(data.decode(), end='')
        c = int(select_word(data, -2, 6))
        z = int(select_word(data, -2, 11))
        
        y = z - c
        x = y + 1
        
        data = f'{x},{y}'
        print(data)
        conn.sendline(data.encode())
        
        data = conn.recvuntil(b'\n')
        print(data.decode(), end='')
        if data.strip().endswith(b'Congratulation! You got the flag!'):
            break
    
    data = conn.recvall()
    print(data.decode(), end='')

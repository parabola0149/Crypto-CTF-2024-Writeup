#!/usr/bin/env python3

from pwn import remote
import subprocess


# host = '00.cr.yp.toc.tf'
host = 'localhost'
port = 17371


def select_word(b, nl, nw):
    return b.splitlines()[nl].split(b' ')[nw]


with remote(host, port) as conn:
    while True:
        data = conn.recvuntil(b', separated by comma: \n')
        print(data.decode(), end='')
        n = int(select_word(data, -1, 3).strip(b','))
        
        result = subprocess.run(['./sturdy_tester', 'bfs01', 'msw', str(n), str(n)], capture_output=True)
        c = int(result.stdout)
        result = subprocess.run(['./sturdy_tester', 'bfs01', 'swm', str(n), str(n)], capture_output=True)
        m = int(result.stdout)
        
        data = f'{m},{c}'
        print(data)
        conn.sendline(data.encode())
        
        data = conn.recvline()
        print(data.decode(), end='')
        if b'CCTF{' in data:
            break

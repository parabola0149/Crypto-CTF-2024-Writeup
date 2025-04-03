#!/usr/bin/env python3

from pwn import remote
import subprocess


# host = '00.cr.yp.toc.tf'
host = 'localhost'
port = 13771


file_1 = 'collision_file_1.bin'
file_2 = 'collision_file_2.bin'

with remote(host, port) as conn:
    REC = []
    STEP = 7
    for _ in range(STEP + 1):
        while True:
            subprocess.run(['./md5_fastcoll', '-o', file_1, file_2], capture_output=True)
            with open(file_1, 'rb') as f:
                d1 = f.read()
            with open(file_2, 'rb') as f:
                d2 = f.read()
            if d1 not in REC and d2 not in REC:
                break
        REC += [d1, d2]
        print(d1.hex(), d2.hex())
        
        conn.recvuntil(b'Please send your first input:  \n')
        conn.sendline(d1.hex().encode())
        
        conn.recvuntil(b'Please send your second input: \n')
        conn.sendline(d2.hex().encode())
    
    data = conn.recvall()
    print(data)

#!/usr/bin/env python3

from pwn import remote
from Crypto.Util.number import isPrime, getPrime, inverse, long_to_bytes
import subprocess
import re


# host = '00.cr.yp.toc.tf'
host = 'localhost'
port = 13371


def select_word(b, nl, nw):
    return b.splitlines()[nl].split(b' ')[nw]


nbit = 1024
e = 65537

while True:
    result = subprocess.run(['./cunningham_chain', f'bits={nbit - 1}', 'length=3', 'kind=1'], capture_output=True)
    match = re.search(rb'bits: ([0-9]+),', result.stdout)
    bits = int(match.group(1))
    print(f'{bits = }')
    if bits != nbit - 1:
        continue
    match = re.search(rb'origin: "([0-9]+)",', result.stdout)
    p_a = int(match.group(1))
    # p_a = 73394155025171296885757898802282475356705770334268142168102878148965878069337481966986951297893863735742954944341586580580799771518333571158595439778641255377047133669911179691390117778164038162046973436922914270185247615066178984447891243651280385648158173061382645195704483369322890382839286123309185880709
    print(f'{p_a = }')
    
    p_b = 2 * p_a + 1
    p_c = 2 * p_b + 1
    assert all(isPrime(i) for i in [p_a, p_b, p_c])
    try:
        d1 = inverse(e, p_a - 1)
        d2 = inverse(e, p_c - 1)
    except ValueError:
        continue
    break


with remote(host, port) as conn:
    conn.recvuntil(b'prime numbers p, q:\n')
    
    p = p_b
    q = getPrime(nbit)
    conn.sendline(f'{p},{q}'.encode())
    data = conn.recvuntil(b'Now, send us the secret string to get the flag: \n')
    c1 = int(select_word(data, 0, -1))
    c2 = int(select_word(data, 1, -1))
    
    m1 = pow(c1, d1, p_a)
    m2 = pow(c2, d2, p_c)
    s = long_to_bytes(m1) + long_to_bytes(m2)
    conn.sendline(s)
    data = conn.recvall()
    print(data)

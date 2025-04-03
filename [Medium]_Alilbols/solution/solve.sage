#!/usr/bin/env sage

from Crypto.Util.number import long_to_bytes


with open('ouput.txt', 'r') as f:
    output = f.read()
h, c = [int(line.split(' = ')[-1]) for line in output.splitlines()]


d = ceil(log(h / 4, 100))
q = 4 * 100 ^ d
B = matrix(ZZ, [
    [h + 1, 1, 0     ],
    [c    , 0, 10 ^ d],
    [q    , 0, 0     ],
])

B = B.LLL()
v = min(B, key=lambda x: x * x)
if v[2] < 0:
    v = -v
m = v[0]

print(long_to_bytes(m))

#!/usr/bin/env sage

from Crypto.Util.number import long_to_bytes


with open('output.txt', 'r') as f:
    output = f.read()
c = output.split(' = ')[-1].strip()


prec = len(c.split('e')[0].replace('.', ''))
prec = floor(log(10 ^ prec, 2))
c = Rational(RealNumber(c))

for b in range(prec, prec + 10):
    for d in divisors(2 * b):
        g = d - 1
        if g <= 0:
            continue
        m = 2 * (c / b + 1 / (g + 1)) - 1
        try:
            m = Integer(m)
        except TypeError:
            continue
        print(f'CCTF{{{long_to_bytes(m).decode()}}}')

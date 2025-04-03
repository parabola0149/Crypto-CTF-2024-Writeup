#!/usr/bin/env sage

from hashlib import sha256


nbit = 4 * 32
C = CyclotomicField(nbit)
F, _ = C.maximal_totally_real_subfield()
_B = nbit >> 2

with open('pubkey_update', 'r') as f:
    output = f.read()
output = [line.strip('[]') for line in output.splitlines()]
pkey = matrix([[F(line[:len(line) // 2]), F(line[-len(line) // 2:])] for line in output])


P.<a, b, c, d> = PolynomialRing(F)
skey = matrix([[a, c], [b, d]])
M = (skey + identity_matrix(2)) * (skey.T - identity_matrix(2)) - pkey
f_list = M.list() + [a * d - b * c - 1]
I = Ideal(f_list)
B = I.groebner_basis()

assert [f.monomials() for f in B] == [[d ^ 2, d, 1], [a, d, 1], [b, d, 1], [c, d, 1]]
f1, f2, f3, f4 = B
sd_list = f1.polynomial(d).change_ring(F).roots(multiplicities=False)

for sd in sd_list:
    [sa], [sb], [sc] = [f.subs(d=sd).polynomial(var).change_ring(F).roots(multiplicities=False) for f, var in zip([f2, f3, f4], [a, b, c])]
    skey = matrix([[sa, sc], [sb, sd]])
    assert (skey + identity_matrix(2)) * (skey.T - identity_matrix(2)) == pkey
    if not all(all(-_B <= j <= _B + 1 for j in i.list()) for i in [sa, sb]):
        continue
    
    flag = str(sum(skey.coefficients()).polynomial()(2^128))
    flag = 'CCTF{' + sha256(flag.encode()).hexdigest() + '}'
    print(flag)

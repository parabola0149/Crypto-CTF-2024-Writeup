#!/usr/bin/env sage

from Crypto.Util.number import long_to_bytes
import itertools
import string


Gp = (2024, 77522466419731346352388161327659294096)
Hp = (187001239996895215821259450553198409012, 158495938026527642884038157170741730943)
Gq = (2024, 92609909821520263487623088269239797003)
Hq = (191534442430060308634251661645421139195, 102273233384427938890774177710170123915)
c  = 15084463560924811262750235394027264639346464192638172940901706702947534963652

nbit = 128


kq = (Hq[0] ^ 3 - Hq[1] ^ 2) * Gq[0] - (Gq[0] ^ 3 - Gq[1] ^ 2) * Hq[0]
[q] = [p for p, _ in factor(kq) if p.nbits() == nbit]
# kq_factors = [5, 5, 11, 23, 521, 1271953, 3018047, 25079843, 5681358923, 28338689773015872412188935645719088445835579, 278451262698064898668334196027031252819]
# assert prod(kq_factors) == kq and all(p.is_prime() for p in kq_factors)
# [q] = [p for p in kq_factors if p.nbits() == nbit]
print(q)
a = (Gq[1] ^ 2 - Gq[0] ^ 3) * inverse_mod(Gq[0], q) % q
kp = Gp[0] ^ 3 + a * Gp[0] - Gp[1] ^ 2
[p] = [p for p, _ in factor(kp) if p.nbits() == nbit]
print(p)

Ep = EllipticCurve(GF(p), [a, 0])
Eq = EllipticCurve(GF(q), [a, 0])
Gp, Hp, Gq, Hq = Ep(Gp), Ep(Hp), Eq(Gq), Eq(Hq)
sp = Hp.log(Gp)
# sp = 50295274560999838770579032062844522666
# assert sp * Gp == Hp
[s] = [i for i in range(sp, p + 1, Gp.order()) if i * Gq == Hq]
print(s)

mp = mod(c, p).nth_root(s, all=True)
mq = mod(c, q).nth_root(s, all=True)
for i in itertools.product(mp, mq):
    m = crt(list(map(ZZ, i)), [p, q])
    flag_content = long_to_bytes(int(m))
    if not all(i in string.printable.encode() for i in flag_content):
        continue
    flag = f'CCTF{{{flag_content.decode()}}}'
    print(flag)

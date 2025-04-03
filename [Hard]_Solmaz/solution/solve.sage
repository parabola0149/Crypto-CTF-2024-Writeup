#!/usr/bin/env sage

from Crypto.Util.number import long_to_bytes
from tqdm import tqdm


P = (1338, 9218578132576071095213927906233283616907115389852794510465118810355739314264)
Q = (3454561753909947353764378180794889923919743476068813953002808647958908878895, 17267599534808803050751150297274989016063324917454246976792837120400888025519)


Px, Py = P
Qx, Qy = Q
r = (Py ^ 2 - Px ^ 3) * Qx - (Qy ^ 2 - Qx ^ 3) * Px

a = 2
while gcd(a, r) != 1:
    a += 1
a = power_mod(a, 4, r)
j = next_prime(2 ^ 31)
for i in tqdm(range(2 ^ 31, 2 ^ 32)):
    if i < j:
        continue
    j = next_prime(j)
    a = power_mod(a, i ^ 2, r)
    g = gcd(a - 1, r)
    if (32 - 1) * 8 + 2 <= g.nbits() < 32 * 8 + 2:
        break
    elif g != 1:
        r //= g
else:
    assert False
p = g
# p = 30126567747372029007183424263223733382328264316268541293679065617875255137317
assert p.is_prime()
print(p)

K = GF(p)
Px, Py = K(Px), K(Py)
c = (Py ^ 2 - Px ^ 3) * Px ^ -1
E = EllipticCurve(K, [c, 0])
P, Q = E(P), E(Q)
m = Q.log(P)
print(f'CCTF{{{long_to_bytes(m).decode()}}}')

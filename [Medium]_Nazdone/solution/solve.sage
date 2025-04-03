#!/usr/bin/env sage

from sage.rings.factorint import factor_trial_division
from Crypto.Util.number import long_to_bytes


n = 301929748923678449872944933611657670834216889867340028357609265175830693931365828840717548752313862343315133541384709574659039910206634528428504034051556622114290811586746168354731258756502196637977942743110508997919976400864419640496894428180120687863921269087080600917900477624095004141559042793509244689248253036809126205146653922738685595903222471152317095497914809983689734189245440774658145462867680027337
c = 104375152140523502741159687674899095271676058870899569351687154311685938980840028326701029233383897490722759532494438442871187152038720886122756131781086198384270569105043114469786514257765392820254951665751573388426239366215033932234329514161827069071792449190823827669673064646681779764841034307000600929149689291216313319444583032339045277433847691961234044840927155960887984372868669401051358701522484473320


divisor_list = []
for c0 in [1, 2]:
    t = n - c0 ^ 3
    t_factor = factor_trial_division(t, limit=2 ^ 24)
    divisor_list += divisors(prod(p ^ e for p, e in t_factor if p.is_prime()))

P.<x> = PolynomialRing(ZZ)
f_factor_list = []
for d in sorted(set(divisor_list)):
    if d <= 1:
        continue
    f = P(n.digits(d))
    f_factor = factor(f)
    try:
        f1, f2, f3 = [i for i, _ in f_factor]
    except ValueError:
        continue
    f_factor_list.append((d, f1, f2, f3))

[(M, f1, f2, f3)] = f_factor_list
p, q, r = [f(M) for f in [f1, f2, f3]]
assert p * q * r == n and all(i.is_prime() for i in [p, q, r])
phi = (p - 1) * (q - 1) * (r - 1)

z1 = max(sum(1 for i in f if i != 0) - 1 for f in [f1, f2, f3])
for z in range(z1, 3 * z1):
    e = M ^ 3 + z - 2
    try:
        d = inverse_mod(e, phi)
    except ZeroDivisionError:
        continue
    m = pow(c, d, n)
    flag = long_to_bytes(int(m))
    if flag.startswith(b'CCTF{'):
        print(flag)

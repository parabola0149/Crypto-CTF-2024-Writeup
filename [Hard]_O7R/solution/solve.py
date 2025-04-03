#!/usr/bin/env python3

from Crypto.Util.number import inverse, long_to_bytes
import itertools
import ast


with open('7seg.txt', 'r') as f:
    output = f.read()
p_7seg, q_7seg, n_7seg, n2_7seg, c_7seg = [ast.literal_eval(line.split(' = ')[-1]) for line in output.splitlines()]

e = 0x10001


table_7seg = ['######-', '-##----', '##-##-#', '####--#', '-##--##', '#-##-##', '#-#####', '###----', '#######', '####-##']
table_7seg = {v: i for i, v in enumerate(table_7seg)}

table_7seg_corrupted = {k: [v] for k, v in table_7seg.items()}
for digit_7seg, digit in table_7seg.items():
    for i in range(len(digit_7seg)):
        if digit_7seg[i] == '-':
            continue
        digit_7seg_2 = list(digit_7seg)
        digit_7seg_2[i] = '-'
        digit_7seg_2 = ''.join(digit_7seg_2)
        table_7seg_corrupted.setdefault(digit_7seg_2, []).append(digit)

p_7seg, q_7seg, n_7seg, n2_7seg, c_7seg = [list(reversed(list_7seg)) for list_7seg in [p_7seg, q_7seg, n_7seg, n2_7seg, c_7seg]]
c = sum(table_7seg[v] * 10 ** i for i, v in enumerate(c_7seg))
param_list = [(0, 0, 0, 0)]
assert len(p_7seg) == len(q_7seg)
for idx in range(len(p_7seg)):
    param_list_2 = []
    modulus = 10 ** (idx + 1)
    for (p_1, q_1, n_1, n2_1), p_2, q_2, n_2, n2_2 in itertools.product(param_list, *[table_7seg_corrupted[list_7seg[idx]] for list_7seg in [p_7seg, q_7seg, n_7seg, n2_7seg]]):
        p, q, n, n2 = [i + j * 10 ** idx for i, j in zip([p_1, q_1, n_1, n2_1], [p_2, q_2, n_2, n2_2])]
        if (p * q) % modulus == n and n ** 2 % modulus == n2:
            param_list_2.append((p, q, n, n2))
    param_list = param_list_2
[(p, q, _, _)] = param_list
n = p * q
n2 = n ** 2
assert all(all((param // 10 ** i) % 10 in table_7seg_corrupted[v] for i, v in enumerate(list_7seg)) for param, list_7seg in zip([p, q, n, n2], [p_7seg, q_7seg, n_7seg, n2_7seg]))

phi = (p - 1) * (q - 1)
d = inverse(e, phi)
m = pow(c, d, n)

print(long_to_bytes(m))

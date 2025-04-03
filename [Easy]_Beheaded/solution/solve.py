#!/usr/bin/env python3

import subprocess
from pathlib import Path


block_size = 16
bits = 16
width, height = 1623, 94

with open('all_flags.enc', 'rb') as f:
    enc_data = f.read()

img_size_1 = (3 * bits // 8) * width * height
img_size_2 = (img_size_1 // block_size + 1) * block_size
size = len(enc_data)
assert size % img_size_2 == 0
header = f'P6\n{width} {height}\n{(1 << bits) - 1}\n'.encode()
background_block = enc_data[:block_size]

dir_path = Path('output')
dir_path.mkdir(exist_ok=True)
for i in range(size // img_size_2):
    img_data = bytearray(header)
    for j in range(0, img_size_2, block_size):
        pos = i * img_size_2 + j
        block = enc_data[pos:pos + block_size]
        img_data += (b'\xff' if block == background_block else b'\x00') * block_size
    file_path = dir_path / f'flag_{i}.png'
    subprocess.run(['convert', '-', file_path], input=img_data)

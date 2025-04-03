#!/usr/bin/env python3

from PIL import Image


x0, y0 = 250, 100
wm, hm = 339.6, 150.52
w1, h1 = 41.41, 111.53
w2, h2 = 15.5, 15.5
pos_7seg = [(1, 0), (2, 1), (2, 3), (1, 4), (0, 3), (0, 1), (1, 2)]

text_7seg = ''
for name in ['p', 'q', 'n', 'n2', 'c']:
    img = Image.open(f'img/{name}.png').convert('L')
    px = img.load()
    width, height = img.size
    
    w, h = round((width - wm) / w1), round((height - hm) / h1)
    list_7seg = []
    for y1 in range(h):
        for x1 in range(w):
            digit_7seg = ''
            for x2, y2 in pos_7seg:
                x = round(x0 + x1 * w1 + x2 * w2)
                y = round(y0 + y1 * h1 + y2 * h2)
                p = px[x, y]
                digit_7seg += '#' if p < 0x80 else '-'
            if digit_7seg != '-' * 7:
                list_7seg.append(digit_7seg)
    text_7seg += f'{name}_7seg = {list_7seg}\n'

with open('7seg.txt', 'w') as f:
    f.write(text_7seg)

# Regenerates erfc_z.in (2,000,000 ErfC arguments) and spouge_grid.in (1,120,152 Spouge (y, m, n) calls).
import random, struct
random.seed(20261027)
def hx(v): return '%016x' % struct.unpack('<Q', struct.pack('<d', v))[0]
with open('erfc_z.in','w') as f:
    for _ in range(1000000): f.write('E %s\n' % hx(random.uniform(0.5, 28.0)))
    for _ in range(500000): f.write('E %s\n' % hx(random.uniform(-6.0, -0.5)))
    for _ in range(500000): f.write('E %s\n' % hx(random.uniform(0.5, 6.0)))
with open('spouge_grid.in','w') as f:
    for m in [30, 57, 99, 120, 150, 233, 300, 451, 600, 1000, 2500]:
        for n in [25, 64, 107, 112, 180, 256, 333, 500, 777, 1200, 5000, 35213]:
            for y in range(600, 60000, 7):
                f.write('S %d %d %d 186063409\n' % (y, m, n))

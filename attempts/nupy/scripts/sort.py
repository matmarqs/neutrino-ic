#!/usr/bin/env python3

# sort data

import sys
import numpy as np
fl = "{:15.5e}"

file = sys.argv[1]

x, y, z = np.loadtxt(file, unpack=True)

a = np.transpose([x, y, z])

L = a.tolist()
L = sorted(L)
N = len(L)

for i in range(N):
    print((3*fl).format(L[i][0], L[i][1], L[i][2]))

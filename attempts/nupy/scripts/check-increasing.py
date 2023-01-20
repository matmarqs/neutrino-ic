#!/usr/bin/env python3

# this script just checks if the data is increasing

import sys
import numpy as np

x, y, z = np.loadtxt(sys.stdin, unpack=True)

for i in range(len(x)-1):
    oldx = x[i]
    newx = x[i+1]
    if newx <= oldx:
        oldy = y[i]
        newy = y[i+1]
        if newy <= oldy:
            print("Failed! Not increasing!")
            print("dado failed: x =", x[i])
            quit()

print("Success! It is increasing!")

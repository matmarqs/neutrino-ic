#!/usr/bin/env python3

import numpy as np
#import math

R, lnNe = np.loadtxt("./nele_bs05op.dat", unpack=True)

num = len(R)

for i in range(num):
    print(" {:.7f}   {:.10e}".format(R[i], 10.0 ** lnNe[i]))

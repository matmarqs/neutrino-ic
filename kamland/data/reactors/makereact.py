#!/usr/bin/env python3

import sys

if len(sys.argv) < 2:
    print(f"Usage: {sys.argv[0]} [input file]", file=sys.stderr)
    exit()

inp = open(sys.argv[1], "r")

name = inp.readline().strip()
react_type = inp.readline().strip()
thermal_cap = inp.readline().strip()

file_name = name.lower().replace(" ", "_") + ".react"
f = open(file_name, "w")

f.write(f"Name                    {name}\n")
f.write(f"ReactorType             {react_type}\n")
f.write(f"ThermalCap[MW_t]        {thermal_cap}\n")
f.write("Year    ElecSupp[GW.h]  RefUnitPow[MW]  TimeOnLine[h]   OpFac[%]    EnergAvFac-An[%]    EnergAvFac-Cu[%]    LF-An[%]    LF-Cu[%]\n")

pos_list = [0, 8, 24, 40, 56, 68, 88, 108, 120]
def spaces(k):
    sp = ""
    for _ in range(k):
        sp += " "
    return sp

for line in inp:
    l = line.split()
    while len(l) < 9:
        l.append("-")
    string = ""
    for n in range(len(pos_list) - 1):
        string += l[n]
        string += spaces(pos_list[n+1] - pos_list[n] - len(l[n]))
    string += l[len(pos_list) - 1]
    f.write(string + '\n')

inp.close()
f.close()

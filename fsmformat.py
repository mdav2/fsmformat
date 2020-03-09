#!/usr/bin/env python
# fsmformat.py
# tries all like atom permutations to find maximum overlap
# between two geometries
# usage: python fsmformat.py 1/geom.xyz 2/geom.xyz
#        python fsmformat.py <path-to-geom1> <path-to-geom2>
# outputs aligned geometries separated by "---" which
# is the input syntax for Q-Chem.
# NB: please check the output geometries! There have been
# some issues with the ordering of the hydrogen atoms.
# Usually the heavies are correctly ordered, and only a couple
# of small changes may be needed.

import re
import sys
import numpy as np
from itertools import permutations

with open(sys.argv[1], 'r') as f:
    xyz1 = f.read()
with open(sys.argv[2], 'r') as f:
    xyz2 = f.read()


def getatoms(astr):
    atoms = re.findall("([A-Z])", astr)
    return atoms


def getpos(astr):
    pos = re.findall("[A-Z]\s+([-.0-9]+)\s+([-.0-9]+)\s+([-.0-9]+)\s", astr)
    return pos


def getall(at, pos, sym):
    lines = []
    O = []
    for idx, val in enumerate(at):
        if val == sym:
            lines.append(idx)
            O.append(pos[idx, :])
    return np.array(lines), np.array(O)


def dist(A):
    return np.sqrt(np.sum(A * A))


def maxco(A, B):
    assert A.shape == B.shape, "Differently sized arrays?"
    ncol = A.shape[1]
    nrow = A.shape[0]
    perms = list(permutations(range(nrow)))
    best = 1E6
    bestperm = None
    for i in perms:
        d = dist(A[i, :] - B)
        if d < best:
            bestperm = i
            best = d

    return bestperm, best


def maxalign(astr1, astr2):
    atoms1 = getatoms(astr1)
    atoms2 = getatoms(astr2)
    pos1 = np.array(getpos(astr1), dtype="float64")
    pos2 = np.array(getpos(astr2), dtype="float64")

    natom = len(atoms1)
    assert natom == len(atoms2), "Not the same number of atoms?"
    output = np.zeros((natom, 3))
    atoms = np.zeros(natom, dtype="<U4")
    for at in ['C', 'H', 'O']:
        orig_C1, Cs1 = getall(atoms1, pos1, at)
        orig_C2, Cs2 = getall(atoms2, pos2, at)
        p, x = maxco(Cs1, Cs2)
        for idx, val in enumerate(orig_C1):
            output[orig_C2[idx]] = Cs1[p, :][idx, :]
            atoms[orig_C2[idx]] = at
    return atoms, output


atoms, output = maxalign(xyz1, xyz2)
ostr = ""
for idx, line in enumerate(atoms):
    ostr = ostr + "{} {} {} {}\n".format(line, *output[idx])

print("made {} align with {}".format(sys.argv[1], sys.argv[2]))
#print("\n")
#print(ostr.strip())
with open("aligned_"+sys.argv[1],'w') as f:
    f.write(str(len(atoms))+"\n\n")
    f.write(ostr.strip())
with open("aligned_"+sys.argv[2],'w') as f:
    f.write(xyz2.strip())
#print('---')
#print(xyz2.strip())

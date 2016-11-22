#!/usr/bin/env python2
#
# MIT License
# 
# Copyright (c) 2016 Anders Steen Christensen
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import sys
import numpy as np

to_deg = 180.0 / np.pi

# Vector perpendicular to (v1-v2) and (v3-v2) centered in v2.
def perp_vector(v1, v2, v3):
    J = np.cross(v2 - v3, v2 - v1) + v2
    J = (J - v2) /(np.sqrt(np.dot(J - v2, J - v2)))
    return J

# Rotate a set of vectors
def rotate(V, J, T, center=None):

    if center is not None:
        V -= center
        
    x = V[0]
    y = V[1]
    z = V[2]
    u = J[0]
    v = J[1]
    w = J[2]
    a = (u*(u*x + v*y + w*z) + (x * (v*v + w*w) - u *(v*y + w*z))*np.cos(T) + np.sqrt(u*u + v*v + w*w)*(-w*y + v*z)*np.sin(T))/(u*u + v*v + w*w)
    b = (v*(u*x + v*y + w*z) + (y * (u*u + w*w) - v *(u*x + w*z))*np.cos(T) + np.sqrt(u*u + v*v + w*w)*( w*x - u*z)*np.sin(T))/(u*u + v*v + w*w)
    c = (w*(u*x + v*y + w*z) + (z * (u*u + v*v) - w *(u*x + v*y))*np.cos(T) + np.sqrt(u*u + v*v + w*w)*(-v*x + u*y)*np.sin(T))/(u*u + v*v + w*w)
    k =  np.array([a, b, c])

    if center is not None:
        k += center

    return k

# Bond angle between three coordinates
def bondangle(a,b,c):
    # In case numpy.dot() returns larger than 1
    # and we cannot take acos() to that number
    acos_out_of_bound = 1.0
    v1 = a - b
    v2 = c - b
    v1 = v1 / np.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)
    v2 = v2 / np.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2)
    dot_product = np.dot(v1,v2)

    if dot_product > acos_out_of_bound:
        dot_product = acos_out_of_bound
    if dot_product < -1.0 * acos_out_of_bound:
        dot_product = -1.0 * acos_out_of_bound

    return np.arccos(dot_product)

class Molecule:
        
    def __init__(self, reac_xyz, prod_xyz, atoms, natoms):

        self.reac_xyz = reac_xyz
        self.prod_xyz = prod_xyz
        self.atoms = atoms
        self.natoms = natoms


def parse_library_file(filename):

    f = open(filename, "r")
    lines = f.readlines()
    f.close()

    natoms = int(lines[0])

    xyz = []

    atoms = []

    for i, line in enumerate(lines[2:2+natoms]):

        tokens = line.split()
        # print "%4i:" % (i+1), line,

        atoms.append(tokens[0])
        x = float(tokens[1])
        y = float(tokens[2])
        z = float(tokens[3])
        xyz.append([x,y,z])


    return atoms, xyz

def parse_xyz(filename):

    f = open(filename, "r")
    lines = f.readlines()
    f.close()

    bonds = []

    natoms = int(lines[0])

    for i, line in enumerate(lines): 

        if "$BONDS" in line:

            tokens = line.split()
            bond_string = "".join(tokens[1:])
            bonds = eval(bond_string)

    reac_xyz = []
    prod_xyz = []

    atoms = []

    # print "Reactant coordinates:"
    for i, line in enumerate(lines[2:2+natoms]):

        tokens = line.split()
        #print "%4i:" % (i+1), line,

        atoms.append(tokens[0])
        x = float(tokens[1])
        y = float(tokens[2])
        z = float(tokens[3])
        reac_xyz.append([x,y,z])

    #print "PRODUCT:"
    for i, line in enumerate(lines[natoms+4:4+2*natoms]):

        tokens = line.split()
        #print "%4i:" % (i+1), line,
        atom = tokens[0]

        if atoms[i] != atom:
            print "Auto-TS error: Mismatch in atoms at position #%s!" % (i+1)
            exit(1)


        x = float(tokens[1])
        y = float(tokens[2])
        z = float(tokens[3])
        prod_xyz.append([x,y,z])



    # print "Found replaceable bonds:", bonds
    # print "System has atoms:", natoms


    mol = Molecule(reac_xyz, prod_xyz, atoms, natoms)

    return bonds, mol


def connect(molecule, mutation, bond, mol_atoms, mut_atoms):

    c = bond[0] - 1
    h = bond[1] - 1

    mol = np.array(molecule)
    mut = np.array(mutation)

    anchor = mol[c]

    translation = mut[0] - anchor

    mut = mut - translation

    for i, x in enumerate(mol):
        if i != h:
            print mol_atoms[i], x[0], x[1], x[2]

    # print mol[h], bond[1]
    # print mol[c], bond[0]
    # print mut[1], "N"
    # print mut[0], "C"

    a = bondangle(mol[h], mol[c], mut[1])

    p = perp_vector(mol[h], mol[c], mut[1])

    for i, x in enumerate(mut):

        mut[i] = rotate(x, p, a, center=mol[c])

        if i != 0:
            print mut_atoms[i], mut[i][0],  mut[i][1], mut[i][2]




if __name__=="__main__":

    input_filename = sys.argv[1]
    input_mutation = sys.argv[2]

    bonds, mol = parse_xyz(input_filename)
    atoms, xyz = parse_library_file(input_mutation)    
    # print mol.reac_xyz
    # print mol.prod_xyz
    # print mol.atoms
    # print mol.natoms
    # print bonds
    # print atoms, xyz
    print 18
    print
    connect(mol.prod_xyz, xyz, bonds[0], mol.atoms, atoms)
    print 18
    print
    connect(mol.reac_xyz, xyz, bonds[0], mol.atoms, atoms)

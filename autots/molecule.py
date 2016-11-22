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

class Molecule:
        
    def __init__(self, filename):

        f = open(filename, "r")
        lines = f.readlines()
        f.close()

        self.bonds = []
        self.reac_xyz = []
        self.prod_xyz = []
        self.atoms = []
        self.natoms = int(lines[0])

        for i, line in enumerate(lines): 

            if "$BONDS" in line:

                tokens = line.split()
                bond_string = "".join(tokens[1:])
                self.bonds = eval(bond_string)

        # print "Reactant coordinates:"
        for i, line in enumerate(lines[2:2+self.natoms]):

            tokens = line.split()
            #print "%4i:" % (i+1), line,

            self.atoms.append(tokens[0])
            x = float(tokens[1])
            y = float(tokens[2])
            z = float(tokens[3])
            self.reac_xyz.append([x,y,z])

        #print "PRODUCT:"
        for i, line in enumerate(lines[self.natoms+4:4+2*self.natoms]):

            tokens = line.split()
            #print "%4i:" % (i+1), line,
            atom = tokens[0]

            if self.atoms[i] != atom:
                print "Auto-TS error: Mismatch in atoms at position #%s!" % (i+1)
                exit(1)


            x = float(tokens[1])
            y = float(tokens[2])
            z = float(tokens[3])
            self.prod_xyz.append([x,y,z])



        # print "Found replaceable bonds:", bonds
        # print "System has atoms:", natoms


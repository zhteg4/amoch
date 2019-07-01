"""
Run LAMMPS (with Dreiding parameters) as a Simulator

Copyright 

LICENSE
"""

import math
from workflow import Simulator, od_ext

class LAMMPS_Dreiding(Simulator):

    def __init__(self, echo_output=True):
        path = 'lammps_dreiding'+od_ext()
        Simulator.__init__(self, path=path, echo_output=echo_output)
        self._lj = []
        self._x6 = []
        self.inputs.add('MINIMIZE_PAIR_STYLE', 'string', \
                        'pair_style        lj/cut 6.0')
        self.inputs.add('MD_PAIR_STYLE', 'string', \
                        'pair_style        buck/coul/long  12.0 12.0')
                              
    def _read_lj(self):
        """Read LJ parameters from data4Lammps output file"""
        with open('LJpaircoeffs.txt', 'r') as f:
            ljstr = f.read()
        e = []
        s = []
        for line in ljstr.split('\n'):
            words = line.split()
            if len(words) > 0:
                e.append(float(words[2]))
                s.append(float(words[3]))
        n = len(e)
        for i in range(n):
            self._lj.append([])
            for j in range(n):
                self._lj[i].append([math.sqrt(e[i]*e[j]), math.sqrt(s[i]*s[j])])

    def _read_x6(self):
        """Read Exp6 parameters from data4Lammps output file"""
        with open('X6paircoeffs.txt', 'r') as f:
            x6str = f.read()
        for line in x6str.split('\n'):
            words = line.split()
            if len(words) > 0:
                t1 = int(words[1])
                t2 = int(words[2])
                A = float(words[3])
                rho = float(words[4])
                C = float(words[5])
                self._x6.append([t1, t2, A, rho, C])
                           
    def _setup(self):
        n = len(self._lj)
        if n == 0:
            self._read_lj()
            n = len(self._lj)
            s = ''
            fmt = 'pair_coeff  {0:d}  {1:d}  {2:g}  {3:g}\n'
            for i in range(n):
                for j in range(i, n):
                    s += fmt.format(i+1,j+1,self._lj[i][j][0],self._lj[i][j][1])
            self.inputs.add('MINIMIZE_PAIR_COEFFS', 'string', s)
        n = len(self._x6)
        if n == 0:
            self._read_x6()
            n = len(self._x6)
            s = ''
            fmt = 'pair_coeff {0:d} {1:d} {2:g} {3:g} {4:g}\n'
            for c in self._x6:
                s += fmt.format(c[0], c[1], c[2], c[3], c[4])
            self.inputs.add('MD_PAIR_COEFFS', 'string', s)


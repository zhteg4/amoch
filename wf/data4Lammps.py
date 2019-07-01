"""
Run data4Lammps as a Simulator

Copyright

LICENSE
"""

import re

from workflow import Simulator, od_ext

_path = '_polymer'

class data4Lammps(Simulator):

    def __init__(self, echo_output=True):
        Simulator.__init__(self, path='data4Lammps'+od_ext(), \
                           echo_output=echo_output)
        self._command = '{}'.format(self.command)
        with open('structure.name', 'w') as f:
            f.write(_path+'\n')

    def _cleanup(self):
        with open('amoch_cell', 'r') as f:
            s = f.read()
        hi = s.strip().split()
        fi = open(_path+'.data', 'r')
        fo = open('polymer.data', 'w')
        while True:
            line = fi.readline()
            if len(line) == 0:
                break
            if 'xlo' in line:
                fo.write('0.0 {:s} xlo xhi\n'.format(hi[0]))
            elif 'ylo' in line:
                fo.write('0.0 {:s} ylo yhi\n'.format(hi[1]))
            elif 'zlo' in line:
                fo.write('0.0 {:s} zlo zhi\n'.format(hi[2]))
            else:
                fo.write(line)
        fi.close()
        fo.close()

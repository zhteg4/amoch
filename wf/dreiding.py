"""
DreidingTyper class

Copyright

LICENSE
"""

import sys

from workflow import Task, od_ext

class DreidingTyper(Task):
    """Apply Dreiding forcefield types to generic AMOCH atom type output"""

    def __init__(self, echo_output=True):
        Task.__init__(self, path='dreiding'+od_ext(), echo_output=echo_output)
        self._atypes = []

    def _write_atom_type_dat(self):
        """Write atom_type.dat with Dreiding types"""
        with open('atom_type.dat.in', 'r') as fi:
            text = fi.read()
        fo = open('atom_type.dat', 'w')
        lines = text.split('\n')
        for i in range(len(lines)):
            self._atypes.append('')
        for line in lines:
            words = line.split()
            if len(words) < 5:
                continue
            ndx = int(words[0])-1
            sym = words[1]
            nbonds = int(words[2])
            resonant = (words[3] == '1')
            Hbond = (words[4] == '1')
            atype = '{:s}'.format(sym)
            if len(sym) == 1:
                atype += '_'
            if sym not in ['H', 'F', 'Cl', 'Br', 'I', 'Na', 'Ca', 'Fe', 'Zn']:
                # Need hybridization
                if resonant:
                    atype += 'R'
                elif Hbond is True:
                    atype += '_HB'
                elif sym == 'O':
                    atype += '{:d}'.format(nbonds+1)
                else:
                    atype += '{:d}'.format(nbonds-1)
            self._atypes[ndx] = atype
            fo.write('{0:d}  {1:s}\n'.format(ndx+1, atype))
        fo.close()

    def _write_bond_type_dat(self):
        """Write bond_type.dat with Dreiding types"""
        with open('bond_type.dat.in', 'r') as fi:
            text = fi.read()
        fo = open('bond_type.dat', 'w') 
        for line in text.split('\n'):
            words = line.split()
            if len(words) < 4:
                continue
            ndx = int(words[0])-1
            ptype1 = int(words[1])-1
            ptype2 = int(words[2])-1
            # order not used
            fmt = '{0:d}  {1:s}  {2:s}\n'
            fo.write(fmt.format(ndx+1, self._atypes[ptype1], self._atypes[ptype2]))
        fo.close()

    def _run(self):
        """Write atom and bond type files with Dreiding types"""
        out = ''
        err = ''
        self._write_atom_type_dat()
        out += 'Wrote atom_type.dat\n'
        self._write_bond_type_dat()
        out += 'Wrote bond_type.dat\n'
        with open('forcefield.name', 'w') as f:
            f.write('DREIDING\n')
        out += 'Wrote forcefield.name\n'
        self.complete = True
        if self.echo_output is True:
            sys.stdout.write(out)
        return True, out, err


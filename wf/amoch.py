"""
Run AMOCH as a Simulator

Copyright

LICENSE
"""

import os.path, re
from glob import glob

from workflow import od_ext, Simulator, uri2path

def extract_int(s, patt):
    """
    Return an int specified as 'patt n' in the string s; return 0 if patt
    is NOT found in s
    """
    m = re.search(patt+'[\s+]([\d]+)', s)
    return int(m.group(1)) if m is not None else 0

class AMOCH(Simulator):

    def __init__(self, echo_output=True):
        Simulator.__init__(self, path='amoch'+od_ext(), echo_output=echo_output)
        self._iter = 0
        with open(self.templates['amoch.in'], 'r') as f:
            s = f.read()
        n = extract_int(s, 'chains')
        m = extract_int(s, 'monomers')
        b = extract_int(s, 'build')
        self._niters = n*m/b if b > 0 else 1

    def _setup(self):
        """Ensure empty restart files exist for the first iteration"""
        if self._iter == 0:
            self.inputs.add('RESTART', 'string', '')
            for pod in self.inputs.od.values():
                if pod['type'] == 'uri':
                    p = uri2path(pod['value'])
                    if os.path.exists(p) is False:
                        open(p, 'w').close()

    def _cleanup(self):
        """
        After the first iteration, set RESTART for subsequent iterations; 
        add the chain restart files to expected inputs; grab the RNG seed
        for future restarts
        """
        if self._iter == 0:
            self.inputs.set_value('RESTART', 'restart .')
        self._iter += 1
        self.complete = (self._iter == self._niters)
        if self.complete is False:
            for f in glob('chain*.dat'):
                self.inputs.add(f, 'uri', 'file://WORKDIR/'+f)

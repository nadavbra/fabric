from __future__ import absolute_import, division, print_function

import sys
import os
from datetime import datetime
from collections import Counter

import numpy as np

from geneffect.util import as_biopython_seq


### Constants ###

ALL_NTS_LIST = list('ACGT')
ALL_NTS_SET = set(ALL_NTS_LIST)


### Project Functions ###

def log(message):
    print('FABRIC|PID-%s [%s]: %s' % (os.getpid(), datetime.now(), message))
    sys.stdout.flush()

    
### General Helper Classes ###

class TimeMeasure(object):
    
    def __init__(self, opening_statement):
        self.opening_statement = opening_statement
        
    def __enter__(self):
        self.start_time = datetime.now()
        log(self.opening_statement)
        
    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.finish_time = datetime.now()
        self.elapsed_time = self.finish_time - self.start_time
        log('Finished after %s.' % self.elapsed_time)
        
        
### Numpy Helper Functions ###

def get_arrays_diff(array1, array2):
    return np.array([value for value, count in (Counter(array1) - Counter(array2)).items() for _ in range(count)])


### Biopython Helper Functions ###

def get_strand_sequence(seq, strand):
    
    seq = as_biopython_seq(seq)
    
    if strand == '+':
        return seq
    else:
        return seq.reverse_complement()
    
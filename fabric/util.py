import sys
import os
from datetime import datetime

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


### Constants ###

ALL_NTS_LIST = list('ACGT')
ALL_NTS_SET = set(ALL_NTS_LIST)


### Project Functions ###

def log(message, end = '\n'):
    print('FABRIC|PID-%s [%s]: %s' % (os.getpid(), datetime.now(), message), end = end)
    sys.stdout.flush()
    
    
### Statistics ###
    
def multipletests_with_nulls(values, method = 'fdr_bh'):
    
    significance = np.full(len(values), False, dtype = bool)
    qvals = np.full(len(values), np.nan, dtype = float)
    mask = pd.notnull(values)
    
    if mask.any():
        significance[np.array(mask)], qvals[np.array(mask)], _, _ = multipletests(values[mask], method = method)
    
    return significance, qvals

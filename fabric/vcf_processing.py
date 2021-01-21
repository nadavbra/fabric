import re
from collections import Counter

import numpy as np
import pandas as pd
        
def extract_variants(vcf_file_handler, only_pass = False, min_AF = 0):
    
    n_filtered_variants = Counter()
    variants = []
    
    for line in vcf_file_handler:
        if line.startswith('##'):
            continue
        elif line.startswith('#'):
            headers = line[1:].strip().split('\t')
        else:
            
            AFs = np.array(list(map(float, _AF_PATTERN.search(line).group(1).split(','))))
            AF_mask = (AFs >= min_AF)
            
            if AF_mask.any():
                
                variant_data = dict(zip(headers, line.strip().split('\t')))
                
                if not only_pass or variant_data['FILTER'] == 'PASS':
                    
                    alts = np.array(variant_data['ALT'].split(','))
                    
                    for alt, AF in zip(alts[AF_mask], AFs[AF_mask]):
                        variants.append([variant_data['ID'], variant_data['CHROM'], int(variant_data['POS']), \
                                variant_data['REF'], alt, AF, float(variant_data['QUAL']), variant_data['FILTER']])
                else:
                    n_filtered_variants['FILTER is not "PASS"'] += 1
            else:
                n_filtered_variants['AF below %f' % min_AF] += 1
    
    variants = pd.DataFrame(variants, columns = ['id', 'chrom', 'pos', 'ref', 'alt', 'AF', 'qual', 'filter'])
    assert (variants['ref'] != variants['alt']).all()
    return variants, n_filtered_variants

_AF_PATTERN = re.compile('AF=(.*?);')

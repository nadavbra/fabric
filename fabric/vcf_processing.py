from __future__ import absolute_import, division, print_function

import re
from collections import Counter

import numpy as np
import pandas as pd
        
from .util import ALL_NTS_SET, log

def extract_variants(vcf_file_handler, only_pass = False, min_AF = 0, verbose = True):
    
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
                    n_filtered_variants['FILTER != "PASS"'] += 1
            else:
                n_filtered_variants['AF below %f' % min_AF] += 1
    
    variants_df = pd.DataFrame(variants, columns = ['id', 'chrom', 'pos', 'ref', 'alt', 'AF', 'qual', 'filter'])
    
    if verbose:
        total_filtered_variants = sum(n_filtered_variants.values())
        total_variants = len(variants_df) + total_filtered_variants
        log('Extracted %d of %d variants. Filtered out %d variants due to: %s.' % (len(variants_df), total_variants, total_filtered_variants, \
                ', '.join(['%s: %d' % entry for entry in n_filtered_variants.items()])))
    
    assert (variants_df['ref'] != variants_df['alt']).all()
    return variants_df
    
def filter_snps(variants_df, verbose = True):
    
    snps_df = variants_df[variants_df['ref'].isin(ALL_NTS_SET) & variants_df['alt'].isin(ALL_NTS_SET)]
    
    if verbose:
        log('Filtered %d SNPs of %d variants (%.2f%%)' % (len(snps_df), len(variants_df), 100 * len(snps_df) / len(variants_df)))
    
    return snps_df
    
def process_snps(snps_df, setup):
    return snps_df.apply(lambda snp_record: setup.geneffect_setup.variant_interpreter.process_snp(snp_record['chrom'], snp_record['pos'], \
            snp_record['ref'], snp_record['alt']), axis = 1)

_AF_PATTERN = re.compile('AF=(.*?);')

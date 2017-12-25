from __future__ import absolute_import, division, print_function

import pandas as pd

from .util import log, TimeMeasure

def load_relevant_mutations(maf_file_path):
    
    all_mutations = pd.read_csv(maf_file_path, delimiter = '\t')
    log('Loaded %d MAF records.' % len(all_mutations))
    
    relevant_mutations = all_mutations[all_mutations['Variant_Type'] == 'SNP']
    log('Filtered %d of %d (%d%%) relevant mutations (SNPs).' % (len(relevant_mutations), len(all_mutations), 100 * len(relevant_mutations) / \
            len(all_mutations)))
    
    return relevant_mutations
    
def interpret_mutations(mutations, setup):

    def _get_mutation_interpretation(mutation_record):
        # It should hold that: Reference_Allele == Tumor_Seq_Allele1 != Tumor_Seq_Allele2 == Allele
        # Tumor_Seq_Allele1 is the ref_nt, and Tumor_Seq_Allele2 is the alt_nt
        snp = setup.geneffect_setup.variant_interpreter.process_snp(_parse_chrom(mutation_record['Chromosome']), mutation_record['Start_Position'], \
                mutation_record['Tumor_Seq_Allele1'], mutation_record['Tumor_Seq_Allele2'])
        return snp

    with TimeMeasure('Interpreting %d mutations...' % len(mutations)):
        mutations['snp_interpretation'] = mutations.apply(_get_mutation_interpretation, axis = 1)

def _parse_chrom(chrom):
    if chrom.lower().startswith('chr'):
        return chrom[3:]
    else:
        return chrom

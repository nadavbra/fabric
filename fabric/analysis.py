from itertools import product

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.api import DescrStatsW, CompareMeans

from .util import ALL_NTS_LIST, log, multipletests_with_nulls
from .gene_score_dist import get_gene_score_distribution
from .gene_testing import test_gene

def analyze_variant_effects(variant_effects, genes, gene_score_models, other_variant_effects_to_compare_against = None):

    varaint_effects_by_gene = variant_effects.groupby('gene_index')
    
    if other_variant_effects_to_compare_against is not None:
        other_variant_effects_by_gene = other_variant_effects_to_compare_against.groupby('gene_index')
        
    results = []

    for i, (gene_index, gene_variant_effects) in enumerate(varaint_effects_by_gene):
        
        log('Analyzing gene %d/%d...' % (i, len(varaint_effects_by_gene)), end = '\r')
        gene = genes.loc[gene_index]
        gene_score_model = gene_score_models[gene_index]
        overall_bg_score_dist = get_gene_score_distribution(gene_score_model, get_nt_substitution_counts(gene_variant_effects))
        missense_bg_score_dist = get_gene_score_distribution(gene_score_model, \
                get_nt_substitution_counts(gene_variant_effects[gene_variant_effects['effect_type'] == 'missense']))
        gene_results = pd.concat([gene, test_gene(gene['cds_len'], gene_variant_effects, overall_bg_score_dist, missense_bg_score_dist)]).rename(gene_index)
        
        if other_variant_effects_to_compare_against is not None:
            if gene_index in other_variant_effects_by_gene:
                
                gene_other_variant_effects = other_variant_effects_by_gene.get_group(gene_index)
                other_bg_score_dist = get_gene_score_distribution(gene_score_model, get_nt_substitution_counts(gene_other_variant_effects))
                        
                reference_z_scores = (gene_variant_effects['effect_scores'] - overall_bg_score_dist.mean()) / overall_bg_score_dist.std()
                other_z_scores = (gene_other_variant_effects['effect_score'] - other_bg_score_dist.mean()) / other_bg_score_dist.std()
                        
                gene_results['diff_z_scores_ci'] = '%f-%f' % CompareMeans(DescrStatsW(reference_z_scores), DescrStatsW(other_z_scores)).tconfint_diff()
                _, gene_results['diff_pval'] = ttest_ind(reference_z_scores, other_z_scores)
            else:
                gene_results['diff_z_scores_ci'] = np.nan
                gene_results['diff_pval'] = np.nan
            
        results.append(gene_results)
        
    results = pd.DataFrame(results).sort_index()
    results['overall_fdr_significance'], results['overall_fdr_qval'] = multipletests_with_nulls(results['overall_pval'])
    results['missense_fdr_significance'], results['missense_fdr_qval'] = multipletests_with_nulls(results['missense_pval'])
    results['types_chi2_fdr_significance'], results['types_chi2_fdr_qval'] = multipletests_with_nulls(results['types_chi2_pval'])
    return results
    
def get_nt_substitution_counts(variants):
    return variants.groupby(['ref', 'alt']).size().reindex(list(product(ALL_NTS_LIST, ALL_NTS_LIST))).fillna(0).astype(int).unstack()

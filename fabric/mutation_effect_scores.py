from __future__ import absolute_import, division, print_function

import numpy as np
import pandas as pd

from .util import log, TimeMeasure

def calc_effect_scores_of_snps(snps_series, setup):
    
    gene_effects, gene_effects_meta_data = _extract_gene_effects(snps_series, setup.uniprot_id_to_gene_index)
    
    with TimeMeasure('Calculating scores for %d gene effects...' % len(gene_effects)):
        gene_effect_scores = calc_effect_scores_of_gene_effects(gene_effects, setup)

    effect_scores_series = pd.Series([dict() for _ in range(len(snps_series))], index = snps_series.index)

    for gene_effect_score, (snp_index, gene_index) in zip(gene_effect_scores, gene_effects_meta_data):
        effect_scores_series[snp_index][gene_index] = gene_effect_score

    return effect_scores_series
    
def calc_effect_scores_of_gene_effects(gene_effects, setup):

    gene_effects = np.array(gene_effects)

    synonymous_mask = np.array([pd.isnull(gene_effect) or gene_effect.is_synonymous() for gene_effect in gene_effects])
    nonsense_mask = np.array([pd.notnull(gene_effect) and gene_effect.is_nonsense() for gene_effect in gene_effects])
    missense_mask = (~nonsense_mask) & (~synonymous_mask)

    scores = np.empty(len(gene_effects))
    scores[nonsense_mask] = 0
    scores[synonymous_mask] = 1

    log('Applying FIRM classifier on %d missense mutations...' % missense_mask.sum())
    # In FIRM, 0-1 scores indicate harmless-harmful scale, while in FABRIC it is the other way around.
    scores[missense_mask] = 1 - setup.firm_classifier.predict_adjusted_proba(gene_effects[missense_mask], thread_pool = setup.thread_pool)

    return scores
    
def load_effect_scores_series(effect_scores_series_file_path):
    return pd.read_csv(effect_scores_series_file_path, index_col = 0, header = None).iloc[:, 0].apply(_parse_int_to_float_dict)

def _extract_gene_effects(snps_series, uniprot_id_to_gene_index):
        
    gene_effects = []
    gene_effects_meta_data = []

    for snp_index, snp in snps_series.iteritems():
        if pd.notnull(snp):
            for gene_effect in snp.cds_gene_effects:
                gene_index = uniprot_id_to_gene_index[gene_effect.affected_gene.uniprot_record.id]
                gene_effects.append(gene_effect)
                gene_effects_meta_data.append((snp_index, gene_index))
                
    return gene_effects, gene_effects_meta_data
    
def _parse_int_to_float_dict(raw_dict):
    
    result = {}
    
    if raw_dict != '{}':
        for raw_entry in raw_dict[1:-1].split(','):
            raw_key, raw_value = raw_entry.split(':')
            result[int(raw_key.strip())] = float(raw_value.strip())
        
    return result

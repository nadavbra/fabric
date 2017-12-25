from __future__ import absolute_import, division, print_function

import os
from collections import defaultdict
import json

import numpy as np
import pandas as pd

from .util import log, TimeMeasure
from .gene_score_dist import get_gene_score_distribution, create_empty_nt_substitution_matrix, GeneScoreModel
from .gene_score_testing import test_gene

def analyze_vcf_genes(setup, gene_bg_scores_dir, snps, gene_indices_to_analyze = None):

    effect_scores_per_gene, overall_nt_substitution_counts_per_gene, missense_nt_substitution_counts_per_gene = \
            get_effect_scores_and_nt_substitution_counts_per_gene(snps)

    gene_results = []
    gene_indices = []
    
    for j, (gene_index, gene_record) in enumerate(setup.gene_dataset.iterrows()):
        if gene_indices_to_analyze is None or gene_index in gene_indices_to_analyze:
            
            log('Analyzing Gene #%d [%d/%d]...' % (gene_index, j + 1, len(setup.gene_dataset)))
            gene = setup.geneffect_setup.genes[gene_index]
            gene_bg_scores_file_path = os.path.join(gene_bg_scores_dir, setup.geneffect_setup._config_setup.ref_genome, \
                    '%s.json' % gene.uniprot_record.id)
            
            observed_scores = np.array(effect_scores_per_gene[gene_index])
            
            if len(observed_scores) == 0:
                log('No observed mutation scores in gene %s. Skipping.' % gene)
                continue
            
            if os.path.exists(gene_bg_scores_file_path):
                with open(gene_bg_scores_file_path, 'r') as f:
                    gene_score_model = GeneScoreModel(json.load(f))
            else:
                log('No background scores file found for gene %s: %s. Skipping.' % (gene, gene_bg_scores_file_path))
                continue
              
            overall_bg_score_dist = get_gene_score_distribution(gene_score_model, overall_nt_substitution_counts_per_gene[gene_index])
            missense_bg_score_dist = get_gene_score_distribution(gene_score_model, missense_nt_substitution_counts_per_gene[gene_index])
            
            test_results = test_gene(len(gene.canonical_cds_isoform), observed_scores, overall_bg_score_dist, missense_bg_score_dist)
            full_results = pd.concat([gene_record, test_results])
            log('Results: %s' % full_results.to_dict())
            
            gene_results.append(full_results)
            gene_indices.append(gene_index)
    
    return pd.DataFrame(gene_results, index = gene_indices)
            
def get_effect_scores_and_nt_substitution_counts_per_gene(snps):
    
    effect_scores_per_gene = defaultdict(list)
    overall_nt_substitution_counts_per_gene = defaultdict(create_empty_nt_substitution_matrix)
    missense_nt_substitution_counts_per_gene = defaultdict(create_empty_nt_substitution_matrix)
        
    with TimeMeasure('Setting effect scores and nt substitution counts per gene...'):
        for _, snp in snps.iterrows():
            for gene_index, score in snp['effect_scores'].items():
                
                effect_scores_per_gene[gene_index].append(score)
                overall_nt_substitution_counts_per_gene[gene_index].loc[snp['ref'], snp['alt']] += 1
                
                if score > 0 and score < 1:
                    missense_nt_substitution_counts_per_gene[gene_index].loc[snp['ref'], snp['alt']] += 1
                    
    return effect_scores_per_gene, overall_nt_substitution_counts_per_gene, missense_nt_substitution_counts_per_gene

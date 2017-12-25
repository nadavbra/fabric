from __future__ import absolute_import, division, print_function

import os
from collections import defaultdict
import json

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.api import DescrStatsW, CompareMeans

from .util import log, TimeMeasure, get_arrays_diff
from .gene_score_dist import get_gene_score_distribution, create_empty_nt_substitution_matrix, GeneScoreModel
from .gene_score_testing import test_gene

def analyze_maf_genes(setup, gene_bg_scores_dir, mutations, only_combined_project = False, analyze_diff = False, gene_indices_to_analyze = None):

    assert not analyze_diff or not only_combined_project
    
    effect_scores_per_gene_per_project, overall_nt_substitution_counts_per_gene_per_project, missense_nt_substitution_counts_per_gene_per_project = \
            get_effect_scores_and_nt_substitution_counts_per_gene_per_project(mutations, only_combined_project = only_combined_project)
    all_projects = set(effect_scores_per_gene_per_project.keys())
    real_projects = all_projects - set(['combined'])
    
    results_per_project = {}
    
    if analyze_diff:
        results_per_project['diff'] = pd.concat([setup.gene_dataset.copy(), pd.DataFrame(index = setup.gene_dataset.index, \
                columns = sorted(real_projects))], axis = 1)
    
    for i, project in enumerate(all_projects):
    
        project_gene_results = []
        project_analyzed_gene_indices = []
    
        for j, (gene_index, gene_record) in enumerate(setup.gene_dataset.iterrows()):
            if gene_indices_to_analyze is None or gene_index in gene_indices_to_analyze:
                
                log('Analyzing project %s [%d/%d]; Gene #%d [%d/%d]...' % (project, i + 1, len(all_projects), \
                        gene_index, j + 1, len(setup.gene_dataset)))
                gene = setup.geneffect_setup.genes[gene_index]
                gene_bg_scores_file_path = os.path.join(gene_bg_scores_dir, setup.geneffect_setup._config_setup.ref_genome, \
                        '%s.json' % gene.uniprot_record.id)
                
                observed_scores = np.array(effect_scores_per_gene_per_project[project][gene_index])
                
                if len(observed_scores) == 0:
                    log('No observed mutation scores in gene %s. Skipping.' % gene)
                    continue
                
                if os.path.exists(gene_bg_scores_file_path):
                    with open(gene_bg_scores_file_path, 'r') as f:
                        gene_score_model = GeneScoreModel(json.load(f))
                else:
                    log('No background scores file found for gene %s: %s. Skipping.' % (gene, gene_bg_scores_file_path))
                    continue
                  
                overall_bg_score_dist = get_gene_score_distribution(gene_score_model, overall_nt_substitution_counts_per_gene_per_project[project][gene_index])
                missense_bg_score_dist = get_gene_score_distribution(gene_score_model, missense_nt_substitution_counts_per_gene_per_project[project][gene_index])
                
                test_results = test_gene(len(gene.canonical_cds_isoform), observed_scores, overall_bg_score_dist, missense_bg_score_dist)
                full_results = pd.concat([gene_record, test_results])
                log('Results: %s' % full_results.to_dict())
                
                project_gene_results.append(full_results)
                project_analyzed_gene_indices.append(gene_index)
                
                if analyze_diff and project in real_projects:
                    
                    other_scores = get_arrays_diff(effect_scores_per_gene_per_project['combined'][gene_index], observed_scores)
                    other_nt_substitution_counts = overall_nt_substitution_counts_per_gene_per_project['combined'][gene_index] - \
                            overall_nt_substitution_counts_per_gene_per_project[project][gene_index]
                    other_bg_score_dist = get_gene_score_distribution(gene_score_model, other_nt_substitution_counts)
                    
                    project_z_scores = (observed_scores - overall_bg_score_dist.mean()) / overall_bg_score_dist.std()
                    other_z_scores = (other_scores - other_bg_score_dist.mean()) / other_bg_score_dist.std()
                    
                    z_scores_diff_ci = CompareMeans(DescrStatsW(project_z_scores), DescrStatsW(other_z_scores)).tconfint_diff()
                    _, z_scores_diff_pval = ttest_ind(project_z_scores, other_z_scores)
                    results_per_project['diff'].loc[gene_index, project] = json.dumps([z_scores_diff_ci, z_scores_diff_pval])
                
        results_per_project[project] = pd.DataFrame(project_gene_results, index = project_analyzed_gene_indices)
          
    results_per_project['diff'] = results_per_project['diff'][pd.notnull(results_per_project['diff']).any(axis = 1)]
    return results_per_project

def get_effect_scores_and_nt_substitution_counts_per_gene_per_project(mutations, only_combined_project = False):
    
    effect_scores_per_gene_per_project = defaultdict(lambda: defaultdict(list))
    overall_nt_substitution_counts_per_gene_per_project = defaultdict(lambda: defaultdict(create_empty_nt_substitution_matrix))
    missense_nt_substitution_counts_per_gene_per_project = defaultdict(lambda: defaultdict(create_empty_nt_substitution_matrix))
        
    with TimeMeasure('Setting effect scores and nt substitution counts per gene per project...'):
        for _, mutation in mutations.iterrows():
        
            tcga_project = None if only_combined_project else mutation['tcga_project']
            ref_nt = mutation['Tumor_Seq_Allele1']
            alt_nt = mutation['Tumor_Seq_Allele2']
            
            for project in ['combined', tcga_project]:
                if not only_combined_project or project == 'combined':
                    for gene_index, score in mutation['effect_scores'].items():
                        
                        effect_scores_per_gene_per_project[project][gene_index].append(score)
                        overall_nt_substitution_counts_per_gene_per_project[project][gene_index].loc[ref_nt, alt_nt] += 1
                        
                        if score > 0 and score < 1:
                            missense_nt_substitution_counts_per_gene_per_project[project][gene_index].loc[ref_nt, alt_nt] += 1
            
    return effect_scores_per_gene_per_project, overall_nt_substitution_counts_per_gene_per_project, missense_nt_substitution_counts_per_gene_per_project

import numpy as np
import pandas as pd
from scipy.stats import chisquare

def test_gene(cds_len, observed_effects, overall_bg_score_dist, missense_bg_score_dist):
    
    n_obs = len(observed_effects)
    mutations_per_nt = n_obs / cds_len
    meta_results = pd.Series([n_obs, mutations_per_nt], index = ['n_obs', 'mutations_per_nt'])
        
    overall_score_results = test_gene_scores(observed_effects, overall_bg_score_dist, False)
    missense_score_results = test_gene_scores(observed_effects, missense_bg_score_dist, True)
    
    observed_type_counts = observed_effects['effect_type'].value_counts().reindex(['synonymous', 'missense', 'nonsense']).fillna(0).astype(int).values
    expected_type_freqs = overall_bg_score_dist.get_type_freqs()
    
    if (expected_type_freqs == 0).any():
        types_chi2_pval = np.nan
    else:
        _, types_chi2_pval = chisquare(observed_type_counts, observed_type_counts.sum() * expected_type_freqs)
    
    type_results = pd.Series([types_chi2_pval, list(observed_type_counts), list(expected_type_freqs)], \
            index = ['types_chi2_pval', 'observed_type_counts', 'expected_type_freqs'])
    
    return pd.concat([meta_results, overall_score_results, missense_score_results, type_results])
            
def test_gene_scores(observed_effects, bg_score_dist, only_missense):
    
    if only_missense:
        observed_scores = observed_effects.loc[observed_effects['effect_type'] == 'missense', 'effect_score']
    else:
        observed_scores = observed_effects['effect_score']
    
    bg_mean = bg_score_dist.mean(only_missense = only_missense)
    bg_std = bg_score_dist.std(only_missense = only_missense)
    
    observed_avg = observed_scores.mean()
    mean_deviation = observed_avg - bg_mean
    z_value = mean_deviation / bg_std
    pval = calc_two_tailed_binning_pval(observed_scores, bg_score_dist, only_missense)
        
    results = pd.Series([bg_mean, bg_std, observed_avg, mean_deviation, z_value, pval], \
            index = ['bg_mean', 'bg_std', 'observed_avg', 'mean_deviation', 'z_value', 'pval'])
            
    if only_missense:
        results.rename(lambda column: 'missense_%s' % column, inplace = True)
    else:
        results.rename(lambda column: 'overall_%s' % column, inplace = True)
            
    return results

def calc_two_tailed_binning_pval(obs_scores, bg_score_dist, only_missense, n_bins = 100):
    
    if bg_score_dist.is_empty():
        return np.nan
    
    if obs_scores.mean() <= bg_score_dist.mean():
        one_tail_func = calc_left_tail_binning_pval
    else:
        one_tail_func = calc_right_tail_binning_pval
        
    return min(1, 2 * one_tail_func(obs_scores, bg_score_dist, only_missense, n_bins = n_bins))
    
def calc_left_tail_binning_pval(obs_scores, bg_score_dist, only_missense, n_bins = 100):
    obs_sum_bin = (n_bins * obs_scores).astype(int).sum()
    bg_bins = bg_score_dist.to_bins(n_bins, only_missense = only_missense)
    bg_sum_bins = iid_sum_bins(bg_bins, len(obs_scores), max_bin = obs_sum_bin)
    return bg_sum_bins[:(obs_sum_bin + 1)].sum()
    
def calc_right_tail_binning_pval(obs_scores, bg_score_dist, only_missense, n_bins = 100):
    reversed_obs_scores = 1 - obs_scores
    reversed_obs_sum_bin = (n_bins * reversed_obs_scores).astype(int).sum()
    reversed_bg_bins = bg_score_dist.to_bins(n_bins, only_missense = only_missense, reversed = True)
    reversed_bg_sum_bins = iid_sum_bins(reversed_bg_bins, len(obs_scores), max_bin = reversed_obs_sum_bin)
    return reversed_bg_sum_bins[:(reversed_obs_sum_bin + 1)].sum()
    
def iid_sum_bins(bins, n, max_bin = None):
    
    n_as_bits = list(map(int, bin(n)[2:][::-1]))
    
    if n_as_bits[0] == 1:
        sum_bins = bins
    else:
        sum_bins = None
        
    for bit in n_as_bits[1:]:
        
        bins = np.convolve(bins, bins)
        
        if max_bin is not None:
            bins = bins[:(max_bin + 1)]
        
        if bit == 1:
            if sum_bins is None:
                sum_bins = bins
            else:
                
                sum_bins = np.convolve(sum_bins, bins)
                
                if max_bin is not None:
                    sum_bins = sum_bins[:(max_bin + 1)]
                
    return sum_bins

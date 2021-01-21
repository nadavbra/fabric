import numpy as np
from scipy.stats import rv_discrete

def get_gene_score_distribution(gene_score_model, nt_substitution_counts):
    nt_substitution_freqs = nt_substitution_counts / nt_substitution_counts.sum().sum()
    return GeneScoreDistribution(gene_score_model, nt_substitution_freqs)
    
class GeneScoreDistribution(object):

    def __init__(self, gene_score_model, substitution_freqs_matrix):
        self.gene_score_model = gene_score_model
        self.substitution_freqs_matrix = substitution_freqs_matrix
        self._decompose_substitution_freqs()
        
    def is_empty(self):
        return self._empty
        
    def get_type_freqs(self):
        if self._empty:
            return np.array(3 * [np.nan])
        else:
            type_freqs_per_substitution = np.array([substitution_score_dist.get_type_freqs() for substitution_score_dist in self._substitution_score_dists])
            return np.dot(self._freqs.reshape(1, -1), type_freqs_per_substitution).flatten()
        
    def sample(self, n, only_missense = False):
        if self._empty:
            return np.array(n * [np.nan])
        else:
            substitution_counts = np.random.multinomial(n, self._freqs)
            samples = np.concatenate([substitution_score_dist.sample(count, only_missense = only_missense) \
                    for count, substitution_score_dist in zip(substitution_counts, self._substitution_score_dists)])
            np.random.shuffle(samples)
            return samples
            
    def mean(self, only_missense = False):
        if self._empty:
            return np.nan
        else:
            return (self._freqs * self._get_substitution_means(only_missense = only_missense)).sum()
        
    def var(self, only_missense = False):
        if self._empty:
            return np.nan
        else:
            
            mean_of_substitution_vars = (self._freqs * self._get_substitution_vars(only_missense = only_missense)).sum()
            substitution_means_dist = rv_discrete(values = (self._get_substitution_means(only_missense = only_missense), self._freqs))
            var_of_substitution_means = substitution_means_dist.var()
            
            # According the the law of total variance.
            return mean_of_substitution_vars + var_of_substitution_means
        
    def std(self, only_missense = False):
        return np.sqrt(self.var(only_missense = only_missense))
        
    def to_bins(self, n_bins, only_missense = False, reversed = False):
        substitution_bins = [substitution_score_dist.to_bins(n_bins, only_missense = only_missense, reversed = reversed) for \
                substitution_score_dist in self._substitution_score_dists]
        return np.dot(self._freqs, substitution_bins)
        
    def _get_substitution_means(self, only_missense):
        return np.array([substitution_score_dist.mean(only_missense = only_missense) for substitution_score_dist in self._substitution_score_dists])
        
    def _get_substitution_vars(self, only_missense):
        return np.array([substitution_score_dist.var(only_missense = only_missense) for substitution_score_dist in self._substitution_score_dists])
        
    def _decompose_substitution_freqs(self):
        
        substitution_freqs = self.substitution_freqs_matrix.stack()
        substitution_freqs = substitution_freqs[substitution_freqs != 0]
        
        if len(substitution_freqs) == 0:
            self._empty = True
        else:
            self._empty = False
            self._substitutions, self._freqs = zip(*substitution_freqs.to_dict().items())
            self._freqs = np.array(self._freqs)
            self._substitution_score_dists = [self.gene_score_model.substitution_score_distributions[substitution] \
                    for substitution in self._substitutions]
    
class GeneScoreModel(object):

    def __init__(self, variants_df):
        self.substitution_score_distributions = {}
        self._parse_data(variants_df)
    
    def _parse_data(self, variants_df):
        for (ref, alt), substitution_variants in variants_df.groupby(['ref', 'alt']):
            n_total = len(substitution_variants)
            n_synonymous = (substitution_variants['effect_type'] == 'synonymous').sum()
            n_nonsense = (substitution_variants['effect_type'] == 'nonsense').sum()
            n_missense = (substitution_variants['effect_type'] == 'missense').sum()
            missense_scores = substitution_variants.loc[substitution_variants['effect_type'] == 'missense', 'effect_score'].values
            self.substitution_score_distributions[(ref, alt)] = SubstitutionScoreDistribution(n_total, n_synonymous, n_nonsense, \
                    n_missense, missense_scores)

class SubstitutionScoreDistribution(object):

    def __init__(self, n_total, n_synonymous, n_nonsense, n_missense, missense_scores):
        assert n_total == n_synonymous + n_nonsense + n_missense
        assert len(missense_scores) == n_missense
        self.f_synonymous = n_synonymous / n_total
        self.f_nonsense = n_nonsense / n_total
        self.f_missense = n_missense / n_total
        self.missense_scores = np.array(missense_scores)
        
    def get_type_freqs(self):
        return np.array([self.f_synonymous, self.f_missense, self.f_nonsense])
        
    def sample(self, n, only_missense = False):
        if only_missense:
            return self._sample_missense(n)
        else:

            rand = np.random.uniform(0, 1, n)
            synonymous_mask = (rand <= self.f_synonymous)
            missense_mask = ((rand > self.f_synonymous) & (rand <= self.f_synonymous + self.f_missense))
            
            samples = np.zeros(n)
            samples[synonymous_mask] = 1
            samples[missense_mask] = self._sample_missense(missense_mask.sum())
            
            return samples
            
    def mean(self, only_missense = False):
        if only_missense:
            return self.missense_scores.mean()
        else:
            return self.f_synonymous + self.f_missense * self.missense_scores.mean()
        
    def var(self, only_missense = False):
        if only_missense:
            return self.missense_scores.var()
        else:
            
            mean_of_type_vars = self.f_missense * self.missense_scores.var()
            type_means_dist = rv_discrete(values = ([1, self.missense_scores.mean(), 0], self.get_type_freqs()))
            var_of_type_means = type_means_dist.var()
            
            # According the the law of total variance.
            return mean_of_type_vars + var_of_type_means
        
    def std(self, only_missense = False):
        return np.sqrt(self.var(only_missense = only_missense))
        
    def to_bins(self, n_bins, only_missense = False, reversed = False):
        
        bins = self._missense_to_bins(n_bins, reversed = reversed)
        
        if not only_missense:
            
            bins *= self.f_missense
            
            if reversed:
                bins[0] += self.f_synonymous
                bins[-1] += self.f_nonsense
            else:
                bins[0] += self.f_nonsense
                bins[-1] += self.f_synonymous
            
        assert np.isclose(bins.sum(), 1)
        return bins
        
    def _sample_missense(self, n):
        return np.random.choice(self.missense_scores, size = n, replace = True)
        
    def _missense_to_bins(self, n_bins, reversed):
        
        scores = 1 - self.missense_scores if reversed else self.missense_scores
        prob_per_score = 1 / len(scores)
        bin_indices = (n_bins * scores).astype(int)
        bins = np.zeros(n_bins + 1)
        
        for bin_index in bin_indices:
            bins[bin_index] += prob_per_score
            
        return bins

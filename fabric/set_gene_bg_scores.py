from __future__ import absolute_import, division, print_function

import os
import json

from geneffect.variant_processing import build_snp_gene_effect

from .util import ALL_NTS_SET, log, get_strand_sequence

def set_gene_bg_scores(setup, gene_bg_scores_dir, override = False):

    gene_bg_scores_subdir = os.path.join(gene_bg_scores_dir, setup.geneffect_setup._config_setup.ref_genome)
    
    if not os.path.exists(gene_bg_scores_subdir):
        os.makedirs(gene_bg_scores_subdir)

    with setup.async_firm_classifier:
        for gene in setup.geneffect_setup.genes:
        
            gene_bg_scores_file_path = os.path.join(gene_bg_scores_subdir, '%s.json' % gene.uniprot_record.id)
        
            if not override and os.path.exists(gene_bg_scores_file_path):
                log('Gene %s already has existing background scores. Skipping.' % gene.uniprot_record.id)
            else:
                log('Processing gene: %s...' % gene)
                _GeneAnalysis(gene, gene_bg_scores_file_path, setup.async_firm_classifier).process()
    
class _GeneAnalysis(object):

    def __init__(self, gene, result_file_path, async_firm_classifier):
        self.gene = gene
        self.result_file_path = result_file_path
        self.substitution_analyses = dict([self._get_substitution_analysis_entry(ref_nt, alt_nt, async_firm_classifier) \
                for ref_nt in ALL_NTS_SET for alt_nt in (ALL_NTS_SET - set([ref_nt]))])
        self.pending_substitution_analyses = set(self.substitution_analyses.keys())
        
    def process(self):
    
        self._create_all_possible_gene_effects()
                    
        for substitution_analysis in self.substitution_analyses.values():
            substitution_analysis.process()
            
    def _get_substitution_analysis_entry(self, ref_nt, alt_nt, async_firm_classifier):
        substitution = (ref_nt, alt_nt)
        substitution_analysis = _SubstitutionAnalysis(ref_nt, alt_nt, async_firm_classifier, self._get_callback_for_substitution(substitution))
        return (substitution, substitution_analysis)
        
    def _get_callback_for_substitution(self, substitution):
        
        def _callback():
            self._report_substitution_analysis_finished(substitution)
            
        return _callback
        
    def _report_substitution_analysis_finished(self, substitution):
        
        self.pending_substitution_analyses.remove(substitution)
        
        if len(self.pending_substitution_analyses) == 0:
            self._save_result()
            
    def _save_result(self):
        with open(self.result_file_path, 'w') as f:
            json.dump(self._get_result(), f)
        
    def _get_result(self):
        return [substitution_analysis.get_result() for substitution_analysis in self.substitution_analyses.values()]
            
    def _create_all_possible_gene_effects(self):
        
        cds_coordinate = 0
    
        for cds_exon in self.gene.canonical_cds_isoform.cds_exons:
            for cds_ref_nt in cds_exon.seq:
                
                cds_coordinate += 1
                positive_strand_ref_nt = str(get_strand_sequence(cds_ref_nt, cds_exon.strand))
                
                for cds_alt_nt in (ALL_NTS_SET - set(cds_ref_nt)):
                    
                    positive_strand_alt_nt = str(get_strand_sequence(cds_alt_nt, cds_exon.strand))
                    substitution = (positive_strand_ref_nt, positive_strand_alt_nt)
                    gene_effect = build_snp_gene_effect(self.gene, cds_exon, cds_coordinate, cds_ref_nt, cds_alt_nt)
                    
                    if gene_effect is not None:
                        self.substitution_analyses[substitution].gene_effects.append(gene_effect)
            
class _SubstitutionAnalysis(object):
    
    def __init__(self, ref_nt, alt_nt, async_firm_classifier, finish_callback):
        self.ref_nt = ref_nt
        self.alt_nt = alt_nt
        self.async_firm_classifier = async_firm_classifier
        self.finish_callback = finish_callback
        self.gene_effects = []
        
    def process(self):
        
        self.n_total = len(self.gene_effects)
        self.n_synonymous = sum([gene_effect.is_synonymous() for gene_effect in self.gene_effects])
        self.n_nonsense = sum([gene_effect.is_nonsense() for gene_effect in self.gene_effects])
        self.missense_gene_effects = [gene_effect for gene_effect in self.gene_effects if gene_effect.is_missense()]
        
        if len(self.missense_gene_effects) == 0:
            self._process_missense_scores_callback([])
        else:
            self.async_firm_classifier.submit_samples(self.missense_gene_effects, self._process_missense_scores_callback)
                    
    def get_result(self):
        return {
            'ref_nt': self.ref_nt,
            'alt_nt': self.alt_nt,
            'n_total': self.n_total,
            'n_synonymous': self.n_synonymous,
            'n_nonsense': self.n_nonsense,
            'n_missense': len(self.missense_gene_effects),
            'missense_scores': self.missense_scores,
        }
    
    def _process_missense_scores_callback(self, raw_missense_scores):
        # In FIRM, 0-1 scores indicate harmless-harmful scale, while in FABRIC it is the other way around.
        self.missense_scores = [1 - score for score in raw_missense_scores]
        self.finish_callback()

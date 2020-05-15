from __future__ import absolute_import, division, print_function

import os

import pandas as pd

import geneffect
import firm

from .util import log

class Setup(object):
    
    def __init__(self, user_specified_ref_genome, gene_dataset_dir = None, thread_pool = None, n_threads = None):
        
        self.geneffect_setup = geneffect.Setup(user_specified_ref_genome)
        firm.setup_uniprot_tracks(self.geneffect_setup)
        self.firm_classifier = firm.load_classifier(self.geneffect_setup)
        self._set_gene_dataset(gene_dataset_dir)
        
        self.thread_pool = thread_pool
        self.n_threads = n_threads
        self._set_async_firm_classifier()
            
    def _set_gene_dataset(self, gene_dataset_dir):
        if gene_dataset_dir is None:
            self.gene_dataset = None
            self.uniprot_id_to_gene_index = None
        else:
            self.gene_dataset = self._get_or_create_gene_dataset(gene_dataset_dir)
            self.uniprot_id_to_gene_index = {uniprot_id: index for index, uniprot_id in self.gene_dataset['uniprot_id'].iteritems()}
            
    def _set_async_firm_classifier(self):
        if self.thread_pool is None:
            self.async_firm_classifier = None
        else:
            self.async_firm_classifier = firm.AsyncClassifier(self.firm_classifier.predict_adjusted_proba, thread_pool = self.thread_pool, \
                    n_threads = self.n_threads)
            
    def _get_or_create_gene_dataset(self, gene_dataset_dir):
        
        gene_dataset = self._create_gene_dataset()
        gene_dataset_csv_file_path = os.path.join(gene_dataset_dir, 'genes_%s.csv' % self.geneffect_setup._config_setup.ref_genome)
        
        if os.path.exists(gene_dataset_csv_file_path):
            assert (gene_dataset['uniprot_id'] == pd.read_csv(gene_dataset_csv_file_path, index_col = 0)['uniprot_id']).all()
            log('The %d genes in the created gene dataset are identical to the ones in the existing CSV file: %s' % (len(gene_dataset), \
                    gene_dataset_csv_file_path))
        else:
            gene_dataset.to_csv(gene_dataset_csv_file_path)
            log('Saved the gene dataset (%d genes) into the CSV file: %s' % (len(gene_dataset), gene_dataset_csv_file_path))
            
        return gene_dataset
        
    def _create_gene_dataset(self):
        return pd.DataFrame(list(map(_create_gene_record, self.geneffect_setup.genes)), columns = ['uniprot_id', 'symbol', 'name', 'refseq_ids', \
                'chr', 'cds_start', 'cds_end'])
        
def _create_gene_record(gene):
    
    cds_coordinates = [coordinate for cds_exon in gene.canonical_cds_isoform.cds_exons for coordinate in \
            [cds_exon.chromosome_start, cds_exon.chromosome_end]]

    return [
        gene.uniprot_record.id,
        gene.symbol,
        gene.name,
        list(gene.refseq_ids),
        gene.canonical_cds_isoform.chromosome,
        min(cds_coordinates),
        max(cds_coordinates),
    ]

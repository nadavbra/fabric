from __future__ import absolute_import, division, print_function

import os
import argparse
import multiprocessing

from .util import log

def get_parser_file_type(parser, must_exist = False):

    def _file_type(path):
    
        path = os.path.expanduser(path)
    
        if must_exist:
            if not os.path.exists(path):
                parser.error('File doesn\'t exist: %s' % path)
            elif not os.path.isfile(path):
                parser.error('Not a file: %s' % path)
            else:
                return path
        else:
        
            dir_path = os.path.dirname(path)
        
            if not os.path.exists(dir_path):
                parser.error('Parent directory doesn\'t exist: %s' % dir_path)
            else:
                return path
    
    return _file_type

def get_parser_directory_type(parser, create_if_not_exists = False):
    
    def _directory_type(path):
    
        path = os.path.expanduser(path)
    
        if not os.path.exists(path):
            if create_if_not_exists:
                if os.path.exists(os.path.dirname(path)):
                    parser.error('Cannot create empty directory (parent directory doesn\'t exist): %s' % path)
                else:
                    os.mkdir(path)
                    return path
            else:
                parser.error('Path doesn\'t exist: %s' % path)
        elif not os.path.isdir(path):
            parser.error('Not a directory: %s' % path)
        else:
            return path
        
    return _directory_type
    
def get_common_args(n_threads = False):
    
    parser = argparse.ArgumentParser(add_help = False)
    parser.add_argument('--reference-genome', metavar = '<GRCh38 or GRCh37/hg19>', dest = 'ref_genome', type = str, required = True, help = \
            'The human reference genome version to use.')
    parser.add_argument('--gene-dataset-dir', metavar = '/path/to/gene_dataset_dir', dest = 'gene_dataset_dir', type = get_parser_directory_type(parser), \
            required = True, help = 'Path to the directory with the gene dataset (a CSV file). If the dataset file doesn\'t exist, ' + \
            'it will automatically be created.')
            
    if n_threads:
        parser.add_argument('--n-threads', dest = 'n_threads', metavar = '<number>', type = int, default = multiprocessing.cpu_count(), help = \
                'The number of threads to use when running FIRM classifier (by default will use all of the machine\'s CPUs).')
                
    return parser
    
def create_thread_pool(args):
    # NOTE: For memory efficiency, it is recommended to create the thread pool before doing any substantial computation.
    # See: https://stackoverflow.com/questions/14749897/python-multiprocessing-memory-usage
    log('Will use %d threads.' % args.n_threads)
    thread_pool = multiprocessing.Pool(args.n_threads)
    return args.n_threads, thread_pool

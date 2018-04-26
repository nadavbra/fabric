from __future__ import absolute_import, division, print_function

from setuptools import setup
        
def readme():
    with open('README.rst', 'r') as f:
        return f.read()

setup(
    name = 'fabric',
    version = '1.0',
    description = 'FABRIC (Functional Alteration Bias Recovery In Coding-regions) is a framework for detecting genes showing functional alteration bias.',
    long_description = readme(),
    url = 'https://github.com/nadavbra/fabric',
    author = 'Nadav Brandes',
    author_email  ='nadav.brandes@mail.huji.ac.il',
    license = 'MIT',
    packages = ['fabric'],
    scripts = [
        'bin/analyze_maf_genes',
        'bin/analyze_vcf_genes',
        'bin/create_vcf_dataset',
        'bin/set_gene_bg_scores',
        'bin/set_maf_gene_effect_scores',
    ],
    install_requires = [
        'numpy',
        'scipy',
        'pandas',
        'biopython',
        'scikit-learn',
        'statsmodels',
        'geneffect==1.1', # https://github.com/nadavbra/geneffect
        'firm', # https://github.com/nadavbra/firm
    ],
)

from setuptools import setup
        
setup(
    name = 'fabric-genetics',
    version = '1.2.0',
    description = 'FABRIC (Functional Alteration Bias Recovery In Coding-regions) is a framework for detecting genes showing functional alteration bias.',
    long_description = 'See: https://github.com/nadavbra/fabric',
    url = 'https://github.com/nadavbra/fabric',
    author = 'Nadav Brandes',
    author_email  ='nadav.brandes@mail.huji.ac.il',
    license = 'MIT',
    packages = ['fabric'],
    scripts = [
        'bin/fabric',
        'bin/vcf_to_csv',
    ],
    install_requires = [
        'numpy',
        'scipy',
        'pandas',
        'statsmodels',
    ],
)

from setuptools import setup, find_packages

VERSION = '0.1'

DESCRIPTION = 'Module to perform instant similarity (iSIM), which allows for comparison of multiple molecules at the same time, yielding the same value as the average pairwise comparisons of other methods.'

setup(
    name='iSIM',
    version=VERSION,
    description=DESCRIPTION,
    url='https://github.com/mqcomplab/iSIM',
    packages=find_packages(),
    install_requires=[
        'numpy', 
        'matplotlib', 
        'pandas', 
        'rdkit', 
        'scipy', 
        'seaborn', 
        'scikit-learn'
    ]
)


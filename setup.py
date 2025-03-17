from setuptools import setup, find_packages

VERSION = '15.0.0'

DESCRIPTION = 'Module to perform instant similarity (iSIM), which allows for comparison of multiple molecules at the same time, yielding the same value as the average pairwise comparisons of other methods.'

setup(
    name='iSIM',
    version=VERSION,
    description=DESCRIPTION,
    url='https://github.com/mqcomplab/iSIM',
    packages=find_packages(),
    include_package_data=True, 
    install_requires=[
        'numpy', 
        'matplotlib', 
        'pandas', 
        'rdkit', #Not available on PyPI. User might need to install separately
        'scipy', 
        'seaborn', 
        'scikit-learn'
    ], 
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
    ],
    zip_safe=False,
)


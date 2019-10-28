#from distutils.core import setup
from setuptools import setup
from distutils.extension import Extension
import sys

# fix problems with pythons terrible import system
import os
file_dir = os.path.dirname(os.path.realpath(__file__))

SRC_DIR = 'deepDegron'

import deepDegron
version = deepDegron.__version__
AUTHOR = 'Collin Tokheim'
EMAIL = 'fake@gmail.com'
URL = 'https://github.com/ctokheim/deepDegron'
DESCRIPTION = 'deepDegron'
PACKAGES = [SRC_DIR,
            ]
setup(name='deepDegron',
        version=version,
        description=DESCRIPTION,
        author=AUTHOR,
        author_email=EMAIL,
        url=URL,
        packages=PACKAGES,
        license='Apache License, Version 2.0',
        install_requires=['numpy', 'scipy', 'pandas', 'pysam'],
        package_data={
        },
        entry_points={
            'console_scripts':[
                'train_cterm_model = deepDegron.train_cterm_nn:cli_main',
                'train_nterm_model = deepDegron.train_nterm_nn:cli_main',
                'deepDegron_test = deepDegron.statistical_test:cli_main',
            ]
        },
        long_description=open('README.md').read(),
        classifiers=['Topic :: Scientific/Engineering :: Bio-Informatics',
                    'Environment :: Console',
                    'Intended Audience :: Developers',
                    'Intended Audience :: Science/Research'],
)

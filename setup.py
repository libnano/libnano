#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
test with:
    python setup.py build_ext --inplace

create wheels with wheel installed:

python setup.py bdist_wheel
python setup.py sdist --formats=gztar


"""

DESCRIPTION = ("Low-level nucleic acid sequence manipulation tools")
with open('README.rst') as fd:
    LONG_DESCRIPTION = fd.read()

DISTNAME = 'libnano'
LICENSE = 'GPLv2'
AUTHORS = "Nick Conway, Ben Pruitt"
EMAIL = "nick.conway@wyss.harvard.edu"
URL = "https://github.com/libnano/libnano"
DOWNLOAD_URL = ''
CLASSIFIERS = [
    'Development Status :: 4 - Beta',
    'Environment :: Console',
    'Intended Audience :: Science/Research',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Cython',
    'Topic :: Scientific/Engineering',
    'License :: OSI Approved :: GNU General Public License v2 (GPLv2)'
]

try:
    from setuptools import (
        setup,
        Extension
    )
except ImportError:
    from distutils.core import (
        setup,
        Extension
    )
from Cython.Build import cythonize
import numpy.distutils.misc_util
import os
import re
import sys
import ast
import shutil

pjoin = os.path.join
rpath = os.path.relpath

# Begin modified code from Flask's version getter
# BSD license
# Copyright (c) 2015 by Armin Ronacher and contributors.
# https://github.com/pallets/flask
# _version_re = re.compile(r'__version__\s+=\s+(.*)')
_version_re = re.compile(r'__version__\: str\s+=\s+(.*)')
with open('libnano/__init__.py', 'rb') as initfile:
    version = str(ast.literal_eval(_version_re.search(
                                   initfile.read().decode('utf-8')).group(1)))
# end Flask derived code

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ platform info ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
if sys.platform == 'win32':
    extra_compile_args = ['']
else:
    extra_compile_args = ['-Wno-unused-function']

# ~~~~~~~~~~~~~~~~~~~~~~~~~~ package / module paths ~~~~~~~~~~~~~~~~~~~~~~~~~ #
PACKAGE_PATH = os.path.abspath(os.path.dirname(__file__))
MODULE_PATH = pjoin(PACKAGE_PATH, 'libnano')
# unused for now
# DATASETS_PATH =     pjoin(MODULE_PATH, 'datasets')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ include dirs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
common_include = ['libnano/src', 'libnano/helpers']

# Non-python files to include in the installation
libnano_files = []

# ~~~~~~~~~~~~~~~~~~~~ configure normal c api extensions ~~~~~~~~~~~~~~~~~~~~ #
normal_extensions = []

# ~~~~~~~~~~~~~~~~~~~~~~~ configure cython extensions ~~~~~~~~~~~~~~~~~~~~~~~ #
cython_extensions = []

# ~~~~~~~~~~~~~~~~~~~ add helpers pyx, pxd, and py filess ~~~~~~~~~~~~~~~~~~~ #
helpers_fp = pjoin(MODULE_PATH, 'helpers')

libnano_files += [rpath(pjoin(helpers_fp, f)) for f in
                  os.listdir(helpers_fp) if ('.py' in f or '.pyx' in f or '.pxd' in f)]

def addExtension(*ext_args, **ext_kwargs):
    global libnano_files
    libnano_files += ext_kwargs['sources']
    cython_extensions.append(Extension(*ext_args, **ext_kwargs))

addExtension(
    'libnano.metric.seqrepeat',
    depends=[],
    sources=['libnano/src/si_seqint.c',
             'libnano/src/sr_seqrepeat.c',
             'libnano/metric/seqrepeat.pyx'],
    include_dirs=common_include + [numpy.get_include()],
    extra_compile_args=extra_compile_args
)

addExtension(
    'libnano.search.seedfinder',
    sources=['libnano/search/seedfinder.pyx'],
    extra_compile_args=extra_compile_args
)

addExtension(
    'libnano.cymem.cymem',
    sources=['libnano/cymem/cymem.pyx'],
    extra_compile_args=extra_compile_args
)

addExtension(
    'libnano.search.seedmatch',
    sources=['libnano/search/seedmatch.pyx',
             'libnano/src/shl_seedhashlist.c',
             'libnano/src/ss_seqstr.c'],
    include_dirs=common_include,
    extra_compile_args=extra_compile_args
)

addExtension(
    'libnano.search.submerpool',
    sources=['libnano/search/submerpool.pyx'],
    include_dirs=common_include,
    extra_compile_args=extra_compile_args
)

addExtension(
    'libnano.seqint',
    depends=[],
    sources=[   'libnano/seqint.pyx',
                'libnano/src/si_seqint.c'],
    include_dirs=common_include,
    extra_compile_args=extra_compile_args
)
libnano_files.append('libnano/seqint.pxd')
libnano_files.append('libnano/datastructures/list_bisect.pxd')

addExtension(
    'libnano.seqstr',
    depends=[],
    sources=['libnano/seqstr.pyx',
             'libnano/src/ss_seqstr.c'],
    include_dirs=common_include + [numpy.get_include()],
    extra_compile_args=extra_compile_args
)
libnano_files.append('libnano/seqstr.pxd')

addExtension(
    'libnano.seqgraph',
    depends=[],
    sources=['libnano/seqgraph.pyx',
             'libnano/src/ss_seqstr.c'],
    include_dirs=common_include + [numpy.get_include()],
    extra_compile_args=extra_compile_args
)

addExtension(
    'libnano.metric.seqscreen',
    depends=[],
    sources=['libnano/src/si_seqint.c',
             'libnano/src/sf_seqscreen.c',
             'libnano/metric/seqscreen.pyx'],
    include_dirs=common_include,
    extra_compile_args=extra_compile_args
)

addExtension(
    'libnano.metric.seqmetric',
    depends=[],
    sources=['libnano/src/sm_seqmetric.c',
             'libnano/metric/seqmetric.pyx'],
    include_dirs=common_include,
    extra_compile_args=extra_compile_args
)

addExtension(
    'libnano.search.restriction',
    depends=[],
    sources=['libnano/search/restriction.pyx'],
    include_dirs=common_include,
    extra_compile_args=extra_compile_args
)
# Restriction dataset(s)
res_data_fp = pjoin(MODULE_PATH, 'datasets')
libnano_files += [rpath(pjoin(root, f), MODULE_PATH) for root, _, files in
                  os.walk(res_data_fp) for f in files if ('.json' in f or '.yaml' in f)]

addExtension(
    'libnano.datastructures.seqrecord.feature',
    depends=[],
    sources=['libnano/datastructures/seqrecord/feature.pyx'],
    include_dirs=common_include,
    extra_compile_args=extra_compile_args
)

addExtension(
    'libnano.datastructures.seqrecord.seqrecordbase',
    depends=[],
    sources=['libnano/datastructures/seqrecord/seqrecordbase.pyx'],
    include_dirs=common_include,
    extra_compile_args=extra_compile_args,
)

addExtension(
    'libnano.simplethermo',
    depends=[],
    sources=['libnano/simplethermo.pyx'],
    include_dirs=common_include,
    extra_compile_args=extra_compile_args,
)

# add header files or extra c files
for path in common_include:
    libnano_files += [rpath(pjoin(path, f)) for f in
                  os.listdir(path) if ('.h' in f or '.c' in f)]

# ~~~~~~~~~~~~~~~~~~~~~ Strip `libnano/` for package data ~~~~~~~~~~~~~~~~~~~~ #
check_libnano = 'libnano/'
lcl = len(check_libnano)
for i, f in enumerate(libnano_files):
    if f.startswith(check_libnano):
        libnano_files[i] = f[lcl:]
# unique the list
libnano_files = list(set(libnano_files))

packages = ['libnano', 'libnano.fileio',
            'libnano.helpers', 'libnano.cymem', 'libnano.search',
            'libnano.datastructures', 'libnano.datastructures.seqrecord',
            'libnano.datasets', 'libnano.metric']

# Commented out by NC 2018.01.05 since we are rolling towards PyPi
script_args = sys.argv[1:]

# ~~~~~~~~~~~~~~~~~~~ remove old built files if specified ~~~~~~~~~~~~~~~~~~~ #
def removeBuiltFiles():
    ''' Remove any existing *.c or *.so files from a previous build process
    '''
    print('Removing previously built files...')
    avoid_dirs = ['src', '\.git', 'build', 'klib', 'scratch']
    ads = '|'.join(['(?:%s)' % d for d in avoid_dirs])
    for root, dirs, files in os.walk(MODULE_PATH):
        if not re.search(ads, root):
            for f in files:
                if f.endswith('.c') or f.endswith('.so') or f.endswith('.pyd'):
                    print('Removing %s' % pjoin(MODULE_PATH, root, f))
                    os.remove(pjoin(MODULE_PATH, root, f))

if '--rmbuilt' in script_args:
    removeBuiltFiles()
    script_args.remove('--rmbuilt')

# ~~~~~~~~~~~~~~~~~~~~~~~ cythonize cython extensions ~~~~~~~~~~~~~~~~~~~~~~~ #
cython_ext_list = cythonize(cython_extensions, include_path=common_include)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ rock and roll ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
install_requires = ['six',
                    'cython',
                    'jinja2',
                    'numpy',
                    'PyYAML',
                    'requests',
                    'primer3-py',
                    'click'
                    ]

setup(
    name=DISTNAME,
    version=version,
    maintainer=AUTHORS,
    packages=packages,
    ext_modules=normal_extensions + cython_ext_list,
    include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs()+common_include,
    package_data={'libnano': libnano_files},
    maintainer_email=EMAIL,
    description=DESCRIPTION,
    license=LICENSE,
    url=URL,
    download_url=DOWNLOAD_URL,
    long_description=LONG_DESCRIPTION,
    classifiers=CLASSIFIERS,
    script_args=script_args,
    zip_safe=False,
    install_requires=install_requires,
    entry_points='''
        [console_scripts]
        geneprober=libnano.scripts.geneprober:cli
    '''
)

#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
=========================================================================
libnano: Low-level nucleic acid sequence manipulation tools
=========================================================================

nucleic acid toolkit

'''

import os
import os.path as op
import re
import sys

DESCRIPTION = 'Low-level nucleic acid sequence manipulation tools'
with open('README.rst') as fd:
    LONG_DESCRIPTION = fd.read()

DISTNAME = 'libnano'
LICENSE = 'GPLv2'
AUTHORS = 'Nick Conway, Ben Pruitt'
EMAIL = 'nick.conway@wyss.harvard.edu'
URL = 'https://github.com/libnano/libnano'
DOWNLOAD_URL = ''
CLASSIFIERS = [
    'Development Status :: 4 - Beta',
    'Environment :: Console',
    'Intended Audience :: Science/Research',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.8',
    'Programming Language :: Python :: 3.9',
    'Programming Language :: Python :: 3.10',
    'Programming Language :: Python :: 3.11',
    'Programming Language :: Cython',
    'Topic :: Scientific/Engineering',
    'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
]


from setuptools import (
    Extension,
    setup,
)

try:
    import numpy as np  # type: ignore
    NUMPY_INCLUDE = [np.get_include()]
except ImportError:
    NUMPY_INCLUDE = [np.get_include()]


def try_cythonize(extension_list, *args, **kwargs):
    '''
    Light cythonize wrapper
    '''
    try:
        from Cython.Build import cythonize  # type: ignore
    except (ImportError, ModuleNotFoundError):
        def cythonize(x, **kwargs):
            return x

    cython_compiler_directives = dict(
        language_level='3',
        c_string_encoding='utf-8',
        # embedsignature=True,
        binding=True,
    )

    return cythonize(
        extension_list,
        compiler_directives=cython_compiler_directives,
    )


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ platform info ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
if sys.platform == 'win32':
    EXTRA_COMPILE_ARGS = ['']
else:
    EXTRA_COMPILE_ARGS = ['-Wno-unused-function']

# ~~~~~~~~~~~~~~~~~~~~~~~~~~ package / module paths ~~~~~~~~~~~~~~~~~~~~~~~~~ #
PACKAGE_PATH = op.abspath(op.dirname(__file__))
MODULE_PATH = op.join(PACKAGE_PATH, 'libnano')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ include dirs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
COMMON_INCLUDE = ['libnano/src', 'libnano/helpers']
# Non-python files to include in the installation
LIBNANO_FILES = []

# ~~~~~~~~~~~~~~~~~~~~~~~ configure cython extensions ~~~~~~~~~~~~~~~~~~~~~~~ #
CYTHON_EXTENSIONS = []

# ~~~~~~~~~~~~~~~~~~~ add helpers pyx, pxd, and py filess ~~~~~~~~~~~~~~~~~~~ #
HELPERS_FP = op.join(MODULE_PATH, 'helpers')

LIBNANO_FILES += [
    op.relpath(op.join(HELPERS_FP, f)) for f in
    os.listdir(HELPERS_FP) if ('.py' in f or '.pyx' in f or '.pxd' in f)
]

NO_NUMPY_CYTHON_MACRO_WARNING = [
    ('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION'),
]


def add_extension(
        *ext_args,
        **ext_kwargs,
):
    global LIBNANO_FILES
    LIBNANO_FILES += ext_kwargs['sources']
    CYTHON_EXTENSIONS.append(
        Extension(*ext_args, **ext_kwargs),
    )


add_extension(
    'libnano.metric.seqrepeat',
    depends=[],
    sources=[
        'libnano/src/si_seqint.c',
        'libnano/src/sr_seqrepeat.c',
        'libnano/metric/seqrepeat.pyx',
    ],
    include_dirs=COMMON_INCLUDE + NUMPY_INCLUDE,
    define_macros=NO_NUMPY_CYTHON_MACRO_WARNING,
    extra_compile_args=EXTRA_COMPILE_ARGS,
)

add_extension(
    'libnano.search.seedfinder',
    sources=['libnano/search/seedfinder.pyx'],
    extra_compile_args=EXTRA_COMPILE_ARGS,
)

add_extension(
    'libnano.cymem.cymem',
    sources=['libnano/cymem/cymem.pyx'],
    extra_compile_args=EXTRA_COMPILE_ARGS,
)

add_extension(
    'libnano.search.seedmatch',
    sources=[
        'libnano/search/seedmatch.pyx',
        'libnano/src/shl_seedhashlist.c',
        'libnano/src/ss_seqstr.c',
    ],
    include_dirs=COMMON_INCLUDE,
    extra_compile_args=EXTRA_COMPILE_ARGS,
)

add_extension(
    'libnano.search.submerpool',
    sources=['libnano/search/submerpool.pyx'],
    include_dirs=COMMON_INCLUDE,
    extra_compile_args=EXTRA_COMPILE_ARGS,
)

add_extension(
    'libnano.seqint',
    depends=[],
    sources=[
        'libnano/seqint.pyx',
        'libnano/src/si_seqint.c',
    ],
    include_dirs=COMMON_INCLUDE,
    extra_compile_args=EXTRA_COMPILE_ARGS,
)
LIBNANO_FILES.append('libnano/seqint.pxd')
LIBNANO_FILES.append('libnano/list_bisect.pxd')

add_extension(
    'libnano.seqstr',
    depends=[],
    sources=[
        'libnano/seqstr.pyx',
        'libnano/src/ss_seqstr.c',
    ],
    include_dirs=COMMON_INCLUDE + NUMPY_INCLUDE,
    extra_compile_args=EXTRA_COMPILE_ARGS,
    define_macros=NO_NUMPY_CYTHON_MACRO_WARNING,
)
LIBNANO_FILES.append('libnano/seqstr.pxd')

add_extension(
    'libnano.seqgraph',
    depends=[],
    sources=[
        'libnano/seqgraph.pyx',
        'libnano/src/ss_seqstr.c',
    ],
    include_dirs=COMMON_INCLUDE + NUMPY_INCLUDE,
    extra_compile_args=EXTRA_COMPILE_ARGS,
    define_macros=NO_NUMPY_CYTHON_MACRO_WARNING,
)

add_extension(
    'libnano.metric.seqscreen',
    depends=[],
    sources=[
        'libnano/src/si_seqint.c',
        'libnano/src/sf_seqscreen.c',
        'libnano/metric/seqscreen.pyx',
    ],
    include_dirs=COMMON_INCLUDE,
    extra_compile_args=EXTRA_COMPILE_ARGS,
)

add_extension(
    'libnano.metric.seqmetric',
    depends=[],
    sources=[
        'libnano/src/sm_seqmetric.c',
        'libnano/metric/seqmetric.pyx',
    ],
    include_dirs=COMMON_INCLUDE,
    extra_compile_args=EXTRA_COMPILE_ARGS,
)

add_extension(
    'libnano.search.restriction',
    depends=[],
    sources=['libnano/search/restriction.pyx'],
    include_dirs=COMMON_INCLUDE,
    extra_compile_args=EXTRA_COMPILE_ARGS,
)
# Restriction dataset(s)
res_data_fp = op.join(MODULE_PATH, 'datasets')
LIBNANO_FILES += [
    op.relpath(op.join(root, f), MODULE_PATH) for root, _, files in
    os.walk(res_data_fp) for f in files if ('.json' in f or '.yaml' in f)
]

add_extension(
    'libnano.seqrecord.feature',
    depends=[],
    sources=['libnano/seqrecord/feature.pyx'],
    include_dirs=COMMON_INCLUDE,
    extra_compile_args=EXTRA_COMPILE_ARGS,
)

add_extension(
    'libnano.seqrecord.seqrecordbase',
    depends=[],
    sources=['libnano/seqrecord/seqrecordbase.pyx'],
    include_dirs=COMMON_INCLUDE,
    extra_compile_args=EXTRA_COMPILE_ARGS,
)

add_extension(
    'libnano.simplethermo',
    depends=[],
    sources=['libnano/simplethermo.pyx'],
    include_dirs=COMMON_INCLUDE,
    extra_compile_args=EXTRA_COMPILE_ARGS,
)


# add header files or extra c files
for path in COMMON_INCLUDE:
    LIBNANO_FILES += [
        op.relpath(op.join(path, f)) for f in
        os.listdir(path) if ('.h' in f or '.c' in f)
    ]

# ~~~~~~~~~~~~~~~~~~~~~ Strip `libnano/` for package data ~~~~~~~~~~~~~~~~~~~~ #
check_libnano = 'libnano/'
lcl = len(check_libnano)
for i, f in enumerate(LIBNANO_FILES):
    if f.startswith(check_libnano):
        LIBNANO_FILES[i] = f[lcl:]

# Unique the list
LIBNANO_FILES = list(set(LIBNANO_FILES))

PACKAGES_LIST = [
    'libnano',
    'libnano.fileio',
    'libnano.helpers',
    'libnano.cymem',
    'libnano.search',
    'libnano.seqrecord',
    'libnano.scripts',
    'libnano.datasets',
    'libnano.metric',
]

# Commented out by NC 2018.01.05 since we are rolling towards PyPi
SCRIPT_ARGS = sys.argv[1:]

# ~~~~~~~~~~~~~~~~~~~ remove old built files if specified ~~~~~~~~~~~~~~~~~~~ #


def remove_built_files():
    '''Remove any existing *.c or *.so files from a previous build process
    '''
    print('Removing previously built files...')
    avoid_dirs = ['src', '.git', 'build', 'klib', 'scratch']
    ads = '|'.join(['(?:%s)' % d for d in avoid_dirs])
    for root, dirs, files in os.walk(MODULE_PATH):
        if not re.search(ads, root):
            for f in files:
                if f.endswith('.c') or f.endswith('.so') or f.endswith('.pyd'):
                    print('Removing %s' % op.join(MODULE_PATH, root, f))
                    os.remove(op.join(MODULE_PATH, root, f))


if '--rmbuilt' in SCRIPT_ARGS:
    remove_built_files()
    SCRIPT_ARGS.remove('--rmbuilt')

# ~~~~~~~~~~~~~~~~~~~~~~~ cythonize cython extensions ~~~~~~~~~~~~~~~~~~~~~~~ #
CYTHON_EXT_LIST = try_cythonize(
    CYTHON_EXTENSIONS,
    include_path=COMMON_INCLUDE,
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ rock and roll ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
INSTALL_REQUIRES = [
    'cython>=0.29.0',
    'numpy>=1.24.0',
    'PyYAML>=6.0',
    'requests>=2.28.1',
    'primer3-py>=1.0.0',
    'click>=8.1.0',
    'ssw-py>=1.0.0',
    'pygments',
    'pandas',
]
if sys.platform == 'win32':
    INSTALL_REQUIRES.append('colorama')

import libnano

setup(
    name=DISTNAME,
    version=libnano.__version__,
    maintainer=AUTHORS,
    packages=PACKAGES_LIST,
    ext_modules=CYTHON_EXT_LIST,
    include_dirs=(
        NUMPY_INCLUDE + COMMON_INCLUDE
    ),
    package_data={'libnano': LIBNANO_FILES},
    maintainer_email=EMAIL,
    description=DESCRIPTION,
    license=LICENSE,
    url=URL,
    download_url=DOWNLOAD_URL,
    long_description=LONG_DESCRIPTION,
    classifiers=CLASSIFIERS,
    script_args=SCRIPT_ARGS,
    zip_safe=False,
    install_requires=INSTALL_REQUIRES,
    entry_points='''
        [console_scripts]
        geneprober=libnano.scripts.geneprober:cli
    ''',
)

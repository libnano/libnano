# libnano

[![Build Status](https://travis-ci.org/libnano/libnano.svg)](https://travis-ci.org/libnano/libnano)

nucleic acid toolkit

Low level Python modules for working with DNA sequences.
The goal is to be fast so most everthing is written in Cython or C.

This is a work in progress. Uses klib `khash.h` for faster C hash tables

## Features:

* `fileio`: genbank and fasta reader/writers, xmfa
* `seqrecord`: less featured replacement for Biopython SeqRecord
* `seqstr`: Reverse Complement and Hamming Distance finders
* `seqgraph`: find cliques of compatible sequences
* `seqint`: convert sequences less than 32 bases to binary
* `seqscreen`: filter sequences based on things like GC content
* `seqsearch`: search sequences for features like restriction sites or
               "submers" (kmers with 0 or more mismatches to a subsequence
               in the target)

## Development Features:

* Build derivative Cython modules with libnano using `libnano.includes()` in
your include_dir in `distutils.core.setup`
* use `cynja` to do templating in Cython
* incorporates a [cymem](https://github.com/syllog1sm/cymem) derivative to help
with memory management.  Adds some features like taking ownership of memory to
help get you in trouble :-)

## Building libnano:

* Run `python setup.py build_ext --inplace` to build in the package dir
* Optional flag "--dev" builds pre-release extensions
* Optional flag "--rmbuilt" removes old *.c and *.so files from previous builds

## git subtree for klib
see [subtree setup](https://blogs.atlassian.com/2013/05/alternatives-to-git-submodule-git-subtree/)

    git subtree add --prefix klib git@github.com:attractivechaos/klib.git master --squash
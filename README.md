# libnano

[![Build Status](https://travis-ci.org/libnano/libnano.svg)](https://travis-ci.org/libnano/libnano)

nucleic acid toolkit

Low level Python modules for working with DNA sequences.<br/>
It was developed to be used as a cog in application specific libraries while avoiding bloat.<br/>
The goal is to be fast so most everthing is written in Cython or C.

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

## Building libnano:

* Run `python setup.py build_ext --inplace` to build in the package dir
* Optional flag "--rmbuilt" removes old *.c and *.so files from previous builds

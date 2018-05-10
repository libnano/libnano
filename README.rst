===========================================================
libnano: Low-level nucleic acid sequence manipulation tools
===========================================================

**libnano** is beta

.. image:: https://secure.travis-ci.org/libnano/libnano.png
  :target: https://travis-ci.org/libnano/libnano
.. image:: https://img.shields.io/pypi/l/libnano.png
  :target: http://www.gnu.org/licenses/gpl-2.0.html
.. image:: https://img.shields.io/pypi/v/libnano.png
  :target: https://pypi.python.org/pypi/libnano

.. _repo: https://github.com/libnano/libnano

nucleic acid toolkit

| Low level Python modules for working with DNA sequences
| It was developed to be used as a cog in application specific libraries while avoiding bloat.<br/>
| The goal is to be fast so most everthing is written in Cython or C.

Features
========

- **fileio**: genbank and fasta reader/writers, xmfa
- **seqrecord**: less featured replacement for Biopython SeqRecord
- **seqstr**: Reverse Complement and Hamming Distance finders
- **seqgraph**: find cliques of compatible sequences
- **seqint**: convert sequences less than 32 bases to binary
- **seqscreen**: filter sequences based on things like GC content
- **seqsearch**: search sequences for features like restriction sites or
               "submers" (kmers with 0 or more mismatches to a subsequence
               in the target)

Installation
============

If you would like to install libnano in your local Python environment
you may use ``pip``::

  $ pip install libnano

or build from source from Github repo_ ::

  $ python setup.py install

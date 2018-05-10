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
| The goal is to be fast so most many things are written in Cython or C.
| This project is under active development, and hopefully documentation will by done soon
| Python 3.6+ only now as the code base is typed for Cython and type-hinted for Python code

Features
========

- **seqstr**: Reverse Complement and Hamming Distance finders
- **seqgraph**: find cliques of compatible sequences
- **seqint**: convert sequences less than 32 bases to binary
- **seqscreen**: filter sequences based on things like GC content
- **ensemblrest**: Use the Ensembl REST API to look up gene and transcript
  information
- **prostr**: work with converting DNA to RNA to AA and back again
- **barcode_tools**: generate DNA barcode sets with length and Hamming distance
  restrictions
- **seqsearch**: search sequences for features like restriction sites or
  "submers" (kmers with 0 or more mismatches to a subsequence in the target)
- **fileio**: genbank and fasta reader/writers, xmfa
- **seqrecord**: less featured replacement for Biopython SeqRecord

Installation
============

If you would like to install libnano in your local Python environment
you may use ``pip``::

  $ pip install libnano

or build from source from Github repo_ ::

  $ python setup.py install

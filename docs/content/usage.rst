Usage
-----

Getting started
~~~~~~~~~~~~~~~

Linkapy is a Python package that is designed to facilitate the integrative analysis of single-cell multi-omics data, where multi-omics means multiple read outs of the same cell.
While an attempt is made to keep Linkapy as general as possible, for now it primarily focuses on data that includes methylation and transcription layers.
If available, transcriptome data should come as one (or more) featureCount tables, that will be just combined into a single matrix.
Methylation data should be provided in the form of 'allcools' files, which are tab-separated files containing methylation information for each cytosine in the genome.
Support for additional formats is planned.

Usage
~~~~~

Linkapy can be used through the command line, or via the API in Python. For the latter, have a look at the `API Reference <../autoapi/index.html>`_.

.. click:: linkapy.CLI:linkapy
   :prog: linkapy CLI
   :nested: full

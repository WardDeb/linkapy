Usage
-----

Getting started
~~~~~~~~~~~~~~~

The starting point for linkapy are processed scNMT-seq data files. 
These include one (or more) featureCount tables, and two types of methylation files ('WGCN' and 'GCHN') representing accessibility and methylation, respectively.
Note that at this moment, linkapy is supposed to be used via the API only.
The required input data can be created from raw fastq-files using the nextflow pipeline provided by the lab of functional epigenetics (KU Leuven, prof. Thienpont).

Creating Matrices
~~~~~~~~~~~~~~~~~

The first step in the analysis is to create matrices from the input data.
This is done by initiating a Parse_scNMT class:

.. code-block:: python

    from linkapy import Parse_scNMT
    parser = Parse_scNMT(METHpath='./', RNApath='./', project_name='scNMT', opath=None, threads=10)
    parser.create_matrices()

Note that the paths will be searched recursively, and look for:

 - RNA : '*gene.tsv'
 - accessibility: '*WCGN*.allc.tsv.gz'
 - methylation: '*GCHN*.allc.tsv.gz'

Upon completion, the matrices will be saved in the specified output directory (opath), as parquet files.

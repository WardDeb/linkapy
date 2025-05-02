Usage
-----

Getting started
~~~~~~~~~~~~~~~

The starting point for linkapy are processed scNMT-seq data files. 
These include one (or more) featureCount tables, and two types of methylation files ('WCGN' and 'GCHN') representing methylation and accessibility, respectively.
Note that at this moment, linkapy can be used over API only (CLI is planned).
The required input data can be created from raw fastq-files using the nextflow pipeline provided by the lab of functional epigenetics (KU Leuven, prof. Thienpont).


TLDR
~~~~

Create muData object from nextflow pipeline output (i.e. featureCounts, allcools).
The methylation data can be parsed in bins (default 10 000 bp long), or several region files can be provided (in bed or bed.gz format).
A qc=True can be added in case you want an aggregated matrix of scaled regions (100bps in total) over all the cells, for QC purposes.
Note that the QC runs serially, and thus can take quite some time depending on the number of cells included, and the number of regions that are covered.

.. code-block:: python

    from linkapy.parsing import Parse_scNMT, Parse_matrices
    parser = Parse_scNMT(
        methpath='./', rnapath='./',
        opath='output_dir', threads=10,
        chromsizes='path_to_chromsizes'
    )
    parser.create_matrices()
    muCreator = Parse_matrices
    muCreator.create_mudata()


.. code-block:: python

    from linkapy.parsing import Parse_scNMT, Parse_matrices
    parser = Parse_scNMT(
        methpath='./', rnapath='./',
        opath='output_dir', threads=10,
        regions = ['bedpath1.bed', 'bvedpath2.bed.gz'],
        qc = True
    )
    parser.create_matrices()
    muCreator = Parse_matrices
    muCreator.create_mudata()

Note that the paths (methpath, rnapath) will be searched recursively, and looks for:

 - RNA : '\*gene.tsv'
 - methylation: '\*WCGN\*.allc.tsv.gz'
 - accessibility: '\*GCHN\*.allc.tsv.gz'

Upon completion, the matrices (arrow format for RNA, mtx format for methylation/accessibility) in the output folder.
Additionaly, a muData object is created there as well, which can be used for further analysis.

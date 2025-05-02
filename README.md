# linkapy
A framework that encompasses a number of techniques applicable to scNMT-seq analysis.

# Installation

Note that a working (rustup)[https://rustup.rs/] installation is required. 
Secondly, if an appropriate python (>= 3.10, < 3.13) and pip version are installed, linkapy can be installed with:

  > pip install .  

Or alternatively, with a working maturin installation:

  > maturin develop --release


# Todo

## core-todo
 - [x] Installeable
 - [x] nextflow output -> anndata
 - [ ] sparse_miss
 - [ ] bins
 - [ ] CNV ~ 'matrix' (cell x regions -> colname)
 - [ ] make modalities optional
 - [ ] exclusion -> blacklist (bed)
 - [ ] annotations
 - [ ] integrations

## meta-todo
 - [x] readthedocs setup
 - [ ] API exposure
 - [ ] CLI
 - [x] actions - lint, action
 - [ ] testdata
 - [ ] test - parse
 - [ ] test - anndata
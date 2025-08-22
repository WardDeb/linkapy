from pathlib import Path
import signal
from rich.console import Console
from rich.logging import RichHandler
from linkapy.linkapy import parse_cools
import logging
import anndata as ad
import mudata as md
import numpy as np
import scipy as sp
import polars as pl
import pandas as pd
from Levenshtein import distance as ls_dist
from difflib import SequenceMatcher
from typing import List

class Linkapy_Parser:
    '''
    Linkapy_Parser mainly functions to create matrices (arrow format for RNA, mtx format for accessibility / methylation)
    from directories containing analyzed scNMT-seq data. Theoretically this could be any type of multi-modal (read: RNA / methylation) data, but the class is written with the scNMT workflow
    from the Thienpont lab (KU Leuven) in mind.
    
    At least one of both items should be provided:
     - methylation_path and/or transcriptome_path
     - regions or chromsizes file (if methylation_path is provided).

    :param str methylation_path: The path to the methylation directory (will be searched recursively!).
    :param str transcriptome_path: The path to the RNA output directory (will be searched recursively!).
    :param str output: The output directory where matrices will be written to. Defaults to current working directory in folder ('linkapy_output').
    :param tuple methylation_pattern: The glob pattern to search methylation path recursively. Defaults to ('GC'). Note that this is a tuple.
    :param tuple transcriptome_pattern: The glob pattern to search transcriptome path recursively. Defaults to ('tsv'). Note that this is a tuple.
    :param bool NOMe: If set, methylation_path will be searched for NOMe-seq data. The methylation path will be searched for patterns ('GCHN', 'WCGN').
    :param int threads: Number of threads to use for parsing. Defaults to 1.
    :param str chromsizes: Path to the chromsizes file for the genome. If set, methylation signal will be aggregated over bins
    :param tuple regions: Path or paths to bed files containing regions to aggregate methylation signal over. Can be gzipped. Note that this is a tuple.
    :param tuple blacklist: Path or paths to bed files containing regions to exclude from the aggregation. Can be gzipped. Note that this is a tuple.
    :param int binsize: Size of the bins to aggregate over. Only relevant if no regions are provided. Defaults to 10000.
    :param str project: Name of the project. Will be treated as a prefix for the output files. Defaults to 'linkapy'.
    '''
    def __init__(
        self, 
        methylation_path=None, 
        transcriptome_path=None, 
        output='linkapy_output',  
        methylation_pattern=('*GC*tsv.gz',),
        methylation_pattern_names=(),
        transcriptome_pattern=('*tsv',), 
        transcriptome_pattern_names=(),
        NOMe=False, 
        threads=1, 
        chromsizes=None, 
        regions=None, 
        blacklist=None, 
        binsize=10000, 
        project='linkapy',
        verbose=False
    ):
        self.output = Path(output)
        self.output.mkdir(parents=True, exist_ok=True)
        self.project = project

        # Set up log
        self.logfile = self.output / f'{self.project}.log'
        self.logger = logging.getLogger()
        rich_handler = RichHandler(rich_tracebacks=True, show_time=False, show_level=True, show_path=False)
        _fmt = logging.Formatter('%(levelname)s - %(asctime)s - %(message)s', datefmt="%H:%M:%S")
        file_handler = logging.FileHandler(self.logfile)
        
        file_handler.setFormatter(_fmt)
        # Set verbosity
        if verbose:
            self.logger.setLevel(logging.DEBUG)
            file_handler.setLevel(logging.DEBUG)
            rich_handler.setLevel(logging.DEBUG)
        else:
            self.logger.setLevel(logging.INFO)
            file_handler.setLevel(logging.INFO)
            rich_handler.setLevel(logging.INFO)
            

        self.logger.addHandler(rich_handler)
        self.logger.addHandler(file_handler)

        console = Console()
        console.rule("[bold green]Linkapy Parser[/bold green]")
        self.logger.info(f"Logging under {self.logfile}")

        # Check parameters
        if not any((methylation_path, transcriptome_path)):
            self.logger.error("No methylation_path or transcriptome_path provided. Exiting.")
            raise ValueError("Missing either transcritpome or methylation path")
        if methylation_path and not any((chromsizes, regions)):
            self.logger.error("Methylation data requires either a chromsizes file or at least one regions file.")
            raise ValueError("Missing regions or chromsizes")
        if chromsizes and regions:
            self.logger.warning("Both chromsizes and regions provided. Chromsizes will be ignored.")
            chromsizes = None
        if methylation_pattern and not methylation_pattern_names and not NOMe:
            self.logger.info("No methylation pattern names provided. The asterisks will be stripped from the patterns to yield labels.")
            methylation_pattern_names = tuple([p.replace('*', '') for p in methylation_pattern])
        if transcriptome_pattern and not transcriptome_pattern_names:
            self.logger.info("No transcriptome pattern names provided. The asterisks will be stripped from the patterns to yield labels.")
            transcriptome_pattern_names = tuple([p.replace('*', '') for p in transcriptome_pattern])
        if NOMe:
            self.logger.info("NOMe flag set. Methylation pattern will be set to ('*GCHN*.tsv.gz', '*WCGN*tsv.gz'), names to ('Acc', 'Meth').")
            methylation_pattern = ('*GCHN*tsv.gz', '*WCGN*.tsv.gz')
            methylation_pattern_names = ('Acc', 'Meth')
        
        # Set up paths
        self.methylation_path = Path(methylation_path) if methylation_path else None
        self.transcriptome_path = Path(transcriptome_path) if transcriptome_path else None
        self.chromsizes = Path(chromsizes) if chromsizes else None
        self.regions = [Path(r) for r in regions] if regions else []
        self.blacklist = [Path(b) for b in blacklist] if blacklist else []
        
        # settings and flags
        self.threads = threads
        self.methylation_pattern = methylation_pattern
        self.methylation_pattern_names = methylation_pattern_names
        self.transcriptome_pattern = transcriptome_pattern
        self.transcriptome_pattern_names = transcriptome_pattern_names
        self.binsize = binsize

        self.logger.debug("Linkapy Parser set up. Parameters:")
        for k, v in self.__dict__.items():
            self.logger.debug(f"{k}: {v}")

        # Validate paths and files.
        self._validate()
        # Discover files to aggregate.
        self._glob()
    
    def _validate(self):
        '''
        Validate the provided paths and parameters.
        '''
        self.logger.info("Validating files and paths.")
        if self.methylation_path and not self.methylation_path.exists():
            self.logger.error(f"Methylation path {self.methylation_path} does not exist.")
            raise FileNotFoundError(f"Methylation path {self.methylation_path} does not exist.")
        if self.transcriptome_path and not self.transcriptome_path.exists():
            self.logger.error(f"Transcriptome path {self.transcriptome_path} does not exist.")
            raise FileNotFoundError(f"Transcriptome path {self.transcriptome_path} does not exist.")
        if self.chromsizes and not self.chromsizes.exists():
            self.logger.error(f"Chromsizes file {self.chromsizes} does not exist.")
            raise FileNotFoundError(f"Chromsizes file {self.chromsizes} does not exist.")
        for r in self.regions:
            if not r.exists():
                self.logger.error(f"Region file {r} does not exist.")
                raise FileNotFoundError(f"Region file {r} does not exist.")
        for b in self.blacklist:
            if not b.exists():
                self.logger.error(f"Blacklist file {b} does not exist.")
                raise FileNotFoundError(f"Blacklist file {b} does not exist.")

        if self.methylation_path:
            if not self.methylation_pattern:
                self.logger.error("No methylation pattern provided. Exiting.")
                raise ValueError("Missing methylation pattern")
            for _ in self.methylation_pattern:
                if '*' not in _:
                    self.logger.warning(f"Methylation pattern {_} doesn't contain an asterisk. Are you sure this is what you want ?")
        if self.transcriptome_path:
            if not self.transcriptome_pattern:
                self.logger.error("No transcriptome pattern provided. Exiting.")
                raise ValueError("Missing transcriptome pattern")
            for _ in self.transcriptome_pattern:
                    if '*' not in _:
                        self.logger.warning(f"Transcriptome pattern {_} doesn't contain an asterisk. Are you sure this is what you want ?")
        
    def _glob(self):
        '''
        Discover files to aggregate over based on the paths and patterns provided.
        '''
        self.logger.info("Globbing files.")
        # If methylation_path is provided, there is at least one pattern (as per validate).
        # The asterisks are stripped form the patterns, and used as keys in a dictionary to keep the globs.
        self.transcriptome_files = {}
        self.methylation_files = {}
        if self.methylation_path:
            for pattern, name in zip(self.methylation_pattern, self.methylation_pattern_names):
                _ = list(self.methylation_path.rglob(pattern))
                assert any(_), f"No files found for pattern \'{pattern}\' in {self.methylation_path}"
                self.methylation_files[name] = _
                self.logger.info(f"Methylation search - pattern \'{pattern}\' - name \'{name}\' = {len(_)} files found.")
        if self.transcriptome_path:
            for pattern in self.transcriptome_pattern:
                _ = list(self.transcriptome_path.rglob(pattern))
                assert any(_), f"No files found for pattern {pattern} in {self.transcriptome_path}"
                self.transcriptome_files[pattern.replace('*', '')] = _
                self.logger.info(f"Transcriptome search - pattern \'{pattern}\' = {len(_)} files found.")

    def parse(self):
        '''
        Parse the globbed files and create the different matrices and their corresponding metadata.
        '''
        self.logger.info("Start parsing files.")
        # RNA
        if self.transcriptome_files:
            (self.output / 'matrices').mkdir(parents=True, exist_ok=True)
            self.logger.info("Parsing RNA files.")
            for pattern, files in self.transcriptome_files.items():
                self.logger.info(f"Parsing RNA pattern \'{pattern}\' - {len(files)} files")
                _prefix = self.output / 'matrices' / f'RNA_{pattern}'
                # Don't rerun if the files already exist.
                _countf = _prefix.with_name(_prefix.name + "_counts.arrow")
                _metaf = _prefix.with_name(_prefix.name + "_meta.arrow")
                if not (_countf.exists() and _metaf.exists()):
                    parse_rna(files, _prefix)
                    self.logger.info(f"RNA pattern \'{pattern}\' written with {_prefix}")
                else:
                    self.logger.info(f"RNA pattern \'{pattern}\' files already exist.")
        
        # Methylation
        if self.methylation_files:
            (self.output / 'matrices').mkdir(parents=True, exist_ok=True)
            self.logger.info("Parsing methylation files.")
            for pattern, files in self.methylation_files.items():
                self.logger.info(f"Parsing methylation pattern \'{pattern}\' - {len(files)} files")
                _prefix = self.output / 'matrices' / f"METH_{pattern}"
                _original_handler = signal.getsignal(signal.SIGINT)
                signal.signal(signal.SIGINT, signal.SIG_DFL)
                _region_labels = [r.name for r in self.regions] if self.regions else []
                parse_cools(
                    [str(i) for i in files],
                    [str(i) for i in self.regions] if self.regions else [],
                    [str(i) for i in self.blacklist] if self.blacklist else [],
                    _region_labels,
                    self.threads,
                    str(_prefix),
                    str(self.chromsizes) if self.chromsizes else 'none',
                    self.binsize if self.binsize else 0,
                )
                signal.signal(signal.SIGINT, _original_handler)

        self.logger.info("Creating MuData object.")
        self.dump_mudata()

    def dump_mudata(self):
        _adatas = []
        _patterns = []
        #if self.transcriptome_files:
        if self.transcriptome_files:
            for pattern in self.transcriptome_files:
                self.logger.info(f"Creating anndata object for \'{pattern}\'")
                _adatas.append(read_rna_to_anndata(self.output / 'matrices' / f'RNA_{pattern}'))
                _patterns.append(f'RNA_{pattern}')
                self.logger.info(f"anndata object for \'{pattern}\' with shape {_adatas[-1].shape}")
        if self.methylation_files:
            for pattern in self.methylation_files:
                self.logger.info(f"Creating anndata object for \'{pattern}\'")
                _adatas.append(read_meth_to_anndata(self.output / 'matrices' / f'METH_{pattern}'))
                _patterns.append(f'METH_{pattern}')
                self.logger.info(f"anndata object for \'{pattern}\' with shape {_adatas[-1].shape}")
        _cells = [_adatas.obs.index.tolist() for _adatas in _adatas]
        self.logger.info(f"{len(_adatas)} anndata objects in total.")
        if len(_cells) > 1:
            self.logger.info("Attempt to match cells across different anndata objects.")
            renamed_obs, rename_df = match_cells(_cells, _patterns)
            if renamed_obs:
                self.logger.info("Matching of cells across anndata objects successfull.")
                rename_df.to_csv(self.output / 'cell_renaming.tsv', sep='\t', index=False)
                self.logger.info(f"Dataframe used to rename cells written to {str(self.output / 'cell_renaming.tsv')}.")
                for new_obs, _ad in zip(renamed_obs, _adatas):
                    _ad.obs.index = new_obs
        self.logger.info("Saving MuData object.")
        md.set_options(pull_on_update=False)
        mudata = md.MuData(
            {
                pattern: _ad for pattern, _ad in zip(_patterns, _adatas)
            }
        )
        mudata.write(self.output / f"{self.project}.h5mu")
        self.logger.info(f"MuData object written to {self.output / f'{self.project}.h5mu'}")

def parse_rna(files, prefix) -> None:
    '''
    Read one or more featureCount files, combine them and write them to a counts and metadata arrow file.
    '''
    metacols = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
    schema = {
        'Geneid': pl.String,
        'Chr': pl.String,
        'Start': pl.UInt32,
        'End': pl.UInt32,
        'Strand': pl.String,
        'Length': pl.UInt32
    }

    metadfs = []
    countdfs = []
    for _f in files:
        df = pl.read_csv(_f, separator='\t', skip_rows=1, has_header=True)
        _schema = schema.copy()
        for sample in df.columns:
            if sample not in _schema:
                _schema[sample] = pl.UInt32
        df = df.select(
            [pl.col(col).cast(_schema[col]) for col in df.columns]
        )
        metadf = df.select(metacols)
        countdf = df.select(pl.exclude(metacols))
        countdf.columns = [c.split('.')[0] for c in countdf.columns]
        metadfs.append(metadf)
        countdfs.append(countdf)
    if len(metadfs) > 1:
        assert all(metadfs[0].equals(df) for df in metadfs[1:])
    # concatenate the countdfs (horizontal)
    countdf = pl.concat(countdfs, how='horizontal')
    countdf.write_ipc(prefix.with_name(prefix.name + "_counts.arrow"), compression='zstd')
    metadfs[0].write_ipc(prefix.with_name(prefix.name + "_meta.arrow"), compression='zstd')


def read_rna_to_anndata(prefix) -> ad.AnnData:
    '''
    From a prefix, read the count matrix, and the metadata, combine them into an AnnData object.
    '''
    _counts = pl.read_ipc(prefix.with_name(prefix.name + "_counts.arrow"), memory_map=False).to_pandas()
    _meta = pl.read_ipc(prefix.with_name(prefix.name + "_meta.arrow"), memory_map=False).to_pandas()
    _meta.index = _meta['Geneid']
    del _meta['Geneid']
    annd = ad.AnnData(
        X=sp.sparse.csr_matrix(_counts.values.T),
        obs=pd.DataFrame(index=list(_counts.columns)),
        var=_meta
    )
    return annd[annd.obs.sort_index().index, :].copy()

def read_meth_to_anndata(prefix) -> ad.AnnData:
    '''
    From a prefix, read the fraction matrices, and their metadata, and combine them into an AnnData object.
    '''
    np.seterr(divide='ignore', invalid='ignore')
    methp = prefix.with_name(prefix.name + ".frac.mtx")
    cellp = prefix.with_name(prefix.name + ".cells.tsv")
    regp = prefix.with_name(prefix.name + ".regions.tsv")
    X = sp.io.mmread(methp).tocsr()
    
    _obs = pl.read_csv(cellp, separator='\t', has_header=False).to_pandas()
    _obs = pd.DataFrame(index=[Path(i).name.split('.')[0] for i in _obs['column_1']])
    _var = pl.read_csv(regp, separator='\t', has_header=True).to_pandas()
    _var.index = _var.index.astype(str)
    annd = ad.AnnData(
        X=X,
        obs=_obs,
        var=_var
    )
    return annd[annd.obs.sort_index().index, :].copy()

def match_cells(_l: List[List[str]], patterns: List[str]) -> tuple[List[List[str]], pd.DataFrame]|tuple[None, None]:
    '''
    Take a list of lists containing putative cell names. Per list, we need a 'best match'.
    This is needed since often an assay or context specific pre- or postfix is used, and we want to match them for the mudata object.
    '''
    a = pd.DataFrame(np.nan, index=range(len(_l[0])), columns=range(len(_l)), dtype="string")
    a[0] = _l[0]
    # Start col_index at 1 since the first column is for the ref cells.
    col_index = 1
    for cell_list in _l[1:]:
        row_index = 0
        for refcell in _l[0]:
            distances = [ls_dist(refcell, cell) for cell in cell_list]
            top_matches = [cell for cell, dist in zip(cell_list, distances) if dist == min(distances)]
            if len(top_matches) > 1 or not top_matches:
                return (None, None)
            a.at[row_index, col_index] = top_matches[0]
            row_index += 1
        col_index += 1
    
    # If any cellname is duplicated inside a column, return None
    for col in a.columns:
        if a[col].duplicated().any():
            return (None, None)
    
    a['common_name'] = a.apply(lambda row: get_common_cellname(row.values), axis=1)
    # Construct a 'renamed' nested list so that the obs can be renamed.
    renl = []
    for pos in range(len(a.columns) - 1):
        _tlis = []
        for _i,r in a.iterrows():
            if pd.isna(r['common_name']):
                _tlis.append(r[patterns[pos]])
            else:
                _tlis.append(r['common_name'])
        renl.append(_tlis)
    a.columns = patterns + ["common_name"]
    return (renl, a)


def get_common_cellname(cellnames: List[str]) -> str:
    ref = cellnames[0]
    for name in cellnames[1:]:
        sm = SequenceMatcher(None, ref, name)
        match = sm.find_longest_match(0, len(ref), 0, len(name))
        ref = ref[match.a: match.a + match.size]
        if not ref:
            return np.nan
    return ref.strip("_.-")

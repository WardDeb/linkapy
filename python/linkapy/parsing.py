from pathlib import Path
import polars as pl
import signal
from rich.logging import RichHandler
from linkapy.linkapy import parse_cools
import logging

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
    :param bool mudata: If set, a MuData object will be created from the matrices and written to the output folder. Defaults to False.
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
        mudata=False, 
        methylation_pattern=('*GC*',), 
        transcriptome_pattern=('*tsv',), 
        NOMe=False, 
        threads=1, 
        chromsizes=None, 
        regions=None, 
        blacklist=None, 
        binsize=10000, 
        project='linkapy',
    ):
        self.output = Path(output)
        self.output.mkdir(parents=True, exist_ok=True)
        self.project = project

        # Set up log
        self.logfile = self.output / f'{self.project}.log'
        self.logger = logging.getLogger()
        self.logger.setLevel(logging.INFO)
        rich_handler = RichHandler(rich_tracebacks=True, show_time=False, show_level=True, show_path=False)
        rich_handler.setLevel(logging.DEBUG)
        _fmt = logging.Formatter('%(levelname)s - %(asctime)s - %(message)s', datefmt="%H:%M:%S")
        file_handler = logging.FileHandler(self.logfile)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(_fmt)

        self.logger.addHandler(rich_handler)
        self.logger.addHandler(file_handler)
        self.logger.info(f"Linkapy Parser - project {self.project}")
        self.logger.info(f"Logging under {self.logfile}")

        # Check parameters
        if not methylation_path and not transcriptome_path:
            self.logger.error("No methylation_path or transcriptome_path provided. Exiting.")
            raise ValueError("Missing either transcritpome or methylation path")
        if methylation_path and (not chromsizes and not regions):
            self.logger.error("Methylation data requires either a chromsizes file or at least one regions file.")
            raise ValueError("Missing regions or chromsizes")
        if chromsizes and regions:
            self.logger.warning("Both chromsizes and regions provided. Chromsizes will be ignored.")
            chromsizes = None
        if NOMe:
            self.logger.info("NOMe flag set. Methylation pattern will be set to ('*GCHN*', '*WCGN*')")
            methylation_pattern = ('*GCHN*', '*WCGN*')
        
        # Set up paths
        self.methylation_path = Path(methylation_path) if methylation_path else None
        self.transcriptome_path = Path(transcriptome_path) if transcriptome_path else None
        self.chromsizes = Path(chromsizes) if chromsizes else None
        self.regions = [Path(r) for r in regions] if regions else []
        self.blacklist = [Path(b) for b in blacklist] if blacklist else []
        
        # settings and flags
        self.threads = threads
        self.mudata = mudata
        self.methylation_pattern = methylation_pattern
        self.transcriptome_pattern = transcriptome_pattern
        self.binsize = binsize

        # Validate paths and files.
        self._validate()
        # Discover files to aggregate.
        self._glob()
    
    def _validate(self):

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
        self.logger.info("Globbing files.")
        # If methylation_path is provided, there is at least one pattern (as per validate).
        # The asterisks are stripped form the patterns, and used as keys in a dictionary to keep the globs.
        
        if self.methylation_path:
            self.methylation_files = {}
            for pattern in self.methylation_pattern:
                _ = list(self.methylation_path.rglob(pattern))
                assert any(_), f"No files found for pattern \'{pattern}\' in {self.methylation_path}"
                self.methylation_files[pattern.replace('*', '')] = _
                self.logger.info(f"Methylation search - pattern \'{pattern}\' = {len(_)} files found.")
        if self.transcriptome_path:
            self.transcriptome_files = {}
            for pattern in self.transcriptome_pattern:
                _ = list(self.transcriptome_path.rglob(pattern))
                assert any(_), f"No files found for pattern {pattern} in {self.transcriptome_path}"
                self.transcriptome_files[pattern.replace('*', '')] = _
                self.logger.info(f"Transcriptome search - pattern \'{pattern}\' = {len(_)} files found.")

    def parse(self):
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
                    _region_labels,
                    self.threads,
                    _prefix.name,
                    str(self.chromsizes) if self.chromsizes else 'none',
                    self.binsize if self.binsize else 0,
                )
                signal.signal(signal.SIGINT, _original_handler)


def parse_rna(files, prefix):
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
        metadfs.append(metadf)
        countdfs.append(countdf)
    if len(metadfs) > 1:
        assert all(metadfs[0].equals(df) for df in metadfs[1:])
    # concatenate the countdfs (horizontal)
    countdf = pl.concat(countdfs, how='horizontal')
    countdf.write_ipc(prefix.with_name(prefix.name + "_counts.arrow"), compression='zstd')
    metadfs[0].write_ipc(prefix.with_name(prefix.name + "_meta.arrow"), compression='zstd')


# class Parse_matrices:
#     '''
#     Parses matrices created previously (with Parse_scNMT) and creates a muon object (written to disk).
#     This is then the starting point to downstream analysis.

#     :param str matrixdir: Directory where the matrices can be found. (required)
#     :param str project: Project name, similar to the one provided in the Parse_scNMT part. Defaults to 'scNMT' (optional)
#     :param str ofile: Name of the output file, defaults to matrixdir / project.h5. (optional)
#     '''

#     def __init__(self, matrixdir, project='scNMT', opath=None):
#         self.matrixdir = Path(matrixdir)
#         self.project = project
#         if not opath:
#             self.opath = self.matrixdir
#         else:
#             self.opath = Path(opath)
#         self.opath.mkdir(parents=True, exist_ok=True)
#         # Initiate a log
#         timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
#         _fmt = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

#         # Logfile
#         log_file = Path(self.opath / f"lap_Parse_matrices_{timestamp}_{uuid.uuid4().hex}.log")
#         logging.basicConfig(level=logging.INFO)
#         self.logger = logging.getLogger(str(log_file))
#         file_handler = logging.FileHandler(log_file)
#         # To file
#         file_handler.setFormatter(_fmt)
#         self.logger.addHandler(file_handler)
#         self.logger.propagate = False

#         _msg(self.logger, "Output directory: " + str(self.opath))
#         self._assert_files()

#     def _assert_files(self):
#         self.rna = self.matrixdir / (self.project + '_rnadf.arrow')
#         self.rnameta = self.matrixdir / (self.project + '_rnameta.arrow')
#         # Accessibility part
#         self.acc = {
#             'cov': self.matrixdir / (self.project + '.acc.cov.mtx'),
#             'meth': self.matrixdir / (self.project + '.acc.meth.mtx'),
#             'site': self.matrixdir / (self.project + '.acc.site.mtx'),
#             'cell': self.matrixdir / (self.project + '.acc.cell.tsv'),
#             'reg': self.matrixdir / (self.project + '.acc.meta.tsv')
#         }
#         self.meth = {
#             'cov': self.matrixdir / (self.project + '.meth.cov.mtx'),
#             'meth': self.matrixdir / (self.project + '.meth.meth.mtx'),
#             'site': self.matrixdir / (self.project + '.meth.site.mtx'),
#             'cell': self.matrixdir / (self.project + '.meth.cell.tsv'),
#             'reg': self.matrixdir / (self.project + '.meth.meta.tsv')
#         }

#         assert self.rna.exists()
#         assert self.rnameta.exists()
#         for _f in self.acc:
#             assert self.acc[_f].exists()
#         for _f in self.meth:
#             assert self.meth[_f].exists()


#     def create_mudata(self, aggtype='fraction'):
#         # Some settings.
#         np.seterr(divide='ignore', invalid='ignore')
#         md.set_options(pull_on_update=False)

#         _msg(self.logger, "-"*100 + "\n" + "Parse_matrices - create_mudata" + "\n" + "-"*100)
#         _msg(self.logger, "Parsing RNA matrices")
#         rnadf = pl.read_ipc(self.rna, memory_map=False).to_pandas()
#         rnameta = pl.read_ipc(self.rnameta, memory_map=False).to_pandas()
#         rnameta.index = rnameta['Geneid']
#         del rnameta['Geneid']
#         rna_adata = anndata.AnnData(
#             X=sp.sparse.csr_matrix(rnadf.values.T),
#             obs=pd.DataFrame(index= [i.replace('_RNA-seq', '').replace("_RNA", "") for i in rnadf.columns]),
#             var=rnameta
#         )
#         _msg(self.logger, "Parsing RNA matrices")
#         _msg(self.logger, f"adata for rna shape = {rna_adata.shape}")
#         _msg(self.logger, "Parsing Accessibility matrices")
#         _m = sp.io.mmread(self.acc['meth']).todense()
#         _c = sp.io.mmread(self.acc['cov']).todense()
#         if aggtype == 'fraction':
#             X = np.zeros_like(_m, dtype=float)
#             X = np.where((_c != 0) & (_m != 0), _m / _c, 0)
#             X = sp.sparse.csr_matrix(X)
#         elif aggtype == 'sum':
#             X = sp.sparse.csr_matrix(_m)
#         _obs = pl.read_csv(self.acc['cell'], separator='\t', has_header=False).to_pandas()
#         cell_names = []
#         for i in _obs['column_1'].to_list():
#             _n = Path(i).name.replace('.GCHN-Both.allc.tsv.gz', '').replace("_NOMe-seq", "").replace("_METH", "")
#             cell_names.append(_n)
#         _obs = pd.DataFrame(index = cell_names)
#         _var = pl.read_csv(self.acc['reg'], separator='\t', has_header=True).to_pandas()
#         # Since there is no check for duplicated values in the regions. We deduplicate here (by name)
#         _var = _var.drop_duplicates(subset='name', keep='first')
#         X = X[:, _var.index]
#         _var = _var.set_index('name')
#         acc_adata = anndata.AnnData(
#             X=X,
#             obs=_obs,
#             var=_var
#         )
#         _msg(self.logger, f"adata for accessibility shape = {acc_adata.shape}")
#         _msg(self.logger, "Parsing Methylation matrices")
#         _m = sp.io.mmread(self.meth['meth']).todense()
#         _c = sp.io.mmread(self.meth['cov']).todense()
#         if aggtype == 'fraction':
#             X = np.zeros_like(_m, dtype=float)
#             X = np.where((_c != 0) & (_m != 0), _m / _c, 0)
#             X = sp.sparse.csr_matrix(X)
#         elif aggtype == 'sum':
#             X = sp.sparse.csr_matrix(_m)
#         _obs = pl.read_csv(self.meth['cell'], separator='\t', has_header=False).to_pandas()
#         cell_names = []
#         for i in _obs['column_1'].to_list():
#             _n = Path(i).name.replace('.WCGN-Both.allc.tsv.gz', '').replace("_NOMe-seq", "").replace("_METH", "")
#             cell_names.append(_n)
#         _obs = pd.DataFrame(index = cell_names)
#         _var = pl.read_csv(self.meth['reg'], separator='\t', has_header=True).to_pandas()
#                 # Since there is no check for duplicated values in the regions. We deduplicate here (by name)
#         _var = _var.drop_duplicates(subset='name', keep='first')
#         X = X[:, _var.index]
#         _var = _var.set_index('name')
#         meth_adata = anndata.AnnData(
#             X=X,
#             obs=_obs,
#             var=_var
#         )
#         _msg(self.logger, f"adata for methylation shape = {meth_adata.shape}")
#         # muData
#         _msg(self.logger, "Creating muData object.")
#         # Take intersection of all obs
#         _msg(self.logger, f"First observations for RNA data = {rna_adata.obs_names[:5]}")
#         _msg(self.logger, f"First observations for ACC data = {acc_adata.obs_names[:5]}")
#         _msg(self.logger, f"First observations for METH data = {meth_adata.obs_names[:5]}")

#         # Sets of obs_names
#         rna_set = set(rna_adata.obs_names)
#         acc_set = set(acc_adata.obs_names)
#         meth_set = set(meth_adata.obs_names)

#         # Cells in all three
#         fincells = list(rna_set & acc_set & meth_set)
#         all_cells = rna_set | acc_set | meth_set
#         dropped_cells = list(all_cells - set(fincells))

#         _msg(self.logger, f"{len(fincells)} surviving cells, {len(dropped_cells)} dropped cells.")
#         if len(fincells) == 0:
#             _msg(self.logger, "No cells in common between RNA, ACC and METH data. Exiting.", lvl='error')
#             sys.exit()
#         if len(dropped_cells) > 0:
#             _msg(self.logger, f"Dropped cells = {dropped_cells}")
#             for _cell in dropped_cells:
#                 _msg(self.logger, f"Cell {_cell} in RNA = {_cell in rna_set}, ACC = {_cell in acc_set}, METH = {_cell in meth_set}")
#         _msg(self.logger, f"Creating object with {len(fincells)} observations.")
#         rna_adata = rna_adata[rna_adata.obs_names.isin(fincells)].copy()
#         acc_adata = acc_adata[acc_adata.obs_names.isin(fincells)].copy()
#         meth_adata = meth_adata[meth_adata.obs_names.isin(fincells)].copy()
#         # Assert variables in acc / meth are equal
#         assert acc_adata.var.index.equals(meth_adata.var.index)
#         _mu = md.MuData(
#             {
#                 "RNA": rna_adata,
#                 "ACC": acc_adata,
#                 "METH": meth_adata
#             }
#         )
#         # Write to disk.
#         _of = self.opath / f"{self.project}.h5mu"
#         _mu.write(_of)
#         _msg(self.logger, f"mudata object written to {_of}")

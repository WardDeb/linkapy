from pathlib import Path
import sys
from rich import print
import logging
from datetime import datetime
import uuid
import polars as pl
import pandas as pd
import numpy as np
import anndata
import scipy as sp
from linkapy.linkapy import parse_cools
import mudata as md
import signal

def _msg(logger, msg, lvl='info'):
    print(msg)
    # For logging sake, remove '-'*100 from the strings.
    msg = msg.replace('-'*100, '')
    if lvl == 'info':
        logger.info(msg)
    elif lvl == 'debug':
        logger.debug(msg)
    elif lvl == 'warning':
        logger.warning(msg)
    else:
        logger.error(msg)
        sys.exit()

class Parse_scNMT:
    '''
    Parse_scNMT mainly functions to create matrices (arrow format for RNA, mtx format for accessibility / methylation)
    from directories containing analyzed scNMT-seq data. Theoretically this could be any type of multi-modal (read: RNA / methyatlion) data, but the class is written with the scNMT workflow
    from the Thienpont lab (KU Leuven) in mind.
    There are two required arguments (methpath and rnapath).
    Note that at least one region should be provided (genes, enhancers, CGI, proms, repeats) or the chromsizes file (for bins).

    :param str methpath: The path to the methylation directory (will be searched recursively!). Searches allcool files with \*WCGN\*allc.tsv.gz and \*GCHN\*.allc.tsv.gz for methylation and accessibility files, respectively. (required)
    :param str rnapath: The path to the RNA output directory (will be searched recursively!). For now looks for \*gene.tsv files (i.e. featureCounts output). Can handle single or multiple files, will be combined. (required)
    :param str project: Name of the project. Defaults to 'scNMT'. Generated matrices will carry the project in their name. (optional)  
    :param str opath: Name of the output directory to store matrices in. Defaults to None, which is the currect working directory. (optional)
    :param str chromsizes: Path to the chromsizes file for the genome. Defaults to None. If set, bins will be included in the accessibility/methylation aggregation. (optional)
    :param int threads: Number of threads to process Allcool files. Defaults to 10. (optional)
    :param str genes: Path to bed file containing genes to aggregate methylation signal over. Can be gzipped. (optional)
    :param str enhancers: Path to bed file containing enhancers to aggregate methylation signal over. Can be gzipped. (optional)
    :param str CGI: Path to bed file containing CGI to aggregate methylation signal over. Can be gzipped. (optional)
    :param str proms: Path to bed file containing proms to aggregate methylation signal over. Can be gzipped. (optional)
    :param str repeats: Path to bed file containing repeats to aggregate methylation signal over. Can be gzipped. (optional)
    '''
    def __init__(self, methpath = './', rnapath = './', project='scNMT', opath=None, chromsizes=None, threads=10, regions = None, qc=False):
        # Initiate a log
        timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        _fmt = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        # Set opath
        if opath:
            self.opath = Path(opath)
            self.opath.mkdir(parents=True, exist_ok=True)
        else:
            self.opath = Path.cwd()

        # Logfile
        log_file = Path(self.opath / f"lap_Parse_SCNMT_{timestamp}_{uuid.uuid4().hex}.log")
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(str(log_file))
        file_handler = logging.FileHandler(log_file)
        # To file
        file_handler.setFormatter(_fmt)
        self.logger.addHandler(file_handler)
        self.logger.propagate = False

        _msg(self.logger, "Output directory: " + str(self.opath))
        # Run QC or not
        self.qc = qc

        # Methylation data
        _m = Path(methpath)
        if not _m.exists():
            _msg(self.logger, f"Error: {methpath} does not exist", lvl='error')
        self.methpath = Path(methpath)

        # RNA data
        _r = Path(rnapath)
        if not _r.exists():
            _msg(self.logger, f"Error: {rnapath} does not exist", lvl='error')
        self.rnapath = Path(rnapath)

        # Parse paths
        self._glob_files()
        self.threads = threads
        self.project = project

        # Regions
        self.chromsizes = chromsizes
        self.regions = []
        self.regionlabels = []
        if regions:
            for _reg in regions:
                if not Path(_reg).exists():
                    _msg(self.logger, f"Error: {_reg} does not exist", lvl='error')
                else:
                    self.regions.append(_reg)
                    self.regionlabels.append(Path(_reg).name.replace('.bed.gz', '').replace('.bed', ''))

        if not self.regions and not self.chromsizes:
            sys.exit("No regions provided, and no chromsizes file provided to construct bins.")
        
    def _glob_files(self):
        # glob for allc.tsv.gz files
        _msg(self.logger, "-"*100 + "\n" + "Parse_scNMT - file globber" + "\n" + "-"*100)
        _msg(self.logger, f"Searching file paths: methylation = [green]{self.methpath}[/green] rna = [green]{self.rnapath}[/green].")
        self.allc_acc_files = list(self.methpath.rglob("*GCHN*.allc.tsv.gz"))
        self.allc_meth_files = list(self.methpath.rglob("*WCGN*.allc.tsv.gz"))
        assert len(self.allc_acc_files) == len(self.allc_meth_files)
        self.rna_files = list(self.rnapath.rglob("*gene.tsv"))
        _msg(self.logger, f"Found {len(self.allc_acc_files)} accessibility files.")
        _msg(self.logger, f"Found {len(self.allc_meth_files)} methylation files.")
        _msg(self.logger, f"Found {len(self.rna_files)} RNA file(s).")

    def _read_rna(self, _f):
        a = pl.read_csv(_f, separator='\t', skip_rows=1, has_header=True)
        a.columns = [i.replace("filtered.", "").replace(".Aligned.sortedByCoord.Processed.out.bam", "") for i in a.columns]
        schema = {
            'Geneid': pl.String,
            'Chr': pl.String,
            'Start': pl.UInt32,
            'End': pl.UInt32,
            'Strand': pl.String,
            'Length': pl.UInt32
        }
        # Fill in rest of counts
        for _s in a.columns:
            if _s not in schema:
                schema[_s] = pl.UInt32
        a = a.select([
            pl.col(col).cast(schema[col]) for col in a.columns
        ])
        metacol = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
        metadf = a.select(metacol)
        cdf = a.select([_c for _c in a.columns if _c not in metacol])
        return (metadf, cdf.lazy())

    def create_matrices(self):
        opath = self.opath
        opath.mkdir(parents=True, exist_ok=True)
        _msg(self.logger, "-"*100 + "\n" + "Parse_scNMT - Parse matrices" + "\n" + "-"*100)

        ## RNA
        _msg(self.logger, "Parsing RNA files.")
        if not Path(opath / (self.project + '_rnadf.arrow')).exists() and not Path(opath / (self.project + '_rnameta.arrow')).exists():    
            # Two situations - one featureCounts.tsv file, or more in wich case we need to merge.
            if len(self.rna_files) == 1:
                metadf, rnadf = self._read_rna(self.rna_files[0])
                rnadf = rnadf.collect()
                rnadf.write_ipc(opath / (self.project + '_rnadf.arrow'), compression='zstd')
                metadf.write_ipc(opath / (self.project + '_rnameta.arrow'), compression='zstd')

            else:
                rnadfs = []
                metadfs = []
                for _f in self.rna_files:
                    metadf, rnadf = self._read_rna(_f)
                    rnadfs.append(rnadf)
                    metadfs.append(metadf)
                # Make sure gene order is ok.
                assert(all(metadfs[0].equals(df) for df in metadfs[1:]))
                rnadf = pl.concat(rnadfs, how="horizontal")
                rnadf = rnadf.collect()
                rnadf.write_ipc(opath / (self.project + '_rnadf.arrow'), compression='zstd')
                metadf.write_ipc(opath / (self.project + '_rnameta.arrow'), compression='zstd')
            _msg(self.logger, f"RNA files written into {opath} 👍")
        else:
            _msg(self.logger, f"RNA files found at {opath} 👍")

        # Accessibility
        accbase = Path(opath / (self.project + '.acc'))
        accfile = Path(opath / (self.project + '.acc.meth.mtx'))
        metafile = Path(opath / (self.project + '.acc.meta.tsv'))
        cellfile = Path(opath / (self.project + '.acc.cell.tsv'))
        original_handler = signal.getsignal(signal.SIGINT)

        if not accfile.exists():
            try:
                signal.signal(signal.SIGINT, signal.SIG_DFL)
                parse_cools(
                    [str(i) for i in self.allc_acc_files],
                    self.regions,
                    self.regionlabels,
                    self.qc,
                    self.threads,
                    str(accbase),
                    str(metafile),
                    str(cellfile)
                )
                _msg(self.logger, f"Acc files written into {opath} 👍")
            finally:
                signal.signal(signal.SIGINT, original_handler)
        else:
            _msg(self.logger, f"Acc files found at {opath} 👍")
        
        # Methylation
        methbase = Path(opath / (self.project + '.meth'))
        methfile = Path(opath / (self.project + '.meth.meth.mtx'))
        metafile = Path(opath / (self.project + '.meth.meta.tsv'))
        cellfile = Path(opath / (self.project + '.meth.cell.tsv'))
        original_handler = signal.getsignal(signal.SIGINT)

        if not methfile.exists():
            try:
                signal.signal(signal.SIGINT, signal.SIG_DFL)
                parse_cools(
                    [str(i) for i in self.allc_meth_files],
                    self.regions,
                    self.regionlabels,
                    self.qc,
                    self.threads,
                    str(methbase),
                    str(metafile),
                    str(cellfile)
                )
                _msg(self.logger, f"Meth files written into {opath} 👍")
            finally:
                signal.signal(signal.SIGINT, original_handler)
        else:
            _msg(self.logger, f"Meth files found at {opath} 👍")


class Parse_matrices:
    '''
    Parses matrices created previously (with Parse_scNMT) and creates a muon object (written to disk).
    This is then the starting point to downstream analysis.

    :param str matrixdir: Directory where the matrices can be found. (required)
    :param str project: Project name, similar to the one provided in the Parse_scNMT part. Defaults to 'scNMT' (optional)
    :param str ofile: Name of the output file, defaults to matrixdir / project.h5. (optional)
    '''

    def __init__(self, matrixdir, project='scNMT', opath=None):
        self.matrixdir = Path(matrixdir)
        self.project = project
        if not opath:
            self.opath = self.matrixdir
        else:
            self.opath = Path(opath)
        self.opath.mkdir(parents=True, exist_ok=True)
        # Initiate a log
        timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        _fmt = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

        # Logfile
        log_file = Path(self.opath / f"lap_Parse_matrices_{timestamp}_{uuid.uuid4().hex}.log")
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(str(log_file))
        file_handler = logging.FileHandler(log_file)
        # To file
        file_handler.setFormatter(_fmt)
        self.logger.addHandler(file_handler)
        self.logger.propagate = False

        _msg(self.logger, "Output directory: " + str(self.opath))
        self._assert_files()

    def _assert_files(self):
        self.rna = self.matrixdir / (self.project + '_rnadf.arrow')
        self.rnameta = self.matrixdir / (self.project + '_rnameta.arrow')
        # Accessibility part
        self.acc = {
            'cov': self.matrixdir / (self.project + '.acc.cov.mtx'),
            'meth': self.matrixdir / (self.project + '.acc.meth.mtx'),
            'site': self.matrixdir / (self.project + '.acc.site.mtx'),
            'cell': self.matrixdir / (self.project + '.acc.cell.tsv'),
            'reg': self.matrixdir / (self.project + '.acc.meta.tsv')
        }
        self.meth = {
            'cov': self.matrixdir / (self.project + '.meth.cov.mtx'),
            'meth': self.matrixdir / (self.project + '.meth.meth.mtx'),
            'site': self.matrixdir / (self.project + '.meth.site.mtx'),
            'cell': self.matrixdir / (self.project + '.meth.cell.tsv'),
            'reg': self.matrixdir / (self.project + '.meth.meta.tsv')
        }

        assert self.rna.exists()
        assert self.rnameta.exists()
        for _f in self.acc:
            assert self.acc[_f].exists()
        for _f in self.meth:
            assert self.meth[_f].exists()


    def create_mudata(self, aggtype='fraction'):
        # Some settings.
        np.seterr(divide='ignore', invalid='ignore')
        md.set_options(pull_on_update=False)

        _msg(self.logger, "-"*100 + "\n" + "Parse_matrices - create_mudata" + "\n" + "-"*100)
        _msg(self.logger, "Parsing RNA matrices")
        rnadf = pl.read_ipc(self.rna, memory_map=False).to_pandas()
        rnameta = pl.read_ipc(self.rnameta, memory_map=False).to_pandas()
        rnameta.index = rnameta['Geneid']
        del rnameta['Geneid']
        rna_adata = anndata.AnnData(
            X=sp.sparse.csr_matrix(rnadf.values.T),
            obs=pd.DataFrame(index= [i.replace('_RNA-seq', '').replace("_RNA", "") for i in rnadf.columns]),
            var=rnameta
        )
        _msg(self.logger, "Parsing RNA matrices")
        _msg(self.logger, f"adata for rna shape = {rna_adata.shape}")
        _msg(self.logger, "Parsing Accessibility matrices")
        _m = sp.io.mmread(self.acc['meth']).todense()
        _c = sp.io.mmread(self.acc['cov']).todense()
        if aggtype == 'fraction':
            X = np.zeros_like(_m, dtype=float)
            X = np.where((_c != 0) & (_m != 0), _m / _c, 0)
            X = sp.sparse.csr_matrix(X)
        elif aggtype == 'sum':
            X = sp.sparse.csr_matrix(_m)
        _obs = pl.read_csv(self.acc['cell'], separator='\t', has_header=False).to_pandas()
        cell_names = []
        for i in _obs['column_1'].to_list():
            _n = Path(i).name.replace('.GCHN-Both.allc.tsv.gz', '').replace("_NOMe-seq", "").replace("_METH", "")
            cell_names.append(_n)
        _obs = pd.DataFrame(index = cell_names)
        _var = pl.read_csv(self.acc['reg'], separator='\t', has_header=True).to_pandas()
        # Since there is no check for duplicated values in the regions. We deduplicate here (by name)
        _var = _var.drop_duplicates(subset='name', keep='first')
        X = X[:, _var.index]
        _var = _var.set_index('name')
        acc_adata = anndata.AnnData(
            X=X,
            obs=_obs,
            var=_var
        )
        _msg(self.logger, f"adata for accessibility shape = {acc_adata.shape}")
        _msg(self.logger, "Parsing Methylation matrices")
        _m = sp.io.mmread(self.meth['meth']).todense()
        _c = sp.io.mmread(self.meth['cov']).todense()
        if aggtype == 'fraction':
            X = np.zeros_like(_m, dtype=float)
            X = np.where((_c != 0) & (_m != 0), _m / _c, 0)
            X = sp.sparse.csr_matrix(X)
        elif aggtype == 'sum':
            X = sp.sparse.csr_matrix(_m)
        _obs = pl.read_csv(self.meth['cell'], separator='\t', has_header=False).to_pandas()
        cell_names = []
        for i in _obs['column_1'].to_list():
            _n = Path(i).name.replace('.WCGN-Both.allc.tsv.gz', '').replace("_NOMe-seq", "").replace("_METH", "")
            cell_names.append(_n)
        _obs = pd.DataFrame(index = cell_names)
        _var = pl.read_csv(self.meth['reg'], separator='\t', has_header=True).to_pandas()
                # Since there is no check for duplicated values in the regions. We deduplicate here (by name)
        _var = _var.drop_duplicates(subset='name', keep='first')
        X = X[:, _var.index]
        _var = _var.set_index('name')
        meth_adata = anndata.AnnData(
            X=X,
            obs=_obs,
            var=_var
        )
        _msg(self.logger, f"adata for methylation shape = {meth_adata.shape}")
        # muData
        _msg(self.logger, "Creating muData object.")
        # Take intersection of all obs
        _msg(self.logger, f"First observations for RNA data = {rna_adata.obs_names[:5]}")
        _msg(self.logger, f"First observations for ACC data = {acc_adata.obs_names[:5]}")
        _msg(self.logger, f"First observations for METH data = {meth_adata.obs_names[:5]}")

        # Sets of obs_names
        rna_set = set(rna_adata.obs_names)
        acc_set = set(acc_adata.obs_names)
        meth_set = set(meth_adata.obs_names)

        # Cells in all three
        fincells = list(rna_set & acc_set & meth_set)
        all_cells = rna_set | acc_set | meth_set
        dropped_cells = list(all_cells - set(fincells))

        _msg(self.logger, f"{len(fincells)} surviving cells, {len(dropped_cells)} dropped cells.")
        if len(fincells) == 0:
            _msg(self.logger, "No cells in common between RNA, ACC and METH data. Exiting.", lvl='error')
            sys.exit()
        if len(dropped_cells) > 0:
            _msg(self.logger, f"Dropped cells = {dropped_cells}")
            for _cell in dropped_cells:
                _msg(self.logger, f"Cell {_cell} in RNA = {_cell in rna_set}, ACC = {_cell in acc_set}, METH = {_cell in meth_set}")
        _msg(self.logger, f"Creating object with {len(fincells)} observations.")
        rna_adata = rna_adata[rna_adata.obs_names.isin(fincells)].copy()
        acc_adata = acc_adata[acc_adata.obs_names.isin(fincells)].copy()
        meth_adata = meth_adata[meth_adata.obs_names.isin(fincells)].copy()
        # Assert variables in acc / meth are equal
        assert acc_adata.var.index.equals(meth_adata.var.index)
        _mu = md.MuData(
            {
                "RNA": rna_adata,
                "ACC": acc_adata,
                "METH": meth_adata
            }
        )
        # Write to disk.
        _of = self.opath / f"{self.project}.h5mu"
        _mu.write(_of)
        _msg(self.logger, f"mudata object written to {_of}")

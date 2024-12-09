from pathlib import Path
import sys
from rich import print
import logging
from datetime import datetime
import uuid
import pandas as pd
import gzip
from concurrent.futures import ProcessPoolExecutor
from itertools import chain
import tempfile
import numpy as np

def get_regions(_f):
    regions = set()
    with gzip.open(_f, 'rb') as f:
        for line in f:
            chrom, pos, strand, context, meth_count, cov, test = line.decode().strip().split('\t')
            _rstr = chrom + '_' + pos + '_' + strand + '_' + context
            regions.add(_rstr)
    return regions

class lap:
    def __init__(self, methpath = './', rnapath = './', threads=10):
        # Initiate a log
        timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        _fmt = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        # Logfile
        log_file = f"lap_{timestamp}_{uuid.uuid4().hex}.log"
        print(f"Log file: {log_file}")
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(log_file)
        file_handler = logging.FileHandler(log_file)
        # To file
        file_handler.setFormatter(_fmt)
        self.logger.addHandler(file_handler)
        self.logger.propagate = False        

        # Methylation data
        _m = Path(methpath)
        if not _m.exists():
            self._msg(f"Error: {methpath} does not exist", lvl='error')
        self.methpath = Path(methpath)

        # RNA data
        _r = Path(rnapath)
        if not _r.exists():
            self._msg(f"Error: {rnapath} does not exist", lvl='error')
        self.rnapath = Path(rnapath)

        # Parse paths
        self._glob_files()
        self.threads = threads
        self.odir = Path('./')

    def _glob_files(self):
        # glob for allc.tsv.gz files
        self._msg(
            f"Searching for methylation files at {self.methpath}, count matrices at {self.rnapath}."
        )
        self.logger.info("Searching for methylation & RNA files.")
        self.allc_acc_files = list(self.methpath.rglob("*WCGN*.allc.tsv.gz"))
        self.allc_meth_files = list(self.methpath.rglob("*GCHN*.allc.tsv.gz"))
        assert len(self.allc_acc_files) == len(self.allc_meth_files)
        self.rna_files = list(self.rnapath.rglob("*gene.tsv"))
        self._msg(
            f"Found {len(self.allc_acc_files)} acc. methylation files, {len(self.allc_meth_files)} meth. methylation files and {len(self.rna_files)} rna files."
        )

    def _msg(self, msg, lvl='info'):
        print(msg)
        if lvl == 'info':
            self.logger.info(msg)
        elif lvl == 'debug':
            self.logger.debug(msg)
        elif lvl == 'warning':
            self.logger.warning(msg)
        else:
            self.logger.error(msg)
            sys.exit()

    def create_matrices(self):
        def read_rna(p):
            a = pd.read_table(p, sep='\t', skiprows=1, header=0, index_col=0)
            a.columns = [_r.replace('filtered.', '').replace('.Aligned.sortedByCoord.Processed.out.bam', '') for _r in a.columns]
            metacols = ['Chr', 'Start' ,'End', 'Strand', 'Length']
            metadf = a[metacols]
            a.drop(columns=metacols, inplace=True)
            return (a, metadf)

        ## RNA
        self._msg("Parsing RNA files.")
        # Two situations - one featureCounts.tsv file, or more in wich case we need to merge.
        if len(self.rna_files) == 1:
            rnadf, metadf = read_rna(self.rna_files[0])
        else:
            rnadfs = []
            metadfs = []
            for _f in self.rna_files:
                rnadf, metadf = read_rna(_f)
                rnadfs.append(rnadf)
                metadfs.append(metadf)
            rnadf = pd.concat(rnadfs, axis=1)
            metadf = pd.concat(metadfs, axis=1)
        
        # Assert that the matrix is not empty
        if rnadf.empty:
            self._msg("Error: RNA count matrix is empty.", lvl='error')
        # Assert that indices for matrix and metadata match.
        assert rnadf.index.equals(metadf.index), self._msg("Error: RNA count matrix and metadata index mismatch.", lvl='error')

        ## Methylation
        ## There are three files per sample, WCGN, HCHN and GCHN. We need only WCGN (methylation) == ACGN & TCGN, and GCHN (accessibility) == GCAN, GCCN, GCTN.
        ## There's probably a more efficient way to do this, though allcools implementation seems a bit filehandle / memory intensive to.
        ## Idea is to loop over all the files, keep a list of regions, and only later parse them properly.
        def parse_allcool(_f, df):
            samplename = _f.name.split('.')[0]
            print(f"Working on sample {samplename}")
            a = pd.read_csv(
                _f,
                compression='gzip',
                header=None,
                sep='\t',
                usecols = [0,1,2,3,4,5], # Chrom, pos, strand, context, meth, cov, (6 == stat test, dropped).
                names = ['chrom', 'pos', 'strand', 'context', 'meth', 'cov']
            )
            a[samplename] = a['meth']/a['cov']
            a.index = a['chrom'].astype(str) + '_' + a['pos'].astype(str) + '_' + a['strand'].astype(str) + '_' + a['context'].astype(str)
            a = a[[samplename]].astype(pd.SparseDtype("float", fill_value=np.nan))
            df = pd.concat([df, a], axis=1)
            print(df.head())
            return(df)

        self._msg(f"Parsing regions for accessibility methylation (with {self.threads} threads)")
        with ProcessPoolExecutor(max_workers=self.threads) as exe:
            regions = list(exe.map(get_regions, self.allc_acc_files))
        regions = set(chain.from_iterable(regions))
        regions = sorted(set(regions), key=lambda x: (x.split('_')[0], int(x.split('_')[1])))
        self._msg(f"{len(regions)} positions found. Creating series for matrix.")
        f_r_pairs = [(i, regions) for i in self.allc_acc_files]
        accdf = pd.DataFrame(index=regions)
        for _f in self.allc_acc_files:
            accdf = parse_allcool(_f, accdf)
        accdf.to_hdf("acc_meth.h5", mode='w')

        self._msg(
            f"Final RNA count matrix dimensions: {rnadf.shape}"
        )
        return(rnadf, metadf)



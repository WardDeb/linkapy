from pathlib import Path
import sys
from rich import print
import logging
from datetime import datetime
import uuid
import pandas as pd

class lap:
    def __init__(self, methpath = './', rnapath = './'):
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

    def _glob_files(self):
        # glob for allc.tsv.gz files
        self._msg(
            f"Searching for methylation files at {self.methpath}, count matrices at {self.rnapath}."
        )
        self.logger.info("Searching for methylation & RNA files.")
        self.allc_files = list(self.methpath.rglob("*.allc.tsv.gz"))
        self.rna_files = list(self.rnapath.rglob("*gene.tsv"))
        self._msg(
            f"Found {len(self.allc_files)} meth and {len(self.rna_files)} rna files."
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
        def read_meth(p):
            # Do something for meth files.
            return pd.DataFrame()

        def read_rna(p):
            a = pd.read_table(p, sep='\t', skiprows=1, header=0, index_col=0)
            a.columns = [_r.replace('filtered.', '').replace('.Aligned.sortedByCoord.Processed.out.bam', '') for _r in a.columns]
            metacols = ['Chr', 'Start' ,'End', 'Strand', 'Length']
            metadf = a[metacols]
            a.drop(columns=metacols, inplace=True)
            return (a, metadf)
        
        ## RNA
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
        
              
        self._msg(
            f"Final RNA count matrix dimensions: {rnadf.shape}"
        )
        return(rnadf, metadf)
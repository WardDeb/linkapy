from pathlib import Path
import sys
from rich import print
import logging
from datetime import datetime
import uuid

class lap:
    def __init__(self, methpath, rnapath):
        # Initiate a log
        timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        log_file = f"lap_{timestamp}_{uuid.uuid4().hex}.log"
        print(f"Log file: {log_file}")
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(log_file)
        self.logger.setLevel(logging.DEBUG)
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(
            logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        )
        self.logger.addHandler(file_handler)

        # Methylation data
        _m = Path(methpath)
        if not _m.exists():
            print(f"Error: {methpath} does not exist")
            sys.exit(1)
        self.methpath = Path(methpath)

        # RNA data
        _r = Path(rnapath)
        if not _r.exists():
            print(f"Error: {rnapath} does not exist")
            sys.exit(1)
        self.rnapath = Path(rnapath)

        # Parse paths
        self._process_files()

    def _process_files(self):
        # glob for allc.tsv.gz files
        self.logger.info("Searching for methylation & RNA files.")
        self.allc_files = list(self.methpath.glob("*.allc.tsv.gz"))
        self.rna_files = list(self.rnapath.glob("*rsem.genes.results.tsv"))
        self.logger.info(f"Found {len(self.allc_files)} meth and {len(self.rna_files)} rna files.")

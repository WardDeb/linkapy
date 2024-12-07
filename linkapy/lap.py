from pathlib import Path
import sys
from rich import print
import logging

class lap:
    def __init__(self, methpath, rnapath):
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
        self.allc_files = list(self.methpath.glob("*.allc.tsv.gz"))
        self.rna_files = list(self.rnapath.glob("*rsem.genes.results.tsv"))

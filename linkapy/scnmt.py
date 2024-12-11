from pathlib import Path
import sys
from rich import print
import logging
from datetime import datetime
import uuid
import polars as pl

class Parse_scNMT:
    def __init__(self, methpath = './', rnapath = './', project='scNMT', threads=10):
        # Initiate a log
        timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        _fmt = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        # Logfile
        log_file = f"lap_Parse_SCNMT_{timestamp}_{uuid.uuid4().hex}.log"
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
        self.project = project

    def _glob_files(self):
        # glob for allc.tsv.gz files
        self._msg("-"*100 + "\n" + "Parse_scNMT - file globber" + "\n" + "-"*100)
        self._msg(f"Searching files: meth = {self.methpath} rna = {self.rnapath}.")
        self.allc_acc_files = list(self.methpath.rglob("*WCGN*.allc.tsv.gz"))
        self.allc_meth_files = list(self.methpath.rglob("*GCHN*.allc.tsv.gz"))
        assert len(self.allc_acc_files) == len(self.allc_meth_files)
        self.rna_files = list(self.rnapath.rglob("*gene.tsv"))
        self._msg(f"Found {len(self.allc_acc_files)} accessibility files.")
        self._msg(f"Found {len(self.allc_meth_files)} methylation files.")
        self._msg(f"Found {len(self.rna_files)} RNA files.")

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

    def create_matrices(self, opath='./'):
        opath = Path(opath)
        opath.mkdir(parents=True, exist_ok=True)
        self._msg("-"*100 + "\n" + "Parse_scNMT - Parse matrices" + "\n" + "-"*100)
        def read_rna(_f):
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

        def parse_allcool(_f, threads):
            samplename = _f.name.split('.')[0]
            print(f"Reading sample {samplename}")
            a = pl.read_csv(
                _f,
                has_header=False,
                separator='\t',
                #usecols = [0,1,2,3,4,5], # Chrom, pos, strand, context, meth, cov, (6 == stat test, dropped).
                new_columns = ['chrom', 'pos', 'strand', 'context', 'meth', 'cov', 'stat'],
                n_threads = threads,
                schema = {
                    'chrom': pl.String,
                    'pos': pl.UInt32,
                    'strand': pl.String,
                    'context': pl.String,
                    'meth': pl.UInt32,
                    'cov': pl.UInt32,
                    'stat': pl.Float32
                }
            )
            a = a.with_columns(
                (pl.col("meth") / pl.col("cov")).alias(samplename).cast(pl.Float32)
            )
            a = a.with_columns(
                pl.concat_str(
                    [
                        pl.col("chrom"),
                        pl.lit("_"),
                        pl.col("pos").cast(pl.String),
                        pl.lit("_"),
                        pl.col("strand"),
                        pl.lit("_"),
                        pl.col("context"),
                    ]
                ).alias("index")
            )
            a = a.select(["index", samplename])
            return(a.lazy())


        ## RNA
        self._msg("Parsing RNA files.")
        # Two situations - one featureCounts.tsv file, or more in wich case we need to merge.
        if len(self.rna_files) == 1:
            metadf, rnadf = read_rna(self.rna_files[0])
            rnadf.write_parquet(
                opath / (self.project + '_rnadf.parquet'),
                compression_level=9
            )
            metadf.write_parquet(
                opat / (self.project + '_rnameta.parquet'),
                compression_level=9
            )

        else:
            rnadfs = []
            metadfs = []
            for _f in self.rna_files:
                metadf, rnadf = read_rna(_f)
                rnadfs.append(rnadf)
                metadfs.append(metadf)
            # Make sure gene order is ok.
            assert(all(metadfs[0].equals(df) for df in metadfs[1:]))
            rnadf = pl.concat(rnadfs, how="horizontal").collect()
            rnadf.write_parquet(
                opath / (self.project + '_rnadf.parquet'),
                compression_level=9
            )
            metadf[0].write_parquet(
                opath / (self.project + '_rnameta.parquet'),
                compression_level=9
            )

        # Accessibility
        # There are three files per sample, WCGN, HCHN and GCHN. We need only WCGN (methylation) == ACGN & TCGN, and GCHN (accessibility) == GCAN, GCCN, GCTN.
        # There's probably a more efficient way to do this, though allcools implementation seems a bit filehandle / memory intensive to.
        if not Path(opath / (self.project + '_acc_df.parquet')).exists():
            self._msg("Parsing accessibility files.")
            accdfs = []
            for _f in self.allc_acc_files:
                accdfs.append(parse_allcool(_f, self.threads))
            
            accdf = accdfs[0]
            for df in accdfs[1:]:
                print("Will try joining treasure now.")
                accdf = accdf.join(df, how="outer", on="index").select(
                    [col for col in accdf.columns if col != 'index_right']
                )
                #accdf = accdf.drop("index_right")
            faccdf = accdf.collect()
            faccdf.write_parquet(
                opath / (self.project + '_acc_df.parquet'),
                compression_level=9
            )

        # Methylation
        if not Path(opath / (self.project + '_meth_df.parquet')).exists():
            self._msg("Parsing methylation files.")
            accdfs = []
            for _f in self.allc_meth_files:
                accdfs.append(parse_allcool(_f, self.threads))

            accdf = accdfs[0]
            for df in accdfs[1:]:
                accdf = accdf.join(df, how="outer", on="index")
                accdf = accdf.drop("index_right")
            faccdf = accdf.collect()
            faccdf.write_parquet(
                opath / (self.project + '_meth_df.parquet'),
                compression_level=9
            )
        
        self._msg(f"Parquet files written to disk: {opath}")



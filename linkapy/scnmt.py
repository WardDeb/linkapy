from pathlib import Path
import sys
from rich import print
import logging
from datetime import datetime
import uuid
import polars as pl
import tempfile

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
        self._msg(f"Searching file paths: methylation = [green]{self.methpath}[/green] rna = [green]{self.rnapath}[/green].")
        self.allc_acc_files = list(self.methpath.rglob("*WCGN*.allc.tsv.gz"))
        self.allc_meth_files = list(self.methpath.rglob("*GCHN*.allc.tsv.gz"))
        assert len(self.allc_acc_files) == len(self.allc_meth_files)
        self.rna_files = list(self.rnapath.rglob("*gene.tsv"))
        self._msg(f"Found {len(self.allc_acc_files)} accessibility files.")
        self._msg(f"Found {len(self.allc_meth_files)} methylation files.")
        self._msg(f"Found {len(self.rna_files)} RNA file(s).")

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

        def parse_allcool(_f, threads, chroms=False, chrom=None):
            if chroms == True:
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
                # Get unique chromosomes
                chroms = sorted(a['chrom'].unique().to_list())
                return chroms
            else:
                samplename = _f.name.split('.')[0]
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
                # subset for chromosome
                
                a = a.filter(pl.col("chrom") == chrom)
                
                a = a.with_columns(
                    (pl.col("meth") / pl.col("cov")).alias("accv").cast(pl.Float32)
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
                a = a.select(["index", "accv"])
                a = a.with_columns(
                    pl.lit(samplename).alias("sample")
                )
                return(a.lazy(), samplename)


        ## RNA
        self._msg("Parsing RNA files.")
        if not Path(opath / (self.project + '_rnadf.parquet')).exists() and not Path(opath / (self.project + '_rnameta.parquet')).exists():    
            # Two situations - one featureCounts.tsv file, or more in wich case we need to merge.
            if len(self.rna_files) == 1:
                metadf, rnadf = read_rna(self.rna_files[0])
                rnadf.sink_parquet(
                    opath / (self.project + '_rnadf.parquet'),
                    compression_level=9
                )
                metadf.write_parquet(
                    opath / (self.project + '_rnameta.parquet'),
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
                rnadf = pl.concat(rnadfs, how="horizontal")
                rnadf.collect().write_parquet(
                    opath / (self.project + '_rnadf.parquet'),
                    compression_level=9
                )
                metadf[0].write_parquet(
                    opath / (self.project + '_rnameta.parquet'),
                    compression_level=9
                )
            self._msg(f"RNA files written into {opath} 👍")
        else:
            self._msg(f"RNA files found at {opath} 👍")

        # Accessibility
        # There are three files per sample, WCGN, HCHN and GCHN. We need only WCGN (methylation) == ACGN & TCGN, and GCHN (accessibility) == GCAN, GCCN, GCTN.
        # There's probably a more efficient way to do this, though allcools implementation seems a bit filehandle / memory intensive to.
        
        # Plug and play the chromosomes to chunks
        # combine chunks

        if not Path(opath / ('.' + self.project + '_accparquets.created')).exists():
            self._msg("No accessibility flag found, Parsing accessibility files.")
            # Parse first file to get all chromosomes
            chroms = parse_allcool(self.allc_acc_files[0], self.threads, chroms=True)
            print(chroms)
            tmpfiles = []
            for _chrom in chroms:
                self._msg(f"Processing chromosome: {_chrom}")
                accdfs = []
                samples = []
                for _f in self.allc_acc_files:
                    lf, sample = parse_allcool(_f, self.threads, chrom=_chrom)
                    accdfs.append(lf)
                    samples.append(sample)
                accdf = pl.concat(accdfs).collect()
                accdf = accdf.pivot(
                    values="accv",
                    index="index",
                    columns="sample"
                )
                _tfile = tempfile.NamedTemporaryFile(delete=True)
                tmpfiles.append(_tfile.name)
                accdf.write_parquet(_tfile, compression_level=9)
            
            # Merge all files
            _dfs = [pl.scan_parquet(_f) for _f in tmpfiles]
            _df_concat = pl.concat(_dfs)

            _df_concat.sink_parquet(
                opath / (self.project + '_acc_df.parquet'),
                compression_level=9
            )

            # Per chromosome, parse and collapse.
            #self._msg("Parsing accessibility files.")
            #accdfs = []
            #samples = []
            #for _f in self.allc_acc_files:
            #    lf, sample = parse_allcool(_f, self.threads)
            #    accdfs.append(lf)
            #    samples.append(sample)
            #accc = pl.concat(accdfs)
            #agf = lambda x: x.mean() # There is nothing to aggregate so just trivial function for first (and only) list val
            #accc.group_by(pl.col("index")).agg(
            #    agf(pl.col("accv").filter(pl.col("sample") == _s)).alias(_s)
            #    for _s in samples
            #).sort("index").sink_parquet(opath / (self.project + '_acc_df.parquet'))




            # accc.write_parquet(opath / (self.project + '_tmp_acc_df.parquet'))
            # _tread = pl.scan_parquet(opath / (self.project + '_tmp_acc_df.parquet'))
            # for c in _tread.collect_streams():
            #     _tpiv = c.pivot(
            #         values="accv",
            #         index="index",
            #         columns="sample"
            #     )
            #     _tpiv.write_parquet(opath / (self.project + '_acc_df.parquet'), compression="zstd", mode="append", compression_level=9)

        # Methylation
        # if not Path(opath / (self.project + '_meth_df.parquet')).exists():
        #     self._msg("Parsing methylation files.")
        #     accdfs = []
        #     for _f in self.allc_meth_files:
        #         accdfs.append(parse_allcool(_f, self.threads))

        #     accdf = accdfs[0]
        #     for df in accdfs[1:]:
        #         accdf = accdf.join(df, how="outer", on="index")
        #         accdf = accdf.drop("index_right")
        #     faccdf = accdf.collect()
        #     faccdf.write_parquet(
        #         opath / (self.project + '_meth_df.parquet'),
        #         compression_level=9
        #     )
        
        self._msg(f"Parquet files written to disk: {opath}")



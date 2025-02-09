from pathlib import Path
import sys
from rich import print
import logging
from datetime import datetime
import uuid
import polars as pl
import tempfile

class Parse_scNMT:
    def __init__(self, methpath = './', rnapath = './', project='scNMT', opath=None, chromsizes=None, threads=10):

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
        log_file = Path(self.opath / f"lap_Parse_SCNMT_{timestamp}_{uuid.uuid4().hex}.log").name
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(log_file)
        file_handler = logging.FileHandler(log_file)
        # To file
        file_handler.setFormatter(_fmt)
        self.logger.addHandler(file_handler)
        self.logger.propagate = False

        self._msg("Output directory: " + str(self.opath))
        # chromsizes, if set.
        self.chromsizes = {}
        self.chunks = []
        self._parse_chromsizes(chromsizes)

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
        self.accchroms = []
        self.methchroms = []

    def _parse_chromsizes(self, chromsizes):
        if not chromsizes:
            self._msg("No chromsizes provided. Data will be parsed per chromosome (could be memory hungry 💀).")
            return
        if not Path(chromsizes).exists():
            self._msg(f"Chromsizes provided but doesn't exist {chromsizes}. Data will be parsed per chromosome (could be memory hungry 💀).")
            return
        with open(chromsizes) as f:
            for line in f:
                chrom, size = line.strip().split('\t')
                # Only retain canonical chromosomes.
                if 'chr' in chrom and chrom != 'chrM':
                    self.chromsizes[chrom] = int(size)
        # Define the chunks
        chunksize = 50000000
        for chrom, size in self.chromsizes.items():
            start = 0
            l_end = min(chunksize, size)
            self.chunks.append([chrom, start, l_end])
            while l_end < size:
                start = l_end
                l_end = min(l_end+chunksize, size)
                self.chunks.append([chrom, start, l_end])

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
        # For logging sake, remove '-'*100 from the strings.
        msg = msg.replace('-'*100, '')
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
        opath = self.opath
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

        def parse_allcool(_f, threads, chrom=None, chroms=False):
            if chroms:
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
                # Check if chrom is a list
                if isinstance(chrom, list):
                    _chr, _start, _end = chrom
                    a = a.filter((pl.col("chrom") == _chr) & (pl.col("pos") >= _start) & (pl.col("pos") <= _end))
                elif chrom:
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
        if not Path(opath / (self.project + '_rnadf.arrow')).exists() and not Path(opath / (self.project + '_rnameta.arrow')).exists():    
            # Two situations - one featureCounts.tsv file, or more in wich case we need to merge.
            if len(self.rna_files) == 1:
                metadf, rnadf = read_rna(self.rna_files[0])
                rnadf = rnadf.collect()
                rnadf.write_ipc(opath / (self.project + '_rnadf.arrow'), compression='zstd')
                metadf.write_ipc(opath / (self.project + '_rnameta.arrow'), compression='zstd')

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
                rnadf = rnadf.collect()
                rnadf.write_ipc(opath / (self.project + '_rnadf.arrow'), compression='zstd')
                metadf.write_ipc(opath / (self.project + '_rnameta.arrow'), compression='zstd')
            self._msg(f"RNA files written into {opath} 👍")
        else:
            self._msg(f"RNA files found at {opath} 👍")

        # Accessibility
        # There are three files per sample, WCGN, HCHN and GCHN. We need only WCGN (methylation) == ACGN & TCGN, and GCHN (accessibility) == GCAN, GCCN, GCTN.
        # There's probably a more efficient way to do this, though allcools implementation seems a bit filehandle / memory intensive to.
        
        if not Path(opath / (self.project + '_acc_df.arrow')).exists():
            self._msg("No accessibility file found, Parsing accessibility files.")
            # Parse first file to get all chromosomes
            self.accchroms = parse_allcool(self.allc_acc_files[0], self.threads, chroms=True)
            self._msg(f"Found {len(self.accchroms)} chromosomes in accessibility files.")
            tmpfiles = []
            if self.chunks:
                for _chrom in self.chunks:
                    self._msg(f"Processing chromosome: {_chrom}")
                    accdfs = []
                    samples = []
                    _fix = 0
                    for _f in self.allc_acc_files:
                        print(f"Processing file {_fix + 1} / {len(self.allc_acc_files)}", end='\r')
                        lf, sample = parse_allcool(_f, self.threads, chrom=_chrom)
                        accdfs.append(lf)
                        samples.append(sample)
                        _fix += 1
                    accdf = pl.concat(accdfs).collect()
                    accdf = accdf.pivot(
                        values="accv",
                        index="index",
                        columns="sample"
                    )
                    for _s in self.allc_acc_files:
                        _s = _s.name.split('.')[0]
                        if _s not in accdf.columns:
                            accdf = accdf.with_columns(pl.lit(float('nan')).cast(pl.Float32).alias(_s))
                    # Fix order.
                    accdf = accdf.select([_s.name.split('.')[0] for _s in self.allc_acc_files])
                    # Assert the number of columns match
                    assert len(accdf.columns) == len(self.allc_acc_files)
                    _tfile = tempfile.NamedTemporaryFile(delete=False)
                    tmpfiles.append(_tfile.name)
                    accdf.write_ipc(_tfile, compression='zstd')
            else:
                for _chrom in self.accchroms:
                    self._msg(f"Processing chromosome: {_chrom}")
                    accdfs = []
                    samples = []
                    _fix = 0
                    for _f in self.allc_acc_files:
                        print(f"Processing file {_fix + 1} / {len(self.allc_acc_files)}", end='\r')
                        lf, sample = parse_allcool(_f, self.threads, chrom=_chrom)
                        accdfs.append(lf)
                        samples.append(sample)
                        _fix += 1
                    accdf = pl.concat(accdfs).collect()
                    accdf = accdf.pivot(
                        values="accv",
                        index="index",
                        columns="sample"
                    )
                    for _s in self.allc_acc_files:
                        _s = _s.name.split('.')[0]
                        if _s not in accdf.columns:
                            accdf = accdf.with_columns(pl.lit(float('nan')).cast(pl.Float32).alias(_s))
                    # Fix order.
                    accdf = accdf.select([_s.name.split('.')[0] for _s in self.allc_acc_files])
                    # Assert the number of columns match
                    assert len(accdf.columns) == len(self.allc_acc_files)
                    _tfile = tempfile.NamedTemporaryFile(delete=False)
                    tmpfiles.append(_tfile.name)
                    accdf.write_ipc(_tfile, compression='zstd')

            
            # Merge all files
            _dfs = [pl.read_ipc(_f, memory_map=False) for _f in tmpfiles]
            _df_concat = pl.concat(_dfs)

            _df_concat.write_ipc(
                opath / (self.project + '_acc_df.arrow'),
                compression='zstd'
            )
            # Remove all the tmpfiles.
            for _f in tmpfiles:
                Path(_f).unlink()
            self._msg(f"Accessibility files written into {opath} 👍")
        else:
            self._msg(f"Accessibility file found in {opath} 👍")
        
        # Methylation
        if not Path(opath / (self.project + '_meth_df.arrow')).exists():
            self._msg("No methylation file found, Parsing methylation files.")
            self.methchroms = parse_allcool(self.allc_meth_files[0], self.threads, chroms=True)
            self._msg(f"Found {len(self.methchroms)} chromosomes in methylation files.")
            tmpfiles = []
            if self.chunks:
                for _chrom in self.chunks:
                    self._msg(f"Processing chromosome: {_chrom}")
                    methdfs = []
                    samples = []
                    _fix = 0
                    for _f in self.allc_meth_files:
                        print(f"Processing file {_fix + 1} / {len(self.allc_meth_files)}", end='\r')
                        lf, sample = parse_allcool(_f, self.threads, chrom=_chrom)
                        methdfs.append(lf)
                        samples.append(sample)
                        _fix += 1
                    methdf = pl.concat(methdfs).collect()
                    methdf = methdf.pivot(
                        values="accv",
                        index="index",
                        columns="sample"
                    )
                    for _s in self.allc_meth_files:
                        _s = _s.name.split('.')[0]
                        if _s not in methdf.columns:
                            methdf = methdf.with_columns(pl.lit(float('nan')).cast(pl.Float32).alias(_s))
                    # Fix order.
                    methdf = methdf.select([_s.name.split('.')[0] for _s in self.allc_meth_files])
                    # Assert the number of columns match
                    assert len(accdf.columns) == len(self.allc_meth_files)
                    _tfile = tempfile.NamedTemporaryFile(delete=False)
                    tmpfiles.append(_tfile.name)
                    methdf.write_ipc(_tfile, compression='zstd')
            else:
                for _chrom in self.methchroms:
                    self._msg(f"Processing chromosome: {_chrom}")
                    methdfs = []
                    samples = []
                    _fix = 0
                    for _f in self.allc_meth_files:
                        print(f"Processing file {_fix + 1} / {len(self.allc_meth_files)}", end='\r')
                        lf, sample = parse_allcool(_f, self.threads, chrom=_chrom)
                        methdfs.append(lf)
                        samples.append(sample)
                        _fix += 1
                    methdf = pl.concat(methdfs).collect()
                    methdf = methdf.pivot(
                        values="accv",
                        index="index",
                        columns="sample"
                    )
                    for _s in self.allc_meth_files:
                        _s = _s.name.split('.')[0]
                        if _s not in methdf.columns:
                            methdf = methdf.with_columns(pl.lit(float('nan')).cast(pl.Float32).alias(_s))
                    # Fix order.
                    methdf = methdf.select([_s.name.split('.')[0] for _s in self.allc_meth_files])
                    # Assert the number of columns match
                    assert len(accdf.columns) == len(self.allc_meth_files)
                    _tfile = tempfile.NamedTemporaryFile(delete=False)
                    tmpfiles.append(_tfile.name)
                    methdf.write_ipc(_tfile, compression='zstd')
            
            # Merge all files
            _dfs = [pl.read_ipc(_f, memory_map=False) for _f in tmpfiles]
            _df_concat = pl.concat(_dfs)

            _df_concat.write_ipc(
                opath / (self.project + '_meth_df.arrow'),
                compression='zstd'
            )
            # Remove all the tmpfiles.
            for _f in tmpfiles:
                Path(_f).unlink()
            self._msg(f"Methylation file written into {opath} 👍")
        else:
            self._msg(f"Methylation file found in {opath} 👍")
        
        self._msg(f"Parquet files written to disk: {opath}")

class Parse_matrices:
    def __init__(self, dir):
        self.matrixdir = Path(dir)
        self._assert_files()

    def _assert_files(self):
        self.rna = self.matrixdir / 'scNMT_rnadf.parquet'
        self.rnameta = self.matrixdir / 'scNMT_rnameta.parquet'
        self.acc = self.matrixdir / 'scNMT_acc_df.parquet'
        self.meth = self.matrixdir / 'scNMT_meth_df.parquet'

        assert self.rna.exists()
        assert self.rnameta.exists()
        assert self.acc.exists()
        assert self.meth.exists()


    def create_muon(self):
        # Read rna parquet file
        #rnadf = pl.read_parquet(self.rna)
        #rnameta = pl.read_parquet(self.rnameta)
        #print(rnadf.head())
        #print(rnameta.head())
        #accdf = pl.read_parquet(self.acc)
        methdf = pl.read_parquet(self.meth)
        #print(accdf.head())
        print(methdf.head())
        
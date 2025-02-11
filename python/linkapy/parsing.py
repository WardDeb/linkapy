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
from scipy.sparse import csr_matrix
from linkapy.linkars import parse_cools

class Parse_scNMT:
    def __init__(self, methpath = './', rnapath = './', project='scNMT', opath=None, chromsizes=None, threads=10, genes=None, enhancers=None, CGI=None, proms=None, repeats=None, binsize=100000):

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

        # Regions
        self.regions = []
        self.regionlabels = []
        if genes:
            self.regions.append(genes)
            self.regionlabels.append('genes')
        if enhancers:
            self.regions.append(enhancers)
            self.regionlabels.append('enhancers')
        if CGI:
            self.regions.append(CGI)
            self.regionlabels.append('CGI')
        if proms:
            self.regions.append(proms)
            self.regionlabels.append('promoters')
        if repeats:
            self.regions.append(repeats)
            self.regionlabels.append('repeats')
        
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

    def _read_rna(_f):
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
        self._msg("-"*100 + "\n" + "Parse_scNMT - Parse matrices" + "\n" + "-"*100)

        ## RNA
        self._msg("Parsing RNA files.")
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
            self._msg(f"RNA files written into {opath} 👍")
        else:
            self._msg(f"RNA files found at {opath} 👍")

        # Accessibility
        accfile = Path(opath / (self.project + 'acc.mtx'))
        metafile = Path(opath / (self.project + 'acc_meta.tsv'))
        if not accfile.exists():
            parse_cools(
                [str(i) for i in self.allc_acc_files],
                self.regions,
                self.regionlabels,
                self.threads,
                str(accfile),
                str(metafile)
            )
        
        # Methylation
        methfile = Path(opath / (self.project + 'meth.mtx'))
        metafile = Path(opath / (self.project + 'meth_meta.tsv'))
        if not methfile.exists():
            parse_cools(
                [str(i) for i in self.allc_meth_files],
                self.regions,
                self.regionlabels,
                self.threads,
                str(methfile),
                str(metafile)
            )
        
        # # # Accessibility
        # # # There are three files per sample, WCGN, HCHN and GCHN. We need only WCGN (methylation) == ACGN & TCGN, and GCHN (accessibility) == GCAN, GCCN, GCTN.
        # # # There's probably a more efficient way to do this, though allcools implementation seems a bit filehandle / memory intensive to.
        # if what == 'all' or what == 'acc':
        #     if not Path(opath / (self.project + '_acc_df_bins.arrow')).exists():
        #         self._msg("Bins accessibility file not found, Parsing accessibility files.")
        #         dfs = []
        #         six = 0
        #         for _f in self.allc_acc_files:
        #             print(f"Sample {six} / {len(self.allc_acc_files)}", end='\r')
        #             dfs.append(parse_allcool(_f, self.binfilters))
        #             six += 1
        #         accdf = pl.concat(dfs, how='horizontal')
        #         accdf.write_ipc(
        #             opath / (self.project + '_acc_df_bins.arrow'),
        #             compression='zstd'
        #         )
        #         del accdf
        #         gc.collect()

        #     if not Path(opath / (self.project + '_acc_df_regs.arrow')).exists():
        #         self._msg("Region accessibility file not found, Parsing accessibility files.")
        #         dfs = []
        #         six = 0
        #         for _f in self.allc_acc_files:
        #             print(f"Sample {six} / {len(self.allc_acc_files)}", end='\r')
        #             dfs.append(parse_allcool(_f, self.filters))
        #             six += 1
        #         accdf = pl.concat(dfs, how='horizontal')
        #         accdf.write_ipc(
        #             opath / (self.project + '_acc_df_regs.arrow'),
        #             compression='zstd'
        #         )
        #         del accdf
        #         gc.collect()
        # if what == 'all' or what == 'meth':
        #     if not Path(opath / (self.project + '_meth_df_bins.arrow')).exists():
        #         self._msg("Bins meth file not found, Parsing methylation files.")
        #         dfs = []
        #         six = 0
        #         for _f in self.allc_meth_files:
        #             print(f"Sample {six} / {len(self.allc_meth_files)}", end='\r')
        #             dfs.append(parse_allcool(_f, self.binfilters))
        #             six += 1
        #         methdf = pl.concat(dfs, how='horizontal')
        #         methdf.write_ipc(
        #             opath / (self.project + '_acc_df_bins.arrow'),
        #             compression='zstd'
        #         )
        #         del methdf
        #         gc.collect()

        #     if not Path(opath / (self.project + '_meth_df_regs.arrow')).exists():
        #         self._msg("Region meth file not found, Parsing methylation files.")
        #         dfs = []
        #         six = 0
        #         for _f in self.allc_meth_files:
        #             print(f"Sample {six} / {len(self.allc_meth_files)}", end='\r')
        #             dfs.append(parse_allcool(_f, self.filters))
        #             six += 1
        #         methdf = pl.concat(dfs, how='horizontal')
        #         methdf.write_ipc(
        #             opath / (self.project + '_acc_df_regs.arrow'),
        #             compression='zstd'
        #         )
        #         del methdf
        #         gc.collect()
        #     # Merge all files
        #     _dfs = [pl.read_ipc(_f, memory_map=False) for _f in tmpfiles]
        #     _df_concat = pl.concat(_dfs)

        #     _df_concat.write_ipc(
        #         opath / (self.project + '_acc_df.arrow'),
        #         compression='zstd'
        #     )
        #     # Remove all the tmpfiles.
        #     for _f in tmpfiles:
        #         Path(_f).unlink()
        #     self._msg(f"Accessibility files written into {opath} 👍")
        # else:
        #     self._msg(f"Accessibility file found in {opath} 👍")
        
        # # Methylation
        # if not Path(opath / (self.project + '_meth_df.arrow')).exists():
        #     self._msg("No methylation file found, Parsing methylation files.")
        #     self.methchroms = parse_allcool(self.allc_meth_files[0], self.threads, chroms=True)
        #     self._msg(f"Found {len(self.methchroms)} chromosomes in methylation files.")
        #     tmpfiles = []
        #     if self.chunks:
        #         for _chrom in self.chunks:
        #             self._msg(f"Processing chromosome: {_chrom}")
        #             methdfs = []
        #             samples = []
        #             _fix = 0
        #             for _f in self.allc_meth_files:
        #                 print(f"Processing file {_fix + 1} / {len(self.allc_meth_files)}", end='\r')
        #                 lf, sample = parse_allcool(_f, self.threads, chrom=_chrom)
        #                 methdfs.append(lf)
        #                 samples.append(sample)
        #                 _fix += 1
        #             methdf = pl.concat(methdfs).collect()
        #             methdf = methdf.pivot(
        #                 values="accv",
        #                 index="index",
        #                 columns="sample"
        #             )
        #             for _s in self.allc_meth_files:
        #                 _s = _s.name.split('.')[0]
        #                 if _s not in methdf.columns:
        #                     methdf = methdf.with_columns(pl.lit(float('nan')).cast(pl.Float32).alias(_s))
        #             # Fix order.
        #             methdf = methdf.select(['index'] + [_s.name.split('.')[0] for _s in self.allc_meth_files])
        #             # Assert the number of columns match
        #             assert len(methdf.columns) == len(self.allc_meth_files) + 1
        #             _tfile = tempfile.NamedTemporaryFile(delete=False)
        #             tmpfiles.append(_tfile.name)
        #             methdf.write_ipc(_tfile, compression='zstd')
        #     else:
        #         for _chrom in self.methchroms:
        #             self._msg(f"Processing chromosome: {_chrom}")
        #             methdfs = []
        #             samples = []
        #             _fix = 0
        #             for _f in self.allc_meth_files:
        #                 print(f"Processing file {_fix + 1} / {len(self.allc_meth_files)}", end='\r')
        #                 lf, sample = parse_allcool(_f, self.threads, chrom=_chrom)
        #                 methdfs.append(lf)
        #                 samples.append(sample)
        #                 _fix += 1
        #             methdf = pl.concat(methdfs).collect()
        #             methdf = methdf.pivot(
        #                 values="accv",
        #                 index="index",
        #                 columns="sample"
        #             )
        #             for _s in self.allc_meth_files:
        #                 _s = _s.name.split('.')[0]
        #                 if _s not in methdf.columns:
        #                     methdf = methdf.with_columns(pl.lit(float('nan')).cast(pl.Float32).alias(_s))
        #             # Fix order.
        #             methdf = methdf.select(['index'] + [_s.name.split('.')[0] for _s in self.allc_meth_files])
        #             # Assert the number of columns match
        #             assert len(methdf.columns) == len(self.allc_meth_files) + 1
        #             _tfile = tempfile.NamedTemporaryFile(delete=False)
        #             tmpfiles.append(_tfile.name)
        #             methdf.write_ipc(_tfile, compression='zstd')
            
        #     # Merge all files
        #     _dfs = [pl.read_ipc(_f, memory_map=False) for _f in tmpfiles]
        #     _df_concat = pl.concat(_dfs)

        #     _df_concat.write_ipc(
        #         opath / (self.project + '_meth_df.arrow'),
        #         compression='zstd'
        #     )
        #     # Remove all the tmpfiles.
        #     for _f in tmpfiles:
        #         Path(_f).unlink()
        #     self._msg(f"Methylation file written into {opath} 👍")
        # else:
        #     self._msg(f"Methylation file found in {opath} 👍")

class Parse_matrices:
    def __init__(self, dir):
        self.matrixdir = Path(dir)
        self._assert_files()

    def _assert_files(self):
        self.rna = self.matrixdir / 'scNMT_rnadf.arrow'
        self.rnameta = self.matrixdir / 'scNMT_rnameta.arrow'
        self.acc = self.matrixdir / 'scNMT_acc_df.arrow'
        self.meth = self.matrixdir / 'scNMT_meth_df.arrow'

        assert self.rna.exists()
        assert self.rnameta.exists()
        assert self.acc.exists()
        assert self.meth.exists()


    def create_muon(self, adir):
        adir = Path(adir)
        adir.mkdir(parents=True, exist_ok=True)
        # RNA part.
        rnadf = pl.read_ipc(self.rna, memory_map=False).to_pandas()
        rnameta = pl.read_ipc(self.rnameta, memory_map=False).to_pandas()
        rnameta.index = rnameta['Geneid']
        del rnameta['Geneid']
        rna_adata = anndata.AnnData(
            X=csr_matrix(rnadf.values),
            var=pd.DataFrame(index=list(rnadf.columns)),
            obs=rnameta
        )
        rna_adata.write_h5ad(
            Path(adir) / 'RNA.h5'
        )

        # Acc part.      
        _samples = list(pl.read_ipc_schema(self.acc).keys())[1::]
        _var = pd.DataFrame(index=_samples)
        X = np.ma.masked_invalid(pl.read_ipc(self.acc, memory_map=False, columns=_samples).to_pandas().values)
        print(X)
        #_obs = pl.read_ipc(self.acc, memory_map=False, columns=['index']).to_pandas()['index'].str.split('_', expand=True)
        #_obs.columns = ['chrom', 'pos', 'strand', 'context']
        #acc_adata = anndata.AnnData(
        #    X=X,
        #    var=_var,
        #    obs=_obs
        #)
        #acc_adata.write_h5ad(
        #    Path(adir) / 'ACC.h5'
        #)
        # # Meth part.
        # _samples = list(pl.read_ipc_schema(self.meth).keys())[1::]
        # _var = pd.DataFrame(index=_samples)
        # X = np.ma.masked_invalid(pl.read_ipc(self.meth, memory_map=False, columns=_samples).to_pandas().values)
        # _obs = pl.read_ipc(self.meth, memory_map=False, columns=['index']).to_pandas()['index'].str.split('_', expand=True)
        # _obs.columns = ['chrom', 'pos', 'strand', 'context']
        # meth_accdata = anndata.AnnData(
        #     X=X,
        #     var=_var,
        #     obs=_obs
        # )
        # meth_accdata.write_h5ad(
        #     Path(adir) / 'METH.h5'
        # )
        

        
        #mameth = np.ma.masked_invalid()
        #print(mameth.values)
        #_methsp = csr_matrix(methdf.values)
        #print(_methsp)
        
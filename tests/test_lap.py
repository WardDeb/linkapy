from linkapy.parsing import Linkapy_Parser

def test_Linkapy_Parser(tmp_path):
    reg = tmp_path / "reg.bed"
    reg.touch()
    cs = tmp_path / "chrom.sizes"
    cs.touch()
    methf = tmp_path / "meth.GC.tsv.gz"
    methf.touch()
    rnaf = tmp_path / "rna.tsv"
    rnaf.touch()

    lp = Linkapy_Parser(
        methylation_path = str(tmp_path),
        transcriptome_path = str(tmp_path),
        output = str(tmp_path),
        methylation_pattern = ('*GC*',),
        transcriptome_pattern = ('*tsv',),
        NOMe = False,
        threads = 10,
        chromsizes = (str(cs),),
        regions = (str(reg),),
        blacklist = (),
        binsize = 10000,
        project = 'linkapy_API'
    )
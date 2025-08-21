from linkapy.parsing import Linkapy_Parser
import mudata as md

def test_parser_chromsizes(tmp_path, bed_files, methylation_files, rna_files):
    lp = Linkapy_Parser(
        methylation_path = str(tmp_path / 'methylation'),
        transcriptome_path = str(tmp_path / 'transcriptome'),
        output = str(tmp_path / 'output'),
        methylation_pattern = ('*WCGN*tsv.gz', '*GCHN*tsv.gz'),
        methylation_pattern_names = ('WCGN', 'GCHN'),
        transcriptome_pattern = ('*tsv',),
        NOMe = False,
        threads = 2,
        chromsizes = str(tmp_path / 'regions' / 'genome.chromsizes'),
        regions = (),
        blacklist = (),
        binsize = 20,
        project = 'chromsizes_test'
    )
    lp.parse()

    # Output checks.
    mo = tmp_path / 'output' / 'chromsizes_test.h5mu'
    assert mo.exists(), "muData object not created."
    mu = md.read(mo)
    print(mu)
    assert mu.shape == (3,52), f"Expected shape (3, 52), got {mu.shape}."
    assert set(mu.obs.index) == set(['cell1', 'cell2', 'cell3']), f"Obs inferral failed, got {mu.obs.index}."
    # Assert values in the muData object.
    mu['METH_WCGN'].X.todense()[0, 0] == 0, f"WCGN cell1 in bin 0-20 should be 0.0. Got {mu['METH_WCGN'].X.todense()[0, 0]}."
    mu['METH_WCGN'].X.todense()[0, 1] == 0.5, f"WCGN cell1 in bin 20-40 should be 0.5. Got {mu['METH_WCGN'].X.todense()[0, 1]}."
    mu['METH_WCGN'].X.todense()[0, 2] == 0.5, f"WCGN cell1 in bin 40-60 should be 0.5. Got {mu['METH_WCGN'].X.todense()[0, 2]}."

def test_parser_regions(tmp_path, bed_files, methylation_files, rna_files):
    lp = Linkapy_Parser(
        methylation_path = str(tmp_path / 'methylation'),
        transcriptome_path = str(tmp_path / 'transcriptome'),
        output = str(tmp_path / 'output'),
        methylation_pattern = ('*WCGN*tsv.gz', '*GCHN*tsv.gz'),
        methylation_pattern_names = ('WCGN', 'GCHN'),
        transcriptome_pattern = ('*tsv',),
        NOMe = False,
        threads = 2,
        chromsizes = None,
        regions = (str(tmp_path / 'regions' / 'gene1.bed'), str(tmp_path / 'regions' / 'gene2.bed')),
        blacklist = (),
        binsize = 20,
        project = 'chromsizes_test'
    )
    lp.parse()

    # Output checks.
    mo = tmp_path / 'output' / 'chromsizes_test.h5mu'
    assert mo.exists(), "muData object not created."
    mu = md.read(mo)
    print(mu)
    assert mu.shape == (3,6), f"Expected shape (3, 6), got {mu.shape}."
    assert set(mu.obs.index) == set(['cell1', 'cell2', 'cell3']), f"Obs inferral failed, got {mu.obs.index}."

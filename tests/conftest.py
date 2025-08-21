import pytest
from pathlib import Path
import gzip

'''
Create fixtures for testing.
The genome here is considered to be 1 chromosome (chr1) with 500 bps in total.
There are 5 methylation (WCGN context) sites:
 - pos 20 TCGG
 - pos 50 ACGT
 - pos 221 ACGA
 - pos 244 ACGA
 - pos 440 ACGA
and 4 methylation (GCHN context) sites:
 - pos 30 GCAA
 - pos 40 GCTC
 - pos 230 GCCC
 - pos 401 GCAG
There are 2 genes in total:
 - gene1: 0-100 bps
 - gene2: 200-250 bps

'''

@pytest.fixture
def methylation_files(tmp_path):
    '''
    simulate allcools.tsv.gz files
    tsv with columns chr, pos, strand, context, meth, cov and sign.
    sign denotes if significantly methylated (1 = no test).
    '''
    methpath = tmp_path / "methylation"
    methpath.mkdir()

    # Cell 1 - avg. meth in gene 1 = 50%
    with gzip.open(methpath / "cell1.WCGN.tsv.gz", "wt") as f:
        f.write('chr1\t20\t+\tTCGG\t25\t50\t1\n')
        f.write('chr1\t50\t+\tACGT\t50\t100\t1\n')

    # Cell 2 - avg. meth in gene 2 = 75%
    with gzip.open(methpath / "cell2.WCGN.tsv.gz", "wt") as f:
        f.write('chr1\t221\t+\tACGA\t75\t100\t1\n')
        f.write('chr1\t244\t+\tACGA\t750\t1000\t1\n')
    
    # Cell 3 - avg. meth in gene 1 = 100%, avg. meth in gene 2 = 0%
    with gzip.open(methpath / "cell3.WCGN.tsv.gz", "wt") as f:
        f.write('chr1\t20\t+\tTCGG\t21\t21\t1\n')
        f.write('chr1\t50\t+\tACGT\t35\t35\t1\n')
        f.write('chr1\t221\t+\tACGA\t0\t100\t1\n')
    
    # Cell 1 - avg. meth in gene 1 = 25%
    with gzip.open(methpath / "cell1.GCHN.tsv.gz", "wt") as f:
        f.write('chr1\t30\t+\tGCAA\t25\t100\t1\n')
        f.write('chr1\t40\t+\tGCTC\t1\t4\t1\n')

    # Cell 2 - avg. meth in gene 2 = 44%
    with gzip.open(methpath / "cell2.GCHN.tsv.gz", "wt") as f:
        f.write('chr1\t230\t+\tGCCC\t44\t100\t1\n')
    
    # Cell 3 - avg. meth in gene 1 = 100%, avg. meth in gene 2 = 0%
    with gzip.open(methpath / "cell3.GCHN.tsv.gz", "wt") as f:
        f.write('chr1\t30\t+\tGCAA\t100\t100\t1\n')
        f.write('chr1\t40\t+\tGCTC\t4\t4\t1\n')
        f.write('chr1\t230\t+\tGCCC\t0\t100\t1\n')

    return methpath

@pytest.fixture
def rna_files(tmp_path):
    '''
    simulate transcriptome.tsv files
    tsv with columns Geneid, Chr, Start, End, Strand, Length
    '''
    rnapath = tmp_path / "transcriptome"
    rnapath.mkdir()

    with open(rnapath / 'count1.tsv', 'w') as f:
        f.write("#Program featureCounts comment line.\n")
        f.write('Geneid\tChr\tStart\tEnd\tStrand\tLength\tcell1_rna\tcell2_rna\n')
        f.write('gene1\tchr1\t0\t100\t+\t101\t50\t100\n')
        f.write('gene2\tchr1\t200\t250\t+\t51\t20\t100\n')
    
    with open(rnapath / 'count2.tsv', 'w') as f:
        f.write("#Program featureCounts comment line.\n")
        f.write('Geneid\tChr\tStart\tEnd\tStrand\tLength\tcell3_rna\n')
        f.write('gene1\tchr1\t0\t100\t+\t101\t33\n')
        f.write('gene2\tchr1\t200\t250\t+\t51\t33\n')
    
    return rnapath

@pytest.fixture
def bed_files(tmp_path):
    '''
    simulate regions.bed files
    '''
    bedpath = tmp_path / "regions"
    bedpath.mkdir()

    with open(bedpath / "gene1.bed", "w") as f:
        f.write('chr1\t0\t100\tgene1\t0\t+\n')

    # gene2: 200-250 bps
    with open(bedpath / "gene2.bed", "w") as f:
        f.write('chr1\t200\t250\tgene2\t0\t+\n')

    # blacklist regions
    with open(bedpath / "blacklist1.bed", 'w') as f:
        f.write('chr1\t40\t60\n')

    with open(bedpath / "blacklist2.bed", 'w') as f:
        f.write('chr1\t235\t250\n')

    with open(bedpath / "genome.chromsizes", 'w') as f:
        f.write('chr1\t500\n')

    return bedpath

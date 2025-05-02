use pyo3::prelude::*;
use pyo3::types::PyList;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use flate2::read::GzDecoder;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::cmp::Ordering;
use sprs::{CsMat, TriMat};
use sprs::io::write_matrix_market;

#[pyfunction]
pub fn parse_cools(
    _coolfiles: Py<PyList>,
    _regions: Py<PyList>,
    _regionlabels: Py<PyList>,
    qc: bool,
    threads: usize,
    obase: &str,
    oregionfile: &str,
    ocellfile: &str,
) -> PyResult<()> {
    let mut coolfiles: Vec<String> = Vec::new();
    let mut regions: Vec<String> = Vec::new();
    let mut regionlabels: Vec<String> = Vec::new();
    Python::with_gil(|py| {
        coolfiles = _coolfiles.extract(py).expect("Failed to retrieve allcoolfiles.");
        regions = _regions.extract(py).expect("Failed to retrieve regions.");
        regionlabels = _regionlabels.extract(py).expect("Failed to retrieve region labels.");
    });
    assert_eq!(regions.len(), regionlabels.len());
    println!("Parsing {} region files.", regions.len());
    // Parse regions.
    let mut parsed_regions: Vec<Region> = Vec::new();
    for (_r, _l) in regions.into_iter().zip(regionlabels.into_iter()) {
        parsed_regions.extend(parse_region(_r, _l));
    }
    parsed_regions.sort_by(|a, b| {
        // First, compare by `chrom`
        let chrom_order = a.chrom.cmp(&b.chrom);
        if chrom_order != Ordering::Equal {
            return chrom_order;
        }

        // If chrom is the same, compare by the first value in `start`
        a.start[0].cmp(&b.start[0])
    });
    println!("Found {} regions.", parsed_regions.len());
    println!("Launching a pool with {} threads to parse allcool files.", threads);
    let pool = ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
    if qc{
        // QC
        // Collapse the cells to keep things lean and mean
        println!("Starting QC.");
        let mut sum_mat: Vec<Vec<f32>> = vec![vec![f32::NAN; 301]; parsed_regions.len()];
        let mut count_mat: Vec<Vec<f32>> = vec![vec![f32::NAN; 301]; parsed_regions.len()];

        let _ = coolfiles
            .iter()
            .for_each(|coolfile| {
                let coolregions = parse_cool(coolfile);
                for (ix, region) in parsed_regions.iter().enumerate() {
                    coolregions
                        .iter()
                        .filter(|x| x.chrom == region.chrom && x.pos >= region.start[0] -100 && x.pos <= (*region.end.last().unwrap() + 100))
                        .for_each(|x| {
                            // We know these have value now. 
                            let frac = x.meth as f32 / x.total as f32;
                            // map pos between 0 and 100
                            let pos;
                            if x.pos >= region.start[0] && x.pos <= *region.end.last().unwrap() {
                                pos = ((x.pos - region.start[0]) as f32 / (*region.end.last().unwrap() - region.start[0]) as f32 * 99.0).round() as usize;
                                if sum_mat[ix][pos+100].is_nan() {
                                    sum_mat[ix][pos+100] = frac;
                                    count_mat[ix][pos+100] = 1.0;
                                } else {
                                    sum_mat[ix][pos+100] += frac;
                                    count_mat[ix][pos+100] += 1.0;
                                }
                            } else {
                                if x.pos < region.start[0] {
                                    pos = (region.start[0] - x.pos) as usize;
                                } else {
                                    pos = (200 + x.pos - region.end.last().unwrap()) as usize;
                                }
                                if sum_mat[ix][pos].is_nan() {
                                    sum_mat[ix][pos] = frac;
                                    count_mat[ix][pos] = 1.0;
                                } else {
                                    sum_mat[ix][pos] += frac;
                                    count_mat[ix][pos] += 1.0;
                                }
                            }
                        });
                }
            });
        let mut mean_mat = TriMat::new((parsed_regions.len(), 301));
        for (ix, row) in sum_mat.iter().enumerate() {
            for (pos, &sum_val) in row.iter().enumerate() {
                let count_val = count_mat[ix][pos];
                if !count_val.is_nan() {
                    mean_mat.add_triplet(ix, pos, sum_val / count_val);
                }
            }
        }

        let oqc = format!("{}.qc.mtx", obase);
        write_matrix_market(oqc, &mean_mat).unwrap();
    }

    {
        // Metrics
        let regvals: Vec<Vec<(f32, f32, f32)>> = pool.install(|| {
            coolfiles
                .par_iter()
                .map(|coolfile| {
                    let coolregions = parse_cool(coolfile);
                    parsed_regions
                        .par_iter()
                        .map(|region| {
                            let (meth_sum, total_sum , sites) = coolregions
                                .iter()
                                .filter(|x| x.chrom == region.chrom && x.pos >= region.start[0] && x.pos <= *region.end.last().unwrap())
                                .fold((f32::NAN, f32::NAN, f32::NAN), |(meth_acc, total_acc, sites), x| {
                                    (
                                        if meth_acc.is_nan() { x.meth as f32 } else { meth_acc + x.meth as f32 },
                                        if total_acc.is_nan() { x.total as f32 } else { total_acc + x.total as f32 },
                                        if sites.is_nan() { 1.0 } else { sites + 1.0 },
                                    )
                                });
                            (meth_sum, total_sum, sites)
                        })
                .collect()
                })
        .collect()
        });
        let (methm, covm, sitem) = tupvec_to_sparse(regvals);
        println!("Finished parsing allcool files.");
        // Create three outfiles from obase
        let ometh = format!("{}.meth.mtx", obase);
        let ocov = format!("{}.cov.mtx", obase);
        let osite = format!("{}.site.mtx", obase);
        write_matrix_market(ometh, &methm).unwrap();
        write_matrix_market(ocov, &covm).unwrap();
        write_matrix_market(osite, &sitem).unwrap();
        println!("Matrices written.");
        println!("Writing metadata to {}.", oregionfile);
        let mut ofile = File::create(oregionfile).unwrap();
        writeln!(ofile, "chrom\tstart\tend\tname\tclass").unwrap();
        for region in parsed_regions {
            writeln!(ofile, "{}\t{}\t{}\t{}\t{}", region.chrom, region.start[0], *region.end.last().unwrap(), region.name, region.class).unwrap();
        }
        let mut ofile = File::create(ocellfile).unwrap();
        for coolfile in coolfiles {
            writeln!(ofile, "{}", coolfile).unwrap();
        }
    }
    Ok(())
}

fn tupvec_to_sparse(dense: Vec<Vec<(f32, f32, f32)>>) -> (CsMat<f32>, CsMat<f32>, CsMat<f32>) {
    let max_row = dense.len();
    let max_col = dense.iter().map(|row| row.len()).max().unwrap_or(0);

    let mut mat1 = TriMat::new((max_row, max_col));
    let mut mat2 = TriMat::new((max_row, max_col));
    let mut mat3 = TriMat::new((max_row, max_col));

    for (i, row) in dense.iter().enumerate() {
        for (j, &(v1, v2, v3)) in row.iter().enumerate() {
            if !v1.is_nan() {
                mat1.add_triplet(i, j, v1);
            }
            if !v2.is_nan() {
                mat2.add_triplet(i, j, v2);
            }
            if !v3.is_nan() {
                mat3.add_triplet(i, j, v3);
            }
        }
    }
    (mat1.to_csr(), mat2.to_csr(), mat3.to_csr())
}

fn parse_cool(_f: &str) -> Vec<CoolRegion> {
    let mut coolregions: Vec<CoolRegion> = Vec::new();
    let reader = BufReader::new(GzDecoder::new(File::open(_f).unwrap()));
    for line in reader.lines() {
        let line = line.unwrap();
        let fields: Vec<&str> = line.split('\t').collect();
        let chrom = fields[0].to_string();
        let pos = fields[1].parse::<u32>().unwrap();
        let meth = fields[4].parse::<u32>().unwrap();
        let total = fields[5].parse::<u32>().unwrap();
        coolregions.push(
            CoolRegion{
                chrom: chrom,
                pos: pos,
                meth: meth,
                total: total,
            }
        );
    }
    coolregions
}

fn parse_region(reg: String, class: String) -> Vec<Region> {
    let mut regions = Vec::new();
    // Get suffix from reg
    let suffix = reg.split('.').last().unwrap();
    // Two options: gz (bed.gz), bed(bed)
    match suffix {
        "gz" => {
            let reader = BufReader::new(GzDecoder::new(File::open(reg).unwrap()));
            let lines = reader.lines();
            for line in lines {
                let line = line.unwrap();
                let fields: Vec<&str> = line.split('\t').collect();
                let chrom = fields[0].to_string();
                let start = fields[1].to_string();
                let end = fields[2].to_string();
                let name: String;
                if fields.len() > 3 {
                    name = fields[3].to_string();
                } else {
                    name = format!("{}:{}-{}", chrom, start, end);
                }
                // check if start, end have commas
                if start.contains(",") {
                    let start: Vec<u32> = start.split(',').map(|x| x.parse::<u32>().unwrap()).collect();
                    let end: Vec<u32> = end.split(',').map(|x| x.parse::<u32>().unwrap()).collect();
                    regions.push(
                        Region{
                            chrom: chrom,
                            start: start,
                            end: end,
                            name: name,
                            class: class.to_string()
                        }
                    );
                } else {
                    let start = start.parse::<u32>().unwrap();
                    let end = end.parse::<u32>().unwrap();
                    regions.push(
                        Region{
                            chrom: chrom,
                            start: vec![start],
                            end: vec![end],
                            name: name,
                            class: class.to_string()
                        }
                    );
                }
            }
            regions
        },
        "bed" => {
            let reader = BufReader::new(File::open(reg).unwrap());
            let lines = reader.lines();
            for line in lines {
                let line = line.unwrap();
                let fields: Vec<&str> = line.split('\t').collect();
                let chrom = fields[0].to_string();
                let start = fields[1].to_string();
                let end = fields[2].to_string();
                let name: String;
                if fields.len() > 3 {
                    name = fields[3].to_string();
                } else {
                    name = format!("{}:{}-{}", chrom, start, end);
                }
                // check if start, end have commas
                if start.contains(",") {
                    let start: Vec<u32> = start.split(',').map(|x| x.parse::<u32>().unwrap()).collect();
                    let end: Vec<u32> = end.split(',').map(|x| x.parse::<u32>().unwrap()).collect();
                    regions.push(
                        Region{
                            chrom: chrom,
                            start: start,
                            end: end,
                            name: name,
                            class: class.to_string()
                        }
                    );
                } else {
                    let start = start.parse::<u32>().unwrap();
                    let end = end.parse::<u32>().unwrap();
                    regions.push(
                        Region{
                            chrom: chrom,
                            start: vec![start],
                            end: vec![end],
                            name: name,
                            class: class.to_string()
                        }
                    );
                }
            }
            regions
        },
        _ => {panic!("File format not supported");}
    }
}

struct Region {
    chrom: String,
    start: Vec<u32>,
    end: Vec<u32>,
    name: String,
    class: String,
}

struct CoolRegion {
    chrom: String,
    pos: u32,
    meth: u32,
    total: u32,
}
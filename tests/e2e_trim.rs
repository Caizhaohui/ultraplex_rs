// IMPORTANT: DO NOT ADD ANY COMMENTS

use std::fs;
use std::io::Write;
use flate2::read::GzDecoder;
use std::io::Read;
use ultraplex_rs::cli::{Args, run};
use ultraplex_rs::align::prefix_match;

#[test]
fn e2e_trim_small_fastq() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("small.fastq");
    let mut f = fs::File::create(&input).unwrap();
    writeln!(f, "@r1").unwrap();
    writeln!(f, "ACGTGGGG").unwrap();
    writeln!(f, "+").unwrap();
    writeln!(f, "IIIIIIII").unwrap();

    let args = Args { inputfastq: input.to_str().unwrap().to_string(), directory: dir.path().to_str().unwrap().to_string(), barcodes: String::new(), outputprefix: "demux".to_string(), nextseq: false, gzip: false, three_prime_only: false, input_2: String::new(), threeprimemismatches: 0, threads: 1, keep_barcode: false, final_min_length: 0, ignore_no_match: false, phredquality: 30 };
    run(args).unwrap();

    let out = dir.path().join("ultraplex_demux_no_match.fastq");
    let content = fs::read_to_string(out).unwrap();
    assert!(content.contains("@r1"));
}

#[test]
fn e2e_test_simple_single_end_named_if_available() {
    let root = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).parent().unwrap().to_path_buf();
    let reads = root.join("tests/test_simple/reads1.fastq.gz");
    let barcodes = root.join("tests/test_simple/barcodes_5_and_3_named.csv");
    if !(reads.exists() && barcodes.exists()) { return; }
    let out_dir = tempfile::tempdir().unwrap();
    let args = Args { inputfastq: reads.to_str().unwrap().to_string(), directory: out_dir.path().to_str().unwrap().to_string(), barcodes: barcodes.to_str().unwrap().to_string(), outputprefix: "single_end_named".to_string(), nextseq: false, gzip: true, three_prime_only: false, input_2: String::new(), threeprimemismatches: 0, threads: 2, keep_barcode: false, final_min_length: 0, ignore_no_match: false, phredquality: 30 };
    run(args).unwrap();
    let files = std::fs::read_dir(out_dir.path()).unwrap().map(|e| e.unwrap().path()).collect::<Vec<_>>();
    assert!(files.iter().any(|p| p.file_name().unwrap().to_string_lossy().contains("ultraplex_single_end_named_")));
}

#[test]
fn e2e_test_simple_length_filter_if_available() {
    let root = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).parent().unwrap().to_path_buf();
    let reads = root.join("tests/test_simple/reads1.fastq.gz");
    let barcodes = root.join("tests/test_simple/barcodes_5_and_3.csv");
    if !(reads.exists() && barcodes.exists()) { return; }
    let out_dir = tempfile::tempdir().unwrap();
    let args = Args { inputfastq: reads.to_str().unwrap().to_string(), directory: out_dir.path().to_str().unwrap().to_string(), barcodes: barcodes.to_str().unwrap().to_string(), outputprefix: "single_end_length60".to_string(), nextseq: false, gzip: true, three_prime_only: false, input_2: String::new(), threeprimemismatches: 0, threads: 2, keep_barcode: false, final_min_length: 60, ignore_no_match: false, phredquality: 30 };
    run(args).unwrap();
    let nm = out_dir.path().join("ultraplex_single_end_length60_no_match.fastq.gz");
    assert!(nm.exists());
}

#[test]
fn e2e_demux_with_barcodes_and_gzip() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("small.fastq");
    let mut f = fs::File::create(&input).unwrap();
    writeln!(f, "@r1").unwrap();
    writeln!(f, "ACGTGGGG").unwrap();
    writeln!(f, "+").unwrap();
    writeln!(f, "IIIIIIII").unwrap();

    let bcsv = dir.path().join("barcodes.csv");
    let mut bf = fs::File::create(&bcsv).unwrap();
    writeln!(bf, "ACGT").unwrap();

    let args = Args { inputfastq: input.to_str().unwrap().to_string(), directory: dir.path().to_str().unwrap().to_string(), barcodes: bcsv.to_str().unwrap().to_string(), outputprefix: "demux".to_string(), nextseq: false, gzip: true, three_prime_only: false, input_2: String::new(), threeprimemismatches: 0, threads: 1, keep_barcode: false, final_min_length: 0, ignore_no_match: false, phredquality: 30 };
    run(args).unwrap();

    let out_bc = dir.path().join("ultraplex_demux_ACGT.fastq.gz");
    assert!(out_bc.exists());
    let mut gz = GzDecoder::new(fs::File::open(out_bc).unwrap());
    let mut content = String::new();
    gz.read_to_string(&mut content).unwrap();
    assert!(content.contains("@r1"));
    let lines: Vec<&str> = content.lines().collect();
    assert!(prefix_match(lines[1].as_bytes(), b"ACGT", 0));
}

#[test]
fn e2e_three_prime_only_with_umi_and_naming() {
    let dir = tempfile::tempdir().unwrap();
    let r1 = dir.path().join("r1.fastq");
    let r2 = dir.path().join("r2.fastq");
    {
        let mut f1 = fs::File::create(&r1).unwrap();
        writeln!(f1, "@r1").unwrap();
        writeln!(f1, "TTTTCCCC").unwrap();
        writeln!(f1, "+").unwrap();
        writeln!(f1, "IIIIIIII").unwrap();
        let mut f2 = fs::File::create(&r2).unwrap();
        writeln!(f2, "@r2").unwrap();
        writeln!(f2, "ACGTNNNN").unwrap();
        writeln!(f2, "+").unwrap();
        writeln!(f2, "IIIIIIII").unwrap();
    }

    let bcsv = dir.path().join("barcodes.csv");
    {
        let mut bf = fs::File::create(&bcsv).unwrap();
        writeln!(bf, "ACGT,NNNN:sampleX").unwrap();
    }

    let args = Args { inputfastq: r2.to_str().unwrap().to_string(), directory: dir.path().to_str().unwrap().to_string(), barcodes: bcsv.to_str().unwrap().to_string(), outputprefix: "demux".to_string(), nextseq: false, gzip: true, three_prime_only: true, input_2: r1.to_str().unwrap().to_string(), threeprimemismatches: 0, threads: 1, keep_barcode: false, final_min_length: 0, ignore_no_match: false, phredquality: 30 };
    run(args).unwrap();

    let out = dir.path().join("ultraplex_demux_sampleX.fastq.gz");
    assert!(out.exists());
    let mut gz = GzDecoder::new(fs::File::open(out).unwrap());
    let mut content = String::new();
    gz.read_to_string(&mut content).unwrap();
    assert!(content.contains("rbc:"));
}
// IMPORTANT: DO NOT ADD ANY COMMENTS

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use csv::ReaderBuilder;
use flate2::write::GzEncoder;
use flate2::Compression;
use needletail::parse_fastx_file;

#[derive(Clone)]
pub struct BarcodeSet {
    pub five_prime: Vec<String>,
    pub three_prime: Vec<String>,
    pub linked: std::collections::HashMap<String, Vec<String>>,
    pub sample_names: std::collections::HashMap<String, String>,
    pub three_p_mismatches: usize,
}

pub fn read_barcodes_csv(path: &str) -> anyhow::Result<BarcodeSet> {
    let mut rdr = ReaderBuilder::new().has_headers(false).from_path(path)?;
    let mut five_prime = Vec::new();
    let mut three_prime = Vec::new();
    let mut linked: HashMap<String, Vec<String>> = HashMap::new();
    let mut sample_names: HashMap<String, String> = HashMap::new();
    let mut fivelength: Option<usize> = None;
    for result in rdr.records() {
        let rec = result?;
        let first = rec.get(0).unwrap_or("").trim().to_uppercase();
        if first.is_empty() { continue; }
        let five_bc = first.split(':').next().unwrap().to_string();
        if fivelength.is_none() { fivelength = Some(five_bc.replace('N', "").len()); } else {
            assert_eq!(fivelength.unwrap(), five_bc.replace('N', "").len());
        }
        five_prime.push(five_bc.clone());
        let mut three_list = Vec::new();
        for i in 1..rec.len() {
            let col = rec.get(i).unwrap_or("").trim();
            if col.is_empty() { continue; }
            let parts: Vec<&str> = col.split(':').collect();
            let bc = parts[0].to_uppercase();
            three_prime.push(bc.clone());
            three_list.push(bc.clone());
            if parts.len() > 1 && !parts[1].is_empty() {
                sample_names.insert(format!("5bc_{}_3bc_{}", five_bc, bc), parts[1].to_string());
            }
        }
        if !three_list.is_empty() { linked.insert(five_bc, three_list); }
    }
    five_prime.sort(); five_prime.dedup();
    three_prime.sort(); three_prime.dedup();
    Ok(BarcodeSet { five_prime, three_prime, linked, sample_names, three_p_mismatches: 0 })
}

pub struct Writers {
    pub default: Box<dyn Write + Send>,
    pub by_barcode: HashMap<String, Box<dyn Write + Send>>,
}

pub fn create_writers(output_dir: &str, prefix: &str, barcodes: &[String], gz: bool) -> anyhow::Result<Writers> {
    let dir = if output_dir.is_empty() { PathBuf::from(".") } else { PathBuf::from(output_dir) };
    if !dir.as_os_str().is_empty() && !dir.exists() { std::fs::create_dir_all(&dir)?; }
    let ext = if gz { "fastq.gz" } else { "fastq" };
    let default_file = File::create(dir.join(format!("ultraplex_{}_no_match.{}", prefix, ext)))?;
    let default: Box<dyn Write + Send> = if gz {
        Box::new(GzEncoder::new(BufWriter::new(default_file), Compression::default()))
    } else {
        Box::new(BufWriter::new(default_file))
    };
    let mut by_barcode = HashMap::new();
    for bc in barcodes {
        let f = File::create(dir.join(format!("ultraplex_{}_{}.{}", prefix, bc, ext)))?;
        let w: Box<dyn Write + Send> = if gz {
            Box::new(GzEncoder::new(BufWriter::new(f), Compression::default()))
        } else {
            Box::new(BufWriter::new(f))
        };
        by_barcode.insert(bc.clone(), w);
    }
    Ok(Writers { default, by_barcode })
}

pub fn write_fastq_record(w: &mut dyn Write, name: &[u8], seq: &[u8], qual: &[u8]) -> anyhow::Result<()> {
    w.write_all(b"@")?;
    w.write_all(name)?;
    w.write_all(b"\n")?;
    w.write_all(seq)?;
    w.write_all(b"\n+\n")?;
    w.write_all(qual)?;
    w.write_all(b"\n")?;
    Ok(())
}

pub fn open_fastx(path: &str) -> anyhow::Result<Box<dyn needletail::FastxReader>> { Ok(parse_fastx_file(path)?) }

pub fn get_writer<'a>(writers: &'a mut Writers, output_dir: &str, prefix: &str, key: &str, gz: bool) -> &'a mut (dyn Write + Send) {
    if !writers.by_barcode.contains_key(key) {
        let dir = if output_dir.is_empty() { PathBuf::from(".") } else { PathBuf::from(output_dir) };
        let ext = if gz { "fastq.gz" } else { "fastq" };
        let f = File::create(dir.join(format!("ultraplex_{}_{}.{}", prefix, key, ext))).expect("create writer");
        let w: Box<dyn Write + Send> = if gz { Box::new(GzEncoder::new(BufWriter::new(f), Compression::default())) } else { Box::new(BufWriter::new(f)) };
        writers.by_barcode.insert(key.to_string(), w);
    }
    writers.by_barcode.get_mut(key).unwrap().as_mut()
}

pub fn rev_comp(seq: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(seq.len());
    for &b in seq.iter().rev() {
        out.push(match b { b'A' => b'T', b'T' => b'A', b'C' => b'G', b'G' => b'C', _ => b'N' });
    }
    out
}

pub fn check_n_position(bcds: &[String]) -> anyhow::Result<()> {
    if bcds.is_empty() { return Ok(()); }
    let ref_b = &bcds[0];
    let positions: Vec<usize> = ref_b.as_bytes().iter().enumerate().filter_map(|(i,&c)| if c!=b'N' { Some(i) } else { None }).collect();
    for b in bcds.iter() {
        let ps: Vec<usize> = b.as_bytes().iter().enumerate().filter_map(|(i,&c)| if c!=b'N' { Some(i) } else { None }).collect();
        if ps != positions { anyhow::bail!("UMI positions not consistent"); }
    }
    Ok(())
}
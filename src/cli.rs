// IMPORTANT: DO NOT ADD ANY COMMENTS

use clap::Parser;
use crate::trim::{quality_trim_index, nextseq_trim_index};
use crate::demux::{read_barcodes_csv, create_writers, write_fastq_record, open_fastx, rev_comp, check_n_position};
use crate::align::{prefix_match, suffix_match, extract_umi_from_suffix};
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use anyhow::Result;

#[derive(Parser, Debug, Clone)]
pub struct Args {
    #[arg(short = 'i', long = "inputfastq", help = "输入 FASTQ 文件路径，支持 .fastq 与 .fastq.gz；推荐 gzip 压缩")]
    pub inputfastq: String,
    #[arg(short = 'd', long = "directory", default_value = "", help = "输出目录，仅用于写出结果文件；不批量读取目录内文件")]
    pub directory: String,
    #[arg(short = 'b', long = "barcodes", default_value = "", help = "条码 CSV。首列 5’ 条码，后续列为链接的 3’ 条码；支持 :样本名")]
    pub barcodes: String,
    #[arg(short = 'o', long = "outputprefix", default_value = "demux", help = "输出前缀，用于命名 ultraplex_<prefix>_<key>.fastq[.gz]")]
    pub outputprefix: String,
    #[arg(long = "nextseq", default_value_t = false, help = "启用 NextSeq 风格的质量修剪（主要针对 3’ 端低质位）")]
    pub nextseq: bool,
    #[arg(long = "gzip", default_value_t = false, help = "以 .fastq.gz 格式写出结果文件")]
    pub gzip: bool,
    #[arg(long = "three_prime_only", default_value_t = false, help = "启用 3’ 条码末端匹配与 UMI 抽取（结合 5’ 前缀）")]
    pub three_prime_only: bool,
    #[arg(short = 'I', long = "input_2", default_value = "", help = "成对测序第二个 FASTQ 路径（预留；后续扩展）")]
    pub input_2: String,
    #[arg(short = 'M', long = "threeprimemismatches", default_value_t = 0, help = "3’ 条码末端匹配允许的错配数")]
    pub threeprimemismatches: usize,
    #[arg(short = 't', long = "threads", default_value_t = 4, help = "并行处理线程数")]
    pub threads: usize,
    #[arg(long = "keep_barcode", default_value_t = false, help = "匹配到 3’ 条码后是否保留条码本体在序列中")]
    pub keep_barcode: bool,
    #[arg(short = 'l', long = "final_min_length", default_value_t = 0, help = "长度过滤阈值，短于该长度的读将跳过写出")]
    pub final_min_length: usize,
    #[arg(long = "ignore_no_match", default_value_t = false, help = "忽略无匹配的读（不写入 no_match 文件）")]
    pub ignore_no_match: bool,
    #[arg(short = 'q', long = "phredquality", default_value_t = 30, help = "质量修剪的 Phred 阈值（默认 30，ASCII 偏移 33）")]
    pub phredquality: i32,
}

pub fn run(args: Args) -> Result<()> {
    let mut out_dir = std::path::PathBuf::from(&args.directory);
    if !args.directory.is_empty() && !args.directory.ends_with('/') {
        out_dir = std::path::PathBuf::from(format!("{}/", args.directory));
    }
    if !out_dir.as_os_str().is_empty() && !out_dir.exists() { std::fs::create_dir_all(&out_dir)?; }

    let mut barcode_set = if !args.barcodes.is_empty() { Some(read_barcodes_csv(&args.barcodes)?) } else { None };
    if let Some(bcs) = &mut barcode_set { bcs.three_p_mismatches = args.threeprimemismatches; }
    if args.three_prime_only {
        if let Some(bcs) = &barcode_set { check_n_position(&bcs.three_prime)?; }
    }
    let mut writers = if let Some(bcs) = &barcode_set {
        create_writers(out_dir.to_str().unwrap_or(""), &args.outputprefix, &bcs.five_prime, args.gzip)?
    } else {
        create_writers(out_dir.to_str().unwrap_or(""), &args.outputprefix, &Vec::new(), args.gzip)?
    };

    let pool = ThreadPoolBuilder::new().num_threads(args.threads).build().unwrap();
    let mut reader = open_fastx(&args.inputfastq)?;
    let mut chunk: Vec<(Vec<u8>, Vec<u8>, Vec<u8>)> = Vec::with_capacity(1024);
    loop {
        match reader.next() {
            Some(Ok(rec)) => {
                let id = rec.id().to_vec();
                let seq = rec.seq().to_vec();
                let qual = rec.qual().map(|q| q.to_vec()).unwrap_or_default();
                chunk.push((id, seq, qual));
                if chunk.len() >= 1024 {
                    let results = pool.install(|| {
                        chunk.par_iter().map(|(id, seq, qual)| {
                            let qual_str = String::from_utf8(qual.clone()).unwrap_or_default();
                            let (start, stop) = if args.nextseq {
                                let idx = nextseq_trim_index(std::str::from_utf8(&seq).unwrap_or(""), &qual_str, args.phredquality, 33);
                                (0, idx)
                            } else { quality_trim_index(&qual_str, args.phredquality, args.phredquality, 33) };
                            let mut out_seq = seq[start..stop].to_vec();
                            let mut out_qual = qual[start..stop].to_vec();
                            let mut head = id.clone();
                            let mut key = String::from("no_match");
                            if let Some(bcs) = &barcode_set {
                                let mut matched = false;
                                if args.three_prime_only {
                                    for (five_bc, three_bcs) in bcs.linked.clone() {
                                        let five_rc = rev_comp(five_bc.as_bytes());
                                        if prefix_match(&out_seq, &five_rc, 0) {
                                            for three_bc in three_bcs {
                                                if suffix_match(&out_seq, three_bc.as_bytes(), bcs.three_p_mismatches) {
                                                    if let Some(umi) = extract_umi_from_suffix(&out_seq, three_bc.as_bytes()) { if !umi.is_empty() { head.extend_from_slice(b"rbc:"); head.extend_from_slice(&umi); } }
                                                    if !args.keep_barcode { let cut = out_seq.len() - three_bc.len(); out_seq = out_seq[..cut].to_vec(); out_qual = out_qual[..cut].to_vec(); }
                                                    let combo = format!("5bc_{}_3bc_{}", five_bc, three_bc);
                                                    key = if let Some(sample) = bcs.sample_names.get(&combo) { sample.clone() } else { combo };
                                                    matched = true; break;
                                                }
                                            }
                                            if matched { break; }
                                        }
                                    }
                                } else {
                                    for bc in &bcs.five_prime { if prefix_match(&out_seq, bc.as_bytes(), 0) { key = bc.clone(); matched = true; break; } }
                                }
                                if !matched && args.ignore_no_match { key = String::from("__skip__"); }
                            }
                            (key, head, out_seq, out_qual)
                        }).collect::<Vec<(String, Vec<u8>, Vec<u8>, Vec<u8>)>>()
                    });
                    for (key, head, out_seq, out_qual) in results {
                        if key == "__skip__" { continue; }
                        if out_seq.len() < args.final_min_length { continue; }
                        if key == "no_match" { write_fastq_record(writers.default.as_mut(), &head, &out_seq, &out_qual)?; } else {
                            let w = crate::demux::get_writer(&mut writers, out_dir.to_str().unwrap_or(""), &args.outputprefix, &key, args.gzip);
                            write_fastq_record(w, &head, &out_seq, &out_qual)?;
                        }
                    }
                    chunk.clear();
                }
            }
            Some(Err(e)) => return Err(e.into()),
            None => {
                if !chunk.is_empty() {
                    let results = pool.install(|| {
                        chunk.par_iter().map(|(id, seq, qual)| {
                            let qual_str = String::from_utf8(qual.clone()).unwrap_or_default();
                            let (start, stop) = if args.nextseq {
                                let idx = nextseq_trim_index(std::str::from_utf8(&seq).unwrap_or(""), &qual_str, args.phredquality, 33);
                                (0, idx)
                            } else { quality_trim_index(&qual_str, args.phredquality, args.phredquality, 33) };
                            let mut out_seq = seq[start..stop].to_vec();
                            let mut out_qual = qual[start..stop].to_vec();
                            let mut head = id.clone();
                            let mut key = String::from("no_match");
                            if let Some(bcs) = &barcode_set {
                                let mut matched = false;
                                if args.three_prime_only {
                                    for (five_bc, three_bcs) in bcs.linked.clone() {
                                        let five_rc = rev_comp(five_bc.as_bytes());
                                        if prefix_match(&out_seq, &five_rc, 0) {
                                            for three_bc in three_bcs {
                                                if suffix_match(&out_seq, three_bc.as_bytes(), bcs.three_p_mismatches) {
                                                    if let Some(umi) = extract_umi_from_suffix(&out_seq, three_bc.as_bytes()) { if !umi.is_empty() { head.extend_from_slice(b"rbc:"); head.extend_from_slice(&umi); } }
                                                    if !args.keep_barcode { let cut = out_seq.len() - three_bc.len(); out_seq = out_seq[..cut].to_vec(); out_qual = out_qual[..cut].to_vec(); }
                                                    let combo = format!("5bc_{}_3bc_{}", five_bc, three_bc);
                                                    key = if let Some(sample) = bcs.sample_names.get(&combo) { sample.clone() } else { combo };
                                                    matched = true; break;
                                                }
                                            }
                                            if matched { break; }
                                        }
                                    }
                                } else {
                                    for bc in &bcs.five_prime { if prefix_match(&out_seq, bc.as_bytes(), 0) { key = bc.clone(); matched = true; break; } }
                                }
                                if !matched && args.ignore_no_match { key = String::from("__skip__"); }
                            }
                            (key, head, out_seq, out_qual)
                        }).collect::<Vec<(String, Vec<u8>, Vec<u8>, Vec<u8>)>>()
                    });
                    for (key, head, out_seq, out_qual) in results {
                        if key == "__skip__" { continue; }
                        if out_seq.len() < args.final_min_length { continue; }
                        if key == "no_match" { write_fastq_record(writers.default.as_mut(), &head, &out_seq, &out_qual)?; } else {
                            let w = crate::demux::get_writer(&mut writers, out_dir.to_str().unwrap_or(""), &args.outputprefix, &key, args.gzip);
                            write_fastq_record(w, &head, &out_seq, &out_qual)?;
                        }
                    }
                }
                break;
            }
        }
    }
    for (_, mut w) in writers.by_barcode.into_iter() { w.flush()?; }
    writers.default.flush()?;
    Ok(())
}

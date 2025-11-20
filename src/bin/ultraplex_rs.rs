use clap::Parser;
use ultraplex_rs::cli::{Args, run};

fn main() -> anyhow::Result<()> { env_logger::init(); let args = Args::parse(); run(args) }
use std::fs;
use std::path::PathBuf;
use std::io::Write;
use clap::CommandFactory;
use ultraplex_rs::cli::Args;

fn main() -> anyhow::Result<()> {
    let mut cmd = Args::command();
    let mut buf = Vec::new();
    cmd.write_help(&mut buf).unwrap();
    let help = String::from_utf8(buf).unwrap();
    let manifest_dir = env!("CARGO_MANIFEST_DIR");
    let mut readme_path = PathBuf::from(manifest_dir);
    readme_path.push("README.md");
    let content = fs::read_to_string(&readme_path)?;
    let start = "<!-- BEGIN:CLI_HELP -->";
    let end = "<!-- END:CLI_HELP -->";
    if let (Some(s), Some(e)) = (content.find(start), content.find(end)) {
        let before = &content[..s + start.len()];
        let after = &content[e..];
        let mut out = String::new();
        out.push_str(before);
        out.push_str("\n````\n");
        out.push_str(&help);
        out.push_str("\n````\n");
        out.push_str(after);
        let mut f = fs::File::create(&readme_path)?;
        f.write_all(out.as_bytes())?;
    }
    Ok(())
}
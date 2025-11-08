//! Reverse (and complement) array-like SAM tags  for reverse alignments.
use std::path::PathBuf;
use std::process;

use anyhow::{Error, Result};
use env_logger::Env;
use structopt::StructOpt;

use revtaglib::revtag;

#[derive(Clone, Debug, StructOpt)]
#[structopt(
    setting = structopt::clap::AppSettings::ColoredHelp,
    setting = structopt::clap::AppSettings::DeriveDisplayOrder,
    rename_all = "kebab-case",
    about
)]
struct Opt {
    /// Input SAM/BAM/CRAM file or stream [default: /dev/stdin]
    #[structopt(short = "i", long = "--input", parse(from_os_str))]
    input: Option<PathBuf>,

    /// Output SAM/BAM/CRAM file or stream [default: /dev/stdout]
    #[structopt(short = "o", long = "--output", parse(from_os_str))]
    output: Option<PathBuf>,

    /// SAM tags with array values to reverse
    #[structopt(long = "--rev")]
    rev: Vec<String>,

    /// SAM tags with array values to reverse complement
    #[structopt(long = "--revcomp")]
    revcomp: Vec<String>,
}

/// Main binary entrypoint.
#[cfg(not(tarpaulin_include))]
fn main() -> Result<(), Error> {
    let env = Env::default().default_filter_or("info");
    let opt = Opt::from_args();

    env_logger::Builder::from_env(env).init();

    // Convert "-" to None for stdin/stdout
    let input = opt.input.and_then(|p| {
        if p.to_str() == Some("-") {
            None
        } else {
            Some(p)
        }
    });

    let output = opt.output.and_then(|p| {
        if p.to_str() == Some("-") {
            None
        } else {
            Some(p)
        }
    });

    match revtag(input, output, opt.rev, opt.revcomp) {
        Ok(exit_code) => process::exit(exit_code),
        Err(except) => panic!("{}", except),
    }
}

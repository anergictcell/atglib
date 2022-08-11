#![doc = include_str!("../README.md")]
#![doc(
    html_logo_url = "https://raw.githubusercontent.com/anergictcell/atg/main/assets/logo_standard.png",
    html_favicon_url = "https://raw.githubusercontent.com/anergictcell/atg/main/assets/favicon.ico"
)]
#![cfg_attr(feature = "with-bench", feature(test))]
#[cfg(all(test, feature = "with-bench"))]
extern crate test;

pub mod bed;
pub mod fasta;
pub mod genepred;
pub mod genepredext;
pub mod gtf;
pub mod models;
pub mod qc;
pub mod refgene;
pub mod tests;
pub mod utils;

use crate::models::TranscriptRead;
use crate::models::Transcripts;
use crate::utils::errors::ReadWriteError;

pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Generic function to read transcript from any possible source
pub fn read_transcripts<R: TranscriptRead>(
    reader: Result<R, ReadWriteError>,
) -> Result<Transcripts, ReadWriteError> {
    match reader {
        Ok(mut r) => r.transcripts(),
        Err(err) => Err(err),
    }
}

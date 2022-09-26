//! Quality Control checks for Transcripts
//!
//! This module provides functions to sanity-check a transcript for various criteria
//!
//! | QC check | Explanation | Non-Coding vs Coding | requires Fasta File |
//! | --- | --- | --- | --- |
//! | Exon | Contains at least one exon | all | no |
//! | Correct CDS Length | The length of the CDS is divisible by 3 | Coding | no |
//! | Correct Start Codon | The CDS starts with `ATG` | Coding | yes |
//! | Correct Stop Codon | The CDS ends with a Stop codon `TAG`, `TAA`, or `TGA` | Coding | yes |
//! | No upstream Start Codon | The 5'UTR does not contain another start codon `ATG` (This test do not make sense biologically. It is totally fine for a transcript to have upstream `ATG` start cordons that are not utilized but the ribosome.) | Coding | yes |
//! | No upstream Stop Codon| The CDS does not contain another in-frame stop-codon | Coding | yes |
//! | No Start codon | The full exon sequence does not contain a start codon `ATG` (Biologically speaking, a non-coding transcript could have `ATG` start codons that are not utilized) | Non-Coding | yes |
//! | Correct Coordinates | The transcript is within the coordinates of the reference genome | all | yes |
//!
//!
//! # QcCheck
//!
//! For interactive use, the most convenient way is to the use the [`QcCheck`] struct to run the full QC suite
//! on a transcript.
//!
//! ```rust
//! use atglib::tests::transcripts::standard_transcript;
//! use atglib::fasta::FastaReader;
//! use atglib::qc::{QcCheck, QcResult};
//! use atglib::models::GeneticCode;
//!
//! let tx = standard_transcript();
//!
//! let mut fasta_reader = FastaReader::from_file("tests/data/small.fasta").unwrap();
//!
//! let code = GeneticCode::default();
//!
//! let qc = QcCheck::new(&tx, &mut fasta_reader, &code);
//! assert_eq!(qc.correct_start_codon(), QcResult::OK)
//! ```
//!
//! # Writer
//!
//! An even easier method is provided via the [`Writer`] struct to run QC checks and save them
//! to an output file right away
//!
//! ```rust
//! use std::io;
//! use atglib::tests;;
//! use atglib::qc::Writer;
//! use atglib::fasta::FastaReader;
//! use atglib::models::TranscriptWrite;
//!
//! let transcripts = vec![tests::transcripts::standard_transcript()];
//!
//! let output = Vec::new(); // substitute this with proper IO (io::stdout())
//! let mut writer = Writer::new(output);
//! writer.fasta_reader(FastaReader::from_file("tests/data/small.fasta").unwrap());
//! writer.write_transcript_vec(&transcripts);
//!
//! let written_output = String::from_utf8(writer.into_inner().unwrap()).unwrap();
//! assert_eq!(written_output, "Test-Gene\tTest-Transcript\tOK\tNOK\tOK\tOK\tOK\tOK\tOK\n");
//! ```
//!
//! # Individual Tests
//!
//! Alternatively, every test can be run in isolation as well. This should only rarely be used
//! as it is less efficient due to a high IO-volume when reading the fasta reference genome
//! sequence for every test separately.
//!
//!
//! # Limitations
//! The QC checks are opinionated and assume the standard start codons `ATG`.
//! If you have use-cases for alternative options, please get in touch with me (open an issue on Github). I'm happy to improve
//! the functionality when needed. I would actually be excited to learn about such non-standard use-cases and implement them.
//!

mod writer;

use crate::fasta::FastaReader;
use crate::models::{CoordinateVector, GeneticCode, Sequence, Transcript};
use crate::utils::errors::FastaError;

pub use writer::Writer;

/// Holds the result of a QC check
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum QcResult {
    /// Test could not be performed (e.g. CDS-length for non-coding transcripts),
    /// so no conclusion could be drawn
    NA,
    /// The test succeeded with an OK results
    OK,
    /// The test failed and gave a NOT OK result
    NOK,
}

impl std::fmt::Display for QcResult {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                QcResult::OK => "OK",
                QcResult::NOK => "NOK",
                QcResult::NA => "N/A",
            }
        )
    }
}

impl From<bool> for QcResult {
    fn from(b: bool) -> Self {
        if b {
            QcResult::OK
        } else {
            QcResult::NOK
        }
    }
}

impl From<Option<bool>> for QcResult {
    fn from(b: Option<bool>) -> Self {
        match b {
            None => QcResult::NA,
            Some(true) => QcResult::OK,
            Some(false) => QcResult::NOK,
        }
    }
}

/// Wrapper struct to run the full QC-check suite on a transcript
///
/// # Examples
///
/// ```rust
/// use atglib::tests::transcripts::standard_transcript;
/// use atglib::fasta::FastaReader;
/// use atglib::qc::{QcCheck, QcResult};
/// use atglib::models::GeneticCode;
///
/// let tx = standard_transcript();
///
/// let mut fasta_reader = FastaReader::from_file("tests/data/small.fasta").unwrap();
///
/// let code = GeneticCode::default();
///
/// let qc = QcCheck::new(&tx, &mut fasta_reader, &code);
/// assert_eq!(qc.correct_start_codon(), QcResult::OK)
/// ```
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct QcCheck {
    exon: QcResult,
    cds_len: QcResult,
    correct_start: QcResult,
    correct_stop: QcResult,
    upstream_stop: QcResult,
    upstream_start: QcResult,
    correct_coordinates: QcResult,
}

impl std::default::Default for QcCheck {
    fn default() -> Self {
        QcCheck {
            exon: QcResult::NA,
            cds_len: QcResult::NA,
            correct_start: QcResult::NA,
            correct_stop: QcResult::NA,
            upstream_start: QcResult::NA,
            upstream_stop: QcResult::NA,
            correct_coordinates: QcResult::NA,
        }
    }
}

impl std::fmt::Display for QcCheck {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.exon,
            self.cds_len,
            self.correct_start,
            self.correct_stop,
            self.upstream_start,
            self.upstream_stop,
            self.correct_coordinates,
        )
    }
}

impl QcCheck {
    /// Runs the full QC-check suite on the transcript
    ///
    /// # Examples
    ///
    /// ```rust
    /// use atglib::tests::transcripts::standard_transcript;
    /// use atglib::fasta::FastaReader;
    /// use atglib::qc::{QcCheck, QcResult};
    /// use atglib::models::GeneticCode;
    ///
    /// let tx = standard_transcript();
    ///
    /// let mut fasta_reader = FastaReader::from_file("tests/data/small.fasta").unwrap();
    /// let code = GeneticCode::default();
    ///
    /// let qc = QcCheck::new(&tx, &mut fasta_reader, &code);
    /// assert_eq!(qc.correct_start_codon(), QcResult::OK)
    /// ```
    pub fn new(transcript: &Transcript, fasta: &mut FastaReader<impl std::io::Read + std::io::Seek>, code: &GeneticCode) -> Self {
        let mut res = QcCheck::default();

        let seq = match Sequence::from_coordinates(
            &transcript.cds_coordinates(),
            &transcript.strand(),
            fasta,
        ) {
            Ok(seq) => seq,
            Err(_) => {
                // if the sequence cannot be generated
                // the transcript coordinates must be wrong
                res.correct_coordinates = QcResult::NOK;
                return res;
            }
        };

        res.correct_coordinates = QcResult::OK;
        res.exon = contains_exon(transcript).into();
        if transcript.is_coding() {
            res.check_cds(transcript, &seq, code);
        };

        let coords = if transcript.is_coding() {
            transcript.utr5_coordinates()
        } else {
            transcript.cds_coordinates()
        };
        res.check_utr(transcript, coords, fasta);

        res
    }

    /// Does the transcript contain at least one exon
    pub fn contains_exon(&self) -> QcResult {
        self.exon
    }

    /// Is the length of the transcript's CDS a multiple of 3,
    /// i.e. has a proper coding frame
    pub fn correct_cds_length(&self) -> QcResult {
        self.cds_len
    }

    /// Is the first codon of the CDS an `ATG`
    pub fn correct_start_codon(&self) -> QcResult {
        self.correct_start
    }

    /// Is the last codon of the CDS a Stop-Codon
    pub fn correct_stop_codon(&self) -> QcResult {
        self.correct_stop
    }

    /// Does the transcript contain unused start codons in the 5'UTR
    /// (or UTR in general for non-coding transcripts)
    pub fn no_upstream_start_codon(&self) -> QcResult {
        self.upstream_start
    }

    /// Does the CDS contain an upstream in-frame stop-codon
    pub fn no_upstream_stop_codon(&self) -> QcResult {
        self.upstream_stop
    }

    /// Does the transcript lie within the reference genome coordinates
    pub fn correct_coordinates(&self) -> QcResult {
        self.correct_coordinates
    }

    fn check_cds(&mut self, transcript: &Transcript, seq: &Sequence, code: &GeneticCode) {
        self.cds_len = correct_cds_length(transcript).into();
        self.correct_start = starts_with_start_codon(seq).into();
        self.correct_stop = ends_with_stop_codon(seq, code).into();

        self.upstream_stop = match early_stop_codon(seq, code) {
            false => QcResult::OK,
            true => QcResult::NOK,
        };
    }

    fn check_utr(
        &mut self,
        transcript: &Transcript,
        utr: CoordinateVector,
        fasta: &mut FastaReader<impl std::io::Read + std::io::Seek>,
    ) {
        match Sequence::from_coordinates(&utr, &transcript.strand(), fasta) {
            Ok(seq) => {
                self.upstream_start = match extra_start_codon(&seq) {
                    false => QcResult::OK,
                    true => QcResult::NOK,
                };
            }
            Err(_) => {
                self.correct_coordinates = QcResult::NOK;
            }
        };
    }

    /// Create an "ideal" QcCheck with all tests passing
    ///
    /// There is no need to use this method, unless you're lazy and want
    /// to have a simple scaffold for unit-tests
    pub fn ideal() -> Self {
        QcCheck {
            exon: QcResult::OK,
            cds_len: QcResult::OK,
            correct_start: QcResult::OK,
            correct_stop: QcResult::OK,
            upstream_start: QcResult::OK,
            upstream_stop: QcResult::OK,
            correct_coordinates: QcResult::OK,
        }
    }
}

/// Returns true if the transcript contains at least one exon
pub fn contains_exon(transcript: &Transcript) -> bool {
    transcript.exon_count() > 0
}

/// Returns `Some(true)` if the transcript is coding and the CDS
/// has a proper reading-frame length, i.e. divisible by 3
///
/// Returns `None` for non-coding transcripts
pub fn correct_cds_length(transcript: &Transcript) -> Option<bool> {
    if !transcript.is_coding() {
        None
    } else {
        let mut len = 0;
        for coords in transcript.cds_coordinates() {
            len += coords.2 - coords.1 + 1
        }
        Some(len % 3 == 0)
    }
}

/// Checks if the transcript is coding and the CDS
/// starts with an `ATG` start codon
///
/// Returns `None` for non-coding transcripts
///
/// This function might return an Error if the Sequence cannot
/// be generated with the provided Fasta Reader. This could be due to
/// IO errors or invalid coordinates of the transcript
pub fn correct_start_codon(
    transcript: &Transcript,
    fasta: &mut FastaReader<impl std::io::Read + std::io::Seek>,
) -> Result<Option<bool>, FastaError> {
    let coords = transcript.start_codon();

    // In the very vast majority of cases, the start codon will be within a single
    // exon, so we optimize for this case first
    let mut seq = match coords.len() {
        0 => return Ok(None),
        1 => fasta.read_sequence(transcript.chrom(), coords[0].0.into(), coords[0].1.into())?,
        _ => {
            let mut tmp_seq = Sequence::with_capacity(3);
            for fragment in coords {
                tmp_seq.append(fasta.read_sequence(
                    transcript.chrom(),
                    fragment.0.into(),
                    fragment.1.into(),
                )?)
            }
            tmp_seq
        }
    };

    if !transcript.forward() {
        seq.reverse_complement();
    }
    Ok(Some(starts_with_start_codon(&seq)))
}

/// Checks if the transcript is coding and the CDS
/// ends with a Stop codon
///
/// Returns `None` for non-coding transcripts
///
/// This function might return an Error if the Sequence cannot
/// be generated with the provided Fasta Reader. This could be due to
/// IO errors or invalid coordinates of the transcript
pub fn correct_stop_codon(
    transcript: &Transcript,
    fasta: &mut FastaReader<impl std::io::Read + std::io::Seek>,
    code: &GeneticCode,
) -> Result<Option<bool>, FastaError> {
    let coords = transcript.stop_codon();

    // In the very vast majority of cases, the stop codon will be within a single
    // exon, so we optimize for this case first
    let mut seq = match coords.len() {
        0 => return Ok(None),
        1 => fasta.read_sequence(transcript.chrom(), coords[0].0.into(), coords[0].1.into())?,
        _ => {
            let mut tmp_seq = Sequence::with_capacity(3);
            for fragment in coords {
                tmp_seq.append(fasta.read_sequence(
                    transcript.chrom(),
                    fragment.0.into(),
                    fragment.1.into(),
                )?)
            }
            tmp_seq
        }
    };

    if !transcript.forward() {
        seq.reverse_complement();
    }
    Ok(Some(ends_with_stop_codon(&seq, code)))
}

/// Checks if the transcript is coding and the CDS
/// does not contain an extra upstream in-frame stop codon
///
/// Returns `None` for non-coding transcripts
///
/// This function might return an Error if the Sequence cannot
/// be generated with the provided Fasta Reader. This could be due to
/// IO errors or invalid coordinates of the transcript
pub fn no_upstream_stop_codon(
    transcript: &Transcript,
    fasta: &mut FastaReader<impl std::io::Read + std::io::Seek>,
    code: &GeneticCode,
) -> Result<Option<bool>, FastaError> {
    if !transcript.is_coding() {
        return Ok(None);
    }
    let seq =
        Sequence::from_coordinates(&transcript.cds_coordinates(), &transcript.strand(), fasta)?;
    if early_stop_codon(&seq, code) {
        Ok(Some(false))
    } else {
        Ok(Some(true))
    }
}

/// Checks if the transcript is coding and the 5'UTR
/// does not contain an start codon
///
/// Returns `None` for non-coding transcripts
///
/// This function might return an Error if the Sequence cannot
/// be generated with the provided Fasta Reader. This could be due to
/// IO errors or invalid coordinates of the transcript
///
/// There is no strict biological impact for this QC check and a `false` result
/// does not indicate that the transcript is incorrect.
pub fn no_upstream_start_codon(
    transcript: &Transcript,
    fasta: &mut FastaReader<impl std::io::Read + std::io::Seek>,
) -> Result<Option<bool>, FastaError> {
    if !transcript.is_coding() {
        return Ok(None);
    }
    let seq =
        Sequence::from_coordinates(&transcript.utr5_coordinates(), &transcript.strand(), fasta)?;
    if extra_start_codon(&seq) {
        Ok(Some(false))
    } else {
        Ok(Some(true))
    }
}

/// Checks if the transcript does not contain an start codon
///
/// This function might return an Error if the Sequence cannot
/// be generated with the provided Fasta Reader. This could be due to
/// IO errors or invalid coordinates of the transcript
///
/// There is no strict biological impact for this QC check and a `false` result
/// does not indicate that the transcript is incorrect.
pub fn no_start_codon(
    transcript: &Transcript,
    fasta: &mut FastaReader<impl std::io::Read + std::io::Seek>,
) -> Result<bool, FastaError> {
    let seq =
        Sequence::from_coordinates(&transcript.exon_coordinates(), &transcript.strand(), fasta)?;
    if extra_start_codon(&seq) {
        Ok(false)
    } else {
        Ok(true)
    }
}

/// checks if any (except for the last) in-frame codon of the sequence
/// is a stop codon
fn early_stop_codon(seq: &Sequence, code: &GeneticCode) -> bool {
    for codon in seq[0..seq.len() - 3].chunks(3) {
        if code.is_stop_codon(codon) {
            return true;
        }
    }
    false
}

fn starts_with_start_codon(cds: &Sequence) -> bool {
    if cds.len() < 3 {
        return false;
    }
    GeneticCode::is_start_codon(&cds[0..3])
}

fn ends_with_stop_codon(cds: &Sequence, code: &GeneticCode) -> bool {
    if cds.len() < 3 {
        return false;
    }
    code.is_stop_codon(&cds[cds.len() - 3..cds.len()])
}

/// checks if any (in or out-of frame) position
/// could be a start codon
fn extra_start_codon(seq: &Sequence) -> bool {
    if seq.len() < 3 {
        return false;
    }
    for idx in 0..seq.len() - 2 {
        let codon = &seq[idx..idx + 3];
        if GeneticCode::is_start_codon(codon) {
            return true;
        }
    }
    false
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::models::Strand;
    use crate::tests::transcripts::standard_transcript;
    use std::str::FromStr;

    #[test]
    fn test_starts_with_start_codon() {
        assert_eq!(
            starts_with_start_codon(&Sequence::from_str("ATGCGATGAT").unwrap()),
            true
        );
        assert_eq!(
            starts_with_start_codon(&Sequence::from_str("ATG").unwrap()),
            true
        );
        assert_eq!(
            starts_with_start_codon(&Sequence::from_str("AT").unwrap()),
            false
        );
        assert_eq!(
            starts_with_start_codon(&Sequence::from_str("TG").unwrap()),
            false
        );
        assert_eq!(
            starts_with_start_codon(&Sequence::from_str("").unwrap()),
            false
        );
        assert_eq!(
            starts_with_start_codon(&Sequence::from_str("CATG").unwrap()),
            false
        );
    }

    #[test]
    fn test_extra_start_codon() {
        assert_eq!(
            extra_start_codon(&Sequence::from_str("ATGCGACGA").unwrap()),
            true
        );
        assert_eq!(extra_start_codon(&Sequence::from_str("ATG").unwrap()), true);
        assert_eq!(extra_start_codon(&Sequence::from_str("AT").unwrap()), false);
        assert_eq!(extra_start_codon(&Sequence::from_str("TG").unwrap()), false);
        assert_eq!(extra_start_codon(&Sequence::from_str("").unwrap()), false);
        assert_eq!(
            extra_start_codon(&Sequence::from_str("CATG").unwrap()),
            true
        );
        assert_eq!(
            extra_start_codon(&Sequence::from_str("CAATG").unwrap()),
            true
        );
    }

    #[test]
    fn test_check_inframe_stop() {
        let code = GeneticCode::default();
        assert_eq!(
            early_stop_codon(&Sequence::from_str("TAGAAA").unwrap(), &code),
            true
        );
        assert_eq!(
            early_stop_codon(&Sequence::from_str("TAG").unwrap(), &code),
            false
        );
        assert_eq!(
            early_stop_codon(&Sequence::from_str("TAGT").unwrap(), &code),
            false
        );
        assert_eq!(
            early_stop_codon(&Sequence::from_str("TAGTA").unwrap(), &code),
            false
        );

        assert_eq!(
            early_stop_codon(&Sequence::from_str("ATGCGATAGTTA").unwrap(), &code),
            true
        );
        assert_eq!(
            early_stop_codon(&Sequence::from_str("ATGCGATAATTA").unwrap(), &code),
            true
        );
        assert_eq!(
            early_stop_codon(&Sequence::from_str("ATGCGATGATTA").unwrap(), &code),
            true
        );
        assert_eq!(
            early_stop_codon(&Sequence::from_str("TAGATGCGATTA").unwrap(), &code),
            true
        );
        assert_eq!(
            early_stop_codon(&Sequence::from_str("TAAATGCGATTA").unwrap(), &code),
            true
        );
        assert_eq!(
            early_stop_codon(&Sequence::from_str("TGAATGCGATTA").unwrap(), &code),
            true
        );
        assert_eq!(
            early_stop_codon(&Sequence::from_str("ATGCGATAG").unwrap(), &code),
            false
        );
        assert_eq!(
            early_stop_codon(&Sequence::from_str("ATGCGATAA").unwrap(), &code),
            false
        );
        assert_eq!(
            early_stop_codon(&Sequence::from_str("ATGCGATGA").unwrap(), &code),
            false
        );

        assert_eq!(
            early_stop_codon(&Sequence::from_str("ATAGATGCGATTA").unwrap(), &code),
            false
        );
        assert_eq!(
            early_stop_codon(&Sequence::from_str("ATAAATGCGATTA").unwrap(), &code),
            false
        );
        assert_eq!(
            early_stop_codon(&Sequence::from_str("ATGAATGCGATTA").unwrap(), &code),
            false
        );

        assert_eq!(
            early_stop_codon(&Sequence::from_str("ATTAGATGCGATTA").unwrap(), &code),
            false
        );
        assert_eq!(
            early_stop_codon(&Sequence::from_str("ATTAAATGCGATTA").unwrap(), &code),
            false
        );
        assert_eq!(
            early_stop_codon(&Sequence::from_str("ATTGAATGCGATTA").unwrap(), &code),
            false
        );
    }

    #[test]
    fn test_incorrect_cds_length() {
        let tx = standard_transcript();

        // Actually, the standard transcript has a wrong CDS length
        assert_eq!(correct_cds_length(&tx), Some(false));
    }

    #[test]
    fn test_correct_cds_length() {
        let mut tx = standard_transcript();
        *tx.exons_mut()[1].cds_start_mut() = Some(23);
        assert_eq!(correct_cds_length(&tx), Some(true));
    }

    #[test]
    fn test_non_coding_cds_length() {
        let mut tx = standard_transcript();
        for ex in tx.exons_mut() {
            *ex.cds_start_mut() = None;
        }
        assert_eq!(correct_cds_length(&tx), None);
    }

    #[test]
    fn test_correct_forward_split_start_codon() {
        let mut fasta = FastaReader::from_file("tests/data/small.fasta").unwrap();
        let tx = standard_transcript();
        assert_eq!(correct_start_codon(&tx, &mut fasta).unwrap(), Some(true));
    }

    #[test]
    fn test_incorrect_forward_split_start_codon() {
        let mut fasta = FastaReader::from_file("tests/data/small.fasta").unwrap();
        let mut tx = standard_transcript();
        *tx.exons_mut()[1].cds_start_mut() = Some(23);
        assert_eq!(correct_start_codon(&tx, &mut fasta).unwrap(), Some(false));

        let mut tx = standard_transcript();
        *tx.exons_mut()[1].cds_start_mut() = Some(24);
        *tx.exons_mut()[2].cds_start_mut() = Some(32);
        assert_eq!(correct_start_codon(&tx, &mut fasta).unwrap(), Some(false));
    }

    #[test]
    fn test_correct_forward_single_start_codon() {
        let mut fasta = FastaReader::from_file("tests/data/small.fasta").unwrap();
        let mut tx = standard_transcript();
        *tx.exons_mut()[1].start_mut() = 20;
        *tx.exons_mut()[1].cds_start_mut() = Some(20);
        assert_eq!(correct_start_codon(&tx, &mut fasta).unwrap(), Some(true));
    }

    #[test]
    fn test_incorrect_forward_single_start_codon() {
        let mut fasta = FastaReader::from_file("tests/data/small.fasta").unwrap();
        let mut tx = standard_transcript();
        *tx.exons_mut()[1].cds_start_mut() = Some(21);
        assert_eq!(correct_start_codon(&tx, &mut fasta).unwrap(), Some(false));
    }

    #[test]
    fn test_correct_reverse_single_start_codon() {
        let mut fasta = FastaReader::from_file("tests/data/small.fasta").unwrap();
        let mut tx = standard_transcript();

        // just to be safe - remove the actual start codon
        *tx.exons_mut()[1].cds_start_mut() = Some(21);

        *tx.strand_mut() = Strand::Minus;
        *tx.exons_mut()[4].end_mut() = 150;
        *tx.exons_mut()[4].cds_end_mut() = Some(132);
        *tx.exons_mut()[4].cds_start_mut() = Some(51);

        assert_eq!(correct_start_codon(&tx, &mut fasta).unwrap(), Some(true));
    }

    #[test]
    fn test_correct_reverse_split_start_codon() {
        let mut fasta = FastaReader::from_file("tests/data/small.fasta").unwrap();
        let mut tx = standard_transcript();

        // just to be safe - remove the actual start codon
        *tx.exons_mut()[1].cds_start_mut() = Some(21);

        *tx.strand_mut() = Strand::Minus;
        *tx.exons_mut()[4].start_mut() = 42;
        *tx.exons_mut()[3].cds_start_mut() = Some(42);
        *tx.exons_mut()[3].cds_end_mut() = Some(42);

        assert_eq!(correct_start_codon(&tx, &mut fasta).unwrap(), Some(true));
    }

    #[test]
    fn test_correct_forward_single_stop_codon() {
        let code = GeneticCode::default();
        let mut fasta = FastaReader::from_file("tests/data/small.fasta").unwrap();
        let mut tx = standard_transcript();

        // TGA
        assert_eq!(
            correct_stop_codon(&tx, &mut fasta, &code).unwrap(),
            Some(true)
        );

        // TAG
        *tx.exons_mut()[2].cds_end_mut() = Some(39);
        *tx.exons_mut()[3].cds_start_mut() = None;
        *tx.exons_mut()[3].cds_end_mut() = None;
        assert_eq!(
            correct_stop_codon(&tx, &mut fasta, &code).unwrap(),
            Some(true)
        );
    }

    #[test]
    fn test_correct_forward_split_stop_codon() {
        let code = GeneticCode::default();
        let mut fasta = FastaReader::from_file("tests/data/small.fasta").unwrap();
        let mut tx = standard_transcript();

        // TAA
        *tx.exons_mut()[2].cds_end_mut() = Some(38);
        *tx.exons_mut()[3].cds_start_mut() = Some(44);
        *tx.exons_mut()[3].cds_end_mut() = Some(44);
        assert_eq!(
            correct_stop_codon(&tx, &mut fasta, &code).unwrap(),
            Some(true)
        );
    }

    #[test]
    fn test_correct_reverse_single_stop_codon() {
        let code = GeneticCode::default();
        let mut fasta = FastaReader::from_file("tests/data/small.fasta").unwrap();
        let mut tx = standard_transcript();

        // TAG ==> CT-A
        *tx.strand_mut() = Strand::Minus;
        *tx.exons_mut()[1].cds_start_mut() = None;
        *tx.exons_mut()[1].cds_end_mut() = None;
        *tx.exons_mut()[2].cds_start_mut() = Some(29);
        *tx.exons_mut()[2].cds_end_mut() = Some(30);
        *tx.exons_mut()[3].cds_start_mut() = Some(38);

        assert_eq!(
            correct_stop_codon(&tx, &mut fasta, &code).unwrap(),
            Some(true)
        );
    }

    #[test]
    fn test_qc_check() {
        let code = GeneticCode::default();
        let mut fasta = FastaReader::from_file("tests/data/small.fasta").unwrap();
        let qc = QcCheck::new(&standard_transcript(), &mut fasta, &code);
        let mut ideal = QcCheck::ideal();
        // the standard transcript does not have a correct CDS-length...
        ideal.cds_len = QcResult::NOK;
        assert_eq!(qc, ideal);
    }
}

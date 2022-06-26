use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::fasta::FastaReader;
use crate::models::{Transcript, TranscriptWrite};
use crate::qc::QcCheck;
use crate::utils::errors::{AtgError, ReadWriteError};

/// Writes [`QcCheck`]s into a `BufWriter`
///
/// # Examples
///
/// ```rust
/// use std::io;
/// use atglib::tests;;
/// use atglib::qc::Writer;
/// use atglib::fasta::FastaReader;
/// use atglib::models::TranscriptWrite;
///
/// let transcripts = vec![tests::transcripts::standard_transcript()];
///
/// let output = Vec::new(); // substitute this with proper IO (io::stdout())
/// let mut writer = Writer::new(output);
/// writer.fasta_reader(FastaReader::from_file("tests/data/small.fasta").unwrap());
/// writer.write_transcript_vec(&transcripts);
///
/// let written_output = String::from_utf8(writer.into_inner().unwrap()).unwrap();
/// assert_eq!(written_output, "Test-Gene\tTest-Transcript\tOK\tNOK\tOK\tOK\tOK\tOK\tOK\n");
/// ```
pub struct Writer<W: std::io::Write> {
    inner: BufWriter<W>,
    fasta_reader: Option<FastaReader<File>>,
}

impl Writer<File> {
    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self, ReadWriteError> {
        match File::create(path.as_ref()) {
            Ok(file) => Ok(Self::new(file)),
            Err(err) => Err(ReadWriteError::new(err)),
        }
    }
}

impl<W: std::io::Write> Writer<W> {
    /// Creates a new generic Writer for any `std::io::Read`` object
    ///
    /// Use this method when you want to write to stdout or
    /// a remote source, e.g. via HTTP
    pub fn new(writer: W) -> Self {
        Writer::from_buf_writer(BufWriter::new(writer))
    }

    /// Constructs a new, empty Writer with the specified capacity.
    pub fn with_capacity(capacity: usize, writer: W) -> Self {
        Writer::from_buf_writer(BufWriter::with_capacity(capacity, writer))
    }

    /// Private constructer method to set default values
    fn from_buf_writer(writer: BufWriter<W>) -> Self {
        Writer {
            inner: writer,
            fasta_reader: None,
        }
    }

    /// Specify a [`FastaReader'](`crate::fasta::FastaReader`) to retrieve
    /// the reference genome sequence.
    ///
    /// You must set a `fasta_reader`, since the `Writer` does not have any
    /// information about the reference genome to use.
    ///
    /// ```rust
    /// use atglib::qc::Writer;
    /// use atglib::fasta::FastaReader;
    ///
    /// let output = Vec::new(); // substitute this with proper IO (io::stdout())
    /// let mut writer = Writer::new(output);
    /// // specify the reference genome fasta file
    /// writer.fasta_reader(FastaReader::from_file("tests/data/small.fasta").unwrap());
    /// ```
    pub fn fasta_reader(&mut self, r: FastaReader<File>) {
        self.fasta_reader = Some(r)
    }

    pub fn flush(&mut self) -> Result<(), AtgError> {
        match self.inner.flush() {
            Ok(res) => Ok(res),
            Err(err) => Err(AtgError::from(err.to_string())),
        }
    }

    pub fn into_inner(self) -> Result<W, AtgError> {
        match self.inner.into_inner() {
            Ok(res) => Ok(res),
            Err(err) => Err(AtgError::from(err.to_string())),
        }
    }

    /// Writes the header row for the tab-separated QC results
    pub fn write_header(&mut self) -> Result<(), std::io::Error> {
        let columns = vec![
            "Gene",
            "transcript",
            "Exon",
            "CDS Length",
            "Correct Start Codon",
            "Correct Stop Codon",
            "No upstream Start Codon",
            "No upstream Stop Codon",
            "Correct Coordinates",
        ];
        self.inner.write_all(columns.join("\t").as_bytes())?;
        self.inner.write_all("\n".as_bytes())
    }
}

impl<W: std::io::Write> TranscriptWrite for Writer<W> {
    /// Writes a single transcript formatted as RefGene with an extra newline
    ///
    /// This method adds an extra newline at the end of the row
    /// to allow writing multiple transcripts continuosly
    fn writeln_single_transcript(&mut self, transcript: &Transcript) -> Result<(), std::io::Error> {
        self.write_single_transcript(transcript)?;
        self.inner.write_all("\n".as_bytes())
    }

    fn write_single_transcript(&mut self, transcript: &Transcript) -> Result<(), std::io::Error> {
        let fasta_reader = match &mut self.fasta_reader {
            None => {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    AtgError::new("No Fasta Reader specified"),
                ))
            }
            Some(r) => r,
        };
        let qc = QcCheck::new(transcript, fasta_reader);
        let tx_name = format!("{}\t{}\t", transcript.gene(), transcript.name());
        self.inner.write_all(tx_name.as_bytes())?;
        self.inner.write_all(qc.to_string().as_bytes())
    }
}

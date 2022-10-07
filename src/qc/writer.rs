use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use log::debug;

use crate::fasta::FastaReader;
use crate::models::{GeneticCode, Transcript, TranscriptWrite};
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
pub struct Writer<W: std::io::Write, R: std::io::Read + std::io::Seek> {
    inner: BufWriter<W>,
    fasta_reader: Option<FastaReader<R>>,
    genetic_code: GeneticCode,
    alternative_genetic_codes: Option<Vec<(String, GeneticCode)>>,
}

impl<R: std::io::Read + std::io::Seek> Writer<File, R> {
    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self, ReadWriteError> {
        match File::create(path.as_ref()) {
            Ok(file) => Ok(Self::new(file)),
            Err(err) => Err(ReadWriteError::new(err)),
        }
    }
}

impl<W: std::io::Write, R:std::io::Read + std::io::Seek> Writer<W, R> {
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
            genetic_code: GeneticCode::default(),
            alternative_genetic_codes: None,
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
    pub fn fasta_reader(&mut self, r: FastaReader<R>) {
        self.fasta_reader = Some(r)
    }

    /// Changes the default genetic code to a custom one
    ///
    /// By default, the standard genetic code is used for translating to AminoAcids
    ///
    pub fn default_genetic_code(&mut self, code: GeneticCode) {
        self.genetic_code = code
    }

    /// Add a custom genetic code that will be used for all transcripts on `chrom`
    ///
    /// For example, when using `QC` for human genomes, you should specify to use
    /// the `vertebrate_mitochondria` genetic code for all transcripts on the mitochondrial
    /// chromsome.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use atglib::qc::Writer;
    /// use atglib::fasta::FastaReader;
    /// use atglib::models::GeneticCode;
    ///
    /// let output = Vec::new(); // substitute this with proper IO (io::stdout())
    /// let mut writer = Writer::new(output);
    /// // specify the reference genome fasta file
    /// writer.fasta_reader(FastaReader::from_file("tests/data/small.fasta").unwrap());
    /// // use vertebrate_mitochondria genetic code for chrM transcripts
    /// writer.add_genetic_code("chrM".to_string(), GeneticCode::vertebrate_mitochondrial());
    ///
    /// let codes = writer.genetic_codes();
    /// assert_eq!(codes[0].0, "*");
    /// assert_eq!(codes[1].0, "chrM");
    ///
    /// // The mitochondrial genetic code contains 4 Stop codons, the standard genetic code has 3:
    /// assert_eq!(codes[0].1.stop_codons().len(), 3);
    /// assert_eq!(codes[1].1.stop_codons().len(), 4);
    /// ```
    ///
    /// # Note
    /// The matching of chromosomes to genetic codes is implemented for use-cases with a maximum of 1 or 2
    /// additional genetic codes. If you have use-cases where you need many different genetic codes
    /// it would be better to split your `Transcripts` beforehand.
    ///
    /// It is optimized for vertebrate genomes (i.e. standard code + vertebrate mitochondrial code)
    pub fn add_genetic_code(&mut self, chrom: String, code: GeneticCode) {
        if self.alternative_genetic_codes.is_none() {
            self.alternative_genetic_codes = Some(Vec::new());
        }
        if let Some(alt_codes) = self.alternative_genetic_codes.as_mut() {
            alt_codes.push((chrom, code));
        }
    }

    /// List all genetic codes that are linked to the QCWriter
    ///
    /// There is no real use for this method, other than for testing
    /// during development
    pub fn genetic_codes(&self) -> Vec<(&str, &GeneticCode)> {
        let mut res = Vec::new();
        res.push(("*", &self.genetic_code));
        if let Some(alt_codes) = &self.alternative_genetic_codes {
            for alt in alt_codes {
                res.push((&alt.0, &alt.1));
            }
        }
        res
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

impl<W: std::io::Write, R: std::io::Read + std::io::Seek> TranscriptWrite for Writer<W, R> {
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

        let mut code = &self.genetic_code;
        if let Some(alts) = &self.alternative_genetic_codes {
            for (chrom, chrom_code) in alts {
                if chrom == transcript.chrom() {
                    debug!(
                        "Using custom genetic code {} for {} on {}",
                        chrom_code,
                        transcript.name(),
                        chrom
                    );
                    code = chrom_code;
                    break;
                }
            }
        }

        let qc = QcCheck::new(transcript, fasta_reader, code);
        let tx_name = format!("{}\t{}\t", transcript.gene(), transcript.name());
        self.inner.write_all(tx_name.as_bytes())?;
        self.inner.write_all(qc.to_string().as_bytes())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::Nucleotide;
    use crate::tests::transcripts::{mito_transcript, standard_transcript};
    #[test]
    fn test_adding_genetic_codes() {
        let output = Vec::new(); // substitute this with proper IO (io::stdout())
        let mut writer = Writer::new(output);
        writer.fasta_reader(FastaReader::from_file("tests/data/small.fasta").unwrap());
        // use vertebrate_mitochondria genetic code for chrM transcripts
        writer.add_genetic_code("chrM".to_string(), GeneticCode::vertebrate_mitochondrial());
        let codes = writer.genetic_codes();
        assert_eq!(codes[0].0, "*");
        assert_eq!(codes[1].0, "chrM");

        assert_eq!(codes[0].1.stop_codons().len(), 3);
        assert_eq!(codes[1].1.stop_codons().len(), 4);

        assert!(codes[0]
            .1
            .is_stop_codon(&[Nucleotide::T, Nucleotide::A, Nucleotide::A]));
        assert!(codes[0]
            .1
            .is_stop_codon(&[Nucleotide::T, Nucleotide::A, Nucleotide::G]));
        assert!(codes[0]
            .1
            .is_stop_codon(&[Nucleotide::T, Nucleotide::G, Nucleotide::A]));

        assert!(codes[1]
            .1
            .is_stop_codon(&[Nucleotide::T, Nucleotide::A, Nucleotide::A]));
        assert!(codes[1]
            .1
            .is_stop_codon(&[Nucleotide::T, Nucleotide::A, Nucleotide::G]));
        assert!(codes[1]
            .1
            .is_stop_codon(&[Nucleotide::A, Nucleotide::G, Nucleotide::A]));
        assert!(codes[1]
            .1
            .is_stop_codon(&[Nucleotide::A, Nucleotide::G, Nucleotide::G]));
    }

    #[test]
    fn test_chromosome_with_alternative_stop_codons() {
        let tx = standard_transcript();
        let mito_tx = mito_transcript();

        // Use only standard genetic code => Fail for mito transcript
        let output = Vec::new();
        let mut writer = Writer::new(output);
        writer.fasta_reader(FastaReader::from_file("tests/data/small.fasta").unwrap());
        writer.writeln_single_transcript(&tx).unwrap();
        writer.writeln_single_transcript(&mito_tx).unwrap();
        let written_output = String::from_utf8(writer.into_inner().unwrap()).unwrap();
        // The Stop Codon is mitochondria specific and the transcript contains an upstream standard stop codon
        //                                                                                                                                     |        |
        //                                                                                                                                     V        V
        assert_eq!(written_output, "Test-Gene\tTest-Transcript\tOK\tNOK\tOK\tOK\tOK\tOK\tOK\nTest-Mito-Gene\tTest-Mito-Transcript\tOK\tOK\tOK\tNOK\tOK\tNOK\tOK\n");

        // use vertebrate_mitochondria genetic code for chrM transcripts
        let output = Vec::new();
        let mut writer = Writer::new(output);
        writer.fasta_reader(FastaReader::from_file("tests/data/small.fasta").unwrap());
        writer.add_genetic_code("chrM".to_string(), GeneticCode::vertebrate_mitochondrial());
        writer.writeln_single_transcript(&tx).unwrap();
        writer.writeln_single_transcript(&mito_tx).unwrap();
        let written_output = String::from_utf8(writer.into_inner().unwrap()).unwrap();
        // That's how it should be if all is defined correctly
        //                                                                                                                                     |       |
        //                                                                                                                                     V       V
        assert_eq!(written_output, "Test-Gene\tTest-Transcript\tOK\tNOK\tOK\tOK\tOK\tOK\tOK\nTest-Mito-Gene\tTest-Mito-Transcript\tOK\tOK\tOK\tOK\tOK\tOK\tOK\n");

        // use vertebrate_mitochondria genetic code for all transcripts => fail for standard transcript
        let output = Vec::new();
        let mut writer = Writer::new(output);
        writer.fasta_reader(FastaReader::from_file("tests/data/small.fasta").unwrap());
        writer.default_genetic_code(GeneticCode::vertebrate_mitochondrial());
        writer.writeln_single_transcript(&tx).unwrap();
        writer.writeln_single_transcript(&mito_tx).unwrap();
        let written_output = String::from_utf8(writer.into_inner().unwrap()).unwrap();
        // The standard transcript uses a Stop codon which is not a stop codon in Mitochondria
        //                                                                   |                                                                   |       |
        //                                                                   V                                                                   V       V
        assert_eq!(written_output, "Test-Gene\tTest-Transcript\tOK\tNOK\tOK\tNOK\tOK\tOK\tOK\nTest-Mito-Gene\tTest-Mito-Transcript\tOK\tOK\tOK\tOK\tOK\tOK\tOK\n");
    }
}

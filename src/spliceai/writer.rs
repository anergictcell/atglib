use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::models::{Transcript, TranscriptWrite, Transcripts};
use crate::utils::errors::AtgError;
use crate::utils::errors::ReadWriteError;
use crate::utils::merge;

const HEADER: &str = "#NAME\tCHROM\tSTRAND\tTX_START\tTX_END\tEXON_START\tEXON_END";

/// Writes [`Transcript`]s into a `BufWriter` to be used as gene annotation
/// by SpliceAI
///
/// # Examples
///
/// ```rust
/// use atglib::tests;;
/// use atglib::spliceai::Writer;
/// use atglib::models::TranscriptWrite;
///
/// let transcripts = vec![tests::transcripts::standard_transcript()];
///
/// let output = Vec::new(); // substitute this with proper IO (io::stdout())
/// let mut writer = Writer::new(output);
/// writer.write_transcript_vec(&transcripts);
///
/// assert_eq!(
/// writer.into_inner().unwrap(), // this is our actual output
/// b"#NAME\tCHROM\tSTRAND\tTX_START\tTX_END\tEXON_START\tEXON_END
/// Test-Gene\tchr1\t+\t11\t55\t11,21,31,41,51,\t15,25,35,45,55,
/// "
/// );
/// ```
pub struct Writer<W: std::io::Write> {
    inner: BufWriter<W>,
    header_written: bool,
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
        Writer {
            inner: BufWriter::new(writer),
            header_written: false,
        }
    }

    pub fn with_capacity(capacity: usize, writer: W) -> Self {
        Writer {
            inner: BufWriter::with_capacity(capacity, writer),
            header_written: false,
        }
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

    pub fn write_header(&mut self) -> Result<(), std::io::Error> {
        if self.header_written {
            return Ok(());
        }
        self.inner.write_all(HEADER.as_bytes())?;
        self.inner.write_all("\n".as_bytes())?;
        self.header_written = true;
        Ok(())
    }
}

impl<W: std::io::Write> TranscriptWrite for Writer<W> {
    /// Writes a single transcript formatted for SpliceAI with an extra newline
    fn writeln_single_transcript(&mut self, transcript: &Transcript) -> Result<(), std::io::Error> {
        self.write_header()?;
        self.write_single_transcript(transcript)?;
        self.inner.write_all("\n".as_bytes())
    }

    fn write_single_transcript(&mut self, transcript: &Transcript) -> Result<(), std::io::Error> {
        self.inner
            .write_all(SpliceAiLine::from(transcript).as_bytes())
    }

    fn write_transcripts(&mut self, transcripts: &Transcripts) -> Result<(), std::io::Error> {
        self.write_header()?;

        for gene in transcripts.genes() {
            let txs = transcripts.by_gene(gene);
            let line = SpliceAiLine::from(txs);
            self.inner.write_all(line.as_bytes())?;
            self.inner.write_all("\n".as_bytes())?;
        }
        Ok(())
    }
}

struct SpliceAiLine {
    line: String,
}

impl SpliceAiLine {
    fn as_bytes(&self) -> &[u8] {
        self.line.as_bytes()
    }
}

impl std::fmt::Display for SpliceAiLine {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.line)
    }
}

impl From<&Transcript> for SpliceAiLine {
    fn from(tx: &Transcript) -> SpliceAiLine {
        SpliceAiLine {
            line: format!(
                "{}\t{}\t{}\t{}\t{}\t{},\t{},",
                tx.gene(),
                tx.chrom(),
                tx.strand(),
                tx.tx_start(),
                tx.tx_end(),
                tx.exons()
                    .iter()
                    .map(|ex| ex.start().to_string())
                    .collect::<Vec<String>>()
                    .join(","),
                tx.exons()
                    .iter()
                    .map(|ex| ex.end().to_string())
                    .collect::<Vec<String>>()
                    .join(",")
            ),
        }
    }
}

impl From<Vec<&Transcript>> for SpliceAiLine {
    fn from(gene: Vec<&Transcript>) -> SpliceAiLine {
        if gene.is_empty() {
            return SpliceAiLine {
                line: "".to_string(),
            };
        }
        let mut min_start = gene[0].tx_start();
        let mut max_end = gene[0].tx_end();
        let mut exons: Vec<(u32, u32)> = Vec::new();
        for tx in &gene {
            if tx.tx_start() < min_start {
                min_start = tx.tx_start()
            }
            if tx.tx_end() > max_end {
                max_end = tx.tx_end()
            }

            for ex in tx.exons() {
                exons.push((ex.start(), ex.end()));
            }
        }
        exons.sort_by(|a, b| a.0.cmp(&b.0));
        let new_exons = merge(&exons);

        SpliceAiLine {
            line: format!(
                "{}\t{}\t{}\t{}\t{}\t{},\t{},",
                gene[0].gene(),
                gene[0].chrom(),
                gene[0].strand(),
                min_start,
                max_end,
                new_exons
                    .iter()
                    .map(|ex| ex.0.to_string())
                    .collect::<Vec<String>>()
                    .join(","),
                new_exons
                    .iter()
                    .map(|ex| ex.1.to_string())
                    .collect::<Vec<String>>()
                    .join(",")
            ),
        }
    }
}

#[cfg(test)]
mod test_spliceailine {
    use super::*;
    use crate::models;
    use crate::tests::transcripts::standard_transcript;

    #[test]
    fn test_single_transcript() {
        let tx = standard_transcript();
        let output = SpliceAiLine::from(&tx);
        assert_eq!(tx.exon_count(), 5);
        assert_eq!(
            output.to_string(),
            "Test-Gene\tchr1\t+\t11\t55\t11,21,31,41,51,\t15,25,35,45,55,".to_string()
        );
    }

    #[test]
    fn test_multiple_transcripts() {
        let mut tx1 = standard_transcript();
        let mut exons = tx1.exons_mut();
        exons.remove(2);
        let mut tx2 = standard_transcript();
        exons = tx2.exons_mut();
        exons.remove(3);
        let output = SpliceAiLine::from(vec![&tx1, &tx2]);
        assert_eq!(tx1.exon_count(), 4);
        assert_eq!(tx2.exon_count(), 4);
        assert_eq!(
            output.to_string(),
            "Test-Gene\tchr1\t+\t11\t55\t11,21,31,41,51,\t15,25,35,45,55,".to_string()
        );
    }

    #[test]
    fn write_multiple_transcripts() {
        let transcripts = vec![standard_transcript(), standard_transcript()];
        let output = Vec::new();
        let mut writer = Writer::new(output);
        writer
            .write_transcript_vec(&transcripts)
            .expect("Error writing into bytevec");
        let written_output = String::from_utf8(writer.into_inner().unwrap()).unwrap();
        let expected_output = [
            "Test-Gene",
            "chr1",
            "+",
            "11",
            "55",
            "11,21,31,41,51,",
            "15,25,35,45,55,",
        ]
        .join("\t")
        .to_string();
        assert_eq!(
            written_output,
            format!(
                "#NAME\tCHROM\tSTRAND\tTX_START\tTX_END\tEXON_START\tEXON_END\n{}\n{}\n",
                expected_output, expected_output
            )
        );
    }

    #[test]
    fn write_transcripts() {
        let mut transcripts = models::Transcripts::with_capacity(3);
        transcripts.push(standard_transcript());

        let mut tx2 = models::TranscriptBuilder::new()
            .name("Test-Transcript_2")
            .chrom("chr1")
            .strand(models::Strand::Plus)
            .gene("Test-Gene")
            .cds_start_stat(models::CdsStat::None)
            .cds_end_stat(models::CdsStat::None)
            .build()
            .unwrap();
        tx2.append_exons(standard_transcript().exons_mut());
        transcripts.push(tx2);

        let mut tx3 = models::TranscriptBuilder::new()
            .name("Test-Transcript_3")
            .chrom("chr1")
            .strand(models::Strand::Plus)
            .gene("Test-Gene_2")
            .cds_start_stat(models::CdsStat::None)
            .cds_end_stat(models::CdsStat::None)
            .build()
            .unwrap();
        tx3.append_exons(standard_transcript().exons_mut());
        transcripts.push(tx3);

        let output = Vec::new();
        let mut writer = Writer::new(output);
        writer
            .write_transcripts(&transcripts)
            .expect("Error writing into bytevec");
        let written_output = String::from_utf8(writer.into_inner().unwrap()).unwrap();
        let expected_output1 = [
            "Test-Gene",
            "chr1",
            "+",
            "11",
            "55",
            "11,21,31,41,51,",
            "15,25,35,45,55,",
        ]
        .join("\t")
        .to_string();
        let expected_output2 = [
            "Test-Gene_2",
            "chr1",
            "+",
            "11",
            "55",
            "11,21,31,41,51,",
            "15,25,35,45,55,",
        ]
        .join("\t")
        .to_string();

        // The `Transcripts` hashmap does not guarantee identical order of iterating
        // the genes. So we must test both possible orders of Test-Gene and Test-Gene2

        let test1 = written_output
            == format!(
                "#NAME\tCHROM\tSTRAND\tTX_START\tTX_END\tEXON_START\tEXON_END\n{}\n{}\n",
                expected_output1, expected_output2
            );

        let test2 = written_output
            == format!(
                "#NAME\tCHROM\tSTRAND\tTX_START\tTX_END\tEXON_START\tEXON_END\n{}\n{}\n",
                expected_output2, expected_output1
            );

        assert!(test1 ^ test2);
    }
}

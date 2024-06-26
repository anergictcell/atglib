use std::collections::HashMap;
use std::convert::TryFrom;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::str::FromStr;

use crate::gtf::{GtfFeature, GtfRecord, GtfRecordsGroup};
use crate::models::{Transcript, TranscriptRead, Transcripts};
use crate::utils::errors::ParseGtfError;
use crate::utils::errors::ReadWriteError;

use super::record::UncheckedGtfRecord;

/// Parses GTF data and creates [`Transcript`]s.
///
/// GTF data can be read from a file, stdin or remote sources
/// All sources are supported that provide a `std::io::Read` implementation.
///
/// # Examples
///
/// ```rust
/// use atglib::gtf::Reader;
/// use atglib::models::TranscriptRead;
///
/// // create a reader from the tests GTF file
/// let reader = Reader::from_file("tests/data/example.gtf");
/// assert_eq!(reader.is_ok(), true);
///
/// // parse the GTF file
/// let transcripts = reader
///     .unwrap()
///     .transcripts()
///     .unwrap();
///
/// assert_eq!(transcripts.len(), 27);
/// ```
pub struct Reader<R> {
    inner: std::io::BufReader<R>,
    line_content: String,
}

impl Reader<File> {
    /// Creates a Reader instance that reads from a File
    ///
    /// Use this method when you want to read from a GTF file
    /// on your local file system
    ///
    /// # Examples
    ///
    /// ```rust
    /// use atglib::gtf::Reader;
    /// use atglib::models::TranscriptRead;
    ///
    /// // create a reader from the tests GTF file
    /// let reader = Reader::from_file("tests/data/example.gtf");
    /// assert_eq!(reader.is_ok(), true);
    ///
    /// // parse the GTF file
    /// let transcripts = reader
    ///     .unwrap()
    ///     .transcripts()
    ///     .unwrap();
    ///
    /// assert_eq!(transcripts.len(), 27);
    /// ```
    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self, ReadWriteError> {
        match File::open(path.as_ref()) {
            Ok(file) => Ok(Self::new(file)),
            Err(err) => Err(ReadWriteError::new(err)),
        }
    }
}

impl<R: std::io::Read> Reader<R> {
    /// creates a new Reader instance from any `std::io::Read` object
    ///
    /// Use this method when you want to read from stdin or from
    /// a remote source, e.g. via HTTP
    pub fn new(reader: R) -> Self {
        Reader {
            inner: BufReader::new(reader),
            line_content: String::with_capacity(200),
        }
    }

    /// Creates a new Reader instance with a known capcity
    ///
    /// Use this when you know the size of your GTF source
    pub fn with_capacity(capacity: usize, reader: R) -> Self {
        Reader {
            inner: BufReader::with_capacity(capacity, reader),
            line_content: String::with_capacity(200),
        }
    }

    /// Returns one line of a GTF file as `GtfRecord`
    ///
    /// This method should rarely be used. GTF files can contain unordered
    /// records and handling lines individually is rarely desired.
    fn line(&mut self) -> Option<Result<GtfRecord, ParseGtfError>> {
        self.line_content.clear();
        match self.inner.read_line(&mut self.line_content) {
            Ok(_) => {}
            Err(x) => {
                return Some(Err(ParseGtfError {
                    message: x.to_string(),
                }))
            }
        }

        if self.line_content.starts_with('#') {
            return self.line();
        }

        if self.line_content.is_empty() {
            return None;
        }

        UncheckedGtfRecord::from_str(&self.line_content)
            .map(|record| {
                if record.standard_feature() {
                    Some(GtfRecord::try_from(record))
                } else {
                    self.line()
                }
            })
            .unwrap_or_else(|err| Some(Err(err)))
    }
}

impl<R: std::io::Read> TranscriptRead for Reader<R> {
    /// Reads in GTF data and returns the final list of `Transcripts`
    ///
    /// # Examples
    ///
    /// ```rust
    /// use atglib::gtf::Reader;
    /// use atglib::models::TranscriptRead;
    ///
    /// // create a reader from the tests GTF file
    /// let reader = Reader::from_file("tests/data/example.gtf");
    ///
    /// // parse the GTF file
    /// let transcripts = reader
    ///     .unwrap()
    ///     .transcripts()
    ///     .unwrap();
    ///
    /// assert_eq!(transcripts.len(), 27);
    /// ```
    fn transcripts(&mut self) -> Result<Transcripts, ReadWriteError> {
        let mut transcript_hashmap: HashMap<String, GtfRecordsGroup> = HashMap::new();
        while let Some(line) = self.line() {
            let gtf_record = line?;
            let transcript = transcript_hashmap
                .entry(gtf_record.transcript().to_string())
                .or_insert_with(|| GtfRecordsGroup::new(gtf_record.transcript()));

            match gtf_record.feature() {
                GtfFeature::Exon => transcript.add_exon(gtf_record),
                GtfFeature::CDS => transcript.add_exon(gtf_record),
                GtfFeature::UTR | GtfFeature::UTR3 | GtfFeature::UTR5 => {
                    transcript.add_exon(gtf_record)
                }
                GtfFeature::StartCodon => {
                    transcript.add_exon(gtf_record);
                }
                GtfFeature::StopCodon => transcript.add_exon(gtf_record),
                _ => {}
            }
        }

        let mut res = Transcripts::with_capacity(transcript_hashmap.len());

        for (_, gtf_transcript) in transcript_hashmap.drain() {
            match Transcript::try_from(gtf_transcript) {
                Ok(transcript) => res.push(transcript),
                Err(err) => {
                    return Err(ReadWriteError::new(format!("Error parsing {err}")));
                }
            }
        }
        Ok(res)
    }
}

#[cfg(test)]
mod test_nm_001385228 {
    use super::*;

    #[test]
    fn test_read() {
        let mut reader = Reader::from_file("tests/data/NM_001385228.1_2.gtf").unwrap();
        let tr = match reader.transcripts() {
            Ok(res) => res,
            _ => panic!("No transcripts could be read"),
        };

        assert_eq!(tr.len(), 1);
        let t = tr.by_name("NM_001385228.1_2")[0];
        assert_eq!(t.exons().len(), 9);
        assert_eq!(t.cds_start().unwrap(), 206105119);
        assert_eq!(t.cds_end().unwrap(), 206135359);

        let start = t.start_codon();
        assert_eq!(start.len(), 1);

        let stop = t.stop_codon();
        assert_eq!(stop.len(), 2);

        assert!(!t.exons()[0].is_coding());
        assert!(!t.exons()[1].is_coding());
        assert!(!t.exons()[2].is_coding());
        assert!(!t.exons()[3].is_coding());
        assert!(!t.exons()[4].is_coding());
        assert!(t.exons()[5].is_coding());
        assert!(t.exons()[6].is_coding());
        assert!(t.exons()[7].is_coding());
        assert!(!t.exons()[8].is_coding());

        assert_eq!(t.exons()[5].cds_start().unwrap(), 206105119);
        assert_eq!(t.exons()[6].cds_start().unwrap(), 206105123);
        assert_eq!(t.exons()[7].cds_start().unwrap(), 206135293);

        assert_eq!(t.exons()[5].cds_end().unwrap(), 206105120);
        assert_eq!(t.exons()[6].cds_end().unwrap(), 206105206);
        assert_eq!(t.exons()[7].cds_end().unwrap(), 206135359);
    }
}

#[cfg(test)]
mod test_nm_201550 {
    use super::*;

    #[test]
    fn test_read() {
        let mut reader = Reader::from_file("tests/data/NM_201550.4.gtf")
            .expect("Something failed reading the GTF file NM_201550.4.gtf");
        let tr = match reader.transcripts() {
            Ok(res) => res,
            _ => panic!("No transcripts could be read"),
        };
        assert_eq!(tr.len(), 1);
        let t = tr.by_name("NM_201550.4")[0];
        assert_eq!(t.exons().len(), 1);
        assert_eq!(t.cds_start().unwrap(), 70003785);
        assert_eq!(t.cds_end().unwrap(), 70004618);

        let start = t.start_codon();
        assert_eq!(start.len(), 1);

        let stop = t.stop_codon();
        assert_eq!(stop.len(), 1);

        assert!(t.exons()[0].is_coding());
    }
}

/*
This one fails due to book-ended exons in the input file
#[cfg(test)]
mod test_nm_001371720 {
    use super::*;

    #[test]
    fn test_read() {
        let mut reader = Reader::from_file("tests/data/NM_001371720.1.gtf")
            .expect("Something failed reading the GTF file NM_001371720.1.gtf");
        let tr = match reader.transcripts() {
            Ok(res) => res,
            _ => panic!("No transcripts could be read")
        };
        assert_eq!(tr.len(), 1);
        let t = tr.by_name("NM_001371720.1")[0];
        assert_eq!(t.exons().len(), 8);
        assert_eq!(t.cds_start().unwrap(), 155158611);
        assert_eq!(t.cds_end().unwrap(), 155162634);

        let start = t.start_codon();
        assert_eq!(start.len(), 1);

        let stop = t.stop_codon();
        assert_eq!(stop.len(), 1);

        assert_eq!(t.exons()[0].is_coding(), true);
        assert_eq!(t.exons()[1].is_coding(), true);
        assert_eq!(t.exons()[2].is_coding(), true);
        assert_eq!(t.exons()[3].is_coding(), true);
        assert_eq!(t.exons()[4].is_coding(), true);
        assert_eq!(t.exons()[5].is_coding(), true);
        assert_eq!(t.exons()[6].is_coding(), true);
        assert_eq!(t.exons()[7].is_coding(), true);
    }
}
*/

#[cfg(test)]
mod test_ighm {
    use super::*;
    use crate::models::CdsStat;

    #[test]
    fn test_read() {
        let mut reader = Reader::from_file("tests/data/id-IGHM.gtf")
            .expect("Something failed reading the GTF file id-IGHM.gtf");
        let tr = match reader.transcripts() {
            Ok(res) => res,
            _ => panic!("No transcripts could be read"),
        };
        assert_eq!(tr.len(), 1);
        let t = tr.by_name("id-IGHM")[0];
        assert_eq!(t.exons().len(), 6);
        assert_eq!(t.cds_start().unwrap(), 106318298);
        assert_eq!(t.cds_end().unwrap(), 106322322);

        let start = t.start_codon();
        assert_eq!(start.len(), 1);

        let stop = t.stop_codon();
        assert_eq!(stop.len(), 1);

        assert!(t.exons()[0].is_coding());
        assert!(t.exons()[1].is_coding());
        assert!(t.exons()[2].is_coding());
        assert!(t.exons()[3].is_coding());
        assert!(t.exons()[4].is_coding());

        assert_eq!(t.cds_start_codon_stat(), CdsStat::Complete);
        assert_eq!(t.cds_stop_codon_stat(), CdsStat::Incomplete);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tests::transcripts;

    #[test]
    fn test_nm_001365057() {
        let transcripts = Reader::from_file("tests/data/NM_001365057.2.gtf")
            .unwrap()
            .transcripts()
            .unwrap();

        assert_eq!(
            transcripts.by_name("NM_001365057.2")[0],
            &transcripts::nm_001365057()
        )
    }

    #[test]
    fn test_nm_001365408() {
        let transcripts = Reader::from_file("tests/data/NM_001365408.1.gtf")
            .unwrap()
            .transcripts()
            .unwrap();

        assert_eq!(
            transcripts.by_name("NM_001365408.1")[0],
            &transcripts::nm_001365408()
        )
    }

    #[test]
    fn test_nm_001371720() {
        let transcripts = Reader::from_file("tests/data/NM_001371720.1.gtf")
            .unwrap()
            .transcripts()
            .unwrap();

        assert_eq!(
            transcripts.by_name("NM_001371720.1")[0],
            &transcripts::nm_001371720(true)
        )
    }

    #[test]
    fn test_nm_201550() {
        let transcripts = Reader::from_file("tests/data/NM_201550.4.gtf")
            .unwrap()
            .transcripts()
            .unwrap();

        assert_eq!(
            transcripts.by_name("NM_201550.4")[0],
            &transcripts::nm_201550()
        )
    }

    #[test]
    fn test_gene_line_returns_none() {
        // no newline at the end
        let transcript = "chr16\tncbiRefSeq.2021-05-17\tgene\t66969419\t66978999\t.\t+\t.\tgene_id \"CES2\"; transcript_id \"NM_001365408.1\"; gene_name \"CES2\";".as_bytes();
        let mut reader = Reader::new(transcript);
        assert!(reader.line().is_none());

        // with newline at the end
        let transcript = "chr16\tncbiRefSeq.2021-05-17\tgene\t66969419\t66978999\t.\t+\t.\tgene_id \"CES2\"; transcript_id \"NM_001365408.1\"; gene_name \"CES2\";\n".as_bytes();
        let mut reader = Reader::new(transcript);
        assert!(reader.line().is_none());

        // no transcript and no gene tag
        let transcript = "chr16\tncbiRefSeq.2021-05-17\tgene\t66969419\t66978999\t.\t+\t.\tgene_name \"CES2\";\n".as_bytes();
        let mut reader = Reader::new(transcript);
        assert!(reader.line().is_none());

        // no transcript tag
        let transcript = "chr16\tncbiRefSeq.2021-05-17\tgene\t66969419\t66978999\t.\t+\t.\tgene_id \"CES2\"; gene_name \"CES2\";\n".as_bytes();
        let mut reader = Reader::new(transcript);
        assert!(reader.line().is_none());

        // no gene tag
        let transcript = "chr16\tncbiRefSeq.2021-05-17\tgene\t66969419\t66978999\t.\t+\t.\ttranscript_id \"NM_001365408.1\"; gene_name \"CES2\";\n".as_bytes();
        let mut reader = Reader::new(transcript);
        assert!(reader.line().is_none());
    }

    #[test]
    fn test_gene_line_is_skipped() {
        let transcript = "chr16\tncbiRefSeq.2021-05-17\tgene\t66969419\t66978999\t.\t+\t.\tgene_id \"CES2\"; transcript_id \"NM_001365408.1\"; gene_name \"CES2\";\n\
        chr1\tncbiRefSeq.2021-05-17\texon\t206100298\t206100445\t.\t-\t.\tgene_id \"SRGAP2B\"; transcript_id \"NM_001385228.1_2\"; exon_number \"9\"; exon_id \"NM_001385228.1_2.9\"; gene_name \"SRGAP2B\";".as_bytes();
        let mut reader = Reader::new(transcript);
        assert!(reader
            .line()
            .expect("The record contains one proper line")
            .is_ok());
        assert!(reader.line().is_none());
    }

    #[test]
    fn test_gene_line_contains_invalid_data() {
        let transcript = "chr16\tncbiRefSeq.2021-05-17\tgene\t66969419\t66978999\t.\tINVALID\t.\tgene_id \"CES2\"; transcript_id \"NM_001365408.1\"; gene_name \"CES2\";\n".as_bytes();
        let mut reader = Reader::new(transcript);
        assert!(reader
            .line()
            .expect("The record contains one proper line")
            .is_err());
        assert!(reader.line().is_none());
    }

    #[test]
    fn test_without_gene_or_transcript_errors() {
        let transcript = "chr16\tncbiRefSeq.2021-05-17\texon\t66969419\t66978999\t.\t+\t.\tgene_id \"CES2\"; transcript_id \"NM_001365408.1\"; gene_name \"CES2\";\n".as_bytes();
        let mut reader = Reader::new(transcript);
        assert!(reader
            .line()
            .expect("The record contains one proper line")
            .is_ok());
        assert!(reader.line().is_none());

        let transcript = "chr16\tncbiRefSeq.2021-05-17\texon\t66969419\t66978999\t.\t+\t.\ttranscript_id \"NM_001365408.1\"; gene_name \"CES2\";\n".as_bytes();
        let mut reader = Reader::new(transcript);
        assert!(reader
            .line()
            .expect("The record contains one proper line")
            .is_err());
        assert!(reader.line().is_none());

        let transcript = "chr16\tncbiRefSeq.2021-05-17\texon\t66969419\t66978999\t.\t+\t.\tgene_id \"CES2\"; gene_name \"CES2\";\n".as_bytes();
        let mut reader = Reader::new(transcript);
        assert!(reader
            .line()
            .expect("The record contains one proper line")
            .is_err());
        assert!(reader.line().is_none());
    }

    #[test]
    fn test_5_utr() {
        let transcript = "chr16\tncbiRefSeq.2021-05-17\t5UTR\t66969419\t66978999\t.\t+\t.\tgene_id \"CES2\"; transcript_id \"NM_001365408.1\"; gene_name \"CES2\";".as_bytes();
        let mut reader = Reader::new(transcript);
        assert_eq!(
            reader
                .line()
                .expect("The record contains one proper line")
                .expect("The data is correct")
                .feature(),
            &GtfFeature::UTR5
        );

        let transcript = "chr16\tncbiRefSeq.2021-05-17\tfive_prime_utr\t66969419\t66978999\t.\t+\t.\tgene_id \"CES2\"; transcript_id \"NM_001365408.1\"; gene_name \"CES2\";".as_bytes();
        let mut reader = Reader::new(transcript);
        assert_eq!(
            reader
                .line()
                .expect("The record contains one proper line")
                .expect("The data is correct")
                .feature(),
            &GtfFeature::UTR5
        );
    }

    #[test]
    fn test_3_utr() {
        let transcript = "chr16\tncbiRefSeq.2021-05-17\t3UTR\t66969419\t66978999\t.\t+\t.\tgene_id \"CES2\"; transcript_id \"NM_001365408.1\"; gene_name \"CES2\";".as_bytes();
        let mut reader = Reader::new(transcript);
        assert_eq!(
            reader
                .line()
                .expect("The record contains one proper line")
                .expect("The data is correct")
                .feature(),
            &GtfFeature::UTR3
        );

        let transcript = "chr16\tncbiRefSeq.2021-05-17\tthree_prime_utr\t66969419\t66978999\t.\t+\t.\tgene_id \"CES2\"; transcript_id \"NM_001365408.1\"; gene_name \"CES2\";".as_bytes();
        let mut reader = Reader::new(transcript);
        assert_eq!(
            reader
                .line()
                .expect("The record contains one proper line")
                .expect("The data is correct")
                .feature(),
            &GtfFeature::UTR3
        );
    }

    #[test]
    fn test_utr() {
        let transcript = "chr16\tncbiRefSeq.2021-05-17\tUTR\t66969419\t66978999\t.\t+\t.\tgene_id \"CES2\"; transcript_id \"NM_001365408.1\"; gene_name \"CES2\";".as_bytes();
        let mut reader = Reader::new(transcript);
        assert_eq!(
            reader
                .line()
                .expect("The record contains one proper line")
                .expect("The data is correct")
                .feature(),
            &GtfFeature::UTR
        );
    }
}

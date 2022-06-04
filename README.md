# ATGlib

_ATGlib_ is a Rust library to work with genomic transcript data. It handles several file formats, such as GTF, GenePred(ext) and Refgene. You can generate bed files, fasta sequences or custom feature sequences.

If you are looking for an actual application, or a command line tool to work with transcripts, GTF files etc, use [ATG](https://crates.io/crates/atg) instead. It is using _ATGlib_ behind the scenes and provides a simple to use interface.


## Documentation
[The library API is mostly documented inline and available on docs.rs](https://docs.rs/atglib)

### Examples

#### Convert GTF to RefGene
```rust
use atglib::gtf::Reader;
use atglib::refgene::Writer;
use atglib::models::{TranscriptRead, TranscriptWrite};

let mut reader = Reader::from_file("tests/data/example.gtf")
    .unwrap_or_else(|_| panic!("Error opening input file."));

let mut writer = Writer::from_file("/dev/null")
    .unwrap_or_else(|_| panic!("Unable to open output file"));

let transcripts = reader.transcripts()
    .unwrap_or_else(|err| panic!("Error parsing GTF: {}", err));

match writer.write_transcripts(&transcripts) {
    Ok(_) => println!("Success"),
    Err(err) => panic!("Error writing RefGene file: {}", err)
};
```


## ToDo / Next tasks
- [ ] Compare transcripts from two different inputs
- [ ] use [Smartstring](https://crates.io/crates/smartstring) or [Smallstr](https://crates.io/crates/smallstr) for gene-symbol, transcript name and chromosome
- [ ] Parallelize input parsing
- [ ] Check if exons can be stored in smaller vec
- [ ] Use std::mem::replace to move out of attributes, e.g. in TranscriptBuilder and remove Copy/Clone traits <https://stackoverflow.com/questions/31307680/how-to-move-one-field-out-of-a-struct-that-implements-drop-trait>

## Known issues
### GTF parsing
- [ ] NM_001371720.1 has two book-ended exons (155160639-155161619 || 155161620-155162101). During input parsing, book-ended features are merged into one exon

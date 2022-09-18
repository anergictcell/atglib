[![Build](https://github.com/anergictcell/atglib/actions/workflows/test_repo.yml/badge.svg)](https://github.com/anergictcell/atglib/actions/workflows/test_repo.yml)
[![crates.io](https://img.shields.io/crates/v/atglib?color=#3fb911)](https://crates.io/crates/atg)
[![doc-rs](https://img.shields.io/docsrs/atglib/latest)](https://docs.rs/atglib/latest/atglib/)

# ATGlib

_ATGlib_ is a Rust library to work with genomic transcript data. It handles several file formats, such as GTF, GenePred(ext) and Refgene. You can generate bed files, fasta sequences or custom feature sequences.

_ATGlib_ is useful for Bioinformaticians and Geneticists. It can be used for all kind of data processing workflows that must handle transcript annotations. The main use case are conversions between various file formats and sanity checks for large scale QC.
It is primarily targeted for human genome, but should work equally well with all other genomes. 

If you are looking for an actual application, or a command line tool to work with transcripts, GTF files etc, use [ATG](https://crates.io/crates/atg) instead. It is using _ATGlib_ behind the scenes and provides a simple to use interface.

## Documentation
[The library API is mostly documented inline and available on docs.rs](https://docs.rs/atglib)

### Examples
For examples on how to use _ATGlib_ on a high-level, you can check the source code of the CLI tool [ATG](https://github.com/anergictcell/atg).


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

#### Work with transcripts directly
```rust
use atglib::tests::transcripts::standard_transcript;
let tx = standard_transcript();

println!(
    "The {} transcript {} is on the {} strand on chromsome {} and has {} exons",
    tx.gene(),
    tx.name(),
    tx.strand(),
    tx.chrom(),
    tx.exon_count()
);

if tx.is_coding() {
    println!(
        "The start and end positions of the CDS are: {}-{}",
        tx.cds_start().unwrap(),
        tx.cds_end().unwrap()
    );
} else {
    println!(
        "The transcript is non-coding and is located at {}:{}-{}",
        tx.chrom(),
        tx.tx_start(),
        tx.tx_end()
    );
}
```


## Limitations
I started _ATGlib_ as a private side project to learn `Rust`. I'm sure there are many parts of the code that are not idiomatic and can be improved. I'd be more than happy to receive feedback and suggestions for improvement. I also encourage everyone, who is interested, to help and contribute the _ATGlib_.

When it comes to functional correctness, I try my best to test the functionality on all levels. _ATGlib_ has a very good test coverage and I also run manual checks before every release with a huge test-set. However, I cannot guarantee that everything is correct, so please use it at your own risk and report any bugs and issues via Github.

## ToDo / Next tasks
- [ ] Add function to compare transcripts from two different inputs
- [ ] use [Smartstring](https://crates.io/crates/smartstring) or [Smallstr](https://crates.io/crates/smallstr) for gene-symbol, transcript name and chromosome
- [ ] Parallelize input parsing
- [ ] Check if exons can be stored in smaller vec
- [ ] Use std::mem::replace to move out of attributes, e.g. in TranscriptBuilder and remove Copy/Clone traits <https://stackoverflow.com/questions/31307680/how-to-move-one-field-out-of-a-struct-that-implements-drop-trait>
- [ ] Change `Codon` to `GenomicCodon`
- [ ] Update error handling and streamling error types

## Known issues
### GTF parsing
- [ ] NM_001371720.1 has two book-ended exons (155160639-155161619 || 155161620-155162101). During input parsing, book-ended features are merged into one exon

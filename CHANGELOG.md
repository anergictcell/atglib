# Changelog

# 0.2.2
- Unpin `serde`, the issue from 0.2.1 was fixed by maintainer
- [Dev] Add first building blocks for benchmarking

# 0.2.1
- Pin serde version (See https://github.com/serde-rs/serde/issues/2538 for context)

# 0.2
- Add GeneticCode to modify the translation table based on the applied genetic code. This change impacts some public functions of the QC-check module.
- Allow generic `Read + Seek` objects for FastaReader. This enables reading directly from S3 or other remote sources.
- Allow `FastaWriter` to write to different files (only one at a time). This means you don't have to initiate a new `FastaWriter` for every output file, but can reuse an existing instance and simply change the output writer.

# 0.1.3
- Add QC-check module to check if transcripts make sense with a given reference genome.

# 0.1.2
- Add ncExon feature-type to the feature sequence output for exons of non-coding transcript.
- Fix `Transcript::utr5_coordinates` to return an empty vector for non-coding transcripts, like `utr3_coordinates`.

# 0.1.1
- Fix the exon numbering to be strand dependent. That means the counting starts at the 5' most exon. Fixes https://github.com/anergictcell/atg/issues/16

# 0.1.0
- Initial version of the standalone library. This is the exact same functionality as in [ATG 0.5.0](https://crates.io/crates/atg/0.5.0).
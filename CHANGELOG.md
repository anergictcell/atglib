# Changelog

# 0.1.2
- Add ncExon feature-type to the feature sequence output for exons of non-coding transcript.
- Fix `Transcript::utr5_coordinates` to return an empty vector for non-coding transcripts, like `utr3_coordinates`.

# 0.1.1
- Fix the exon numbering to be strand dependent. That means the counting starts at the 5' most exon. Fixes https://github.com/anergictcell/atg/issues/16

# 0.1.0
- Initial version of the standalone library. This is the exact same functionality as in [ATG 0.5.0](https://crates.io/crates/atg/0.5.0).
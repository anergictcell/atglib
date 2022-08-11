//! Write to a custom format useful for SpliceAi splice predictions
//!
//! [SpliceAI](https://github.com/Illumina/SpliceAI) requires a gene annotation file.
//! The file format is very similar to `GenePred`, but has only one record per gene and
//! is thus not directly compatible with `GenePred`.
//! Even though it is possible to run SpliceAI when the annotation file is grouped by
//! transcripts, this [Writer](`crate::spliceai::Writer`) contains functionality to
//! provide a gene-based annotation file.
//!
//! The annotation file is created by merging all transcripts of a gene into a single
//! transcript. The merging considers all exons of all transcripts. See
//! [merge](`crate::utils::merge`) for details how exons are merged.
//!
//! ```text
//! Transcript 1   ----XXXXXX----XXXX------------XXXXXX---XXXX
//! Transcript 2   ----XXXXXX------------XXXX----XXXXXX---XXXX
//! Transcript 3   -------XXX------------------XXXXX------XXXX
//! Merged         ----OOOOOO----OOOO----OOOO--OOOOOOOO---OOOO
//! ```


mod writer;

pub use crate::spliceai::writer::Writer;

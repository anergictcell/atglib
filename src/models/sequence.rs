use core::str::FromStr;
use std::convert::TryFrom;
use std::fmt;
use std::fs::File;

use crate::fasta::FastaReader;
use crate::models::transcript::CoordinateVector;
use crate::utils::errors::{AtgError, FastaError};

// UTF-8 encoding of all nucleotides
const UPPERCASE_A: u8 = 0x41;
const UPPERCASE_C: u8 = 0x43;
const UPPERCASE_G: u8 = 0x47;
const UPPERCASE_T: u8 = 0x54;
const UPPERCASE_N: u8 = 0x4e;
const LOWERCASE_A: u8 = 0x61;
const LOWERCASE_C: u8 = 0x63;
const LOWERCASE_G: u8 = 0x67;
const LOWERCASE_T: u8 = 0x74;
const LOWERCASE_N: u8 = 0x64;

const LF: u8 = 0xa;
const CR: u8 = 0xd;

/// Nucleotide is a single DNA nucleotide (A C G T N)
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Nucleotide {
    A,
    C,
    G,
    T,
    N,
}

impl Nucleotide {
    /// Crates a `Nucleotide` from a character
    pub fn new(c: &char) -> Result<Self, AtgError> {
        match c {
            'a' | 'A' => Ok(Self::A),
            'c' | 'C' => Ok(Self::C),
            'g' | 'G' => Ok(Self::G),
            't' | 'T' => Ok(Self::T),
            'n' | 'N' => Ok(Self::N),
            _ => Err(AtgError::new("Invalid nucleotide")),
        }
    }

    /// Returns the complementary nucleotide
    pub fn complement(&self) -> Self {
        match self {
            Self::A => Self::T,
            Self::C => Self::G,
            Self::G => Self::C,
            Self::T => Self::A,
            Self::N => Self::N,
        }
    }

    // Returns the UTF-8 encoding of the Nucleotide string representation
    pub fn to_bytes(self) -> u8 {
        match self {
            Self::A => UPPERCASE_A,
            Self::C => UPPERCASE_C,
            Self::G => UPPERCASE_G,
            Self::T => UPPERCASE_T,
            Self::N => UPPERCASE_N,
        }
    }

    // Returns an &str of the Nucleotide
    pub fn to_str(self) -> &'static str {
        match self {
            Self::A => "A",
            Self::C => "C",
            Self::G => "G",
            Self::T => "T",
            Self::N => "N",
        }
    }
}

impl FromStr for Nucleotide {
    type Err = AtgError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "a" | "A" => Ok(Self::A),
            "c" | "C" => Ok(Self::C),
            "g" | "G" => Ok(Self::G),
            "t" | "T" => Ok(Self::T),
            "n" | "N" => Ok(Self::N),
            _ => Err(AtgError::new("Invalid nucleotide")),
        }
    }
}

impl fmt::Display for Nucleotide {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::A => "A",
                Self::C => "C",
                Self::G => "G",
                Self::T => "T",
                Self::N => "N",
            }
        )
    }
}

impl TryFrom<&char> for Nucleotide {
    type Error = AtgError;
    fn try_from(c: &char) -> Result<Self, Self::Error> {
        match c {
            'a' | 'A' => Ok(Self::A),
            'c' | 'C' => Ok(Self::C),
            'g' | 'G' => Ok(Self::G),
            't' | 'T' => Ok(Self::T),
            'n' | 'N' => Ok(Self::N),
            '\n' | '\r' => Err(AtgError::new("newline")),
            _ => panic!("invalid nucleotide {}", c),
        }
    }
}

impl TryFrom<&u8> for Nucleotide {
    type Error = AtgError;
    fn try_from(b: &u8) -> Result<Nucleotide, AtgError> {
        match b {
            &LOWERCASE_A | &UPPERCASE_A => Ok(Self::A),
            &LOWERCASE_C | &UPPERCASE_C => Ok(Self::C),
            &LOWERCASE_G | &UPPERCASE_G => Ok(Self::G),
            &LOWERCASE_T | &UPPERCASE_T => Ok(Self::T),
            &LOWERCASE_N | &UPPERCASE_N => Ok(Self::N),
            &LF | &CR => Err(AtgError::new("newline")),
            _ => panic!("invalid nucleotide {}", b),
        }
    }
}

impl From<&Nucleotide> for char {
    fn from(n: &Nucleotide) -> Self {
        match n {
            Nucleotide::A => 'A',
            Nucleotide::C => 'C',
            Nucleotide::G => 'G',
            Nucleotide::T => 'T',
            Nucleotide::N => 'N',
        }
    }
}

impl From<&Nucleotide> for u8 {
    fn from(n: &Nucleotide) -> u8 {
        n.to_bytes()
    }
}

/// A DNA sequence consisting of Nucleotides.
///
/// It provides some utility methods, like
///[`reverse_complement`](`Sequence::reverse_complement`)
pub struct Sequence {
    sequence: Vec<Nucleotide>,
}

impl FromStr for Sequence {
    type Err = AtgError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut sequence: Vec<Nucleotide> = vec![];
        for c in s.chars() {
            sequence.push(Nucleotide::new(&c)?)
        }
        Ok(Self { sequence })
    }
}

impl fmt::Display for Sequence {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut s = String::with_capacity(self.len());
        for n in &self.sequence {
            s.push(n.into())
        }
        write!(f, "{}", s)
    }
}

impl Default for Sequence {
    fn default() -> Self {
        Self::new()
    }
}

impl Sequence {
    /// Creates a new sequence
    ///
    /// # Examples
    ///
    /// ```rust
    /// use atglib::models::Sequence;
    ///
    /// let seq = Sequence::new();
    /// assert_eq!(seq.len(), 0)
    /// ```
    pub fn new() -> Self {
        Sequence {
            sequence: Vec::new(),
        }
    }

    /// Creates a new sequence with the specified capacity
    ///
    /// Use this method if you know in advance the final size of the Sequence.
    /// It creates an empty Sequence, but one with an initial buffer that can
    /// hold capacity Nucleotides.
    ///
    /// It is important to note that although the returned Sequence has the capacity specified,
    /// the Sequence will have a zero length
    ///
    /// # Examples
    ///
    /// ```rust
    /// use atglib::models::Sequence;
    ///
    /// let mut seq = Sequence::with_capacity(5);
    /// assert_eq!(seq.len(), 0);
    ///
    /// // this will not re-allocate memory, since all nucleotides fit into capacity
    /// for c in vec!['A', 'C', 'G', 'T', 'N'] {
    ///     seq.push_char(&c).unwrap();
    /// }
    /// assert_eq!(seq.len(), 5);
    /// ```
    ///
    pub fn with_capacity(capacity: usize) -> Self {
        Sequence {
            sequence: Vec::with_capacity(capacity),
        }
    }

    /// Creates a new `Sequence` from a raw bytes nucleotide sequence, ignoring newlines
    ///
    /// The `len` value is not required to be correct, it helps with allocating the right
    /// amount of memory for the Sequence. The length of the byte vector can be misleading
    /// because it might contain newlines and other whitespace that is not added to the
    /// Sequence.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use atglib::models::Sequence;
    ///
    /// let seq = Sequence::from_raw_bytes("AC".as_bytes(), 2).unwrap();
    /// assert_eq!(seq.len(), 2);
    /// let seq = Sequence::from_raw_bytes("A\nC\r\nGT".as_bytes(), 2).unwrap();
    /// assert_eq!(seq.len(), 4);
    /// ```
    ///
    pub fn from_raw_bytes(bytes: &[u8], len: usize) -> Result<Self, AtgError> {
        let mut seq = Self::with_capacity(len);
        for b in bytes {
            if let Ok(n) = Nucleotide::try_from(b) {
                seq.push(n)?
            }
        }
        Ok(seq)
    }

    /// Creates a new Sequence from genomic coordinates
    ///
    /// Use this method if you want to have the sequenced of gapped features,
    /// e.g. exons of a Transcript.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use atglib::fasta::FastaReader;
    /// use atglib::models::{Sequence, Strand};
    /// use atglib::tests::transcripts::standard_transcript;
    ///
    /// let tx = standard_transcript();
    /// let coordinates = tx.cds_coordinates();
    ///
    /// let mut fasta_reader = FastaReader::from_file("tests/data/small.fasta").unwrap();
    ///
    /// let seq = Sequence::from_coordinates(&coordinates, &Strand::Plus, &mut fasta_reader).unwrap();
    /// assert_eq!(seq.len(), 11);
    /// ```
    pub fn from_coordinates(
        coordinates: &CoordinateVector,
        strand: &crate::models::Strand,
        fasta_reader: &mut FastaReader<File>,
    ) -> Result<Sequence, FastaError> {
        let capacity: u32 = coordinates.iter().map(|x| x.2 - x.1 + 1).sum();
        let mut seq = Sequence::with_capacity(capacity as usize);

        for segment in coordinates {
            seq.append(fasta_reader.read_sequence(segment.0, segment.1.into(), segment.2.into())?)
        }
        if strand == &crate::models::Strand::Minus {
            seq.reverse_complement()
        }
        Ok(seq)
    }

    /// Returns the length of the Sequence
    /// # Examples
    ///
    /// ```rust
    /// use atglib::models::Sequence;
    ///
    /// let seq = Sequence::from_raw_bytes("AC".as_bytes(), 2).unwrap();
    /// assert_eq!(seq.len(), 2);
    /// ```
    ///
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    /// Returns true if the Sequence contains no Nucleotides.
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    /// Appends a `char` as Nucleotide to the back of a collection.
    ///
    /// # Examples
    /// ```rust
    /// use atglib::models::{Nucleotide, Sequence};
    ///
    /// let mut seq = Sequence::from_raw_bytes("AC".as_bytes(), 2).unwrap();
    /// assert_eq!(seq.to_string(), "AC".to_string());
    ///
    /// seq.push_char(&'T').unwrap();
    /// assert_eq!(seq.to_string(), "ACT".to_string());
    pub fn push_char(&mut self, c: &char) -> Result<(), AtgError> {
        self.sequence.push(Nucleotide::try_from(c)?);
        Ok(())
    }

    /// Appends a Nucleotide to the back of a collection.
    ///
    /// # Examples
    /// ```rust
    /// use atglib::models::{Nucleotide, Sequence};
    ///
    /// let mut seq = Sequence::from_raw_bytes("AC".as_bytes(), 2).unwrap();
    /// assert_eq!(seq.to_string(), "AC".to_string());
    ///
    /// seq.push(Nucleotide::new(&'T').unwrap()).unwrap();
    /// assert_eq!(seq.to_string(), "ACT".to_string());
    /// ```
    pub fn push(&mut self, n: Nucleotide) -> Result<(), AtgError> {
        self.sequence.push(n);
        Ok(())
    }

    /// Clears the sequence, removing all nucleotides.
    ///
    /// Note that this method has no effect on the allocated capacity.
    ///
    /// # Examples
    /// ```rust
    /// use atglib::models::{Sequence};
    ///
    /// let mut seq = Sequence::from_raw_bytes("AC".as_bytes(), 2).unwrap();
    /// seq.clear();
    /// assert_eq!(seq.len(), 0);
    /// ```
    pub fn clear(&mut self) {
        self.sequence.clear()
    }

    /// Moves all the elements of `other` into `Self`, leaving `other` empty.
    ///
    /// # Examples
    /// ```rust
    /// use atglib::models::{Nucleotide, Sequence};
    ///
    /// let mut seq = Sequence::from_raw_bytes("AC".as_bytes(), 2).unwrap();
    /// assert_eq!(seq.to_string(), "AC".to_string());
    ///
    /// let seq_2 = Sequence::from_raw_bytes("GT".as_bytes(), 2).unwrap();
    /// seq.append(seq_2);
    /// assert_eq!(seq.to_string(), "ACGT".to_string());
    pub fn append(&mut self, other: Sequence) {
        self.sequence.append(&mut other.into_inner())
    }

    /// Unwraps the Sequence, returning the underlying Vector of [`Nucleotide`]s
    fn into_inner(self) -> Vec<Nucleotide> {
        self.sequence
    }

    /// Changes `Self` to the complementary sequence
    ///
    /// # Examples
    /// ```rust
    /// use atglib::models::Sequence;
    ///
    /// let mut seq = Sequence::from_raw_bytes("AC".as_bytes(), 2).unwrap();
    /// assert_eq!(seq.to_string(), "AC".to_string());
    ///
    /// seq.complement();
    /// assert_eq!(seq.to_string(), "TG".to_string());
    /// ```
    pub fn complement(&mut self) {
        for n in &mut self.sequence {
            *n = n.complement();
        }
    }

    /// Reverses the `Sequence`, in place
    ///
    /// # Examples
    /// ```rust
    /// use atglib::models::Sequence;
    ///
    /// let mut seq = Sequence::from_raw_bytes("AC".as_bytes(), 2).unwrap();
    /// assert_eq!(seq.to_string(), "AC".to_string());
    ///
    /// seq.reverse();
    /// assert_eq!(seq.to_string(), "CA".to_string());
    /// ```
    pub fn reverse(&mut self) {
        self.sequence.reverse()
    }

    /// Changes `Self` into the reverse complement sequence
    ///
    /// # Examples
    /// ```rust
    /// use atglib::models::Sequence;
    ///
    /// let mut seq = Sequence::from_raw_bytes("AC".as_bytes(), 2).unwrap();
    /// assert_eq!(seq.to_string(), "AC".to_string());
    ///
    /// seq.reverse_complement();
    /// assert_eq!(seq.to_string(), "GT".to_string());
    /// ```
    pub fn reverse_complement(&mut self) {
        self.reverse();
        self.complement();
    }

    /// Returns the Sequence as a byte array of UTF-8 encoded nucleotides
    ///
    /// # Examples
    /// ```rust
    /// use atglib::models::Sequence;
    ///
    /// let seq = Sequence::from_raw_bytes("AC".as_bytes(), 2).unwrap();
    /// assert_eq!(seq.to_bytes(), [0x41, 0x43]);
    /// ```
    pub fn to_bytes(&self) -> Vec<u8> {
        self.sequence.iter().map(|n| n.to_bytes()).collect()
    }

    /// Writes the sequence into a target String
    ///
    /// The function does not make any assumption about the state of the
    /// target string. The client is expected to clear the string beforehand
    ///
    /// # Examples
    /// ```rust
    /// use atglib::models::Sequence;
    ///
    /// let seq = Sequence::from_raw_bytes("AC".as_bytes(), 2).unwrap();
    /// let mut my_string = String::new();
    /// seq.write_into_string(& mut my_string);
    /// assert_eq!(my_string, "AC");
    /// assert_eq!(my_string.len(), 2);
    /// ```
    pub fn write_into_string(&self, target: &mut String) {
        for c in &self.sequence {
            target.push(c.into())
        }
    }

    /// Returns an iterator over chunk_size [`crate::models::Nucleotide`]s at a time,
    /// starting at the beginning of the Sequence.
    ///
    /// The chunks are slices and do not overlap. If chunk_size does not divide
    /// the length of the slice, then the last chunk will not have length chunk_size.
    ///
    /// # Panics
    /// Panics if chunk_size is 0.
    ///
    /// # Examples
    /// ```rust
    /// use atglib::models::{Nucleotide, Sequence};
    ///
    /// let seq = Sequence::from_raw_bytes("ATGCTA".as_bytes(), 2).unwrap();
    /// let mut iter = seq.chunks(3);
    /// assert_eq!(iter.next().unwrap(), &[Nucleotide::A, Nucleotide::T, Nucleotide::G]);
    /// assert_eq!(iter.next().unwrap(), &[Nucleotide::C, Nucleotide::T, Nucleotide::A]);
    /// ```
    pub fn chunks(&self, chunk_size: usize) -> std::slice::Chunks<'_, Nucleotide> {
        self.sequence.chunks(chunk_size)
    }

    /// Returns the leftmost position of other in self
    ///
    /// # Panics
    /// Panics if other is empty
    ///
    /// # Examples
    /// ```rust
    /// use atglib::models::{Nucleotide, Sequence};
    ///
    /// let seq = Sequence::from_raw_bytes("ATGCTA".as_bytes(), 2).unwrap();
    /// assert_eq!(seq.position([Nucleotide::T, Nucleotide::A]), Some(4));
    /// ```
    pub fn position<T>(&self, other: T) -> Option<usize>
    where
        T: AsRef<[Nucleotide]>,
    {
        let query = other.as_ref();
        assert!(
            !query.is_empty(),
            "empty sequence was passed to Sequence::position"
        );
        for (i, nuc) in self.sequence[0..self.sequence.len() - query.len() + 1]
            .iter()
            .enumerate()
        {
            if nuc == &query[0] {
                let subsequence = &self.sequence[i..i + query.len()];
                if subsequence == query {
                    return Some(i);
                }
            }
        }
        None
    }

    /// Returns true if other is a subsequence of self
    ///
    /// # Panics
    /// Panics if other is empty
    ///
    /// # Examples
    /// ```rust
    /// use atglib::models::{Nucleotide, Sequence};
    ///
    /// let seq = Sequence::from_raw_bytes("ATGCTA".as_bytes(), 2).unwrap();
    /// assert_eq!(seq.contains([Nucleotide::T, Nucleotide::A]), true);
    /// assert_eq!(seq.contains([Nucleotide::A, Nucleotide::A]), false);
    /// ```
    pub fn contains<T>(&self, other: T) -> bool
    where
        T: AsRef<[Nucleotide]>,
    {
        self.position(other).is_some()
    }

    /// # Examples
    /// ```rust
    /// use atglib::models::{Nucleotide, Sequence};
    ///
    /// let seq = Sequence::from_raw_bytes("ATGCTA".as_bytes(), 2).unwrap();
    /// let seq_2 = Sequence::from_raw_bytes("ATGCTA".as_bytes(), 2).unwrap();
    /// let seq_3 = Sequence::from_raw_bytes("ATG".as_bytes(), 2).unwrap();
    ///
    /// assert_eq!(seq.equals(&seq_2), true);
    /// assert_eq!(seq.equals(&seq_3), false);
    ///
    /// assert_eq!(seq_3.equals([Nucleotide::A, Nucleotide::T, Nucleotide::G]), true);
    /// ```
    pub fn equals<T>(&self, other: T) -> bool
    where
        T: AsRef<[Nucleotide]>,
    {
        let query = other.as_ref();
        if query.len() != self.len() {
            return false;
        }
        for (a, b) in query.iter().zip(self.sequence.iter()) {
            if a != b {
                return false;
            }
        }
        true
    }
}

impl AsRef<[Nucleotide]> for Sequence {
    fn as_ref(&self) -> &[Nucleotide] {
        &self.sequence
    }
}

impl From<Sequence> for Vec<u8> {
    fn from(s: Sequence) -> Vec<u8> {
        s.to_bytes()
    }
}

/// implementing slice indexing operations for Sequence
/// so that seq[1..3] operations are possible.
impl<Idx> std::ops::Index<Idx> for Sequence
where
    Idx: std::slice::SliceIndex<[Nucleotide]>,
{
    type Output = Idx::Output;

    fn index(&self, idx: Idx) -> &Self::Output {
        &self.sequence[idx]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_create_sequence() {
        let s = "ATCGACGATCGATCGATGAGCGATCGACGATCGCGCTATCGCTA";
        let seq = Sequence::from_str(&s).unwrap();

        assert_eq!(seq.len(), 44);
        assert_eq!(seq.to_string(), s.to_string())
    }

    #[test]
    fn test_chunks() {
        let s = "ATCGACGATCGATCGATGAGCGATCGACGATCGCGCTATCGCTA";
        let seq = Sequence::from_str(&s).unwrap();

        let mut iter = seq.chunks(3);
        assert_eq!(
            iter.next().unwrap(),
            &[Nucleotide::A, Nucleotide::T, Nucleotide::C]
        );
        assert_eq!(
            iter.next().unwrap(),
            &[Nucleotide::G, Nucleotide::A, Nucleotide::C]
        );
    }

    #[test]
    fn test_contains() {
        let s = "ATGCGA";
        let seq = Sequence::from_str(&s).unwrap();

        assert_eq!(seq.contains(vec![Nucleotide::A]), true);
        assert_eq!(seq.contains(vec![Nucleotide::C]), true);
        assert_eq!(seq.contains(vec![Nucleotide::G]), true);
        assert_eq!(seq.contains(vec![Nucleotide::T]), true);
        assert_eq!(seq.contains(vec![Nucleotide::N]), false);
        assert_eq!(seq.contains(vec![Nucleotide::A, Nucleotide::T]), true);
        assert_eq!(seq.contains(vec![Nucleotide::T, Nucleotide::G]), true);
        assert_eq!(seq.contains(vec![Nucleotide::A, Nucleotide::C]), false);
        assert_eq!(seq.contains(vec![Nucleotide::G, Nucleotide::T]), false);
        assert_eq!(seq.contains(vec![Nucleotide::G, Nucleotide::A]), true);
        assert_eq!(
            seq.contains(vec![Nucleotide::C, Nucleotide::G, Nucleotide::A]),
            true
        );
        assert_eq!(
            seq.contains(vec![Nucleotide::G, Nucleotide::A, Nucleotide::C]),
            false
        );
    }

    #[test]
    #[should_panic]
    fn test_contains_fails() {
        let s = "ATGCGA";
        let seq = Sequence::from_str(&s).unwrap();
        assert_eq!(seq.contains(vec![]), true);
    }
}

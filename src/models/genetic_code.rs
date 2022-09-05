// The genetic code is stored in an array with 64 elements (one for each triplet)
// There is one standard implementation and several variants. Variants only overwrite
// the specific differences they have compared to the standard code.
// The lookup occurs through defined calculations:
// A => positions 0-15
// C => posisions 16-31
// G => positions 32-47
// T => positions 48-63
// https://www.ncbi.nlm.nih.gov/IEB/ToolBox/SDKDOCS/SEQFEAT.HTML

use std::{convert::TryFrom, fmt::Display};

pub use crate::models::sequence::Nucleotide;
use crate::utils::errors::AtgError;

type GeneticCodeLookup = [AminoAcid; 64];

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AminoAcid {
    Ter,
    A,
    C,
    D,
    E,
    F,
    G,
    H,
    I,
    K,
    L,
    M,
    N,
    P,
    Q,
    R,
    S,
    T,
    V,
    W,
    Y,
}

impl AminoAcid {
    pub fn single_letter(&self) -> char {
        match self {
            AminoAcid::A => 'A',
            AminoAcid::C => 'C',
            AminoAcid::D => 'D',
            AminoAcid::E => 'E',
            AminoAcid::F => 'F',
            AminoAcid::G => 'G',
            AminoAcid::H => 'H',
            AminoAcid::I => 'I',
            AminoAcid::K => 'K',
            AminoAcid::L => 'L',
            AminoAcid::M => 'M',
            AminoAcid::N => 'N',
            AminoAcid::P => 'P',
            AminoAcid::Q => 'Q',
            AminoAcid::R => 'R',
            AminoAcid::S => 'S',
            AminoAcid::T => 'T',
            AminoAcid::V => 'V',
            AminoAcid::W => 'W',
            AminoAcid::Y => 'Y',
            AminoAcid::Ter => '*',
        }
    }
}

impl AsRef<str> for AminoAcid {
    fn as_ref(&self) -> &str {
        match self {
            AminoAcid::A => "Ala",
            AminoAcid::C => "Cys",
            AminoAcid::D => "Asp",
            AminoAcid::E => "Glu",
            AminoAcid::F => "Phe",
            AminoAcid::G => "Gly",
            AminoAcid::H => "His",
            AminoAcid::I => "Ile",
            AminoAcid::K => "Lys",
            AminoAcid::L => "Leu",
            AminoAcid::M => "Met",
            AminoAcid::N => "Asn",
            AminoAcid::P => "Pro",
            AminoAcid::Q => "Gln",
            AminoAcid::R => "Arg",
            AminoAcid::S => "Ser",
            AminoAcid::T => "Thr",
            AminoAcid::V => "Val",
            AminoAcid::W => "Trp",
            AminoAcid::Y => "Tyr",
            AminoAcid::Ter => "Ter",
        }
    }
}

impl Display for AminoAcid {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.as_ref())
    }
}

impl TryFrom<char> for AminoAcid {
    type Error = AtgError;
    fn try_from(s: char) -> Result<AminoAcid, AtgError> {
        match s {
            'A' => Ok(AminoAcid::A),
            'C' => Ok(AminoAcid::C),
            'D' => Ok(AminoAcid::D),
            'E' => Ok(AminoAcid::E),
            'F' => Ok(AminoAcid::F),
            'G' => Ok(AminoAcid::G),
            'H' => Ok(AminoAcid::H),
            'I' => Ok(AminoAcid::I),
            'K' => Ok(AminoAcid::K),
            'L' => Ok(AminoAcid::L),
            'M' => Ok(AminoAcid::M),
            'N' => Ok(AminoAcid::N),
            'P' => Ok(AminoAcid::P),
            'Q' => Ok(AminoAcid::Q),
            'R' => Ok(AminoAcid::R),
            'S' => Ok(AminoAcid::S),
            'T' => Ok(AminoAcid::T),
            'V' => Ok(AminoAcid::V),
            'W' => Ok(AminoAcid::W),
            'Y' => Ok(AminoAcid::Y),
            '*' => Ok(AminoAcid::Ter),
            _ => Err(AtgError::new("Invalid amino acid")),
        }
    }
}

impl TryFrom<&str> for AminoAcid {
    type Error = AtgError;
    fn try_from(s: &str) -> Result<AminoAcid, AtgError> {
        match s {
            "A" | "Ala" => Ok(AminoAcid::A),
            "C" | "Cys" => Ok(AminoAcid::C),
            "D" | "Asp" => Ok(AminoAcid::D),
            "E" | "Glu" => Ok(AminoAcid::E),
            "F" | "Phe" => Ok(AminoAcid::F),
            "G" | "Gly" => Ok(AminoAcid::G),
            "H" | "His" => Ok(AminoAcid::H),
            "I" | "Ile" => Ok(AminoAcid::I),
            "K" | "Lys" => Ok(AminoAcid::K),
            "L" | "Leu" => Ok(AminoAcid::L),
            "M" | "Met" => Ok(AminoAcid::M),
            "N" | "Asn" => Ok(AminoAcid::N),
            "P" | "Pro" => Ok(AminoAcid::P),
            "Q" | "Gln" => Ok(AminoAcid::Q),
            "R" | "Arg" => Ok(AminoAcid::R),
            "S" | "Ser" => Ok(AminoAcid::S),
            "T" | "Thr" => Ok(AminoAcid::T),
            "V" | "Val" => Ok(AminoAcid::V),
            "W" | "Trp" => Ok(AminoAcid::W),
            "Y" | "Tyr" => Ok(AminoAcid::Y),
            "Ter" | "*" | "Stop" => Ok(AminoAcid::Ter),
            _ => Err(AtgError::new("Invalid amino acid")),
        }
    }
}

/// The genetic code is basically a lookup table from DNA codons to AminoAcids
pub struct GeneticCode {
    code: GeneticCodeLookup,
}

impl Default for GeneticCode {
    /// Crates the [standard genetic code](https://en.wikipedia.org/wiki/DNA_codon_table)
    fn default() -> GeneticCode {
        GeneticCode {
            code: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
                .chars()
                .map(|c| c.try_into().unwrap()) // cannot fail
                .collect::<Vec<AminoAcid>>()
                .try_into()
                .unwrap(), // cannot fail
        }
    }
}

impl GeneticCode {
    /// Creates the genetic code of [vertrebrate mitochondria](https://en.wikipedia.org/wiki/Vertebrate_mitochondrial_code)
    pub fn vertebrate_mitochondria() -> GeneticCode {
        GeneticCode {
            code: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG"
                .chars()
                .map(|c| c.try_into().unwrap()) // cannot fail
                .collect::<Vec<AminoAcid>>()
                .try_into()
                .unwrap(), // cannot fail
        }
    }

    /// Translates a codon into an AminoAcid
    ///
    /// If the codon contains an `N` nucleotide, the method will return `AtgError`
    pub fn translate(&self, codon: &[Nucleotide; 3]) -> Result<AminoAcid, AtgError> {
        let first_offset = codon[0].as_ncbi_int()? * 16;
        let second_offset = codon[1].as_ncbi_int()? * 4;
        let third_offset = codon[2].as_ncbi_int()?;
        Ok(self.code[first_offset + second_offset + third_offset])
    }

    /// Returns a vector of all codons that code for the AminoAcid
    pub fn reverse_lookup(&self, aa: &AminoAcid) -> Vec<[Nucleotide; 3]> {
        self.code
            .iter()
            .enumerate()
            .filter(|(_, elmt)| *elmt == aa)
            .map(|(idx, _)| {
                // `idx` cannot be larger than 63 (len of `self.code`): `pos1` must be 0 <= 3
                let (pos1, remainder) = (idx / 16, idx % 16);
                // `remainder` cannot be larger than 15: `pos2` must be 0 <= 3
                let (pos2, pos3) = (remainder / 4, remainder % 4);
                [
                    Nucleotide::from(pos1), // cannot fail, due to 0 <= pos <= 3 conditions, described above
                    Nucleotide::from(pos2), // cannot fail, due to 0 <= pos <= 3 conditions, described above
                    Nucleotide::from(pos3), // cannot fail, due to 0 <= pos <= 3 conditions, described above
                ]
            })
            .collect::<Vec<[Nucleotide; 3]>>()
    }

    /// Returns all possible Stop codons
    pub fn stop_codons(&self) -> Vec<[Nucleotide; 3]> {
        self.reverse_lookup(&AminoAcid::Ter)
    }

    /// Returns true if the codon is a stop codon
    ///
    /// This method can only handle codons with exaccly 3 nucleotides. All other input will return `false`
    pub fn is_stop_codon(&self, codon: &[Nucleotide]) -> bool {
        if codon.len() != 3 {
            return false;
        }

        matches!(
            self.translate(&[codon[0], codon[1], codon[2]]),
            Ok(AminoAcid::Ter)
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_transcript() {
        let code = GeneticCode::default();
        assert_eq!(
            code.translate(&[Nucleotide::A, Nucleotide::T, Nucleotide::G])
                .unwrap()
                .as_ref(),
            "Met"
        );
        assert_eq!(
            code.translate(&[Nucleotide::T, Nucleotide::T, Nucleotide::T])
                .unwrap()
                .as_ref(),
            "Phe"
        );
        assert_eq!(
            code.translate(&[Nucleotide::T, Nucleotide::A, Nucleotide::A])
                .unwrap()
                .as_ref(),
            "Ter"
        );
        assert_eq!(
            code.translate(&[Nucleotide::T, Nucleotide::A, Nucleotide::G])
                .unwrap()
                .as_ref(),
            "Ter"
        );
    }

    #[test]
    fn test_reverse_lookup() {
        let code = GeneticCode::default();
        for c in "ACDEFGHIKLMNPQRSTVWY*".chars() {
            let aa = AminoAcid::try_from(c).unwrap();
            for nt in code.reverse_lookup(&aa) {
                assert_eq!(aa, code.translate(&nt).unwrap());
            }
        }
    }

    #[test]
    fn test_stop_codon() {
        let code = GeneticCode::default();
        let ter = code.stop_codons();
        assert_eq!(ter.len(), 3);
        assert!(ter.contains(&[Nucleotide::T, Nucleotide::A, Nucleotide::A]));
        assert!(ter.contains(&[Nucleotide::T, Nucleotide::A, Nucleotide::G]));
        assert!(ter.contains(&[Nucleotide::T, Nucleotide::G, Nucleotide::A]));

        assert!(!ter.contains(&[Nucleotide::A, Nucleotide::G, Nucleotide::A]));
        assert!(!ter.contains(&[Nucleotide::A, Nucleotide::G, Nucleotide::G]));
    }

    #[test]
    fn test_mito_stop_codon() {
        let mito_code = GeneticCode::vertebrate_mitochondria();
        let mito_ter = mito_code.stop_codons();
        assert_eq!(mito_ter.len(), 4);
        assert!(mito_ter.contains(&[Nucleotide::T, Nucleotide::A, Nucleotide::A]));
        assert!(mito_ter.contains(&[Nucleotide::T, Nucleotide::A, Nucleotide::G]));
        assert!(mito_ter.contains(&[Nucleotide::A, Nucleotide::G, Nucleotide::A]));
        assert!(mito_ter.contains(&[Nucleotide::A, Nucleotide::G, Nucleotide::G]));

        assert!(!mito_ter.contains(&[Nucleotide::T, Nucleotide::G, Nucleotide::A]));
    }
}

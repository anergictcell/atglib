use std::fmt;

use crate::models::sequence::Nucleotide;
use crate::utils::errors::AtgError;
use crate::models::AminoAcid;


/// The genetic code lookup table is an Array with 64 amino acids
/// The lookup occurs through defined calculations:
/// T => positions 0-15
/// C => posisions 16-31
/// A => positions 32-47
/// G => positions 48-63
/// https://www.ncbi.nlm.nih.gov/IEB/ToolBox/SDKDOCS/SEQFEAT.HTML
type GeneticCodeLookup = [AminoAcid; 64];


/// The genetic code is basically a lookup table from DNA codons to AminoAcids
///
/// # Note
/// The genetic code does not support alternative start codons. *ATGlib* considers only `ATG` as start codon.
///
/// # Examples
/// ```
/// use atglib::models::{AminoAcid, GeneticCode, Nucleotide};
/// let code = GeneticCode::default();
/// assert_eq!(
///     code.translate(&[Nucleotide::A, Nucleotide::T, Nucleotide::G])
///         .unwrap(),
///     AminoAcid::M
/// );
/// ```
#[derive(Debug)]
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

impl PartialEq for GeneticCode {
    fn eq(&self, other: &Self) -> bool {
        self.code == other.code
    }
}
impl Eq for GeneticCode {}


impl fmt::Display for GeneticCode {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.code.iter().map(|aa| aa.single_letter()).collect::<String>())
    }
}

impl GeneticCode {
    /// Creates a new, custom genetic code
    ///
    /// The aa_table must be the amino acid translation for each possible codon, in the same order
    /// as the [NCBI genetic code tables](https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/gc.prt)
    ///
    /// The genetic code translation is provided as an array with 64 elements (one for each codon)
    /// # Examples
    /// ```
    /// use atglib::models::{AminoAcid, GeneticCode, Nucleotide};
    /// // use the genetic code table for `Yeast Mitochondrial`
    /// let code = GeneticCode::new("FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG").unwrap();
    /// assert_eq!(
    ///     code.translate(&[Nucleotide::T, Nucleotide::G, Nucleotide::A])
    ///         .unwrap(),
    ///     AminoAcid::W
    /// );
    /// ```
    pub fn new(aa_table: &str) -> Result<GeneticCode, AtgError> {
        if aa_table.len() != 64 {
            return Err(AtgError::new("aa_table has wrong length. 64 amino acids are required"))
        }

        let allowed_chars = "*ACDEFGHIKLMNPQRSTVWY";

        for letter in aa_table.chars() {
            if !allowed_chars.contains(letter) {
                return Err(AtgError::new(format!("Invalid amino acid {}", letter)))
            }
        }
        Ok(GeneticCode {
            code: aa_table
                .chars()
                .map(|c| c.try_into().unwrap()) // cannot fail
                .collect::<Vec<AminoAcid>>()
                .try_into()
                .unwrap(), // cannot fail
        })
    }

    /// Creates the genetic code of [vertrebrate mitochondria](https://en.wikipedia.org/wiki/Vertebrate_mitochondrial_code)
    ///
    /// # Examples
    /// ```
    /// use atglib::models::{AminoAcid, GeneticCode, Nucleotide};
    /// let code = GeneticCode::vertebrate_mitochondria();
    /// assert_eq!(
    ///     code.translate(&[Nucleotide::A, Nucleotide::G, Nucleotide::A])
    ///         .unwrap(),
    ///     AminoAcid::Ter
    /// );
    /// ```
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

    /// Tries to create a genetic code from either a known code or from the provided lookup
    ///
    /// This method can be used as a helper function for user-facing applications where users
    /// can provide either the name of the code to use (e.g. `vertebrate_mitochondria`) or provide the
    /// actual lookup table (e.g. `FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG`)
    ///
    /// # Examples
    /// ```
    /// use atglib::models::{AminoAcid, GeneticCode, Nucleotide};
    /// 
    /// assert_eq!(
    ///     GeneticCode::vertebrate_mitochondria(),
    ///     GeneticCode::guess("vertebrate_mitochondria").unwrap()
    /// );
    ///
    /// assert_eq!(
    ///     GeneticCode::vertebrate_mitochondria(),
    ///     GeneticCode::guess("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG").unwrap()
    /// );
    /// ```
    pub fn guess(code: &str) -> Result<GeneticCode, AtgError> {
        match code {
            "standard" | "default" => Ok(GeneticCode::default()),
            "vertebrate_mitochondria" => Ok(GeneticCode::vertebrate_mitochondria()),
            _ => match GeneticCode::new(code) {
                Ok(code) => Ok(code),
                Err(_) => Err(AtgError::new("Genetic code not known or invalid"))
            }
        }
    }


    /// Translates a provided codon into an AminoAcid
    ///
    /// If the codon contains an `N` nucleotide, the method will return `AtgError`
    ///
    /// # Examples
    /// ```
    /// use atglib::models::{AminoAcid, GeneticCode, Nucleotide};
    /// let code = GeneticCode::default();
    /// let aa = code.translate(&[Nucleotide::A, Nucleotide::T, Nucleotide::G]).unwrap();
    /// assert_eq!(
    ///     aa,
    ///     AminoAcid::M
    /// );
    /// ```
    pub fn translate(&self, codon: &[Nucleotide; 3]) -> Result<AminoAcid, AtgError> {
        let first_offset = codon[0].as_ncbi_int()? * 16;
        let second_offset = codon[1].as_ncbi_int()? * 4;
        let third_offset = codon[2].as_ncbi_int()?;
        Ok(self.code[first_offset + second_offset + third_offset])
    }

    /// Returns a vector of all codons that code for the provided AminoAcid
    ///
    /// # Examples
    /// ```
    /// use atglib::models::{AminoAcid, GeneticCode, Nucleotide};
    /// let code = GeneticCode::default();
    /// let codons = code.reverse_lookup(&AminoAcid::M);
    /// assert_eq!(
    ///     codons[0],
    ///     [Nucleotide::A, Nucleotide::T, Nucleotide::G]
    /// );
    /// ```
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
    ///
    /// # Examples
    /// ```
    /// use atglib::models::{AminoAcid, GeneticCode, Nucleotide};
    /// let code = GeneticCode::default();
    /// let codons = code.stop_codons();
    /// assert!(codons.contains(&[Nucleotide::T, Nucleotide::A, Nucleotide::A]));
    /// ```
    pub fn stop_codons(&self) -> Vec<[Nucleotide; 3]> {
        self.reverse_lookup(&AminoAcid::Ter)
    }

    /// Returns true if the provided codon is a stop codon
    ///
    /// This method can only handle codons with exaccly 3 nucleotides. All other input will return `false`
    ///
    /// # Examples
    /// ```
    /// use atglib::models::{AminoAcid, GeneticCode, Nucleotide};
    /// let code = GeneticCode::default();
    /// assert!(code.is_stop_codon(&[Nucleotide::T, Nucleotide::A, Nucleotide::A]));
    /// ```
    pub fn is_stop_codon(&self, codon: &[Nucleotide]) -> bool {
        if codon.len() != 3 {
            return false;
        }

        matches!(
            self.translate(&[codon[0], codon[1], codon[2]]),
            Ok(AminoAcid::Ter)
        )
    }


    /// Returns `true` if the codon is a start codon
    ///
    /// This method considers only the canonical `ATG` start codon and does
    /// not include non-standard start codons
    pub fn is_start_codon(codon: &[Nucleotide]) -> bool {
        codon == [Nucleotide::A, Nucleotide::T, Nucleotide::G]
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

    #[test]
    fn test_start_codon() {
        let seq = vec![
            Nucleotide::C,
            Nucleotide::A,
            Nucleotide::T,
            Nucleotide::G,
            Nucleotide::T
        ];
        assert!(GeneticCode::is_start_codon(&[Nucleotide::A, Nucleotide::T, Nucleotide::G]));
        assert!(GeneticCode::is_start_codon(&seq[1..4]));
        assert!(!GeneticCode::is_start_codon(&seq[0..3]));
        assert!(!GeneticCode::is_start_codon(&seq[1..5]));
        assert!(!GeneticCode::is_start_codon(&seq));
    }
}

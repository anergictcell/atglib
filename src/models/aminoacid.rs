use core::fmt::Display;

use crate::utils::errors::AtgError;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
/// Enum of all standard amino acids
///
/// Use the single-letter code to create `AminoAcid` instances, except for the Termination codon, wich uses `Ter`
///
/// # Examples
///
/// ```
/// use atglib::models::AminoAcid;
///
/// // derive from `slice`
/// assert_eq!(AminoAcid::M, AminoAcid::try_from("Met").unwrap());
/// assert_eq!(AminoAcid::M, AminoAcid::try_from("M").unwrap());
///
/// // derive from `char`
/// assert_eq!(AminoAcid::Ter, AminoAcid::try_from('*').unwrap());
///
/// // the `Display` and `as_ref` traits return the three letter code
/// assert_eq!(format!("{}", AminoAcid::A), "Ala".to_string());
///
/// assert_eq!(AminoAcid::A.as_ref(), "Ala".to_string());
///
/// // Get the single-letter code
/// assert_eq!(AminoAcid::A.single_letter(), 'A');
/// ```
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
    /// Returns the single-letter code of the `AminoAcid`
    ///
    /// # Examples
    /// ```
    /// use atglib::models::AminoAcid;
    /// let aa = AminoAcid::M;
    /// assert_eq!(aa.single_letter(), 'M');
    /// ```
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

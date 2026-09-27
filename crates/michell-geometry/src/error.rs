use std::fmt;

/// Errors reported by this crate.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Error {
    /// The B-spline data does not describe a valid clamped tensor-product surface.
    InvalidSpline(String),
    /// The surface is a valid spline but not a valid hull half-breadth function.
    InvalidGeometry(String),
    /// Physically meaningless conditions (non-positive speed, density, ...).
    InvalidConditions(String),
    /// Ill-conditioned or inconsistent input data (e.g. a least-squares fit
    /// with too few samples).
    InvalidInput(String),
    /// A file could not be parsed.
    Parse(String),
    /// The file was parsed but uses features this crate does not support.
    Unsupported(String),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::InvalidSpline(msg) => write!(f, "invalid B-spline: {msg}"),
            Error::InvalidGeometry(msg) => write!(f, "invalid hull geometry: {msg}"),
            Error::InvalidConditions(msg) => write!(f, "invalid conditions: {msg}"),
            Error::InvalidInput(msg) => write!(f, "invalid input: {msg}"),
            Error::Parse(msg) => write!(f, "parse error: {msg}"),
            Error::Unsupported(msg) => write!(f, "unsupported: {msg}"),
        }
    }
}

impl std::error::Error for Error {}

pub type Result<T> = std::result::Result<T, Error>;

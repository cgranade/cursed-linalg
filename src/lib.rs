#![cfg_attr(all(), doc = include_str!("../README.md"))]
#![warn(missing_docs)]
#![warn(rustdoc::missing_doc_code_examples)]

mod decompositions;
pub use decompositions::*;

mod inv;
pub use inv::*;

mod expm;
pub use expm::*;

mod pow;
pub use pow::*;

mod sp_func;
pub use sp_func::*;

pub mod error;

// Don't expose array_ext publicly. These are traits we need internally, but
// aren't stable enough to be part of the public API.
mod array_ext;

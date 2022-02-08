#![cfg_attr(all(), doc = include_str!("../README.md"))]
#![warn(missing_docs)]
#![warn(rustdoc::missing_doc_code_examples)]

mod decompositions;
pub use decompositions::*;

mod inv;
pub use inv::*;

mod expm;
pub use expm::*;

pub mod error;

mod array_ext;

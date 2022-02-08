It's not BLASed, it's cursed!

# Introduction

This crate implements common linear algebra functionality in Rust directly in terms of [`ndarray`] data structures, without a dependency on BLAS or LAPACK. Implementing in pure-Rust makes it easier to target linear algebra applications without requiring any additional shared libraries.

# Status

- This crate is under development, and may be missing essential features, and functions may hit `todo!` or `unimplemented!` panics.
- Performance and numerical stability improvements may be needed for your application.
- Unit, integration, and performance tests are being developed, and may not have full coverage.

# Licensed

[`cursed-linalg`] is licensed under the MIT license.

# Acknowledgments

- Portions of this library were ported from [MathNet.Numerics](https://github.com/mathnet/mathnet-numerics) under the MIT license.
- [`cauchy`] is used to abstract over different representations of complex scalars.
- [`thiserror`] is used to create error enums.
- [`miette`] is used to provide nice diagnostics for errors.

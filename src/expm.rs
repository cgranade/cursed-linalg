use std::{ops::Neg, process::Output, fmt::Error};

use cauchy::Scalar;
use ndarray::{Data, RawData, Array2, ArrayBase, RawArrayView, Ix2, OwnedRepr, ScalarOperand};

use crate::{approximate_factorial, error::LinalgError, Identity, Inv, MatrixPower};

// NB: We use the notation of
//     https://www.cs.jhu.edu/~misha/ReadingSeminar/Papers/Moler03.pdf.
//
//     In particular:
//
//         e^{A} ≈ R_pq(A) ≔ [D_pq(A)]⁻¹ N_pq(A)
//
//     where
//
//         N_pq(A) ≔ Σ_j=0^p [(p + q - j)! p!] / [(p + q)! j! (p - j)!] A^j
//
//     and
//
//         D_pq(A) ≔ Σ_j=0^q [(p + q - j)! p!] / [(p + q)! j! (q - j)!] (-A)^j.

/// Types that support the matrix exponential $e^{A}$.
pub trait Expm {
    /// Errors that can result from the [`expm`] method.
    ///
    /// [`expm`]: Expm::expm
    type Error;
    /// Type used to represent $e^{A}$ when calling [`expm`].
    ///
    /// [`expm`]: Expm::expm
    type Output;

    /// Returns the matrix exponential $e^{A}$ of a matrix $A$.
    fn expm(&self) -> Result<Self::Output, Self::Error>;
}

impl<A, S, E> Expm for ArrayBase<S, Ix2>
where
    S: Data + RawData<Elem = A>,
    A: Scalar + ScalarOperand,
    E: std::convert::From<LinalgError>,
    ArrayBase<S, Ix2>: Inv<Error = E> + Identity,
// TODO: Allow simpler type constraints.
// impl<T, E> Expm for T
// where
//     T: Identity + Inv<Error = E> + MatrixPower<Error = E>,
{
    type Error = E;
    type Output = ArrayBase<OwnedRepr<A>, Ix2>;

    fn expm(&self) -> Result<Self::Output, Self::Error> {
        // TODO: allow generalizing p, q
        let p = 16;
        let q = 16;

        let n = pade_n::<A, S, E>(self, p, q)?;
        let d_inv = pade_d::<A, S, E>(self, p, q)?.inv()?;

        Ok(d_inv.dot(&n))
    }
}

// TODO: Allow simpler type constraints.
// fn pade_d<'a, T: Identity + MatrixPower<Error = E>, E>(mtx: &'a T, p: u32, q: u32) -> Result<T, E>
// where &'a T: Neg<Output = T> {
fn pade_d<A: Scalar + ScalarOperand, S: Data + RawData<Elem = A>, E>(mtx: &ArrayBase<S, Ix2>, p: u32, q: u32) -> Result<ArrayBase<OwnedRepr<A>, Ix2>, E>
where
    ArrayBase<S, Ix2>: Inv<Error = E> + Identity,
    E: std::convert::From<LinalgError>,
{
    // D_pq(A) ≔ Σ_j=0^q [(p + q - j)! p!] / [(p + q)! j! (q - j)!] (-A)^j
    let mut result = mtx.eye_like().to_owned();
    for j in 0..q + 1 {
        let coeff =
            (approximate_factorial::<_, f64>(p + q - j) * approximate_factorial::<_, f64>(p)) /
            (approximate_factorial::<_, f64>(p + q) * approximate_factorial::<_, f64>(j) * approximate_factorial::<_, f64>(q - j));
        let coeff = A::from(coeff).unwrap();
        let mtx_pow = (-mtx.to_owned()).matrix_power(j.into())?;
        let term = mtx_pow * coeff;
        result = result + term;
    }
    Ok(result)
}

fn pade_n<A: Scalar + ScalarOperand, S: Data + RawData<Elem = A>, E>(mtx: &ArrayBase<S, Ix2>, p: u32, q: u32) -> Result<ArrayBase<OwnedRepr<A>, Ix2>, E>
where
    ArrayBase<S, Ix2>: Inv<Error = E> + Identity,
    E: std::convert::From<LinalgError>,
{
    // N_pq(A) ≔ Σ_j=0^p [(p + q - j)! p!] / [(p + q)! j! (p - j)!] A^j
    // [(p + q - 0)! p!] / [(p + q)! (p)!]
    let mut result = mtx.eye_like().to_owned();
    for j in 0..p + 1 {
        let coeff =
            (approximate_factorial::<_, f64>(p + q - j) * approximate_factorial::<_, f64>(p)) /
            (approximate_factorial::<_, f64>(p + q) * approximate_factorial::<_, f64>(j) * approximate_factorial::<_, f64>(p - j));
        let coeff = A::from(coeff).unwrap();
        let mtx_pow = mtx.to_owned().matrix_power(j.into())?;
        let term = mtx_pow * coeff;
        result = result + term;
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use cauchy::c64;
    use num_traits::{Zero, One};
    use ndarray::array;

    use crate::Expm;
    use crate::error::LinalgError;

    #[test]
    fn expm_works_for_rz() -> Result<(), LinalgError> {
        let hz = c64::new(0.0, -1.0) * array![
            [c64::one(), c64::zero()],
            [c64::zero(), -c64::one()]
        ];
        let u = hz.expm();
        // FIXME: returns the wrong answer, so we print while we diagnose.
        println!("{:?}", u);
        panic!("Test currently fails.");
    }
}

use cauchy::Scalar;
use ndarray::{Array1, Array2, Data, OwnedRepr, RawData};

use crate::{LUDecomposable, LUDecomposition};

/// Types that support the matrix inverse $A^{-1}$.
pub trait Inv {
    /// Errors that can result from the [`inv`] method.
    ///
    /// [`inv`]: Inv::inv
    type Error;
    /// Type used to represent $A^{-1}$ when calling [`inv`].
    ///
    /// [`inv`]: Inv::inv
    type Output;

    /// Computes the inverse of a given matrix.
    fn inv(&self) -> Result<Self::Output, Self::Error>;
}

impl<M, A, S, E> Inv for M
where
    M: LUDecomposable<Elem = A, Repr = S, OwnedRepr = OwnedRepr<A>, Error = E>,
    S: Data + RawData<Elem = A>,
    A: Scalar,
    // TODO: Allow decomposables and decompositions to have different errors.
    <M as LUDecomposable>::Output: LUDecomposition<A, OwnedRepr<A>, Error = E>,
{
    type Error = M::Error;
    type Output =
        <<M as LUDecomposable>::Output as LUDecomposition<A, OwnedRepr<A>>>::MatrixSolution;

    fn inv(&self) -> Result<Self::Output, Self::Error> {
        // TODO: Check determinant here.
        let lu = self.lu()?;
        let id_mtx = Array2::from_diag(&Array1::ones(lu.order()));
        lu.solve_matrix(&id_mtx)
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;
    use cauchy::c64;
    use ndarray::array;
    use num_traits::Zero;

    use crate::{error::LinalgError, Inv};

    #[test]
    fn inv_works_f64() -> Result<(), LinalgError> {
        let mtx = array![[6.0, 18.0, 3.0], [2.0, 12.0, 1.0], [4.0, 15.0, 3.0]];
        // TODO: Actually write the test!
        let inv = mtx.inv()?;
        println!("{:?}", inv);
        Ok(())
    }

    #[test]
    fn inv_works_c64() -> Result<(), LinalgError> {
        let mtx = array![
            [c64::zero(), c64::new(0.0, 1.0), c64::new(2.0, 0.0)],
            [c64::new(0.0, 3.0), c64::new(4.0, 0.0), c64::new(0.0, 5.0)],
            [c64::new(6.0, 0.0), c64::new(0.0, 7.0), c64::new(8.0, 0.0)]
        ];
        let inv = mtx.inv()?;
        let expected_inv = array![
            [
                c64::new(-0.69791667, 0.),
                c64::new(0., -0.0625),
                c64::new(0.13541667, 0.)
            ],
            [
                c64::new(0., -0.0625),
                c64::new(0.125, 0.),
                c64::new(0., -0.0625)
            ],
            [
                c64::new(0.46875, 0.),
                c64::new(0., -0.0625),
                c64::new(-0.03125, 0.)
            ],
        ];

        for (actual, expected) in inv.iter().zip(expected_inv.iter()) {
            assert_abs_diff_eq!(actual, expected, epsilon = 1e-6);
        }
        Ok(())
    }
}

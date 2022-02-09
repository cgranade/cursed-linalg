use cauchy::Scalar;
use ndarray::{linalg::Dot, ArrayBase, Data, Ix2, OwnedRepr, RawData, Array2, Array1, Array, Dim};
use num_traits::Zero;

use crate::{error::LinalgError, Inv, LUDecomposable};

// TODO: Move trait into new mod.
// TODO: Consider dropping in favor of ndarray::eye; this trait is mostly
//       useful for associating ndarray::eye with a type so that it can be
//       more easily called from other traits.
pub trait Identity: Sized {
    fn eye(size: usize) -> Self;
    fn eye_like(&self) -> Self;
}

impl<A> Identity for Array<A, Ix2>
where
    A: Scalar + Zero + Clone,
{
    fn eye(size: usize) -> Self {
        Self::eye(size)
    }

    fn eye_like(&self) -> Self {
        <Array2<A> as Identity>::eye(std::cmp::min(
            self.shape()[0],
            self.shape()[1]
        ))
    }
}

/// Types that support matrix powers $A^{n}$.
pub trait MatrixPower: Sized {
    /// Errors that can result from the [`matrix_power`] method.
    ///
    /// [`matrix_power`]: MatrixPower::matrix_power
    type Error;

    fn matrix_power(&self, pow: i64) -> Result<Self, Self::Error>;
}

impl<M> MatrixPower for M
where
    M: Inv<Output = Self> + Dot<Self, Output = Self> + Identity,
{
    type Error = M::Error;

    fn matrix_power(&self, pow: i64) -> Result<Self, Self::Error> {
        if pow < 0 {
            return self.inv()?.matrix_power(-pow);
        }

        let mut result = self.eye_like();

        // TODO: use binary exponentiation; the simple for loop is just to
        //       quickly get something that can be tested.
        for _ in 0..pow {
            result = result.dot(self);
        }

        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;
    use ndarray::{array, Array2};

    use crate::{MatrixPower, error::LinalgError};

    #[test]
    fn pauli_x_squares_to_ident() -> Result<(), LinalgError> {
        let x: Array2<f64> = array![
            [0.0, 1.0],
            [1.0, 0.0]
        ];
        let actual = x.matrix_power(2)?;
        let expected = array![
            [1.0, 0.0],
            [0.0, 1.0]
        ];
        for (actual, expected) in actual.iter().zip(expected.iter()) {
            assert_abs_diff_eq!(actual, expected, epsilon = 1e-6);
        }
        Ok(())
    }
}
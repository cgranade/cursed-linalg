use ndarray::{
    Array, ArrayBase, Axis, Data, DataMut, Dimension, Ix2, OwnedRepr, RawData, RemoveAxis,
};
use num_traits::Zero;

use crate::error::LinalgError;

pub(crate) trait ShapeExt {
    fn require_square(&self) -> Result<usize, LinalgError>;
}

impl<S, D> ShapeExt for ArrayBase<S, D>
where
    S: RawData,
    D: Dimension,
{
    fn require_square(&self) -> Result<usize, LinalgError> {
        let dim = self.raw_dim();
        if dim[0] != dim[1] {
            Err(LinalgError::NotSquare(dim[0], dim[1]))
        } else {
            Ok(dim[0])
        }
    }
}

pub trait MatrixExt
where
    Self: Sized,
{
    type Output;
    fn lower_triangular(&self) -> Self::Output;
    fn upper_triangular(&self) -> Self::Output;
}
impl<S, A> MatrixExt for ArrayBase<S, Ix2>
where
    S: Data<Elem = A>,
    A: Zero + Clone,
{
    type Output = ArrayBase<OwnedRepr<A>, Ix2>;

    /// Returns a copy of the matrix with all elements above the diagonal set to zero.
    fn lower_triangular(&self) -> Self::Output {
        Array::from_shape_fn(self.raw_dim(), |(i, j)| {
            if i >= j {
                self[(i, j)].clone()
            } else {
                A::zero()
            }
        })
    }

    fn upper_triangular(&self) -> Self::Output {
        Array::from_shape_fn(self.raw_dim(), |(i, j)| {
            if i < j {
                self[(i, j)].clone()
            } else {
                A::zero()
            }
        })
    }
}

pub(crate) trait RemoveAxisExt {
    fn swap_index_axis(&mut self, axis: Axis, idx_source: usize, idx_dest: usize);
}
impl<S, D> RemoveAxisExt for ArrayBase<S, D>
where
    S: Data + DataMut,
    <S as RawData>::Elem: Clone,
    D: Dimension + RemoveAxis,
{
    fn swap_index_axis(&mut self, axis: Axis, idx_source: usize, idx_dest: usize) {
        // TODO: Avoid copying both; this is needed currently to avoid mutably
        //       borrowing something which has an outstanding immutable
        //       borrow.
        let source = self.index_axis(axis, idx_source).clone().to_owned();
        let dest = self.index_axis(axis, idx_dest).clone().to_owned();
        self.index_axis_mut(axis, idx_source).assign(&dest);
        self.index_axis_mut(axis, idx_dest).assign(&source);
    }
}

#[cfg(test)]
mod tests {
    use super::MatrixExt;
    use ndarray::{array, Array2};

    #[test]
    fn upper_triangular_works() {
        let mtx: Array2<f64> = array![[6.0, 18.0, 3.0], [2.0, 12.0, 1.0], [4.0, 15.0, 3.0]];
        let upper = mtx.upper_triangular();
        assert_eq!(
            upper,
            array![[0.0, 18.0, 3.0], [0.0, 0.0, 1.0], [0.0, 0.0, 0.0]]
        );
    }

    #[test]
    fn lower_triangular_works() {
        let mtx: Array2<f64> = array![[6.0, 18.0, 3.0], [2.0, 12.0, 1.0], [4.0, 15.0, 3.0]];
        let lower = mtx.lower_triangular();
        assert_eq!(
            lower,
            array![[6.0, 0.0, 0.0], [2.0, 12.0, 0.0], [4.0, 15.0, 3.0]]
        );
    }
}

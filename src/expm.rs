use crate::Inv;

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

impl<T> Expm for T
where
    T: Inv,
{
    type Error = <Self as Inv>::Error;
    type Output = Self;

    fn expm(&self) -> Result<Self::Output, Self::Error> {
        // TODO: allow generalizing p, q.
        todo!()
    }
}

// fn pade_dinv()

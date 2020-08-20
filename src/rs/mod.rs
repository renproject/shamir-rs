use secp256k1::scalar::Scalar;
use std::cmp;

pub mod poly;

use poly::{EEAState, Poly};

pub fn decode<'a, I>(points: I, k: usize) -> Option<(Poly, Option<Vec<Scalar>>)>
where
    I: Iterator<Item = (&'a Scalar, &'a Scalar)> + ExactSizeIterator + Clone,
{
    let indices = points.clone().map(|(index, _)| index);
    let precompute = Precompute::new(indices);
    decode_with_precompute(&precompute, points, k)
}

pub struct Precompute {
    g0: Poly,
    basis: Vec<Poly>,
}

impl Precompute {
    pub fn new<'a, I>(indices: I) -> Self
    where
        I: Iterator<Item = &'a Scalar> + ExactSizeIterator + Clone,
    {
        let g0 = precompute_poly_g0(indices.clone());
        let basis = poly::lagrange_basis(indices);
        Precompute { g0, basis }
    }
}

#[allow(clippy::many_single_char_names)]
pub fn decode_with_precompute<'a, I>(
    precompute: &Precompute,
    points: I,
    k: usize,
) -> Option<(Poly, Option<Vec<Scalar>>)>
where
    I: Iterator<Item = (&'a Scalar, &'a Scalar)> + ExactSizeIterator + Clone,
{
    let n = points.len();
    if k > n {
        panic!("k is too large for the given codeword");
    }
    let values = points.clone().map(|(_, value)| value);

    let g1 = Poly::interpolate_using_basis(precompute.basis.iter().zip(values));
    let mut eea_state = EEAState::new(&precompute.g0, &g1);
    while eea_state.rnext.degree() >= (n + k) / 2 {
        eea_state.step();
    }

    let (g, v) = (&eea_state.rnext, &eea_state.tnext);
    let mut f1 = Poly::with_capacity(cmp::max(g.degree(), v.degree()) + 1);
    let mut r = Poly::with_capacity(v.degree() + 1);
    f1.divide(&mut r, g, v);

    if r.is_zero() && f1.degree() < k {
        let indices = points.map(|(index, _)| index);
        let mut err_indices = Vec::new();
        for index in indices {
            if v.evaluate_at(index).is_zero() {
                err_indices.push(*index);
            }
        }
        Some((
            f1,
            if err_indices.is_empty() {
                None
            } else {
                Some(err_indices)
            },
        ))
    } else {
        None
    }
}

pub fn precompute_poly_g0<'a, I>(indices: I) -> Poly
where
    I: Iterator<Item = &'a Scalar> + ExactSizeIterator,
{
    let mut linear_term = Poly::new_normalised_linear(&Scalar::default());
    let mut g0 = Poly::with_capacity(indices.len());
    g0.set_to_one();

    for index in indices {
        linear_term[0].negate(index);
        g0.mul_assign(&linear_term);
    }

    g0
}

#[cfg(test)]
mod tests {
    use super::*;
    use secp256k1::scalar;

    #[test]
    fn decoding_with_errors() {
        const N: usize = 10;
        const K: usize = 3;
        let mut indices = [Scalar::default(); N];
        scalar::randomise_scalars_using_thread_rng(&mut indices);
        let precompute = Precompute::new(indices.iter());
        let mut poly = Poly::with_capacity(K);
        poly.randomise_using_thread_rng(K - 1);

        let mut codeword = [Scalar::default(); N];
        for (index, value) in indices.iter().zip(&mut codeword) {
            *value = poly.evaluate_at(index);
        }

        // No errors.
        let res = decode_with_precompute(&precompute, indices.iter().zip(codeword.iter()), K);
        assert!(res.is_some());
        let (res_poly, errors) = res.unwrap();
        assert_eq!(res_poly, poly);
        assert!(errors.is_none());

        // The number of errors is small enough that the polynomial can be reconstructed. The
        // indices at which there were errors should also be returned.
        let threshold = (N - K) / 2;
        for i in 1..=threshold {
            codeword[i - 1].randomise_using_thread_rng();
            let res = decode_with_precompute(&precompute, indices.iter().zip(codeword.iter()), K);
            assert!(res.is_some());
            let (res_poly, errors) = res.unwrap();
            assert_eq!(res_poly, poly);
            assert!(errors.is_some());
            assert_eq!(errors.unwrap().as_slice(), &indices[..i]);
        }

        // If there are too many errors, reconstruction of the polynomial is not possible. It is
        // possible that there can be so many errors that the codeword now represents a different
        // polynomial, but since we are randomly changing the elements in the code word, the
        // probability that this would happen is negligible.
        for i in threshold + 1..N {
            codeword[i - 1].randomise_using_thread_rng();
            let res = decode_with_precompute(&precompute, indices.iter().zip(codeword.iter()), K);
            assert!(res.is_none());
        }
    }
}

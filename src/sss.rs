use crate::rs::poly;
use secp256k1::scalar::Scalar;

#[derive(Copy, Clone, Debug, Default, Eq, PartialEq)]
pub struct Share {
    pub index: Scalar,
    pub value: Scalar,
}

impl Share {
    pub fn add(&mut self, a: &Share, b: &Share) {
        if a.index != b.index {
            panic!("cannot add shares with different indices")
        }
        self.index = a.index;
        self.value.add_mut(&a.value, &b.value);
    }

    pub fn add_assign(&mut self, a: &Share) {
        if self.index != a.index {
            panic!("cannot add shares with different indices")
        }
        self.value.add_assign_mut(&a.value);
    }

    pub fn scale(&mut self, share: &Share, scalar: &Scalar) {
        self.index = share.index;
        self.value.mul_mut(&share.value, scalar);
    }

    pub fn scale_assign(&mut self, scalar: &Scalar) {
        self.value.mul_assign_mut(scalar);
    }
}

pub fn share_secret(indices: &[Scalar], secret: &Scalar, k: usize) -> Vec<Share> {
    let mut shares = Vec::with_capacity(indices.len());
    shares.resize(indices.len(), Share::default());
    share_secret_in_place(shares.as_mut_slice(), indices, secret, k);
    shares
}

pub fn share_secret_in_place(
    dst_shares: &mut [Share],
    indices: &[Scalar],
    secret: &Scalar,
    k: usize,
) {
    let n = indices.len();
    if n < k {
        panic!("invalid inputs")
    }
    let mut coeff = Scalar::new_random_using_thread_rng();
    for (i, _) in indices.iter().enumerate() {
        dst_shares[i].value = coeff;
    }
    for _ in (1..k - 1).rev() {
        coeff.randomise_using_thread_rng();
        for (j, index) in indices.iter().enumerate() {
            dst_shares[j].value.mul_assign_mut(index);
            dst_shares[j].value.add_assign_mut(&coeff);
        }
    }
    for (i, index) in indices.iter().enumerate() {
        dst_shares[i].index = *index;
        dst_shares[i].value.mul_assign_mut(index);
        dst_shares[i].value.add_assign_mut(secret);
    }
}

pub fn share_secret_and_get_coeffs_in_place(
    dst_shares: &mut [Share],
    coeff_slice: &mut [Scalar],
    indices: &[Scalar],
    secret: &Scalar,
) {
    let n = indices.len();
    let k = coeff_slice.len();
    if n < k {
        panic!("invalid inputs")
    }
    coeff_slice[0] = *secret;
    coeff_slice
        .iter_mut()
        .skip(1)
        .for_each(|c| c.randomise_using_thread_rng());
    for (i, index) in indices.iter().enumerate() {
        dst_shares[i].index = *index;
        poly::eval_scalar_slice_in_place(&mut dst_shares[i].value, coeff_slice.iter(), index);
    }
}

pub fn interpolate_shares_at_zero<'a, I>(shares: I) -> Scalar
where
    I: Iterator<Item = &'a Share> + Clone,
{
    let mut result = Scalar::default();
    interpolate_shares_at_zero_in_place(&mut result, shares);
    result
}

pub fn interpolate_shares_at_zero_in_place<'a, I>(dst: &mut Scalar, shares: I)
where
    I: Iterator<Item = &'a Share> + Clone,
{
    let mut numerator = Scalar::default();
    let mut denominator = Scalar::default();
    let mut tmp = Scalar::default();
    dst.clear();
    for Share { index: i, value } in shares.clone() {
        eval_lagrange_basis_at_zero_in_place(
            &mut tmp,
            i,
            shares.clone().map(|Share { index, .. }| index),
            &mut numerator,
            &mut denominator,
        );
        tmp.mul_assign_mut(value);
        dst.add_assign_mut(&tmp);
    }
}

pub(crate) fn eval_lagrange_basis_at_zero_in_place<'a, I>(
    eval: &mut Scalar,
    i: &Scalar,
    indices: I,
    numerator: &mut Scalar,
    denominator: &mut Scalar,
) where
    I: Iterator<Item = &'a Scalar>,
{
    numerator.set_u64(1);
    denominator.set_u64(1);
    for j in indices {
        if i == j {
            continue;
        }
        numerator.mul_assign_mut(j);
        eval.sub_mut(j, i);
        denominator.mul_assign_mut(eval)
    }
    eval.divide_mut(numerator, denominator);
}

pub fn shares_are_k_consistent(shares: &[Share], k: usize) -> bool {
    if shares.len() < k {
        panic!("not enough shares for given threshold")
    }
    let secret = interpolate_shares_at_zero(shares[..k].iter());
    shares_are_k_consistent_with_secret(shares, &secret, k)
}

pub fn shares_are_k_consistent_with_secret(shares: &[Share], secret: &Scalar, k: usize) -> bool {
    if shares.len() < k {
        panic!("not enough shares for given threshold")
    }
    
    let mut reconstructed_secret = Scalar::default();
    
    if shares.len() == k {
        interpolate_shares_at_zero_in_place(&mut reconstructed_secret, shares.iter());
        if reconstructed_secret != *secret {
            return false;
        }
    }
    
    for i in 0..(shares.len() - k) {
        interpolate_shares_at_zero_in_place(&mut reconstructed_secret, shares[i..i + k].iter());
        if reconstructed_secret != *secret {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use secp256k1::scalar;

    #[test]
    fn poly_eval_at_zero() {
        let mut coeffs = [Scalar::default(); 10];
        scalar::randomise_scalars_using_thread_rng(&mut coeffs);
        let eval = poly::eval_scalar_slice(&coeffs, &Scalar::zero());
        assert_eq!(eval, coeffs[0]);
    }

    #[test]
    fn poly_eval_at_one() {
        let mut coeffs = [Scalar::default(); 10];
        scalar::randomise_scalars_using_thread_rng(&mut coeffs);
        let eval = poly::eval_scalar_slice(&coeffs, &Scalar::one());
        let actual = coeffs.iter().fold(Scalar::zero(), |mut acc, c| {
            acc.add_assign_mut(c);
            acc
        });
        assert_eq!(eval, actual);
    }

    #[test]
    fn shares_eq_to_poly_eval() {
        const N: usize = 10;
        const K: usize = 5;

        let mut indices = [Scalar::default(); N];
        scalar::randomise_scalars_using_thread_rng(&mut indices);

        let secret = Scalar::new_random_using_thread_rng();
        let mut shares = [Share::default(); N];
        let mut coeffs = [Scalar::default(); K];
        share_secret_and_get_coeffs_in_place(&mut shares, &mut coeffs, &indices, &secret);

        for Share { index, value } in &shares {
            let eval = poly::eval_scalar_slice(&coeffs, index);
            assert_eq!(eval, *value)
        }
    }

    #[test]
    fn sss_share_with_ceoffs_and_open() {
        const N: usize = 10;
        const K: usize = 5;

        let mut indices = [Scalar::default(); N];
        scalar::randomise_scalars_using_thread_rng(&mut indices);

        let secret = Scalar::new_random_using_thread_rng();
        let mut shares = [Share::default(); N];
        let mut coeffs = [Scalar::default(); K];
        share_secret_and_get_coeffs_in_place(&mut shares, &mut coeffs, &indices, &secret);
        let reconstructed = interpolate_shares_at_zero(shares.iter());
        assert_eq!(reconstructed, secret);
    }

    #[test]
    fn sss_share_and_open() {
        let mut indices = [Scalar::default(); 10];
        scalar::randomise_scalars_using_thread_rng(&mut indices);
        let k = 5;

        let secret = Scalar::new_random_using_thread_rng();
        let shares = share_secret(&indices, &secret, k);
        let reconstructed = interpolate_shares_at_zero(shares.iter());
        assert_eq!(reconstructed, secret);
    }

    #[test]
    fn sss_shares_are_k_consistent() {
        let mut indices = [Scalar::default(); 10];
        scalar::randomise_scalars_using_thread_rng(&mut indices);
        let k = 5;

        let secret = Scalar::new_random_using_thread_rng();
        let shares = share_secret(&indices, &secret, k);
        assert!(shares_are_k_consistent(&shares, k));
    }

    #[test]
    fn sss_share_addition() {
        const N: usize = 10;
        let mut indices = [Scalar::default(); N];
        scalar::randomise_scalars_using_thread_rng(&mut indices);
        let k = 5;
        let secret1 = Scalar::new_random_using_thread_rng();
        let secret2 = Scalar::new_random_using_thread_rng();
        let shares1 = share_secret(&indices, &secret1, k);
        let shares2 = share_secret(&indices, &secret2, k);

        let mut sum_shares = [Share::default(); N];
        for i in 0..N {
            sum_shares[i].add(&shares1[i], &shares2[i]);
        }
        let mut sum_secret = Scalar::default();
        sum_secret.add_mut(&secret1, &secret2);
        assert!(shares_are_k_consistent_with_secret(
            &sum_shares,
            &sum_secret,
            k
        ))
    }

    #[test]
    fn sss_share_scaling() {
        const N: usize = 10;
        let mut indices = [Scalar::default(); N];
        scalar::randomise_scalars_using_thread_rng(&mut indices);
        let k = 5;
        let secret = Scalar::new_random_using_thread_rng();
        let scale = Scalar::new_random_using_thread_rng();
        let shares = share_secret(&indices, &secret, k);

        let mut scaled_shares = [Share::default(); N];
        for i in 0..N {
            scaled_shares[i].scale(&shares[i], &scale);
        }
        let mut scaled_secret = Scalar::default();
        scaled_secret.mul_mut(&secret, &scale);
        assert!(shares_are_k_consistent_with_secret(
            &scaled_shares,
            &scaled_secret,
            k
        ))
    }
}

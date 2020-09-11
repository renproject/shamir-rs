use secp256k1::group::Gej;
use secp256k1::scalar::Scalar;
use std::ops::{Deref, DerefMut};

use crate::ped;
use crate::sss::{self, Share};

#[derive(Copy, Clone, Debug, Default, Eq, PartialEq)]
pub struct VShare {
    pub share: Share,
    pub decommitment: Scalar,
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct SharingCommitment(Vec<Gej>);

impl SharingCommitment {
    pub fn with_capacity(n: usize) -> Self {
        SharingCommitment(Vec::with_capacity(n))
    }

    pub fn default_with_len(n: usize) -> Self {
        let mut v = Vec::with_capacity(n);
        v.resize_with(n, Gej::default);
        SharingCommitment(v)
    }

    pub fn new_from_vec(v: Vec<Gej>) -> Self {
        SharingCommitment(v)
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn push(&mut self, com: Gej) {
        self.0.push(com);
    }

    pub fn add_assign_mut(&mut self, other: &SharingCommitment) {
        assert_eq!(self.len(), other.len());
        for (com, other_com) in self.0.iter_mut().zip(other.0.iter()) {
            com.add_assign(other_com);
        }
    }

    pub fn scale_assign_mut(&mut self, scale: &Scalar) {
        self.0.iter_mut().for_each(|c| c.scalar_mul_assign(scale));
    }
}

impl Deref for SharingCommitment {
    type Target = [Gej];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for SharingCommitment {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

pub struct VSharing {
    pub vshares: Vec<VShare>,
    pub commitment: SharingCommitment,
}

impl VShare {
    pub fn add_mut(&mut self, a: &VShare, b: &VShare) {
        self.share.add(&a.share, &b.share);
        self.decommitment.add_mut(&a.decommitment, &b.decommitment);
    }

    pub fn add_assign_mut(&mut self, a: &VShare) {
        self.share.add_assign(&a.share);
        self.decommitment.add_assign_mut(&a.decommitment);
    }

    pub fn scale_mut(&mut self, share: &VShare, scalar: &Scalar) {
        self.share.scale(&share.share, scalar);
        self.decommitment.mul_mut(&share.decommitment, scalar);
    }

    pub fn scale_assign_mut(&mut self, scalar: &Scalar) {
        self.share.scale_assign(scalar);
        self.decommitment.mul_assign_mut(scalar);
    }
}

impl<'a> Into<&'a Share> for &'a VShare {
    fn into(self) -> &'a Share {
        &self.share
    }
}

pub fn vshare_secret(
    h: &Gej,
    indices: &[Scalar],
    secret: &Scalar,
    k: usize,
) -> (Vec<VShare>, SharingCommitment) {
    vshare_secret_and_decommitment(
        h,
        indices,
        secret,
        &Scalar::new_random_using_thread_rng(),
        k,
    )
}

pub fn vshare_secret_and_decommitment(
    h: &Gej,
    indices: &[Scalar],
    secret: &Scalar,
    decommitment: &Scalar,
    k: usize,
) -> (Vec<VShare>, SharingCommitment) {
    let n = indices.len();
    let mut vshares = Vec::with_capacity(n);
    let mut sharing_commitment = Vec::with_capacity(k);
    vshares.resize_with(n, VShare::default);
    sharing_commitment.resize_with(k, Gej::default);
    vshare_secret_and_decommitment_in_place(
        &mut vshares,
        &mut sharing_commitment,
        h,
        indices,
        secret,
        decommitment,
    );
    (vshares, SharingCommitment(sharing_commitment))
}

pub fn vshare_secret_in_place(
    dst_vshares: &mut [VShare],
    dst_sharing_commitment: &mut [Gej],
    h: &Gej,
    indices: &[Scalar],
    secret: &Scalar,
) {
    vshare_secret_and_decommitment_in_place(
        dst_vshares,
        dst_sharing_commitment,
        h,
        indices,
        secret,
        &Scalar::new_random_using_thread_rng(),
    )
}

pub fn vshare_secret_and_decommitment_in_place(
    dst_vshares: &mut [VShare],
    dst_sharing_commitment: &mut [Gej],
    h: &Gej,
    indices: &[Scalar],
    secret: &Scalar,
    decommitment: &Scalar,
) {
    let n = indices.len();
    let k = dst_sharing_commitment.len();
    if n < k || k < 2 {
        panic!("invalid inputs")
    }
    let mut fcoeff = Scalar::new_random_using_thread_rng();
    let mut gcoeff = Scalar::new_random_using_thread_rng();
    ped::ped_commit_in_place(&mut dst_sharing_commitment[k - 1], h, &fcoeff, &gcoeff);
    for (i, _) in indices.iter().enumerate() {
        dst_vshares[i].share.value = fcoeff;
        dst_vshares[i].decommitment = gcoeff;
    }
    for i in (1..k - 1).rev() {
        fcoeff.randomise_using_thread_rng();
        gcoeff.randomise_using_thread_rng();
        ped::ped_commit_in_place(&mut dst_sharing_commitment[i], h, &fcoeff, &gcoeff);
        for (j, index) in indices.iter().enumerate() {
            dst_vshares[j].share.value.mul_assign_mut(index);
            dst_vshares[j].share.value.add_assign_mut(&fcoeff);
            dst_vshares[j].decommitment.mul_assign_mut(index);
            dst_vshares[j].decommitment.add_assign_mut(&gcoeff);
        }
    }
    ped::ped_commit_in_place(&mut dst_sharing_commitment[0], h, secret, decommitment);
    for (i, index) in indices.iter().enumerate() {
        dst_vshares[i].share.index = *index;
        dst_vshares[i].share.value.mul_assign_mut(index);
        dst_vshares[i].share.value.add_assign_mut(secret);
        dst_vshares[i].decommitment.mul_assign_mut(index);
        dst_vshares[i].decommitment.add_assign_mut(decommitment);
    }
}

pub fn vshare_is_valid(vshare: &VShare, sharing_commitment: &[Gej], h: &Gej) -> bool {
    let expected = poly_eval_gej_slice_in_exponent(sharing_commitment, &vshare.share.index);
    let actual = ped::ped_commit(h, &vshare.share.value, &vshare.decommitment);
    expected == actual
}

pub fn interpolate_shares_at_zero<'a, I>(vshares: I) -> (Scalar, Scalar)
where
    I: Iterator<Item = &'a VShare> + Clone,
{
    let mut secret = Scalar::default();
    let mut decommitment = Scalar::default();
    interpolate_shares_at_zero_in_place(&mut secret, &mut decommitment, vshares);
    (secret, decommitment)
}

pub fn interpolate_shares_at_zero_in_place<'a, I>(
    secret_dst: &mut Scalar,
    decommitment_dst: &mut Scalar,
    shares: I,
) where
    I: Iterator<Item = &'a VShare> + Clone,
{
    let mut numerator = Scalar::default();
    let mut denominator = Scalar::default();
    let mut tmp1 = Scalar::default();
    let mut tmp2 = Scalar::default();
    secret_dst.clear();
    decommitment_dst.clear();
    for VShare {
        share: Share { index: i, value },
        decommitment,
    } in shares.clone()
    {
        sss::eval_lagrange_basis_at_zero_in_place(
            &mut tmp1,
            i,
            shares.clone().map(|vs| &vs.share.index),
            &mut numerator,
            &mut denominator,
        );
        tmp2.mul_mut(&tmp1, value);
        secret_dst.add_assign_mut(&tmp2);
        tmp2.mul_mut(&tmp1, decommitment);
        decommitment_dst.add_assign_mut(&tmp2);
    }
}

pub fn poly_eval_gej_slice_in_exponent(coeffs: &[Gej], point: &Scalar) -> Gej {
    let mut eval = Gej::default();
    poly_eval_gej_slice_in_exponent_in_place(&mut eval, coeffs, point);
    eval
}

pub fn poly_eval_gej_slice_in_exponent_in_place(dst: &mut Gej, coeffs: &[Gej], point: &Scalar) {
    *dst = coeffs.last().copied().unwrap_or_default();
    coeffs.iter().rev().skip(1).for_each(|c| {
        dst.scalar_mul_assign(point);
        dst.add_assign(c);
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    use secp256k1::{
        group::Gej,
        scalar::{self, Scalar},
    };

    #[test]
    fn verifiable_sharing_has_valid_shares() {
        const N: usize = 10;
        const K: usize = 5;
        let h = Gej::new_random_using_thread_rng();
        let mut indices = [Scalar::default(); N];
        scalar::randomise_scalars_using_thread_rng(&mut indices);
        let mut shares = [VShare::default(); N];
        let mut sharing_commitment = [Gej::default(); K];

        let secret = Scalar::new_random_using_thread_rng();
        vshare_secret_in_place(&mut shares, &mut sharing_commitment, &h, &indices, &secret);
        for share in &shares {
            assert!(vshare_is_valid(share, &sharing_commitment, &h))
        }
    }

    #[test]
    fn shares_with_modified_values_are_invalid() {
        const N: usize = 10;
        const K: usize = 5;
        let h = Gej::new_random_using_thread_rng();
        let mut indices = [Scalar::default(); N];
        scalar::randomise_scalars_using_thread_rng(&mut indices);
        let mut shares = [VShare::default(); N];
        let mut sharing_commitment = [Gej::default(); K];

        let secret = Scalar::new_random_using_thread_rng();
        vshare_secret_in_place(&mut shares, &mut sharing_commitment, &h, &indices, &secret);
        for share in &mut shares {
            share.share.value.randomise_using_thread_rng();
            assert!(!vshare_is_valid(share, &sharing_commitment, &h))
        }
    }

    #[test]
    fn shares_with_modified_indices_are_invalid() {
        const N: usize = 10;
        const K: usize = 5;
        let h = Gej::new_random_using_thread_rng();
        let mut indices = [Scalar::default(); N];
        scalar::randomise_scalars_using_thread_rng(&mut indices);
        let mut shares = [VShare::default(); N];
        let mut sharing_commitment = [Gej::default(); K];

        let secret = Scalar::new_random_using_thread_rng();
        vshare_secret_in_place(&mut shares, &mut sharing_commitment, &h, &indices, &secret);
        for share in &mut shares {
            share.share.index.randomise_using_thread_rng();
            assert!(!vshare_is_valid(share, &sharing_commitment, &h))
        }
    }

    #[test]
    fn shares_with_modified_decommitments_are_invalid() {
        const N: usize = 10;
        const K: usize = 5;
        let h = Gej::new_random_using_thread_rng();
        let mut indices = [Scalar::default(); N];
        scalar::randomise_scalars_using_thread_rng(&mut indices);
        let mut shares = [VShare::default(); N];
        let mut sharing_commitment = [Gej::default(); K];

        let secret = Scalar::new_random_using_thread_rng();
        vshare_secret_in_place(&mut shares, &mut sharing_commitment, &h, &indices, &secret);
        for share in &mut shares {
            share.decommitment.randomise_using_thread_rng();
            assert!(!vshare_is_valid(share, &sharing_commitment, &h))
        }
    }
}

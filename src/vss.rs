use crate::ped::{self, PedCommitment};
use crate::sss::{self, Share};
use secp256k1::group::Gej;
use secp256k1::scalar::Scalar;

#[derive(Copy, Clone, Default)]
pub struct VShare {
    pub share: Share,
    pub decommitment: Scalar,
}

pub type SharingCommitment = Vec<PedCommitment>;

impl VShare {
    pub fn add(&mut self, a: &VShare, b: &VShare) {
        self.share.add(&a.share, &b.share);
        self.decommitment.add(&a.decommitment, &b.decommitment);
    }

    pub fn add_assign(&mut self, a: &VShare) {
        self.share.add_assign(&a.share);
        self.decommitment.add_assign(&a.decommitment);
    }

    pub fn scale(&mut self, share: &VShare, scalar: &Scalar) {
        self.share.scale(&share.share, scalar);
        self.decommitment.mul(&share.decommitment, scalar);
    }

    pub fn scale_assign(&mut self, scalar: &Scalar) {
        self.share.scale_assign(scalar);
        self.decommitment.mul_assign(scalar);
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
    let n = indices.len();
    let mut vshares = Vec::with_capacity(n);
    let mut sharing_commitment = Vec::with_capacity(k);
    vshares.resize_with(n, VShare::default);
    sharing_commitment.resize_with(k, PedCommitment::default);
    vshare_secret_in_place(&mut vshares, &mut sharing_commitment, h, indices, secret);
    (vshares, sharing_commitment)
}

pub fn vshare_secret_in_place(
    dst_vshares: &mut [VShare],
    dst_sharing_commitment: &mut [PedCommitment],
    h: &Gej,
    indices: &[Scalar],
    secret: &Scalar,
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
            dst_vshares[j].share.value.mul_assign(index);
            dst_vshares[j].share.value.add_assign(&fcoeff);
            dst_vshares[j].decommitment.mul_assign(index);
            dst_vshares[j].decommitment.add_assign(&gcoeff);
        }
    }
    gcoeff.randomise_using_thread_rng();
    ped::ped_commit_in_place(&mut dst_sharing_commitment[0], h, secret, &gcoeff);
    for (i, index) in indices.iter().enumerate() {
        dst_vshares[i].share.index = *index;
        dst_vshares[i].share.value.mul_assign(index);
        dst_vshares[i].share.value.add_assign(secret);
        dst_vshares[i].decommitment.mul_assign(index);
        dst_vshares[i].decommitment.add_assign(&gcoeff);
    }
}

pub fn vshare_is_valid(vshare: &VShare, sharing_commitment: &[PedCommitment], h: &Gej) -> bool {
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
        tmp2.mul(&tmp1, value);
        secret_dst.add_assign(&tmp2);
        tmp2.mul(&tmp1, decommitment);
        decommitment_dst.add_assign(&tmp2);
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
        let mut sharing_commitment = [PedCommitment::default(); K];

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
        let mut sharing_commitment = [PedCommitment::default(); K];

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
        let mut sharing_commitment = [PedCommitment::default(); K];

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
        let mut sharing_commitment = [PedCommitment::default(); K];

        let secret = Scalar::new_random_using_thread_rng();
        vshare_secret_in_place(&mut shares, &mut sharing_commitment, &h, &indices, &secret);
        for share in &mut shares {
            share.decommitment.randomise_using_thread_rng();
            assert!(!vshare_is_valid(share, &sharing_commitment, &h))
        }
    }
}

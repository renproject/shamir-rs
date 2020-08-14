#![feature(test)]

extern crate test;

use secp256k1::{
    group::Gej,
    scalar::{self, Scalar},
};
use shamir::{
    ped::PedCommitment,
    sss,
    vss::{self, VShare},
};
use test::Bencher;

#[bench]
fn bench_vshare_secret(b: &mut Bencher) {
    const N: usize = 100;
    const K: usize = 33;
    let h = Gej::new_random_using_thread_rng();

    let mut indices = [Scalar::default(); N];
    scalar::randomise_scalars_using_thread_rng(&mut indices);
    let mut shares = [VShare::default(); N];
    let mut sharing_commitment = [PedCommitment::default(); K];
    let secret = Scalar::new_random_using_thread_rng();
    b.iter(|| {
        vss::vshare_secret_in_place(&mut shares, &mut sharing_commitment, &h, &indices, &secret)
    });
}

#[bench]
fn bench_reconstruct_vshared_secret(b: &mut Bencher) {
    const N: usize = 100;
    const K: usize = 33;
    let h = Gej::new_random_using_thread_rng();

    let mut indices = [Scalar::default(); N];
    scalar::randomise_scalars_using_thread_rng(&mut indices);
    let mut shares = [VShare::default(); N];
    let mut sharing_commitment = [PedCommitment::default(); K];
    let secret = Scalar::new_random_using_thread_rng();
    vss::vshare_secret_in_place(&mut shares, &mut sharing_commitment, &h, &indices, &secret);
    let mut reconstructed = Scalar::default();
    b.iter(|| sss::interpolate_shares_at_zero_in_place(&mut reconstructed, &shares[..K]));
}

#[bench]
fn bench_verify_share(b: &mut Bencher) {
    const N: usize = 100;
    const K: usize = 33;
    let h = Gej::new_random_using_thread_rng();

    let mut indices = [Scalar::default(); N];
    scalar::randomise_scalars_using_thread_rng(&mut indices);
    let mut shares = [VShare::default(); N];
    let mut sharing_commitment = [PedCommitment::default(); K];
    let secret = Scalar::new_random_using_thread_rng();
    vss::vshare_secret_in_place(&mut shares, &mut sharing_commitment, &h, &indices, &secret);
    b.iter(|| vss::vshare_is_valid(&shares[0], &sharing_commitment, &h));
}

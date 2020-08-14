#![feature(test)]

extern crate test;

use secp256k1::scalar::{self, Scalar};
use shamir::sss::{self, Share};
use test::Bencher;

#[bench]
fn bench_share_secret(b: &mut Bencher) {
    const N: usize = 100;
    let k = 33;

    let mut indices = [Scalar::default(); N];
    scalar::randomise_scalars_using_thread_rng(&mut indices);
    let mut shares = [Share::default(); N];
    let secret = Scalar::new_random_using_thread_rng();
    b.iter(|| sss::share_secret_in_place(&mut shares, &indices, &secret, k));
}

#[bench]
fn bench_reconstruct_secret(b: &mut Bencher) {
    const N: usize = 100;
    let k = 33;

    let mut indices = [Scalar::default(); N];
    scalar::randomise_scalars_using_thread_rng(&mut indices);
    let mut shares = [Share::default(); N];
    let secret = Scalar::new_random_using_thread_rng();
    sss::share_secret_in_place(&mut shares, &indices, &secret, k);
    let mut reconstructed = Scalar::default();
    b.iter(|| sss::interpolate_shares_at_zero_in_place(&mut reconstructed, &shares[..k]));
}

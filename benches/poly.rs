#![feature(test)]

extern crate test;

use secp256k1::scalar::{self, Scalar};
use shamir::rs::poly::{self, Poly};
use test::Bencher;

#[bench]
fn bench_poly_add_same_degree(b: &mut Bencher) {
    let mut poly1 = Poly::with_capacity(100);
    let mut poly2 = Poly::with_capacity(100);
    let mut poly_sum = Poly::with_capacity(100);

    poly1.randomise_using_thread_rng(99);
    poly2.randomise_using_thread_rng(99);
    b.iter(|| poly_sum.add(&poly1, &poly2));
}

#[bench]
fn bench_poly_mul_same_degree(b: &mut Bencher) {
    let mut poly1 = Poly::with_capacity(50);
    let mut poly2 = Poly::with_capacity(50);
    let mut poly_prod = Poly::with_capacity(100);

    poly1.randomise_using_thread_rng(49);
    poly2.randomise_using_thread_rng(49);
    b.iter(|| poly_prod.mul(&poly1, &poly2));
}

#[bench]
fn bench_poly_interpolate(b: &mut Bencher) {
    let mut indices = [Scalar::default(); 100];
    let mut values = [Scalar::default(); 100];
    scalar::randomise_scalars_using_thread_rng(&mut indices);
    scalar::randomise_scalars_using_thread_rng(&mut values);
    let basis = poly::lagrange_basis(indices.iter());
    let mut interp = Poly::with_capacity(100);

    b.iter(|| interp.interpolate_using_basis_in_place(basis.iter().zip(values.iter())));
}

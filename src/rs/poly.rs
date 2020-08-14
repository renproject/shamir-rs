use secp256k1::scalar::Scalar;
use std::cmp;

pub struct Poly {
    coeffs: Vec<Scalar>,
}

impl Poly {
    pub fn with_capacity(n: usize) -> Self {
        let mut poly = Poly {
            coeffs: Vec::with_capacity(n),
        };
        poly.coeffs.push(Scalar::new_zero());
        poly
    }

    pub fn randomise_using_thread_rng(&mut self, degree: usize) {
        self.coeffs.clear();
        self.coeffs
            .resize_with(degree + 1, Scalar::new_random_using_thread_rng)
    }

    pub fn set(&mut self, a: &Poly) {
        self.coeffs.resize_with(a.coeffs.len(), Scalar::default);
        self.coeffs.copy_from_slice(&a.coeffs[..]);
    }

    pub fn degree(&self) -> usize {
        self.coeffs.len() - 1
    }

    pub fn leading_coefficient(&self) -> &Scalar {
        self.coeffs.last().unwrap()
    }

    pub fn is_zero(&self) -> bool {
        self.coeffs.len() == 1 && self.coeffs[0].is_zero()
    }

    pub fn zero(&mut self) {
        self.coeffs.clear();
        self.coeffs.push(Scalar::new_zero());
    }

    pub fn interpolate(points: &[(Scalar, Scalar)]) -> Self {
        let mut interp = Poly::with_capacity(points.len());
        let mut basis_poly = Poly::with_capacity(points.len());

        for (i, y) in points {
            lagrange_basis_poly_in_place(&mut basis_poly, points.iter().map(|(x, _)| x), i);
            interp.add_scaled_assign(&basis_poly, y);
        }
        interp
    }

    pub fn interpolate_using_basis<'a>(
        basis_and_values: impl Iterator<Item = (&'a Poly, &'a Scalar)> + ExactSizeIterator,
    ) -> Self {
        let mut interp = Poly::with_capacity(basis_and_values.len());
        interp.interpolate_using_basis_in_place(basis_and_values);
        interp
    }

    pub fn interpolate_using_basis_in_place<'a>(
        &mut self,
        basis_and_values: impl Iterator<Item = (&'a Poly, &'a Scalar)> + ExactSizeIterator,
    ) {
        self.zero();
        self.coeffs.reserve_exact(basis_and_values.len() - 1);
        basis_and_values.for_each(|(b, v)| self.add_scaled_assign(b, v))
    }

    pub fn evaluate_at(&self, point: &Scalar) -> Scalar {
        let mut eval = self.coeffs.last().copied().unwrap_or_default();
        for c in self.coeffs.iter().rev().skip(1) {
            eval.mul_assign(point);
            eval.add_assign(c);
        }
        eval
    }

    pub fn scale(&mut self, a: &Poly, scale: &Scalar) {
        self.coeffs.resize_with(a.coeffs.len(), Scalar::default);
        self.coeffs
            .iter_mut()
            .enumerate()
            .for_each(|(i, c)| c.mul(&a.coeffs[i], scale));
    }

    pub fn scale_assign(&mut self, scale: &Scalar) {
        self.coeffs.iter_mut().for_each(|c| c.mul_assign(scale));
    }

    pub fn add(&mut self, a: &Poly, b: &Poly) {
        let (shorter, longer, slen, llen) = if a.degree() < b.degree() {
            (a, b, a.coeffs.len(), b.coeffs.len())
        } else {
            (b, a, b.coeffs.len(), a.coeffs.len())
        };
        self.coeffs.clear();
        self.coeffs.resize_with(llen, Scalar::default);
        for i in 0..llen {
            self.coeffs[i].add(&shorter.coeffs[i], &longer.coeffs[i]);
        }
        self.coeffs[slen..].copy_from_slice(&longer.coeffs[slen..]);
    }

    pub fn add_assign(&mut self, a: &Poly) {
        let min_len = cmp::min(self.coeffs.len(), a.coeffs.len());
        for i in 0..min_len {
            self.coeffs[i].add_assign(&a.coeffs[i]);
        }
        if self.coeffs.len() < a.coeffs.len() {
            let a_tail = &a.coeffs[self.coeffs.len()..];
            self.coeffs.extend_from_slice(a_tail);
        }
    }

    pub fn add_scaled_assign(&mut self, a: &Poly, scale: &Scalar) {
        let mut prod = Scalar::default();
        self.coeffs.resize_with(a.coeffs.len(), Scalar::default);
        for i in 0..a.coeffs.len() {
            prod.mul(&a.coeffs[i], scale);
            self.coeffs[i].add_assign(&prod);
        }
    }

    pub fn mul(&mut self, a: &Poly, b: &Poly) {
        let prod_degree = a.degree() + b.degree();
        self.coeffs.clear();
        self.coeffs.resize_with(prod_degree + 1, Scalar::default);
        let a_len = a.coeffs.len();
        let b_len = b.coeffs.len();
        let mut prod = Scalar::default();
        for i in 0..a_len {
            for j in 0..b_len {
                prod.mul(&a.coeffs[i], &b.coeffs[j]);
                self.coeffs[i + j].add_assign(&prod);
            }
        }
    }

    pub fn mul_assign(&mut self, a: &Poly) {
        let prod_degree = self.degree() + a.degree();
        let self_len = self.coeffs.len() as isize;
        let a_len = a.coeffs.len() as isize;
        self.coeffs.resize_with(prod_degree + 1, Scalar::default);
        let mut prod = Scalar::default();
        for i in (0..prod_degree as isize + 1).rev() {
            let (a_lower, a_upper) = (cmp::max(0, i + 1 - self_len), cmp::min(a_len, i + 1));
            let (self_lower, self_upper) = (cmp::max(0, i + 1 - a_len), cmp::min(self_len, i + 1));
            assert!(a_upper - a_lower == self_upper - self_lower);
            prod.mul(
                &a.coeffs[a_lower as usize],
                &self.coeffs[(self_upper - 1) as usize],
            );
            self.coeffs[i as usize] = prod;
            for j in 1..(a_upper - a_lower) {
                prod.mul(
                    &a.coeffs[(a_lower + j) as usize],
                    &self.coeffs[(self_upper - 1 - j) as usize],
                );
                self.coeffs[i as usize].add_assign(&prod);
            }
        }
    }

    pub fn divide(&mut self, r: &mut Poly, a: &Poly, b: &Poly) {
        if b.is_zero() {
            panic!("cannot divide by zero")
        }

        self.coeffs.clear();
        r.set(a);
        if a.degree() < b.degree() {
            self.coeffs.push(Scalar::new_zero());
            return;
        }
        self.coeffs
            .resize_with(a.degree() - b.degree() + 1, Scalar::default);
        let mut t = Scalar::default();
        let mut prod = Scalar::default();
        let mut b_coeff_inv = Scalar::default();
        b_coeff_inv.inverse(&b.coeffs.last().unwrap());

        while !r.is_zero() && r.degree() >= b.degree() {
            let degree_diff = r.degree() - b.degree();
            t.mul(r.leading_coefficient(), &b_coeff_inv);
            self.coeffs[degree_diff] = t;
            t.negate_assign();
            for (i, c) in b.coeffs.iter().enumerate() {
                prod.mul(&t, &c);
                r.coeffs[degree_diff + i].add_assign(&prod);
            }
            let _ = r.coeffs.pop();
        }
    }
}

impl PartialEq for Poly {
    fn eq(&self, other: &Self) -> bool {
        if self.coeffs.len() != other.coeffs.len() {
            return false;
        }
        for (self_coeff, other_coeff) in self.coeffs.iter().zip(other.coeffs.iter()) {
            if self_coeff != other_coeff {
                return false;
            }
        }
        true
    }
}

pub fn lagrange_basis<'a, I>(indices: I) -> Vec<Poly>
where
    I: Iterator<Item = &'a Scalar> + ExactSizeIterator + Clone,
{
    let l = indices.len();
    let mut basis = Vec::with_capacity(l);
    basis.resize_with(l, || Poly::with_capacity(l));

    for (i, basis_poly) in indices.clone().zip(basis.iter_mut()) {
        lagrange_basis_poly_in_place(basis_poly, indices.clone(), i);
    }
    basis
}

pub fn lagrange_basis_poly_in_place<'a>(
    dst: &mut Poly,
    indices: impl Iterator<Item = &'a Scalar>,
    index: &Scalar,
) {
    dst.coeffs.clear();
    dst.coeffs.push(Scalar::new_one());
    let mut prod_term = Poly::with_capacity(2);
    prod_term.coeffs.push(Scalar::new_one());
    let mut tmp = Scalar::default();
    let mut denominator = Scalar::new_one();
    for i in indices {
        if i == index {
            continue;
        }
        prod_term.coeffs[0].negate(i);
        dst.mul_assign(&prod_term);
        tmp.sub(index, i);
        denominator.mul_assign(&tmp)
    }
    denominator.inverse_assign();
    dst.scale_assign(&denominator);
}

#[cfg(test)]
mod tests {
    use super::*;
    use secp256k1::scalar;

    #[test]
    fn multiplication_and_division_are_inverses() {
        let mut a = Poly::with_capacity(10);
        let mut b = Poly::with_capacity(10);
        let mut q = Poly::with_capacity(10);
        let mut r = Poly::with_capacity(10);
        let mut prod = Poly::with_capacity(20);
        let mut reconstructed = Poly::with_capacity(20);

        a.randomise_using_thread_rng(9);
        b.randomise_using_thread_rng(9);

        prod.mul(&a, &b);
        q.divide(&mut r, &prod, &b);
        reconstructed.mul(&q, &b);
        reconstructed.add_assign(&r);
        assert!(reconstructed == prod);
    }

    #[test]
    fn mul_equals_mul_assign() {
        let alen = 7;
        let blen = 10;
        let mut a = Poly::with_capacity(alen);
        let mut b = Poly::with_capacity(blen);
        let mut prod = Poly::with_capacity(alen + blen);
        let mut prod_assign = Poly::with_capacity(alen + blen);

        a.randomise_using_thread_rng(alen - 1);
        b.randomise_using_thread_rng(blen - 1);
        prod_assign.set(&a);

        prod.mul(&a, &b);
        prod_assign.mul_assign(&b);
        assert!(prod == prod_assign);
    }

    #[test]
    fn interpolation_matches_given_points() {
        let mut points = [(Scalar::default(), Scalar::default()); 10];
        for (x, y) in &mut points {
            x.randomise_using_thread_rng();
            y.randomise_using_thread_rng();
        }
        let interp = Poly::interpolate(&points);
        for (x, y) in &points {
            assert!(interp.evaluate_at(x) == *y);
        }
    }

    #[test]
    fn lagrange_basis_is_correct() {
        let mut indices = [Scalar::default(); 10];
        scalar::randomise_scalars_using_thread_rng(&mut indices);
        let basis = lagrange_basis(indices.iter());
        for (i, basis_poly) in basis.iter().enumerate() {
            for (j, index) in indices.iter().enumerate() {
                if i == j {
                    assert!(basis_poly.evaluate_at(index).is_one());
                } else {
                    assert!(basis_poly.evaluate_at(index).is_zero());
                }
            }
        }
    }

    #[test]
    fn interpolation_matches_given_points_using_basis() {
        let mut points = [(Scalar::default(), Scalar::default()); 10];
        for (x, y) in &mut points {
            x.randomise_using_thread_rng();
            y.randomise_using_thread_rng();
        }
        let basis = lagrange_basis(points.iter().map(|(x, _)| x));
        let interp = Poly::interpolate_using_basis(basis.iter().zip(points.iter().map(|(_, y)| y)));
        for (x, y) in &points {
            assert!(interp.evaluate_at(x) == *y);
        }
    }
}

use secp256k1::group::Gej;
use secp256k1::scalar::Scalar;

pub fn ped_commit(h: &Gej, value: &Scalar, decommitment: &Scalar) -> Gej {
    let mut commitment = Gej::default();
    ped_commit_in_place(&mut commitment, h, value, decommitment);
    commitment
}

pub fn ped_commit_in_place(dst: &mut Gej, h: &Gej, value: &Scalar, decommitment: &Scalar) {
    let mut hpow = Gej::default();
    hpow.scalar_mul(h, decommitment);
    dst.scalar_base_mul(value);
    dst.add_assign(&hpow);
}

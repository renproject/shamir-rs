use secp256k1::group::Gej;
use secp256k1::scalar::Scalar;

pub type PedCommitment = Gej;

pub fn ped_commit(h: &Gej, value: &Scalar, decommitment: &Scalar) -> PedCommitment {
    let mut commitment = PedCommitment::default();
    ped_commit_in_place(&mut commitment, h, value, decommitment);
    commitment
}

pub fn ped_commit_in_place(
    dst: &mut PedCommitment,
    h: &Gej,
    value: &Scalar,
    decommitment: &Scalar,
) {
    let mut hpow = Gej::default();
    hpow.scalar_mul(h, decommitment);
    dst.scalar_base_mul(value);
    dst.add_assign(&hpow);
}

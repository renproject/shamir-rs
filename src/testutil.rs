use crate::vss::SharingCommitment;
use secp256k1::group::Gej;

pub fn random_commitment_using_thread_rng(k: usize) -> SharingCommitment {
    let mut commitment = SharingCommitment::with_capacity(k);
    for _ in 0..k {
        commitment.push(Gej::new_random_using_thread_rng());
    }
    commitment
}

pub fn zero_commitment(k: usize) -> SharingCommitment {
    let mut commitment = SharingCommitment::with_capacity(k);
    for _ in 0..k {
        commitment.push(Gej::infinity());
    }
    commitment
}

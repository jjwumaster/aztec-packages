use dep::std::hash::{pedersen_with_separator, sha256};
use crate::constants_gen::GENERATOR_INDEX__SIGNATURE_PAYLOAD;

fn bytes_to_sha256_in_a_field<N>(bytes_to_hash: [u8; N]) -> Field {
    let sha256_hashed = sha256(bytes_to_hash);

    // Convert it to a field element
    let mut v = 1;
    let mut high = 0 as Field;
    let mut low = 0 as Field;

    for i in 0..16 {
        high = high + (sha256_hashed[15 - i] as Field) * v;
        low = low + (sha256_hashed[16 + 15 - i] as Field) * v;
        v = v * 256;
    }

    // Abuse that a % p + b % p = (a + b) % p and that low < p
    let hash_in_a_field = low + high * v;

    hash_in_a_field
}

fn compute_message_hash<N>(args: [Field; N]) -> Field {
  // @todo @lherskind We should probably use a separate generator for this, 
  //                  to avoid any potential collisions with payloads.
  pedersen_with_separator(args, GENERATOR_INDEX__SIGNATURE_PAYLOAD)[0]
}
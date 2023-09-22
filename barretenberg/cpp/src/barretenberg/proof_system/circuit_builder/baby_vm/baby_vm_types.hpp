#pragma once

namespace proof_system::baby_vm {

template <typename FF> struct VMOperation {
    bool add = false;
    bool mul = false;
    bool eq = false;
    bool reset = false;
    FF base_point = 0;
    FF mul_scalar_full = 0;
};

} // namespace proof_system::baby_vm
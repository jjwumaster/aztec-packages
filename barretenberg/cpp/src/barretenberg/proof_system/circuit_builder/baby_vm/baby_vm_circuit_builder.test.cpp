#include "baby_vm_circuit_builder.hpp"
#include <gtest/gtest.h>

using namespace barretenberg;

namespace {
auto& engine = numeric::random::get_debug_engine();
}

namespace proof_system::baby_vm_circuit_builder_tests {

TEST(BabyVMCircuitBuilderTests, BaseCase)
{
    using Flavor = ;
    using FF = typename Flavor::FF;

    proof_system::ECCVMCircuitBuilder<Flavor> trace;
    FF a = Fr::random_element();
    FF b = Fr::random_element();
    Fr x = Fr::random_element();
    Fr y = Fr::random_element();

    typename G1::element expected_1 = a + b*x ;

    trace.add_accumulate(a);
    trace.mul_accumulate(b, x);
    trace.eq_and_reset(expected_1);

    bool result = circuit.check_circuit();
    EXPECT_EQ(result, true);
}

} // namespace proof_system::baby_vm_circuit_builder_tests
#pragma once

#include "./baby_vm_types.hpp"
#include "barretenberg/ecc/curves/bn254/fr.hpp"
#include "barretenberg/ecc/curves/grumpkin/grumpkin.hpp"
#include "barretenberg/honk/flavor/ecc_vm.hpp"
#include "barretenberg/honk/proof_system/lookup_library.hpp"
#include "barretenberg/honk/proof_system/permutation_library.hpp"
#include "barretenberg/proof_system/relations/relation_parameters.hpp"

namespace proof_system {

template <typename Flavor> class BabyVMCircuitBuilder {
  public:
    using FF = typename Flavor::FF;

    static constexpr size_t NUM_POLYNOMIALS = Flavor::NUM_ALL_ENTITIES;
    static constexpr size_t NUM_WIRES = Flavor::NUM_WIRES;

    using VMOperation = baby_vm::VMOperation<FF>;
    std::vector<VMOperation> vm_operations;
    using Polynomial = barretenberg::Polynomial<FF>;
    using RawPolynomials = typename Flavor::RawPolynomials;

    void add_accumulate(const FF& to_add)
    {
        vm_operations.emplace_back(VMOperation{
            .add = true,
            .mul = false,
            .eq = false,
            .reset = false,
            .base_point = to_add,
        });
    }

    void mul_accumulate(const FF& to_mul, const FF& scalar)
    {
        vm_operations.emplace_back(VMOperation{
            .add = false,
            .mul = true,
            .eq = false,
            .reset = false,
            .base_point = to_mul,
            .mul_scalar_full = scalar,
        });
    }

    void eq_and_reset(const FF& expected)
    {
        vm_operations.emplace_back(VMOperation{
            .add = false,
            .mul = false,
            .eq = true,
            .reset = true,
            .base_point = expected,
            .mul_scalar_full = 0,
        });
    }

    RawPolynomials compute_full_polynomials()
    {
        RawPolynomials polynomials;
        const auto transcript_state =
            ECCVMTranscriptBuilder<Flavor>::compute_transcript_state(vm_operations, get_number_of_muls());

        for (size_t i = 0; i < transcript_state.size(); ++i) {
            rows.transcript_accumulator_empty[i] = transcript_state[i].accumulator_empty;
            rows.transcript_add[i] = transcript_state[i].q_add;
            rows.transcript_mul[i] = transcript_state[i].q_mul;
            rows.transcript_eq[i] = transcript_state[i].q_eq;
            rows.transcript_reset_accumulator[i] = transcript_state[i].q_reset_accumulator;
            rows.transcript_msm_transition[i] = transcript_state[i].msm_transition;
            rows.transcript_pc[i] = transcript_state[i].pc;
            rows.transcript_msm_count[i] = transcript_state[i].msm_count;
            rows.transcript_x[i] = transcript_state[i].base_x;
            rows.transcript_y[i] = transcript_state[i].base_y;
            rows.transcript_z1[i] = transcript_state[i].z1;
            rows.transcript_z2[i] = transcript_state[i].z2;
            rows.transcript_z1zero[i] = transcript_state[i].z1_zero;
            rows.transcript_z2zero[i] = transcript_state[i].z2_zero;
            rows.transcript_op[i] = transcript_state[i].opcode;
            rows.transcript_accumulator_x[i] = transcript_state[i].accumulator_x;
            rows.transcript_accumulator_y[i] = transcript_state[i].accumulator_y;
            rows.transcript_msm_x[i] = transcript_state[i].msm_output_x;
            rows.transcript_msm_y[i] = transcript_state[i].msm_output_y;
            rows.transcript_collision_check[i] = transcript_state[i].collision_check;
        }

        return polynomials;
    }

    bool check_circuit()
    {

        proof_system::RelationParameters<typename Flavor::FF> params{
            .eta = 0,
            .beta = beta,
            .gamma = gamma,
            .public_input_delta = 0,
            .lookup_grand_product_delta = 0,
            .beta_sqr = beta_sqr,
            .beta_cube = beta_cube,
            .eccvm_set_permutation_delta = eccvm_set_permutation_delta,
        };

        auto rows = compute_full_polynomials();
        const size_t num_rows = rows[0].size();

        const auto evaluate_relation = [&]<typename Relation>(const std::string& relation_name) {
            auto relation = Relation();
            typename Relation::RelationValues result;
            for (auto& r : result) {
                r = 0;
            }
            constexpr size_t NUM_SUBRELATIONS = result.size();

            for (size_t i = 0; i < num_rows; ++i) {
                typename Flavor::RowPolynomials row;
                for (size_t j = 0; j < NUM_POLYNOMIALS; ++j) {
                    row[j] = rows[j][i];
                }
                relation.add_full_relation_value_contribution(result, row, params, 1);

                bool x = true;
                for (size_t j = 0; j < NUM_SUBRELATIONS; ++j) {
                    if (result[j] != 0) {
                        info("Relation ", relation_name, ", subrelation index ", j, " failed at row ", i);
                        x = false;
                    }
                }
                if (!x) {
                    return false;
                }
            }
            return true;
        };

        bool result = true;
        result = result && evaluate_relation.template operator()<BabyVMAdditionRelation<FF>>("BabyVMAdditionRelation");
        result = result && evaluate_relation.template operator()<BabyVMMultiplicationRelation<FF>>(
                               "BabyVMMultiplicationRelation");

        return result;
    }

    [[nodiscard]] size_t get_num_gates() const { return 0 }

    [[nodiscard]] size_t get_circuit_subgroup_size(const size_t num_rows) const
    {

        const auto num_rows_log2 = static_cast<size_t>(numeric::get_msb64(num_rows));
        size_t num_rows_pow2 = 1UL << (num_rows_log2 + (1UL << num_rows_log2 == num_rows ? 0 : 1));
        return num_rows_pow2;
    }
};
} // namespace proof_system

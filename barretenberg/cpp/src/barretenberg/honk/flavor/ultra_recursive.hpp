#pragma once
#include "barretenberg/ecc/curves/bn254/g1.hpp"
#include "barretenberg/honk/pcs/commitment_key.hpp"
#include "barretenberg/honk/pcs/kzg/kzg.hpp"
#include "barretenberg/polynomials/barycentric.hpp"
#include "barretenberg/polynomials/univariate.hpp"

#include "barretenberg/honk/flavor/ultra.hpp"
#include "barretenberg/honk/transcript/transcript.hpp"
#include "barretenberg/polynomials/evaluation_domain.hpp"
#include "barretenberg/polynomials/polynomial.hpp"
#include "barretenberg/proof_system/circuit_builder/ultra_circuit_builder.hpp"
#include "barretenberg/proof_system/flavor/flavor.hpp"
#include "barretenberg/proof_system/relations/auxiliary_relation.hpp"
#include "barretenberg/proof_system/relations/elliptic_relation.hpp"
#include "barretenberg/proof_system/relations/gen_perm_sort_relation.hpp"
#include "barretenberg/proof_system/relations/lookup_relation.hpp"
#include "barretenberg/proof_system/relations/permutation_relation.hpp"
#include "barretenberg/proof_system/relations/ultra_arithmetic_relation.hpp"
#include "barretenberg/srs/factories/crs_factory.hpp"

#include <array>
#include <concepts>
#include <span>
#include <string>
#include <type_traits>
#include <vector>

#include "barretenberg/stdlib/primitives/curves/bn254.hpp"
#include "barretenberg/stdlib/primitives/field/field.hpp"

namespace proof_system::honk::flavor {

/**
 * @brief The recursive counterpart to the "native" Ultra flavor.
 * @details This flavor can be used to instantiate a recursive Ultra Honk verifier for a proof created using the
 * conventional Ultra flavor. It is similar in structure to its native counterpart with two main differences: 1) the
 * curve types are stdlib types (e.g. field_t instead of field) and 2) it does not specify any Prover related types
 * (e.g. Polynomial, ExtendedEdges, etc.) since we do not emulate prover computation in circuits, i.e. it only makes
 * sense to instantiate a Verifier with this flavor.
 *
 * @note Unlike conventional flavors, "recursive" flavors are templated by a builder (much like native vs stdlib types).
 * This is because the flavor itself determines the details of the underlying verifier algorithm (i.e. the set of
 * relations), while the Builder determines the arithmetization of that algorithm into a circuit.
 *
 * @tparam BuilderType Determines the arithmetization of the verifier circuit defined based on this flavor.
 */
template <typename BuilderType> class UltraRecursive_ {
  public:
    using CircuitBuilder = BuilderType; // Determines arithmetization of circuit instantiated with this flavor
    using Curve = plonk::stdlib::bn254<CircuitBuilder>;
    using GroupElement = typename Curve::Element;
    using Commitment = typename Curve::Element;
    using CommitmentHandle = typename Curve::Element;
    using FF = typename Curve::ScalarField;

    // Note(luke): Eventually this may not be needed at all
    using VerifierCommitmentKey = pcs::VerifierCommitmentKey<Curve>;

    static constexpr size_t NUM_WIRES = flavor::Ultra::NUM_WIRES;
    // The number of multivariate polynomials on which a sumcheck prover sumcheck operates (including shifts). We often
    // need containers of this size to hold related data, so we choose a name more agnostic than `NUM_POLYNOMIALS`.
    // Note: this number does not include the individual sorted list polynomials.
    static constexpr size_t NUM_ALL_ENTITIES = 43;
    // The number of polynomials precomputed to describe a circuit and to aid a prover in constructing a satisfying
    // assignment of witnesses. We again choose a neutral name.
    static constexpr size_t NUM_PRECOMPUTED_ENTITIES = 25;
    // The total number of witness entities not including shifts.
    static constexpr size_t NUM_WITNESS_ENTITIES = 11;

    // define the tuple of Relations that comprise the Sumcheck relation
    using Relations = std::tuple<proof_system::UltraArithmeticRelation<FF>,
                                 proof_system::UltraPermutationRelation<FF>,
                                 proof_system::LookupRelation<FF>,
                                 proof_system::GenPermSortRelation<FF>,
                                 proof_system::EllipticRelation<FF>,
                                 proof_system::AuxiliaryRelation<FF>>;

    static constexpr size_t MAX_RELATION_LENGTH = get_max_relation_length<Relations>();

    // MAX_RANDOM_RELATION_LENGTH = algebraic degree of sumcheck relation *after* multiplying by the `pow_zeta` random
    // polynomial e.g. For \sum(x) [A(x) * B(x) + C(x)] * PowZeta(X), relation length = 2 and random relation length = 3
    static constexpr size_t MAX_RANDOM_RELATION_LENGTH = MAX_RELATION_LENGTH + 1;
    static constexpr size_t NUM_RELATIONS = std::tuple_size<Relations>::value;

    // define the container for storing the univariate contribution from each relation in Sumcheck
    using TupleOfTuplesOfUnivariates = decltype(create_relation_univariates_container<FF, Relations>());
    using TupleOfArraysOfValues = decltype(create_relation_values_container<FF, Relations>());

  private:
    template <typename DataType, typename HandleType>
    /**
     * @brief A base class labelling precomputed entities and (ordered) subsets of interest.
     * @details Used to build the proving key and verification key.
     */
    class PrecomputedEntities : public PrecomputedEntities_<DataType, HandleType, NUM_PRECOMPUTED_ENTITIES> {
      public:
        DataType& q_m = std::get<0>(this->_data);
        DataType& q_c = std::get<1>(this->_data);
        DataType& q_l = std::get<2>(this->_data);
        DataType& q_r = std::get<3>(this->_data);
        DataType& q_o = std::get<4>(this->_data);
        DataType& q_4 = std::get<5>(this->_data);
        DataType& q_arith = std::get<6>(this->_data);
        DataType& q_sort = std::get<7>(this->_data);
        DataType& q_elliptic = std::get<8>(this->_data);
        DataType& q_aux = std::get<9>(this->_data);
        DataType& q_lookup = std::get<10>(this->_data);
        DataType& sigma_1 = std::get<11>(this->_data);
        DataType& sigma_2 = std::get<12>(this->_data);
        DataType& sigma_3 = std::get<13>(this->_data);
        DataType& sigma_4 = std::get<14>(this->_data);
        DataType& id_1 = std::get<15>(this->_data);
        DataType& id_2 = std::get<16>(this->_data);
        DataType& id_3 = std::get<17>(this->_data);
        DataType& id_4 = std::get<18>(this->_data);
        DataType& table_1 = std::get<19>(this->_data);
        DataType& table_2 = std::get<20>(this->_data);
        DataType& table_3 = std::get<21>(this->_data);
        DataType& table_4 = std::get<22>(this->_data);
        DataType& lagrange_first = std::get<23>(this->_data);
        DataType& lagrange_last = std::get<24>(this->_data);

        std::vector<HandleType> get_selectors() override
        {
            return { q_m, q_c, q_l, q_r, q_o, q_4, q_arith, q_sort, q_elliptic, q_aux, q_lookup };
        };
        std::vector<HandleType> get_sigma_polynomials() override { return { sigma_1, sigma_2, sigma_3, sigma_4 }; };
        std::vector<HandleType> get_id_polynomials() override { return { id_1, id_2, id_3, id_4 }; };

        std::vector<HandleType> get_table_polynomials() { return { table_1, table_2, table_3, table_4 }; };
    };

    /**
     * @brief Container for all witness polynomials used/constructed by the prover.
     * @details Shifts are not included here since they do not occupy their own memory.
     */
    template <typename DataType, typename HandleType>
    class WitnessEntities : public WitnessEntities_<DataType, HandleType, NUM_WITNESS_ENTITIES> {
      public:
        DataType& w_l = std::get<0>(this->_data);
        DataType& w_r = std::get<1>(this->_data);
        DataType& w_o = std::get<2>(this->_data);
        DataType& w_4 = std::get<3>(this->_data);
        DataType& sorted_1 = std::get<4>(this->_data);
        DataType& sorted_2 = std::get<5>(this->_data);
        DataType& sorted_3 = std::get<6>(this->_data);
        DataType& sorted_4 = std::get<7>(this->_data);
        DataType& sorted_accum = std::get<8>(this->_data);
        DataType& z_perm = std::get<9>(this->_data);
        DataType& z_lookup = std::get<10>(this->_data);

        std::vector<HandleType> get_wires() override { return { w_l, w_r, w_o, w_4 }; };
        // The sorted concatenations of table and witness data needed for plookup.
        std::vector<HandleType> get_sorted_polynomials() { return { sorted_1, sorted_2, sorted_3, sorted_4 }; };
    };

    /**
     * @brief A base class labelling all entities (for instance, all of the polynomials used by the prover during
     * sumcheck) in this Honk variant along with particular subsets of interest
     * @details Used to build containers for: the prover's polynomial during sumcheck; the sumcheck's folded
     * polynomials; the univariates consturcted during during sumcheck; the evaluations produced by sumcheck.
     *
     * Symbolically we have: AllEntities = PrecomputedEntities + WitnessEntities + "ShiftedEntities". It could be
     * implemented as such, but we have this now.
     */
    template <typename DataType, typename HandleType>
    class AllEntities : public AllEntities_<DataType, HandleType, NUM_ALL_ENTITIES> {
      public:
        DataType& q_c = std::get<0>(this->_data);
        DataType& q_l = std::get<1>(this->_data);
        DataType& q_r = std::get<2>(this->_data);
        DataType& q_o = std::get<3>(this->_data);
        DataType& q_4 = std::get<4>(this->_data);
        DataType& q_m = std::get<5>(this->_data);
        DataType& q_arith = std::get<6>(this->_data);
        DataType& q_sort = std::get<7>(this->_data);
        DataType& q_elliptic = std::get<8>(this->_data);
        DataType& q_aux = std::get<9>(this->_data);
        DataType& q_lookup = std::get<10>(this->_data);
        DataType& sigma_1 = std::get<11>(this->_data);
        DataType& sigma_2 = std::get<12>(this->_data);
        DataType& sigma_3 = std::get<13>(this->_data);
        DataType& sigma_4 = std::get<14>(this->_data);
        DataType& id_1 = std::get<15>(this->_data);
        DataType& id_2 = std::get<16>(this->_data);
        DataType& id_3 = std::get<17>(this->_data);
        DataType& id_4 = std::get<18>(this->_data);
        DataType& table_1 = std::get<19>(this->_data);
        DataType& table_2 = std::get<20>(this->_data);
        DataType& table_3 = std::get<21>(this->_data);
        DataType& table_4 = std::get<22>(this->_data);
        DataType& lagrange_first = std::get<23>(this->_data);
        DataType& lagrange_last = std::get<24>(this->_data);
        DataType& w_l = std::get<25>(this->_data);
        DataType& w_r = std::get<26>(this->_data);
        DataType& w_o = std::get<27>(this->_data);
        DataType& w_4 = std::get<28>(this->_data);
        DataType& sorted_accum = std::get<29>(this->_data);
        DataType& z_perm = std::get<30>(this->_data);
        DataType& z_lookup = std::get<31>(this->_data);
        DataType& table_1_shift = std::get<32>(this->_data);
        DataType& table_2_shift = std::get<33>(this->_data);
        DataType& table_3_shift = std::get<34>(this->_data);
        DataType& table_4_shift = std::get<35>(this->_data);
        DataType& w_l_shift = std::get<36>(this->_data);
        DataType& w_r_shift = std::get<37>(this->_data);
        DataType& w_o_shift = std::get<38>(this->_data);
        DataType& w_4_shift = std::get<39>(this->_data);
        DataType& sorted_accum_shift = std::get<40>(this->_data);
        DataType& z_perm_shift = std::get<41>(this->_data);
        DataType& z_lookup_shift = std::get<42>(this->_data);

        std::vector<HandleType> get_wires() override { return { w_l, w_r, w_o, w_4 }; };
        // Gemini-specific getters.
        std::vector<HandleType> get_unshifted() override
        {
            return { q_c,           q_l,   q_r,      q_o,     q_4,     q_m,          q_arith, q_sort,
                     q_elliptic,    q_aux, q_lookup, sigma_1, sigma_2, sigma_3,      sigma_4, id_1,
                     id_2,          id_3,  id_4,     table_1, table_2, table_3,      table_4, lagrange_first,
                     lagrange_last, w_l,   w_r,      w_o,     w_4,     sorted_accum, z_perm,  z_lookup

            };
        };
        std::vector<HandleType> get_to_be_shifted() override
        {
            return { table_1, table_2, table_3, table_4, w_l, w_r, w_o, w_4, sorted_accum, z_perm, z_lookup };
        };
        std::vector<HandleType> get_shifted() override
        {
            return { table_1_shift, table_2_shift, table_3_shift,      table_4_shift, w_l_shift,     w_r_shift,
                     w_o_shift,     w_4_shift,     sorted_accum_shift, z_perm_shift,  z_lookup_shift };
        };

        AllEntities() = default;

        AllEntities(const AllEntities& other)
            : AllEntities_<DataType, HandleType, NUM_ALL_ENTITIES>(other){};

        AllEntities(AllEntities&& other)
            : AllEntities_<DataType, HandleType, NUM_ALL_ENTITIES>(other){};

        AllEntities& operator=(const AllEntities& other)
        {
            if (this == &other) {
                return *this;
            }
            AllEntities_<DataType, HandleType, NUM_ALL_ENTITIES>::operator=(other);
            return *this;
        }

        AllEntities& operator=(AllEntities&& other)
        {
            AllEntities_<DataType, HandleType, NUM_ALL_ENTITIES>::operator=(other);
            return *this;
        }

        ~AllEntities() = default;
    };

  public:
    /**
     * @brief The verification key is responsible for storing the the commitments to the precomputed (non-witnessk)
     * polynomials used by the verifier.
     *
     * @note Note the discrepancy with what sort of data is stored here vs in the proving key. We may want to resolve
     * that, and split out separate PrecomputedPolynomials/Commitments data for clarity but also for portability of our
     * circuits.
     */
    class VerificationKey : public VerificationKey_<PrecomputedEntities<Commitment, CommitmentHandle>> {
      public:
        /**
         * @brief Construct a new Verification Key with stdlib types from a provided native verification key
         *
         * @param builder
         * @param native_key Native verification key from which to extract the precomputed commitments
         */
        VerificationKey(CircuitBuilder* builder, auto native_key)
            : VerificationKey_<PrecomputedEntities<Commitment, CommitmentHandle>>(native_key->circuit_size,
                                                                                  native_key->num_public_inputs)
        {
            this->q_m = Commitment::from_witness(builder, native_key->q_m);
            this->q_l = Commitment::from_witness(builder, native_key->q_l);
            this->q_r = Commitment::from_witness(builder, native_key->q_r);
            this->q_o = Commitment::from_witness(builder, native_key->q_o);
            this->q_4 = Commitment::from_witness(builder, native_key->q_4);
            this->q_c = Commitment::from_witness(builder, native_key->q_c);
            this->q_arith = Commitment::from_witness(builder, native_key->q_arith);
            this->q_sort = Commitment::from_witness(builder, native_key->q_sort);
            this->q_elliptic = Commitment::from_witness(builder, native_key->q_elliptic);
            this->q_aux = Commitment::from_witness(builder, native_key->q_aux);
            this->q_lookup = Commitment::from_witness(builder, native_key->q_lookup);
            this->sigma_1 = Commitment::from_witness(builder, native_key->sigma_1);
            this->sigma_2 = Commitment::from_witness(builder, native_key->sigma_2);
            this->sigma_3 = Commitment::from_witness(builder, native_key->sigma_3);
            this->sigma_4 = Commitment::from_witness(builder, native_key->sigma_4);
            this->id_1 = Commitment::from_witness(builder, native_key->id_1);
            this->id_2 = Commitment::from_witness(builder, native_key->id_2);
            this->id_3 = Commitment::from_witness(builder, native_key->id_3);
            this->id_4 = Commitment::from_witness(builder, native_key->id_4);
            this->table_1 = Commitment::from_witness(builder, native_key->table_1);
            this->table_2 = Commitment::from_witness(builder, native_key->table_2);
            this->table_3 = Commitment::from_witness(builder, native_key->table_3);
            this->table_4 = Commitment::from_witness(builder, native_key->table_4);
            this->lagrange_first = Commitment::from_witness(builder, native_key->lagrange_first);
            this->lagrange_last = Commitment::from_witness(builder, native_key->lagrange_last);
        };
    };

    /**
     * @brief A field element for each entity of the flavor.
     */
    class AllValues : public AllEntities<FF, FF> {
      public:
        using Base = AllEntities<FF, FF>;
        using Base::Base;
        AllValues(std::array<FF, NUM_ALL_ENTITIES> _data_in) { this->_data = _data_in; }
    };

    /**
     * @brief A container for commitment labels.
     * @note It's debatable whether this should inherit from AllEntities. since most entries are not strictly needed. It
     * has, however, been useful during debugging to have these labels available.
     *
     */
    class CommitmentLabels : public AllEntities<std::string, std::string> {
      public:
        CommitmentLabels()
        {
            this->w_l = "W_L";
            this->w_r = "W_R";
            this->w_o = "W_O";
            this->w_4 = "W_4";
            this->z_perm = "Z_PERM";
            this->z_lookup = "Z_LOOKUP";
            this->sorted_accum = "SORTED_ACCUM";

            // The ones beginning with "__" are only used for debugging
            this->q_c = "__Q_C";
            this->q_l = "__Q_L";
            this->q_r = "__Q_R";
            this->q_o = "__Q_O";
            this->q_4 = "__Q_4";
            this->q_m = "__Q_M";
            this->q_arith = "__Q_ARITH";
            this->q_sort = "__Q_SORT";
            this->q_elliptic = "__Q_ELLIPTIC";
            this->q_aux = "__Q_AUX";
            this->q_lookup = "__Q_LOOKUP";
            this->sigma_1 = "__SIGMA_1";
            this->sigma_2 = "__SIGMA_2";
            this->sigma_3 = "__SIGMA_3";
            this->sigma_4 = "__SIGMA_4";
            this->id_1 = "__ID_1";
            this->id_2 = "__ID_2";
            this->id_3 = "__ID_3";
            this->id_4 = "__ID_4";
            this->table_1 = "__TABLE_1";
            this->table_2 = "__TABLE_2";
            this->table_3 = "__TABLE_3";
            this->table_4 = "__TABLE_4";
            this->lagrange_first = "__LAGRANGE_FIRST";
            this->lagrange_last = "__LAGRANGE_LAST";
        };
    };

    class VerifierCommitments : public AllEntities<Commitment, CommitmentHandle> {
      public:
        VerifierCommitments(std::shared_ptr<VerificationKey> verification_key)
        {
            this->q_m = verification_key->q_m;
            this->q_l = verification_key->q_l;
            this->q_r = verification_key->q_r;
            this->q_o = verification_key->q_o;
            this->q_4 = verification_key->q_4;
            this->q_c = verification_key->q_c;
            this->q_arith = verification_key->q_arith;
            this->q_sort = verification_key->q_sort;
            this->q_elliptic = verification_key->q_elliptic;
            this->q_aux = verification_key->q_aux;
            this->q_lookup = verification_key->q_lookup;
            this->sigma_1 = verification_key->sigma_1;
            this->sigma_2 = verification_key->sigma_2;
            this->sigma_3 = verification_key->sigma_3;
            this->sigma_4 = verification_key->sigma_4;
            this->id_1 = verification_key->id_1;
            this->id_2 = verification_key->id_2;
            this->id_3 = verification_key->id_3;
            this->id_4 = verification_key->id_4;
            this->table_1 = verification_key->table_1;
            this->table_2 = verification_key->table_2;
            this->table_3 = verification_key->table_3;
            this->table_4 = verification_key->table_4;
            this->lagrange_first = verification_key->lagrange_first;
            this->lagrange_last = verification_key->lagrange_last;
        }
    };
};

} // namespace proof_system::honk::flavor

// A contract used for testing a random hodgepodge of small features from simulator and end-to-end tests.
contract Test {
    use dep::aztec::{
        abi,
        abi::PrivateContextInputs,
        context::PrivateContext,
        oracle::{
            get_public_key::get_public_key,
            context::get_portal_address,
            rand::rand
        },
        log::emit_unencrypted_log,
        types::vec::BoundedVec,
        constants_gen::EMPTY_NULLIFIED_COMMITMENT,
    };

    #[event]
    struct ExampleEvent {
        value: Field,
    }

    #[aztec(private)]
    fn constructor() {}

    #[aztec(private)]
    fn getPublicKey(
        address: Field,
    ) -> [Field; 2]{
        let pub_key = get_public_key(address);
        
        [pub_key.x, pub_key.y]
    }

    // Get the portal contract address through an oracle call
    #[aztec(private)]
    fn getPortalContractAddress(
        aztec_address: Field
    ) -> Field {
        get_portal_address(aztec_address)
    }

    // Get the address of the l1 portal for this contract (taken from the input context)
    #[aztec(private)]
    fn getThisPortalAddress() -> Field {
        context.this_portal_address()
    }

    // Get the address of this contract (taken from the input context)
    #[aztec(private)]
    fn getThisAddress() -> Field {
        context.this_address()
    }

    // Test codegen for Aztec.nr interfaces
    // See yarn-project/acir-simulator/src/client/private_execution.test.ts 'nested calls through autogenerated interface'
    // Note; this function is deliberately NOT annotated with #[aztec(private)] due to its use in tests
    fn testCodeGen(
        inputs: PrivateContextInputs,
        aField: Field,
        aBool: bool,
        aNumber: u32,
        anArray: [Field; 2],
        aStruct: DummyNote,
        aDeepStruct: DeepStruct,
    ) -> distinct pub abi::PrivateCircuitPublicInputs {
        let mut args: BoundedVec<Field, 17> = BoundedVec::new(0);
        args.push(aField);
        args.push(aBool as Field);
        args.push(aNumber as Field);
        args.push_array(anArray);
        args.push(aStruct.amount);
        args.push(aStruct.secretHash);
        args.push(aDeepStruct.aField);
        args.push(aDeepStruct.aBool as Field);
        args.push(aDeepStruct.aNote.amount);
        args.push(aDeepStruct.aNote.secretHash);
        for note in aDeepStruct.manyNotes {
            args.push(note.amount);
            args.push(note.secretHash);            
        }
        let args_hash = abi::hash_args(args.storage);
        let mut context = PrivateContext::new(inputs, args_hash);
        context.return_values.push(args_hash);
        context.finish()
    }

    // Purely exists for testing
    #[aztec(public)]
    fn createL2ToL1MessagePublic(
        amount: Field,
        secretHash: Field,
    ) {
        // Create a commitment to the amount
        let note = DummyNote::new(amount, secretHash);
        
        // Public oracle call to emit new commitment.
        context.message_portal(note.get_commitment());
        
    }

    // Purely exists for testing
    #[aztec(public)]
    fn createNullifierPublic(
        amount: Field,
        secretHash: Field,
    ) {
        // Create a commitment to the amount
        let note = DummyNote::new(amount, secretHash);
        
        // Public oracle call to emit new commitment.
        context.push_new_nullifier(note.get_commitment(), EMPTY_NULLIFIED_COMMITMENT);
    }

    // Forcefully emits a nullifier (for testing purposes)
    #[aztec(private)]
    fn emit_nullifier(nullifier: Field) {
        context.push_new_nullifier(nullifier, EMPTY_NULLIFIED_COMMITMENT);
    }

    // docs:start:is-time-equal
    #[aztec(public)]
    fn isTimeEqual(
        time: Field,
    ) -> Field {
        assert(context.timestamp() == time);
        time
    }
    // docs:end:is-time-equal

    // docs:start:emit_unencrypted_log
    #[aztec(public)]
    fn emit_unencrypted(
        value: Field
    ) -> Field {
        emit_unencrypted_log(&mut context, value);
        0
    }
    // docs:end:emit_unencrypted_log

    // Purely exists for testing
    unconstrained fn getRandom(
        kindaSeed: Field
    ) -> Field {
        kindaSeed * rand()
    }

    struct DummyNote {
        amount: Field,
        secretHash: Field
    }

    impl DummyNote {
        fn new(amount: Field, secretHash: Field) -> pub Self {
            Self {
                amount: amount,
                secretHash: secretHash
            }
        }

        fn get_commitment(self) -> pub Field {
            dep::std::hash::pedersen([self.amount, self.secretHash])[0]
        }
    }

    struct DeepStruct {
        aField: Field,
        aBool: bool,
        aNote: DummyNote,
        manyNotes: [DummyNote; 3],
    }
}

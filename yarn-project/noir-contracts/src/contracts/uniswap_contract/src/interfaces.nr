use dep::aztec::{
  context::{ PrivateContext, PublicContext, Context },
  selector::compute_selector,
  types::address::AztecAddress,
};

struct Token {
  address: Field,
}

impl Token {
  pub fn at(address: Field) -> Self {
    Self { address }
  }

  pub fn transfer_public(self: Self, context: PublicContext, from: Field, to: Field, amount: Field, nonce: Field) {
    let _transfer_return_values = context.call_public_function(
      self.address,
      compute_selector("transfer_public((Field),(Field),Field,Field)"),
      [from, to, amount, nonce]
    );
  }

  pub fn unshield(self: Self, context: &mut PrivateContext, from: Field, to: Field,  amount: Field, nonce: Field) {
    let _return_values = context.call_private_function(
      self.address,
      compute_selector("unshield((Field),(Field),Field,Field)"),
      [from, to, amount, nonce]
    );
  }
}

struct TokenBridge {
  address: Field,
}

impl TokenBridge {
  pub fn at(address: Field) -> Self {
    Self { address }
  }

  pub fn token(self: Self, context: PublicContext) -> AztecAddress {
    let return_values = context.call_public_function(self.address, compute_selector("get_token()"), []);
    AztecAddress::new(return_values[0])
  }

  pub fn claim_public(self: Self, context: PublicContext, to: Field, amount: Field, canceller: Field, msg_key: Field, secret: Field) {
    let _return_values = context.call_public_function(
      self.address,
      compute_selector("claim_public((Field),Field,Field,Field,Field)"),
      [to, amount, canceller, msg_key, secret]
    );
  }

  pub fn exit_to_l1_public(self: Self, context: PublicContext, recipient: Field, amount: Field, callerOnL1: Field, nonce: Field) {
    let _return_values = context.call_public_function(
      self.address,
      compute_selector("exit_to_l1_public((Field),Field,(Field),Field)"),
      [recipient, amount, callerOnL1, nonce]
    );
  }
}

use dep::aztec::note::note_interface::NoteInterface;
use dep::aztec::note::note_header::NoteHeader;
use dep::aztec::note::utils::compute_unique_siloed_note_hash;
use dep::aztec::oracle::get_secret_key::get_secret_key;

global PUBLIC_KEY_NOTE_LEN: Field = 3;

// Stores a public key composed of two fields
// TODO: Do we need to include a nonce, in case we want to read/nullify/recreate with the same pubkey value?
struct PublicKeyNote {
    x: Field,
    y: Field,
    owner: Field, // We store the owner address only to get the secret key to compute the nullifier
    header: NoteHeader,
}

impl PublicKeyNote {
    pub fn new(x: Field, y: Field, owner: Field) -> Self {
        PublicKeyNote {
            x,
            y,
            owner,
            header: NoteHeader::empty(),
        }
    }

    // serialize the note as 3 fields
    pub fn serialize(self) -> [Field; PUBLIC_KEY_NOTE_LEN] {
        [self.x, self.y, self.owner]
    }

    pub fn compute_nullifier(self) -> Field {
        let unique_siloed_note_hash = compute_unique_siloed_note_hash(PublicKeyNoteMethods, self);
        let secret = get_secret_key(self.owner);
        // TODO(#1205) Should use a non-zero generator index.
        dep::std::hash::pedersen([
            unique_siloed_note_hash,
            secret.low,
            secret.high,
        ])[0]
    }

    pub fn set_header(&mut self, header: NoteHeader) {
        self.header = header;
    }
}

fn deserialize(preimage: [Field; PUBLIC_KEY_NOTE_LEN]) -> PublicKeyNote {
    PublicKeyNote {
        x: preimage[0],
        y: preimage[1],
        owner: preimage[2],
        header: NoteHeader::empty(),
    }
}

fn serialize(note: PublicKeyNote) -> [Field; PUBLIC_KEY_NOTE_LEN] {
    note.serialize()
}

fn compute_note_hash(note: PublicKeyNote) -> Field {
    // TODO(#1205) Should use a non-zero generator index.
    dep::std::hash::pedersen(note.serialize())[0]
}

fn compute_nullifier(note: PublicKeyNote) -> Field {
    note.compute_nullifier()
}

fn get_header(note: PublicKeyNote) -> NoteHeader {
    note.header
}

fn set_header(note: &mut PublicKeyNote, header: NoteHeader) {
    note.set_header(header);
}

global PublicKeyNoteMethods = NoteInterface {
    deserialize,
    serialize,
    compute_note_hash,
    compute_nullifier,
    get_header,
    set_header,
};

//! `boatmath-web` — a browser front end for the `michell` tools.
//!
//! Upload a hull and see it the way the physics sees it: cut into sections,
//! then load it into cases and queue studies on them. The computations are
//! [`boatmath`]'s, shared with the CLI; this crate keeps them in a store
//! ([`store`]) and runs a queue of them ([`worker`]).

pub use boatmath::*;

pub mod store;
pub mod worker;

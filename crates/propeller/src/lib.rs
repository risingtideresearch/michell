//! Propellers and motors: the best Wageningen B-series propeller for a
//! thrust at a speed ([`bseries`]), and the motors that can drive it, each
//! at its best reduction and operating point ([`motor`]).
//!
//! Ported from [propopt](https://github.com/risingtideresearch/propopt)'s web
//! cores, `web/propcore.js` and `web/motorcore.js`, with its motor database
//! vendored in `data/`. The golden tests hold the port to the original's
//! answers (`tests/golden/make.js` makes them).

// `!(x > 0.0)` is used throughout, as in the original, to reject NaN too.
#![allow(clippy::neg_cmp_op_on_partial_ord)]

pub mod bseries;
pub mod motor;

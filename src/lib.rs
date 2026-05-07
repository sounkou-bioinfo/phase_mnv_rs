//! Library core for `phase_tools-rs`.
//!
//! The command-line binaries in this package are thin frontends over this
//! library. The public module layout is being stabilized incrementally; prefer
//! adding reusable genomics kernels here before exposing new binary-only logic.

pub mod assembly;
pub mod io;

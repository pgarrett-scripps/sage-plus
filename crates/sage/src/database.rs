//! Search-space definition, digestion and the fragment index.

mod builder;
mod digest;
mod estimate;
mod fragment_index;
mod limit;
mod mass_offset;
mod query;

pub use builder::*;
pub use digest::*;
pub use estimate::*;
pub use fragment_index::*;
pub use limit::Limit;
pub use mass_offset::*;
pub use query::*;

#[cfg(test)]
#[path = "../tests/unit/database.rs"]
mod test;

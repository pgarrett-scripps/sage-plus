//! Values a modification declares for some or all of its sites.
//!
//! A modification field of this type accepts either a list, which applies to
//! every site of the modification, or a map from site to list:
//!
//! ```json
//! "Phospho": {"mass": 79.966331, "sites": ["S", "T", "Y"], "immonium_ions": {"Y": [216.0420]}}
//! "Acetyl":  {"mass": 42.010565, "sites": ["K"], "immonium_ions": [126.0913]}
//! ```
//!
//! Map keys must be sites that the modification declares in `sites`, spelled
//! the same way; a site the map does not list gets nothing.

use std::collections::BTreeMap;

use serde::{
    de::{self, MapAccess, SeqAccess, Visitor},
    Deserialize, Deserializer, Serialize,
};

/// A list for every site, or a map from declared site to list.
#[derive(Clone, Debug, PartialEq, Serialize, schemars::JsonSchema)]
#[serde(untagged)]
pub enum SiteMap<T> {
    /// Applies to every site of the modification.
    All(Vec<T>),
    /// Applies to the listed sites only; keys are sites from `sites`.
    Sites(BTreeMap<String, Vec<T>>),
}

impl<T> Default for SiteMap<T> {
    fn default() -> Self {
        Self::All(Vec::new())
    }
}

impl<T> From<Vec<T>> for SiteMap<T> {
    fn from(values: Vec<T>) -> Self {
        Self::All(values)
    }
}

impl<T> SiteMap<T> {
    /// True when no site gets a value.
    pub fn is_empty(&self) -> bool {
        match self {
            Self::All(values) => values.is_empty(),
            Self::Sites(sites) => sites.values().all(Vec::is_empty),
        }
    }

    /// Values for `site`, spelled as in the modification's `sites`.
    pub fn for_site(&self, site: &str) -> &[T] {
        match self {
            Self::All(values) => values,
            Self::Sites(sites) => sites.get(site).map_or(&[], Vec::as_slice),
        }
    }

    /// Every value, for every site.
    pub fn values(&self) -> impl Iterator<Item = &T> {
        let (all, sites) = match self {
            Self::All(values) => (Some(values), None),
            Self::Sites(sites) => (None, Some(sites)),
        };
        all.into_iter()
            .flatten()
            .chain(sites.into_iter().flat_map(|sites| sites.values().flatten()))
    }

    /// Check that every map key is one of `sites`. The error names the field,
    /// the modification and the declared sites.
    pub fn validate_sites<S: AsRef<str>>(
        &self,
        field: &str,
        modification: &str,
        sites: &[S],
    ) -> Result<(), String> {
        let Self::Sites(map) = self else {
            return Ok(());
        };
        for key in map.keys() {
            if !sites.iter().any(|site| site.as_ref() == key) {
                let declared = sites
                    .iter()
                    .map(|site| format!("`{}`", site.as_ref()))
                    .collect::<Vec<_>>()
                    .join(", ");
                return Err(format!(
                    "`{field}` key `{key}` of modification `{modification}` is not one of its sites ({declared}); use a declared site or a list for all sites"
                ));
            }
        }
        Ok(())
    }
}

impl<'de, T: Deserialize<'de>> Deserialize<'de> for SiteMap<T> {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        struct SiteMapVisitor<T>(std::marker::PhantomData<T>);

        impl<'de, T: Deserialize<'de>> Visitor<'de> for SiteMapVisitor<T> {
            type Value = SiteMap<T>;

            fn expecting(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                formatter.write_str("a list for every site, or a map from site to list")
            }

            fn visit_seq<A: SeqAccess<'de>>(self, seq: A) -> Result<Self::Value, A::Error> {
                Vec::deserialize(de::value::SeqAccessDeserializer::new(seq)).map(SiteMap::All)
            }

            fn visit_map<A: MapAccess<'de>>(self, map: A) -> Result<Self::Value, A::Error> {
                BTreeMap::deserialize(de::value::MapAccessDeserializer::new(map))
                    .map(SiteMap::Sites)
            }
        }

        deserializer.deserialize_any(SiteMapVisitor(std::marker::PhantomData))
    }
}

#[cfg(test)]
#[path = "../tests/unit/site_map.rs"]
mod test;

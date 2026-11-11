//! A count or length that may be `"unlimited"`.
//!
//! Moved out of `database.rs` unchanged; `database::Limit` is the public path.

use serde::{Deserialize, Serialize};
use std::borrow::Cow;

/// A count or length that is a number or `"unlimited"`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Limit<T> {
    Limited(T),
    Unlimited,
}

impl<T: Copy> Limit<T> {
    /// The number, or `None` when unlimited.
    pub fn limited(self) -> Option<T> {
        match self {
            Limit::Limited(value) => Some(value),
            Limit::Unlimited => None,
        }
    }
}

impl<T> From<T> for Limit<T> {
    fn from(value: T) -> Self {
        Limit::Limited(value)
    }
}

impl<T: std::fmt::Display> std::fmt::Display for Limit<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Limit::Limited(value) => value.fmt(f),
            Limit::Unlimited => f.write_str("unlimited"),
        }
    }
}

impl<T: Serialize> Serialize for Limit<T> {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        match self {
            Limit::Limited(value) => value.serialize(serializer),
            Limit::Unlimited => serializer.serialize_str("unlimited"),
        }
    }
}

impl<'de, T: TryFrom<u64>> Deserialize<'de> for Limit<T> {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        struct Visitor<T>(std::marker::PhantomData<T>);
        impl<T: TryFrom<u64>> serde::de::Visitor<'_> for Visitor<T> {
            type Value = Limit<T>;

            fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                f.write_str("a nonnegative integer or \"unlimited\"")
            }

            fn visit_u64<E: serde::de::Error>(self, value: u64) -> Result<Self::Value, E> {
                T::try_from(value)
                    .map(Limit::Limited)
                    .map_err(|_| E::custom(format!("{value} is too large; use \"unlimited\"")))
            }

            fn visit_i64<E: serde::de::Error>(self, value: i64) -> Result<Self::Value, E> {
                match u64::try_from(value) {
                    Ok(value) => self.visit_u64(value),
                    Err(_) => Err(E::invalid_value(
                        serde::de::Unexpected::Signed(value),
                        &self,
                    )),
                }
            }

            fn visit_str<E: serde::de::Error>(self, value: &str) -> Result<Self::Value, E> {
                match value {
                    "unlimited" => Ok(Limit::Unlimited),
                    _ => Err(E::invalid_value(serde::de::Unexpected::Str(value), &self)),
                }
            }
        }
        deserializer.deserialize_any(Visitor(std::marker::PhantomData))
    }
}

impl<T: schemars::JsonSchema> schemars::JsonSchema for Limit<T> {
    fn inline_schema() -> bool {
        true
    }

    fn schema_name() -> Cow<'static, str> {
        format!("Limit_{}", T::schema_name()).into()
    }

    fn json_schema(generator: &mut schemars::SchemaGenerator) -> schemars::Schema {
        schemars::json_schema!({
            "anyOf": [generator.subschema_for::<T>(), {"const": "unlimited"}]
        })
    }
}

# Consumer processing adapter

Source: MannLabs/timsrust, crates.io timsrust-tdf 0.6.6, Apache-2.0. The published crate records upstream commit `80e235a331057e0ef922eb3c784b7a631c9b9666`, path `crates/timsrust-tdf`. The Apache-2.0 license is retained in LICENSE. This adapter is not publishable.

This consumer-owned adapter preserves upstream DDA and DIA spectrum processing, window splitting, smoothing, centroiding and range-based converters. Its SQL metadata facade uses timsrust-reader MetadataDatabase and its raw frame reader uses timsrust-reader TdfReader. No private reader source is copied here.

The existing generic reader contract uses shared references. A mutex serializes raw frame decoding, while spectrum processing retains its existing parallelism. Reader errors remain in the error source chain. miniTDF and unrelated formats keep their original readers.

The private reader owns bounds, decompression and read-only database validation. Processing policy remains with the consumer. The original compression modules remain solely for public error-type compatibility and are not constructed by the adapter.

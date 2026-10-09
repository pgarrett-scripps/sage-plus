# Third-party notices

Sage Plus reads Bruker TDF files with `sage-plus-tdf`, which derives its frame
decompression from MannLabs TimsRust 0.6.6. That crate is licensed
`Apache-2.0 AND MIT` and carries the upstream MIT text, the Apache 2.0 text and its
source attribution (`LICENSE`, `LICENSE-APACHE` and `NOTICE.md`) in its own
repository.
The Sage Plus license remains in the root `LICENSE` file.

## Unimod

Sage Plus compiles modification names and masses from the Unimod database
(https://www.unimod.org) into the executable. The data, in
`crates/sage/data/unimodifications.json`, is a JSON extract of Unimod's
`unimod.xml` (data version 17:06:2025 11:31) and carries Unimod's notice:
"Copyright (C) 2002-2006 Unimod; this information may be copied, distributed
and/or modified under certain conditions, but it comes WITHOUT ANY WARRANTY;
see the accompanying Design Science License for more details." The extract
remains under the Design Science License. Its source data is that JSON file in
the Sage Plus source distribution.

The notice, a description of the conversion, and the full Design Science
License text are in `crates/sage/data/LICENSE-unimod.txt`. Executable archives
include it as `LICENSE-unimod.txt`, and containers include it at
`/app/licenses/unimod.txt`.

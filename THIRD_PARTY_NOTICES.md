# Third-party notices

Sage Plus includes a patched copy of MannLabs `filemanager` 0.6.6. Its crate
metadata declares Apache-2.0, while the pinned repository root contains an MIT
license. The source distribution retains that upstream copyright and MIT text
unchanged in `vendor/filemanager/LICENSE` and includes the standard Apache 2.0
text in `vendor/filemanager/LICENSE-APACHE`. The same directory records the source
revision, original file hashes, and patch description.

Executable archives include both texts as `LICENSE-filemanager` and
`LICENSE-filemanager-APACHE`. Containers include them at
`/app/licenses/filemanager.txt` and `/app/licenses/filemanager-APACHE.txt`.
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

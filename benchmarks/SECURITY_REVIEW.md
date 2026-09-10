Hardening dependency review

Reviewed September 10, 2026 with cargo-audit 0.22.2 against a freshly downloaded RustSec database.
The scan examines the complete lockfile. Raw before and after JSON reports are retained with the
local hardening evidence. The scheduled CI scan fetches current advisories on each run.

**Resolved in the candidate**

| Dependency | Reviewed version | Candidate version | Advisory |
|---|---|---|---|
| crossbeam-epoch | 0.9.18 | 0.9.20 | RUSTSEC-2026-0204 |
| h2 | 0.4.13 | 0.4.16 | RUSTSEC-2026-0258 |
| anyhow | 1.0.102 | 1.0.103 | RUSTSEC-2026-0190, unsoundness |
| memmap2 | 0.9.10 | 0.9.11 | RUSTSEC-2026-0186, unsoundness |
| direct quick-xml | 0.31.0 | 0.41.0 | RUSTSEC-2026-0194 and RUSTSEC-2026-0195 |
| timsrust family | 0.6.4 | 0.6.5 | Replace the yanked direct release |

The mzML reader now uses the patched direct XML dependency. Workspace tests, native mzMLb tests,
minimal-feature tests, strict Clippy, and Rust 1.88 checks pass with these versions. These checks
include the committed real Bruker TDF fixture but do not exercise every remote storage service.

**Beta.3 dependency disposition**

The September 10 online audit passes with zero vulnerabilities, zero unsoundness
warnings, and no yanked dependency warnings. Two unmaintained dependencies remain:
`instant` and `paste`, associated with the native HDF5 stack. They are maintenance
warnings and do not fail the publication gate. No advisory is suppressed.

The published filemanager 0.6.5 and current upstream source still require
object_store 0.11, which brings quick-xml 0.37.5 and two availability advisories,
RUSTSEC-2026-0194 and RUSTSEC-2026-0195. The ordinary CLI, API runner, and MCP reject
remote Bruker URLs, which limits exposure there. Direct library path handling has a
wider surface. The release removes the old dependency instead of relying on a
claim that every affected path is unreachable.

The root Cargo patch uses the published filemanager 0.6.5 source from timsrust
commit `5e09572fa3f1aace86e6ba64244a01fd0850fbf3`. It updates object_store to 0.14.1,
imports the new extension trait, converts byte ranges to the new `u64` API, and
checks object-size conversion. The crate declares Apache-2.0, but its pinned
repository root has an MIT license. Both texts, the original copyright, original
file hashes, compatibility tests, and removal plan are retained in `vendor/filemanager`.
Archives and containers include the dependency license and third-party notice.

The patch tests exercise upload, metadata, byte ranges, listing, buffered file
upload, and download. A loopback HTTP S3 fixture checks XML listing and range
responses. All 47 upstream documentation tests also pass. This does not certify
all providers' authentication. Workspace tests separately cover the committed
real Bruker TDF fixture.
Remove the patch after a published dependency chain resolves the parser findings,
then repeat compatibility, workspace, and audit checks.

The update also removes the legacy rustls-pemfile warning. Updating chacha20 from
0.10.1 to 0.10.2 removes its yanked-release warning. The fresh online scan includes
the registry check because an offline advisory scan alone previously missed it.
Raw JSON is retained in the local beta.3 evidence bundle. Current release status
and final verification are recorded in [the beta.3 checklist](BETA3_RELEASE.md).

Reproduce the scan:

```shell
cargo install cargo-audit --version 0.22.2 --locked
cargo audit --file Cargo.lock --deny unsound
```

Read the current advisories in the [RustSec database](https://rustsec.org/advisories/).
An audit with no findings would still not establish that all dependencies are vulnerability free.

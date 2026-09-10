# Beta.3 release checklist

Prepared September 10, 2026 for `v0.1.0-beta.3` on `codex/release-beta3`.
This is a hardening release. It does not change scoring defaults or establish
broader scientific calibration.

## Local gates

| Gate | Status |
|---|---|
| Workspace version and changelog | beta.3, matching tag check passed |
| Online dependency audit | Passed, zero vulnerabilities, unsoundness findings, or yanked dependencies |
| Maintenance warnings | `instant` and `paste` remain tracked |
| Default and minimal workspace tests | 399 passed in each configuration |
| Patched filemanager compatibility | 2 new tests and 47 upstream documentation tests passed |
| Strict Clippy and Rust 1.88 compatibility check | Passed |
| Python benchmark tooling | 16 tests passed |
| Formatting and schema synchronization | Passed |
| Workspace line coverage | 81.70%, above the 80% gate |
| Optimized native archive and smoke search | Final verification pending |
| Four paired representative workloads | Final verification pending |
| Bounded entrapment comparison | September 4 evidence retained, final candidate verification pending |

Workspace tests include the committed real Bruker TDF and Thermo RAW fixtures,
mzMLb, MCP worker isolation, and the new hardening regressions. The storage patch
also exercises S3 XML and HTTP ranges through a loopback fixture. These checks do
not certify every cloud provider's authentication or acquisition workflow.

## Hosted gates

- [ ] Release preparation committed and available for review.
- [ ] Required Rust CI checks pass on the release preparation commit.
- [ ] Windows and macOS smoke tests pass.
- [ ] Manual release workflow builds all seven archives and the AMD64 container.
- [ ] Native archive downloaded, checksum verified, and smoke-tested.
- [ ] Maintainer merges the reviewed preparation into `main`.
- [ ] Annotated `v0.1.0-beta.3` tag points to the verified release source.
- [ ] Tag workflow publishes the prerelease, checksums, and versioned container.

GitHub currently authenticates as `afk-sapien` with read-only access to
`pgarrett-scripps/sage-plus`. Publishing requires a maintainer session or a
maintainer to complete merge, workflow dispatch, and tagging. A local archive
does not substitute for the hosted platform and publication gates.

## Release scope and evidence

The release fixes configured batching, gzip completion, event ordering, input
validation, cancellation boundaries, and MCP persistence. It introduces explicit
local overwrite, run-summary schema 9, verified benchmark manifests, and CI gates.
The pinned storage patch removes the remaining vulnerable transitive XML parser.

The [September 4 results](HARDENING_RESULTS.md) retain the original hashes,
timing investigations, and 20-seed entrapment evidence. The experiment uses one
HEK file and an incomplete 2,000-protein target FASTA. It has low power near 1%,
and estimates at higher cutoffs remain above nominal in both builds. Independent
study, protein-group, localization, transfer, and library calibration remain in
the [scientific protocol](SCIENTIFIC_PROTOCOL.md).

Reproduce the local experiments with [the hardening commands](HARDENING_USAGE.md).

See the [security review](SECURITY_REVIEW.md), [changelog](../CHANGELOG.md), and
[maintainer release procedure](../RELEASING.md) for the dependency disposition,
compatibility changes, and publication commands.

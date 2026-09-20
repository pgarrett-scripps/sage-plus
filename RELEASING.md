# Releasing Sage Plus

Sage Plus releases are built from annotated `v*` tags. The release workflow validates the tag,
runs the workspace tests, builds every supported executable archive, generates SHA-256 checksums,
builds the Linux AMD64 container, and only then publishes the GitHub release.

The workflow publishes executable archives to GitHub Releases and a container to
`ghcr.io/pgarrett-scripps/sage-plus`. It does not publish to crates.io, Bioconda, or Homebrew.

Release binaries are compiled on runners with the same CPU architecture as their targets. The
musl binaries are compiled inside the matching architecture of the official Rust Alpine image.
This is required because the bundled HDF5 configuration used by mzMLb support executes target
probes while it builds and therefore cannot use a conventional cross compiler. Native builds set
the supported CMake policy floor to 3.5 so bundled libraries with older declarations also configure
under CMake 4. The macOS builds predefine `fdopen` to itself while compiling the bundled zlib 1.2.11
source. This avoids an obsolete classic Mac compatibility branch that Xcode 16.3 and newer expose
through `TARGET_OS_MAC`.

## One-time repository setup

Keep the repository's default Actions token permission read-only. The release workflow grants
`contents: write` and `packages: write` only to the jobs that need them.

Protect `main` and require these Rust workflow checks before merging:

- `Format and release build`
- `Coverage`
- `Test Rust 1.88.0`
- `Test Rust 1.97.1`

Do not require the release workflow on ordinary branches. It is intended for manual packaging
checks and release tags.

## Account and destination check

Before remote changes, inspect `gh auth status` and the Git push URL. Use an
already-configured account with access to the existing repository. For the owner
account in this workspace:

```shell
gh auth switch --hostname github.com --user pgarrett-scripps
gh api user --jq .login
git remote get-url --push origin
```

The destination must be `pgarrett-scripps/sage-plus`. Do not create a personal fork
as a workaround for an inactive maintainer account.

## Prepare a release

The beta.3 security dependency gate passes, as recorded in
[`benchmarks/SECURITY_REVIEW.md`](benchmarks/SECURITY_REVIEW.md). Track the remaining hosted
validation and publication gates for the current release in
[`benchmarks/BETA6_RELEASE.md`](benchmarks/BETA6_RELEASE.md), and the previous release's in
[`benchmarks/BETA4_RELEASE.md`](benchmarks/BETA4_RELEASE.md).
The workflow audits the complete lockfile and runs the storage patch compatibility tests before
packaging. Archives include analytical schemas, the changelog, and third-party notices and licenses.

1. Set `[workspace.package].version` in `Cargo.toml`. All Sage Plus crates inherit this version.
2. Add a matching `## [vX.Y.Z]` or prerelease section to `CHANGELOG.md`, leaving a new empty
   `## [Unreleased]` section above it.
3. Update release-specific documentation, then regenerate and verify the lockfile:

   ```shell
   cargo check --workspace
   bash scripts/check-release-version.sh
   cargo fmt --all -- --check
   cargo test --workspace --locked
   cargo build --release --workspace --locked
   ```

4. Open a release preparation pull request against `main`. If it comes from a fork, a maintainer
   must approve the pending workflow runs in Actions. Wait for every required Rust and dependency
   check to pass, review the evidence, and merge the preparation into `main`.
5. Run `Release Sage Plus` manually from the Actions page. A manual run builds and retains all
   archives and validates the Docker build, but does not publish a release or container.
6. Download the `release-dist` artifact and inspect at least the archive for the maintainer's
   native platform.

## Publish

Create and push exactly one annotated tag after the preparation commit is on `main`:

```shell
git switch main
git pull --ff-only origin main
git tag -a v0.1.0-beta.6 -m "Sage Plus v0.1.0-beta.6"
git push origin v0.1.0-beta.6
```

The tag starts the release workflow. Prerelease identifiers such as `-beta.1` cause GitHub to mark
the release as a prerelease. Stable versions are marked as the latest release and also update the
container's `latest` tag.

If any job fails, fix the cause and publish a new version. Never move or overwrite a tag that users
may already have fetched.

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

## Repository settings

Keep the repository's default Actions token permission read-only. The release workflow grants
`contents: write` and `packages: write` only to the jobs that need them.

Branch protection on `main` blocks force-pushes and deletion but does not require status checks,
so release branches are merged locally and pushed. The Rust workflow runs on every push to `main`
and on pull requests; check its result after pushing rather than waiting on it before merging.

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

The workflow audits the complete lockfile and runs the storage patch compatibility tests before
packaging. Archives include analytical schemas, the changelog, and third-party notices and licenses.

1. On a release branch, set `[workspace.package].version` in `Cargo.toml`. All Sage Plus crates
   inherit this version.
2. Move the `## [Unreleased]` entries in `CHANGELOG.md` into a matching `## [vX.Y.Z]` or
   prerelease section, leaving a new empty `## [Unreleased]` section above it.
3. Update release-specific documentation, such as version examples in `README.md` and `DOCS.md`.
4. Merge the release branch into `main` locally and run the checks there:

   ```shell
   git switch main
   git pull --ff-only origin main
   git merge --no-ff release-branch
   bash scripts/check-release-version.sh
   cargo fmt --all -- --check
   cargo clippy --workspace --all-targets --all-features --locked -- -D warnings
   cargo test --workspace --locked
   cargo build --release --workspace --locked
   ```

5. Push `main`. The Rust workflow runs on the push. Optionally, run `Release Sage Plus` manually
   from the Actions page: a manual run builds and retains all archives and validates the Docker
   build, but does not publish a release or container.

## Publish

Create and push exactly one annotated tag on the release commit on `main`:

```shell
git tag -a vX.Y.Z-beta.N -m "Sage Plus vX.Y.Z-beta.N"
git push origin vX.Y.Z-beta.N
```

The tag starts the release workflow. Prerelease identifiers such as `-beta.1` cause GitHub to mark
the release as a prerelease. Stable versions are marked as the latest release and also update the
container's `latest` tag.

If any job fails, fix the cause and publish a new version. Never move or overwrite a tag that users
may already have fetched. Release checklists from earlier betas are kept in
[`benchmarks/archive/`](benchmarks/archive/).

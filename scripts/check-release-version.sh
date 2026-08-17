#!/usr/bin/env bash
set -euo pipefail

repository_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repository_root"

release_version="$(
  sed -nE 's/^version = "([^"]+)"$/\1/p' Cargo.toml | head -n 1
)"

if [[ -z "$release_version" ]]; then
  echo "Could not read [workspace.package].version from Cargo.toml" >&2
  exit 1
fi

crate_manifests=(
  crates/sage/Cargo.toml
  crates/sage-cli/Cargo.toml
  crates/sage-cloudpath/Cargo.toml
  crates/sage-mcp/Cargo.toml
)

for manifest in "${crate_manifests[@]}"; do
  if ! rg --quiet --fixed-strings 'version.workspace = true' "$manifest"; then
    echo "$manifest must inherit the workspace release version" >&2
    exit 1
  fi
  if ! rg --quiet --fixed-strings 'publish = false' "$manifest"; then
    echo "$manifest must set publish = false for GitHub-only releases" >&2
    exit 1
  fi
done

if [[ $# -gt 0 ]]; then
  release_tag="$1"
  if [[ ! "$release_tag" =~ ^v[0-9]+\.[0-9]+\.[0-9]+(-[0-9A-Za-z.-]+)?$ ]]; then
    echo "Release tag must look like v1.2.3 or v1.2.3-beta.1: $release_tag" >&2
    exit 1
  fi

  tag_version="${release_tag#v}"
  if [[ "$tag_version" != "$release_version" ]]; then
    echo "Release tag $release_tag does not match workspace version $release_version" >&2
    exit 1
  fi

  if ! rg --quiet --fixed-strings "## [$release_tag]" CHANGELOG.md; then
    echo "CHANGELOG.md does not contain a section for $release_tag" >&2
    exit 1
  fi

  unreleased_content="$(
    awk '
      /^## \[Unreleased\]$/ { in_unreleased = 1; next }
      in_unreleased && /^## / { exit }
      in_unreleased && NF { print }
    ' CHANGELOG.md
  )"
  if [[ -n "$unreleased_content" ]]; then
    echo "CHANGELOG.md still contains entries under [Unreleased]" >&2
    exit 1
  fi
fi

printf '%s\n' "$release_version"

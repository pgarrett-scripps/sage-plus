#!/usr/bin/env bash
# Check that the container the release workflow builds can actually be built.
#
# The release workflow does not build the image from the repository: it stages a
# minimal context in `.release/docker` and builds from that. A `COPY` added to
# the Dockerfile for a path the staging step does not copy therefore passes
# every local build and fails only in the release run, after the tag exists.
# That is what happened to v0.1.0-beta.20 (`COPY assets/contaminants/` against a
# staging step that copied Dockerfile, THIRD_PARTY_NOTICES.md and licenses
# only).
#
# This script replays the workflow's own staging commands -- read out of the
# workflow file, so the two cannot drift -- and checks that every `COPY` source
# in the Dockerfile exists in the result.
#
#   scripts/check-release-container.sh                 # staging check only
#   scripts/check-release-container.sh --build PATH    # also docker build, with
#                                                      # PATH as the musl `sage`
#
# The staging check needs nothing but the repository and runs in a second, which
# is why CI runs it. `--build` needs Docker and a statically linked binary; use
# it before tagging, with the musl artifact of the release run.
set -euo pipefail

repository_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repository_root"

workflow=.github/workflows/release.yml
binary=""
if [[ ${1:-} == --build ]]; then
  binary=${2:?--build needs the path to a statically linked sage binary}
  [[ -x $binary ]] || { echo "not executable: $binary" >&2; exit 1; }
fi

# The staging step, verbatim from the workflow.
staging=$(awk '
  /name: Prepare minimal Docker context/ { inside = 1; next }
  inside && /^      - name:/ { exit }
  inside { print }
' "$workflow")
[[ -n $staging ]] || { echo "could not find the staging step in $workflow" >&2; exit 1; }

context=.release/docker
archive=sage-plus-check-x86_64-unknown-linux-musl
rm -rf "$context" .release/input
mkdir -p ".release/input/$archive"
if [[ -n $binary ]]; then
  # install, not cp: cp onto an existing file keeps that file's mode, and a
  # staging-only run leaves a non-executable placeholder here.
  install -m 755 "$binary" ".release/input/$archive/sage"
else
  : > ".release/input/$archive/sage"   # a placeholder: the staging check only needs the path
fi

# Replay every mkdir and cp the workflow runs, with the release tag's archive
# name replaced by the one above and the `tar` of the downloaded artifact
# skipped (the artifact is what `--build` supplies).
while IFS= read -r line; do
  line=${line#"${line%%[![:space:]]*}"}
  case $line in
    mkdir*|cp*) ;;
    *) continue ;;
  esac
  line=${line//\$\{\{ needs.version.outputs.release_tag \}\}/check}
  line=${line//sage-plus-check-x86_64-unknown-linux-musl/$archive}
  eval "archive=$archive; $line"
done <<< "$staging"

missing=0
while IFS= read -r source; do
  [[ -n $source ]] || continue
  if [[ ! -e "$context/$source" ]]; then
    echo "Dockerfile copies '$source', which the release workflow does not stage" >&2
    missing=1
  fi
done < <(awk '/^COPY /{ print $2 }' Dockerfile)
if (( missing )); then
  echo "Add it to the 'Prepare minimal Docker context' step of $workflow" >&2
  exit 1
fi
echo "every Dockerfile COPY source is staged by $workflow"

if [[ -z $binary ]]; then
  exit 0
fi

docker build -t sage-plus-release-check:local "$context"
echo "--- sage --version inside the image"
docker run --rm --network none sage-plus-release-check:local sage --version
echo "--- the bundled contaminant FASTA is readable inside the image"
docker run --rm --network none sage-plus-release-check:local \
  head -c 120 /app/assets/contaminants/universal-contaminants.fasta
echo

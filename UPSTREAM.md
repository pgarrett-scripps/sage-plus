# Upstream relationship

Sage Plus is an independently maintained downstream distribution of
[Sage](https://github.com/lazear/sage). The Git history, MIT license, authorship, and citation
metadata from Sage are intentionally preserved. Sage Plus is not an official Sage release.

The `main` branch is the stable integration branch for Sage Plus. Upstream Sage changes are
merged into it and tested against the full Sage Plus feature set.

Upstream `master` is merged into `main` up to `4a051b2` ("chore: update ci", 2026-10-22), so the
merge base is the upstream tip as of that commit. Check the current state with
`git merge-base main upstream/master` and `git log main..upstream/master`.

Where upstream and Sage Plus solve the same problem, the merge keeps the Sage Plus
implementation and does not carry the upstream one:

| Upstream | Sage Plus keeps |
| --- | --- |
| named PTMs (`44be82a`) | named modifications with explicit sites |
| radical z ions, `Kind::ZRadical` (`cb9934f`) | `z_dot` ions, `Kind::ZDot` |
| per-mod variable modification limits (`2c9922e`) | per-modification `max_count`/`max_total_count` (see [DOCS.md](DOCS.md#static-and-variable-behavior)) |
| timsrust 0.6.7, OpenTFRaw | `sage-plus-tdf` and `sage-plus-raw` from crates.io |
| telemetry and its URL change (`a5eaa01`) | telemetry removed |

Upstream fixes that Sage Plus lacked are taken in full: the LFQ integration-window and alignment
fixes (`75a5300`, `63f10e7`) and total PSM ordering (`df16dec`). One upstream change is
deliberately deferred rather than rejected: `44be82a` also replaces the string-keyed target/decoy
competition in `crates/sage/src/fdr.rs` with a structural peptide key, which depends on upstream's
`ModId` type and is its own task.

To merge upstream:

```shell
git remote add upstream https://github.com/lazear/sage.git
git fetch upstream
git switch main
git merge upstream/master
cargo test --workspace
```

Changes intended for upstream Sage should be developed as narrow branches based on the lowest
clean upstream dependency possible. The Sage Plus integration branch should never be merged into
an upstream pull-request branch.

When publishing scientific work that uses Sage Plus, cite the original Sage paper listed in
`CITATION.cff` and describe the Sage Plus version or commit used.

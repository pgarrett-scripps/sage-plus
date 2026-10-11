# Upstream relationship

Sage Plus is an independently maintained downstream distribution of
[Sage](https://github.com/lazear/sage). The Git history, MIT license, authorship, and citation
metadata from Sage are intentionally preserved. Sage Plus is not an official Sage release.

The `main` branch is the stable integration branch for Sage Plus. Upstream Sage changes are
merged into it and tested against the full Sage Plus feature set.

The current merge base with upstream `master` is `d74024d` ("fix: track charges on deisotoped
ions (#220)", committed 2026-06-05). Later upstream commits are not merged, including per-mod
variable modification limits (`2c9922e`, 2026-08-16), a telemetry URL change (`a5eaa01`), and a
CI benchmarking dispatch (`f929506`). Sage Plus has its own per-modification `max_count` and
`max_total_count` settings instead (see [DOCS.md](DOCS.md#static-and-variable-behavior)). Check the
current state with `git merge-base main upstream/master` and
`git log main..upstream/master`.

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

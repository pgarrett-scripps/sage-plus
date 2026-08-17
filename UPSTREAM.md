# Upstream relationship

Sage Plus is an independently maintained downstream distribution of
[Sage](https://github.com/lazear/sage). The Git history, MIT license, authorship, and citation
metadata from Sage are intentionally preserved. Sage Plus is not an official Sage release.

The `main` branch is the stable integration branch for Sage Plus. Upstream Sage changes are
merged into it periodically and tested against the full Sage Plus feature set:

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

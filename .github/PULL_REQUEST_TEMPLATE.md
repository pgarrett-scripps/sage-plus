## What this changes

<!-- One or two sentences. Link the issue it closes, if there is one. -->

## Why it cannot live upstream

<!--
Sage Plus is a downstream distribution of https://github.com/lazear/sage and keeps its diff
against upstream small. If this change would also suit upstream Sage, send it there instead.
If it cannot (it depends on a Sage Plus feature, reader or output), say why in one line and
delete this comment. For a pure bug fix or a docs change, write "n/a".
-->

## Checklist

- [ ] `cargo fmt`
- [ ] `cargo clippy --workspace --all-targets --all-features --locked -- -D warnings` is clean
      (this is the exact command CI runs)
- [ ] `cargo test --workspace` passes
- [ ] Tests added for new behaviour, in `crates/<crate>/tests/unit/` for unit tests or
      `crates/<crate>/tests/` for integration tests — not inside a production source file
- [ ] Workspace line coverage stays at or above 80% (`cargo llvm-cov --workspace --summary-only`)
- [ ] Public API, CLI flags, configuration keys and output column names are unchanged, or the
      change is described below
- [ ] Every new configuration option defaults to the previous behaviour, so existing
      configurations produce the same search
- [ ] No `[patch.crates-io]` section or path override for `sage-plus-tdf` / `sage-plus-raw` is
      committed
- [ ] `UPSTREAM.md` updated if this adds a place where Sage Plus and upstream solve the same
      problem differently
- [ ] `CHANGELOG.md` updated under `[Unreleased]`, without a release heading or a version bump

## Effect on search results

<!--
Does this change the output of an existing search? If yes, say which fixture or dataset you
compared before and after, and what moved. If no, say "none" and how you checked.
-->

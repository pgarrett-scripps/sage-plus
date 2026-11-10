# Security policy

## Reporting a vulnerability

Report security vulnerabilities privately, through GitHub's security advisories:

**[Open a draft security advisory](https://github.com/pgarrett-scripps/sage-plus/security/advisories/new)**

Please do not open a public issue, and do not post the details in a pull request or discussion,
until a fix is available.

Include, as far as you can:

- what an attacker can do, and what they need (a crafted input file, a configuration, network
  access);
- the input that triggers it — a crafted mzML, mzMLb, Bruker `.d`/TDF, Thermo `.raw` or spectral
  library file, or the smallest excerpt of one that still reproduces the problem;
- the Sage Plus version (`sage --version`) or commit hash, and the operating system;
- the output, panic message or backtrace.

You will get an acknowledgement in the advisory thread. Sage Plus is maintained as research
software by a small group, so there is no guaranteed response time; if a report sees no reply
after two weeks, comment on the advisory.

## Supported versions

Fixes land on `main` and go out in the next release. Older releases are not patched. Sage Plus
is a downstream distribution of [Sage](https://github.com/lazear/sage) and is not an official
Sage release.

If a vulnerability reproduces in unmodified upstream Sage, it affects upstream users too; say so
in the advisory and we will coordinate with the upstream maintainers rather than disclose it
here first.

## Scope

In scope: the `sage` binary and the four workspace crates (`sage-core`, `sage-cli`,
`sage-cloudpath`, `sage-dia`), including the parsers they use for untrusted input files, and the
cloud object-store paths.

Out of scope: vulnerabilities in dependencies with no exploitable path from Sage Plus — report
those to the dependency, and open a normal issue here if the lockfile needs to move. Advisories
that only `cargo audit` surfaces, with no reachable call, are handled as ordinary dependency
updates in public issues and pull requests.

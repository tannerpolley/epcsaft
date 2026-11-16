# Private release maintenance

`tools/refresh-main-wheel` distributes private development wheels and requires
a private repository. Public source releases do not use that tool.

The public repository is the release source. Commit the candidate, then use
`git archive --format=tar --prefix=epcsaft-<version>/ <commit>` with an explicit
list of tracked paths, excluding all `AGENTS.md` and `CLAUDE.md` files,
`.agents/`, `.claude/`, `.serena/`, `analyses/`, `validation/`, `build/` and
`docs/scientific/latex/`. Compress with `gzip -n` and retain the exact bytes.
Do not use the removed `tools/release-snapshot`. Ignored local files are never
inputs to `git archive`; verify that no CC BY-NC packet or analysis input is
tracked or exported. Check each exported packet's `license_or_use_basis` and
run `data/tools/validate_packets.py` in the extracted tree. Packet locks and the
parameter index already describe the public tree; do not regenerate them unless
validation demonstrates a mismatch. Check the archive for secrets, HSL binaries,
journal PDFs and personal paths before publication.

Retain the archive, candidate revision, membership/exclusion report and checks
under the task's `build/issue-<number>/`. Build and test the conda recipe against
that local archive with its SHA-256, outside the home directory. Prepare the
corresponding-source set with the archive, exact dependency sources, Boost
headers, notices, patch and recipe, and verify its inner and outer checksums.
The committed recipe keeps the previous published archive hash until the new
asset exists; the local test recipe and source-set recipe use the candidate hash.
Tags, uploads, releases and Windows probe branches require release authority.

After merges, `tools/refresh-main-wheel` compares retained previous/candidate
prepared environments before packaging main. Select the retained baseline with
`EPCSAFT_PREVIOUS_PREPARED`; the candidate defaults to main’s prepared build tree.
Set `EPCSAFT_NUMERICAL_EVIDENCE` to space-separated existing evidence paths.
Wheels are `auditwheel`-repaired `manylinux_2_39` wheels carrying Ipopt, MUMPS and
Ceres; before retaining one it runs `tools/check-wheel-container.sh`, which refuses
any HSL member and runs `verify-installed-wheel.py` and the Python suite in a
digest-pinned `python:3.13-slim` container without those host libraries.
It retains the wheel in `build/wheels/main-<short-commit>/`, verifies a draft's
uploaded bytes and full commit tag, then publishes the private prerelease
`main-<full-commit>` with `latest=false`; a different identity stops without
replacement. `build/wheels/main-current/` is updated after confirmed publication.
The [downstream guide](../../downstream/README.md) owns retries, authenticated
download and file pins.

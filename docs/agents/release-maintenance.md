# Private release maintenance

`tools/refresh-main-wheel` distributes private development wheels and requires
a private repository. Public source releases do not use that tool.

Maintainers prepare a local snapshot with `python3 tools/release-snapshot` after
committing the candidate. It writes a new tree and a JSON manifest under
`build/release/`, excludes agent configuration, the manuscript (`docs/scientific/latex/`), ignored files and `build/`, and replaces files larger
than 50 MiB with SHA-256 notes. It also excludes every packet whose own
`packet.toml` `license_or_use_basis` is CC BY-NC, and the CC BY-NC analysis
inputs named in the tool, and lists both in the manifest; a packet without a
readable basis stops the snapshot. In the exported tree only, it drops the
excluded packets' locks and regenerates `data/docs/parameter-index.md`, so the
shipped packet validator passes; `data/tools/restore_feos.sh` restores them. The manifest records the HSL and journal-PDF
path checks and personal-path counts. Complete the secret scan before any
publication; snapshot generation does not publish or tag a release.

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
`main-<full-commit>` with `latest=false`. Retries reuse the retained wheel and
resume missing assets in the same draft. An interruption before an upload can
resume; a failed transfer may leave an incomplete `starter` asset. Incomplete,
unequal or unreadable existing assets stop for owner diagnosis; the tool never
deletes or replaces them. A different identity also stops without replacement.
`build/wheels/main-current/` is updated after confirmed publication. The
[downstream guide](../../downstream/README.md) owns authenticated download and file pins.

# Downstream application workspaces

The three official applications consume the Engine. Each application's own thread
owns its repository; the Engine asks for application work or evidence by opening an
issue in that repository and edits no application repository. Engine agents may read
application source in independent clones under this directory. The owner's
2026-09-21 instruction leaves the main local application checkouts outside
`downstream/` unchanged.

| Application | GitHub repository | Input and migration requirements |
|---|---|---|
| Amine-Thermodynamics | [repository](https://github.com/tannerpolley/Amine-Thermodynamics) | `docs/coordination/greenfield-migration.md` |
| MEA-Absorption-Column | [repository](https://github.com/tannerpolley/MEA-Absorption-Column) | `docs/coordination/greenfield-migration.md` |
| Lithium_Extraction | [repository](https://github.com/tannerpolley/Lithium_Extraction) | `docs/coordination/greenfield-migration.md` |

The Engine tracks only this guide, not application source, Git repositories,
environments or build caches. There are no submodules or runtime sibling-source
imports. A fresh Engine worktree therefore contains the guide; when reading needs
application source, clone its GitHub repository into `downstream/<application>`.

The original `codex/greenfield-migration-draft` clones contain retained local edits.
Preserve them and their `copy-provenance.json` record. Their `origin` remotes point
to the main local application checkouts; do not fetch, pull or push those remotes.
Their `github` remote is the GitHub repository; fetch it to read current application `main`.
A fresh direct GitHub clone can use its normal `origin` remote. Do not reset a
retained draft to make it match a newer base.

Read the selected application's `AGENTS.md`, `docs/scientific/CONTEXT.md` and
migration requirements before work. GitHub issues own the live scope, authorization,
dependencies and execution queue; the Engine's
[scientific context](../docs/scientific/CONTEXT.md) owns shared responsibilities.
An input-requirements checkpoint does not authorize numerical-kernel implementation,
application adoption or qualification.

All numerical work consumes one explicitly identified, non-editable Engine wheel.
The private GitHub prerelease `main-<full-40-character-commit>` supplies the wheel,
its `<wheel-name>.sha256` checksum and `release.json`. Each published release,
commit tag and local per-commit copy is retained forever. `tools/refresh-main-wheel`
is the only publisher; an application selects a qualified commit under its own
adoption authority.

Publisher retries reuse the retained wheel and resume missing assets in the same
draft after an interruption before an upload. A failed transfer can leave an
incomplete `starter` asset ([GitHub upload behavior](https://docs.github.com/en/rest/releases/assets#upload-a-release-asset)).
Incomplete, unequal or unreadable existing assets stop for owner diagnosis;
the publisher never deletes or replaces assets.

Authenticate with `gh auth login`, then download into the application's persistent
wheel directory before `uv sync --locked`. Replace the commit below with the
selected full revision. Keep credentials out of dependency URLs and logs.

```bash
commit=<full-40-character-commit>
tag="main-$commit"
wheel_dir="build/engine-wheels/$tag"
mkdir -p "$wheel_dir"
gh release download "$tag" --repo tannerpolley/ePC-SAFT --dir "$wheel_dir" \
    --pattern 'epcsaft-*.whl' --pattern '*.sha256' --pattern release.json
(cd "$wheel_dir" && sha256sum --check --strict ./*.sha256)
gh api "repos/tannerpolley/ePC-SAFT/releases/tags/$tag" \
    --jq '.assets[] | select(.name | endswith(".whl")) | .digest'
gh api "repos/tannerpolley/ePC-SAFT/git/ref/tags/$tag" --jq .object.sha
```

The asset `digest` GitHub reports must equal `sha256:` plus the checksum just
checked, and the returned tag target must equal the selected full commit. Check that
`release.json` records that commit, the downloaded wheel name and its SHA-256.
It also records the Python/platform tags, toolchain identity and paths to existing
numerical evidence; the record itself makes no physical-validation claim. A
failed checksum or mismatched identity stops adoption. Download retries use a new
empty directory and compare bytes before retaining a file; never overwrite an
existing pin with different bytes.

Pin the exact saved filename in the application's `pyproject.toml`, for example:

```toml
[project]
dependencies = ["epcsaft"]

[tool.uv.sources]
epcsaft = { path = "build/engine-wheels/main-<full-commit>/epcsaft-<version>-cp313-cp313-manylinux_2_39_x86_64.whl" }
```

Wheels from #576 on carry their Ipopt, MUMPS and Ceres libraries (`epcsaft.libs/`)
and need glibc 2.39 or newer; HSL is never in the wheel. They are for the owner's
own use only: sharing one needs a licence review first. Earlier `linux_x86_64`
wheels need those libraries installed on the host.

Each application's own thread edits its dependency declaration, regenerates its
lock with `uv lock`, and runs its existing acceptance checks. Subsequent setup
obtains and verifies this exact file before `uv sync --locked`; it does not compile
the Engine. Keep the wheel after installation because `epcsaft.provenance()` reads
and hashes its `file:` installation source. Record the release tag, full commit,
local path and SHA-256 in the application's existing numerical evidence, with
its inputs and claim limits. Never pin `~/.cache/epcsaft` or an Engine task's
scratch directory.

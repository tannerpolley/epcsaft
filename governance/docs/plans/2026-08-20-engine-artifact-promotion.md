# Engine artifact promotion after protected candidate completion

Status: completed after task
`019fdeb0-c2fc-7cc1-8d5e-1ba8166126ac` reported completion and released its
candidate path.

This is a Governance/build cutover only. It changes no Engine equations,
application chemistry, or MEA-Thermodynamics files.

## Subjects

The verified current-main subject is the unified `epcsaft` wheel built from
Engine commit `6d1c434a67edfdad1f29f624311969a499ad83c1`:

| Field | Value |
| --- | --- |
| Distribution | `epcsaft` |
| Version | `0.2.0.dev0` |
| Filename | `epcsaft-0.2.0.dev0-cp313-cp313-linux_x86_64.whl` |
| Bytes | `4659285` |
| SHA-256 | `69466464ae323af3b2702ccf56ca00de67d389cd4b51f87b605cc5d6fce54a21` |
| Engine tree | `f03504ac53b66f1c44eba376da98cf0f539be3de` |

The pre-promotion active subject was the `21ed202d...` wheel from Engine
commit `c712fe9730faef0732341ccdabc4d799fef2d953`. Its immutable object
remains at
`objects/sha256/21ed202df4aa5f8a516aa6400a3a09fe9b5e1d267935f6c5411b92e8faf9d049/`.

The formerly protected non-active candidate was SHA-256
`67fef817d5c48115c40c0cb34d4883add8a5fc31cea4ab558ad57cd101d7dd33` at the
task-owned candidate path under the external artifact store. The task reported
completion before disposition; the bytes are now retained only in the recorded
compressed archive.

## Promotion gate and transaction

After the protected task completes:

1. Read its completed status and re-hash the candidate before disposition.
   If the path is missing or its bytes differ, stop and preserve the current
   active subject.
2. Create a new compressed candidate-disposition archive under
   `artifacts/archive/`, with separate `wheels/` members for the superseded
   active `21ed...` wheel and the protected `67f...` candidate. Verify both
   member hashes before removing the uncompressed candidate path. Record the
   archive hash, byte count, member paths, hashes, distributions, and
   `superseded-replay` retention in `manifests/artifact-dispositions.json`.
3. Verify the current-main wheel hash and copy it into a temporary sibling of
   the sole `active/epcsaft/` alias. Replace the alias with one atomic rename;
   never mutate the old object and never create a second active alias.
4. Update `manifests/artifact-store-policy.json` to bind the admitted
   distribution to commit `6d1c434...` and SHA-256 `69466464...`, then run
   `python governance/tools/artifact_store.py normalize`. This updates both
   the external and tracked one-entry indexes.
5. Require all of the following before installation is reopened: artifact
   store `check`, artifact store `resolve --distribution epcsaft`, management
   preflight, the proportional installed-wheel Engine gate, and a clean
   repository cleanup audit. If any check fails, stop and restore the active
   alias from its preserved `21ed...` object before changing further records.

The candidate archive is provenance only and is never an installation source.
The old active object is preserved for replay provenance; no fallback channel
or alternate active alias is admitted.

## Completed promotion record

The candidate and pre-promotion alias were retained in
`archive/candidate-wheel-dispositions-20260820.tar.zst` (SHA-256
`a7f92963174fbdf9e21bf30b0044aa444f30be100baf94328ceafe5807497171`,
9,192,185 bytes) with these exact members:

```text
wheels/epcsaft-active-20260820/21ed202d/epcsaft-0.2.0.dev0-cp313-cp313-linux_x86_64.whl
wheels/epcsaft-m3-m7-candidate-20260820/67fef817/epcsaft-0.2.0.dev0-cp313-cp313-linux_x86_64.whl
```

The active alias and the new immutable object both hash to
`69466464ae323af3b2702ccf56ca00de67d389cd4b51f87b605cc5d6fce54a21`, and
`artifact_store.py normalize`, `check`, `resolve --distribution epcsaft`,
the installed-wheel verifier, and Management preflight passed.

The focused Bülow replay evaluated all four retained reactive rows exactly:
statuses and solver statuses matched, all compared compositions/densities/
pressures/ARDs matched bit-for-bit, and every numeric difference was `0.0`.
It used unchanged defaults and did not rerun the 24-row neutral sulfolane-CO2
regression or modify retained evidence.

## Post-promotion installed-wheel replay

After the artifact-store and Engine gates pass, run the focused installed-wheel
Bülow 2021 sour-gas reactive replay for its four retained MDEA/CO2 states. The
replay must use the promoted non-editable wheel, the parameter document's
H+/MDEAH+ and OH-/HCO3-/CO3-- species, and its default
`exclude-same-sign-ion-pairs` behavior. The parameter document contains no
explicit same-sign `k_ij` records, so this replay exercises the ionic-pair
exclusion path directly. Compare the four outputs with the retained evidence
under unchanged tolerances and record exact differences if any; preserve the
retained evidence and do not rewrite it merely to refresh hashes. Do not rerun
the 24-row neutral sulfolane-CO2 regression for this promotion because it does
not exercise ionic pair exclusion.

## Downstream lock handoff

The downstream Engine lock must identify the promoted wheel by exact bytes, not
by version alone. Values derived from the verified current-main installed wheel
are:

```text
distribution = epcsaft
version = 0.2.0.dev0
engine_commit = 6d1c434a67edfdad1f29f624311969a499ad83c1
engine_tree = f03504ac53b66f1c44eba376da98cf0f539be3de
wheel_filename = epcsaft-0.2.0.dev0-cp313-cp313-linux_x86_64.whl
wheel_sha256 = 69466464ae323af3b2702ccf56ca00de67d389cd4b51f87b605cc5d6fce54a21
public_api_signature_identity = sha256:6c6ad454b5080780841cabfe843307860e35bcd0a1ae90b9fbd21eede78ae309
public_header_sha256.epcsaft = fef47a5e86f2f1a6e7026b63d227cda116439284b69f8bba2d00589c36c18f61
public_header_sha256.epcsaft.equilibrium = 93c1605bcf94064c3c96f96de0103b25f49102fec8c620804de9a3d7e8b954c4
public_header_sha256.epcsaft.regression = 171527d7a0b360b26130af45174a68d5fc07f532bbaee5c42b835d65fcc59b94
```

The Governance object path after promotion is
`artifacts/objects/sha256/69466464ae323af3b2702ccf56ca00de67d389cd4b51f87b605cc5d6fce54a21/epcsaft-0.2.0.dev0-cp313-cp313-linux_x86_64.whl`.

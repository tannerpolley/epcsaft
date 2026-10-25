# Khudaida historical/current state difference: diagnosed

**Frozen record — not rerunnable as a current-main calculation.** Record snapshot:
commit `ca719012c31e8b41edcb98ae82a388117b58d343`. Original wheel identities, source
inputs and retained numerical evidence keep their recorded provenance. The completed
association-volume diagnosis compares historical states on wheel `684b213a…`; later
corrected-model qualification belongs to the settler leaf.

The current Engine's fixed-feed two-liquid states for the eight Khudaida
Figure 2 feeds (293.15 K, 100 kPa, water + ethanol + isobutanol + Na+ + Cl-)
differ from the retained historical HELD2 rows by up to `0.00807` in mole
fraction, `0.00537 mol` in phase amount and `0.00430` in `G/RT`. The difference is
caused by the model, not the solver. Packet v5/v6 changed three unlike-pair
association volumes.

## Cause

Parameter bundle 2, used for the historical rows, stored the unlike-pair
association volumes precomputed at 293.15 K with the cited Gross–Sadowski
2002 Eq. 3 (Wolbach–Sandler) rule, which includes the size factor

```text
kappa_ij = sqrt(kappa_i kappa_j) * [ sqrt(sigma_i sigma_j) / ((sigma_i + sigma_j)/2) ]^3
```

The packet-v5/v6 translation replaced those numbers with the Engine rule
`arithmetic-energy-geometric-volume`, which is `sqrt(kappa_i kappa_j)` with no
size factor, while keeping the Gross–Sadowski Eqs. 2–3 provenance. The energy
rule (arithmetic mean) is unchanged. Water sigma is the Table 5 correlation at
293.15 K (`2.7990 Å`).

| pair | bundle 2 (historical) | rule with size factor | current v6 rule | change |
|---|---:|---:|---:|---:|
| water/ethanol | 0.03798098 | 0.03798098 | 0.03821013 | +0.60% |
| water/isobutanol | 0.01299624 | 0.01299624 | 0.01327787 | +2.17% |
| ethanol/isobutanol | 0.01119849 | 0.01119849 | 0.01125192 | +0.48% |

## Evidence

[`scripts/diagnose.py`](scripts/diagnose.py) evaluates each historical
state under the current Engine and re-solves the same feed. The
[per-point rows](results/diagnosis.csv) and [summary](results/summary.json)
use wheel `684b213a…`. Maxima over the eight feeds:

| quantity | current v6 model | v6 with size-factor volumes |
|---|---:|---:|
| max \|Δμ/RT\| between phases, neutral species, at historical state | `0.0332` | `4.1e-8` |
| max \|Δ(μ_Na+μ_Cl)/RT\| between phases, at historical state | `3.3e-4` | `3.4e-7` |
| G/RT of historical state minus solved state | `+2.7e-5` to `+3.2e-5` | `1.2e-10` |
| solved state: max \|Δμ/RT\|, neutral / NaCl | `4.7e-10` / `9.3e-10` | `4.8e-10` / `9.6e-10` |
| solved vs historical: mole fraction | `0.00807` | `1.2e-9` |
| solved vs historical: phase amount (mol) | `0.00537` | `8.7e-10` |
| solved vs historical: G/RT | `0.00430` | `4.9e-12` |

The historical rows close the material balance to `1.5e-13`. Under the current
v6 model the historical state is not an equilibrium: chemical potentials differ
between the phases by about 0.03 RT, and its Gibbs energy is higher than the
current solution at every feed. The current solution is the equilibrium of the
v6 model. With the bundle-2 association volumes restored, the current Engine
reproduces the historical HELD2 states to the historical solver's own
stationarity level. It also reproduces the reported total `G/RT` to `5e-12`.
The EOS, reference convention, IDAES mapping and equilibrium algorithm are
therefore not causes. The former `1e-4` lineage gate had no uncertainty
derivation and measured this parameter change, not a solver regression.

Which volume convention is correct is a model-input question. Khudaida et al.
cite Gross–Sadowski 2002 Eqs. 2–3, and the repository's pair-parameter research
note records that the Wolbach–Sandler volume rule carries the size factor. The
historical state is also closer to the data. In the previous version of this
record (commit `4b990757`), the grand formula-basis AAD against Khudaida Table 3
was `0.01105` historical versus `0.01284` current. Against the digitized
organic branch of the paper's own ePC-SAFT curve it was `0.00238` versus
`0.00698`.

## Resolution

Issue #155 found that every cited source uses the size factor (Gross–Sadowski
2002 Eq. 3). Packet v7 names the Engine rule `wolbach-sandler`, which applies
the factor with the state water `sigma(T)`. The [settler reference](../2026-khudaida-idaes-settler/README.md)
is requalified on v7 and reproduces these historical states. The packet-v6
input of this diagnosis is `data/packets/khudaida-2026-figure-2-electrolyte-lle/6`;
#341 removed the byte-identical leaf copy.

## Frozen diagnosis limits

Claim limits: declared two-liquid local branches at these eight feeds only;
no global-stability or predictive-accuracy claim.

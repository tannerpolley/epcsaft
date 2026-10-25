# IAPWS one-state imposed-constant solution verification

Current scope is the native liquid solution of 2H2O = H3O+ + OH− from one mole
water at 298.15 K and 100000 Pa with the source constant imposed on the common
molality convention. This is numerical solution verification. The general
IAPWS correlation and independent solved activity product are not verified;
source molalities, activities and ln Q remain unavailable pending Engine #793.

Run from the repository root with this checkout's own editable Engine:

```bash
env -u EPCSAFT_HSL_LIBRARY nice -n 10 .venv/bin/python analyses/reproduction/2019-iapws/scripts/reproduce.py
cd analyses/reproduction/2019-iapws
QUARTO_PYTHON=../../../.venv/bin/python nice -n 10 quarto render notebook.qmd
mkdir -p ../../../build/issue-789/2019-iapws/render
mv notebook.html ../../../build/issue-789/2019-iapws/render/
```

The producer compiles a small native reader under `build/issue-789/`, solves
through the Engine, evaluates the raw native rows at the returned coordinates,
and contracts the three identified Chemical rows. It does not reconstruct
activities or solve a caller extent equation. `results/comparison.csv` records
every numerical check with its tolerance and evidence label;
`results/diagnostics.json` retains the actual state, species amounts, raw and
classified rows, reference status and cause. Charge and ion equality coincide
for this single reaction; the shared 1e-12 mol check also enforces the accepted
1e-10 mol charge ceiling. Reference terminal change uses 5e-5 as a local
convergence estimate. `notebook.qmd` reads these files without running the model.

The source case and `inputs/water_self_ionization_iapws_r11_07_2019.json` retain
release R11-07(2019), Eqs. 1–4/Tables 1–3, source constant and mass provenance.
The selected `held-cameretti-sadowski-2008/4` declaration replaces only the
admitted water molar mass by 0.018015268 kg/mol; original packet mass and both
identities remain in `results/run-identity.json`. No packet is edited. Ion
masses retain their packet values; no independent exact physical mass claim
follows. Source density is not an EOS target. The older
`results/water-self-ionization.json` remains historical under its own
fingerprint. No pH, caloric or global-stability claim is made. Independent
delivered review remains open.

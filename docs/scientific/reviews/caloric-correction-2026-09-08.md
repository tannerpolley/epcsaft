# Residual internal-energy correction

2026-09-08 — correction verified.

The owner authorized correcting the missing temperature factor immediately.
The state result now follows `residual_internal_energy` in `engine/docs/equations.md`:
$u^{\mathrm{res}}=-RT^2a_T$. Adjacent entropy, enthalpy and heat-capacity
expressions were inspected against the canonical identities without changing them.
The existing mixture checks additionally verify $u^{\mathrm{res}}=h^{\mathrm{res}}-RT(Z-1)$
and $u^{\mathrm{res}}=RTa^{\mathrm{res}}+Ts^{\mathrm{res}}$.

This bounded correction does not qualify M1 or its prerequisites.

Validation: `python3 tools/run-numerical-evidence.py F1` passed; the complete
native `model_checks` executable passed, including neutral, associating,
electrolyte and nine-species mixture identities. `bash tools/check-all.sh`
passed against the rebuilt non-editable wheel. No tolerance was changed.

The F1 state value and analytic reference both equal
-6854.895753759838 J/mol; the enthalpy identity agrees within its existing bound.
All four temperature-ladder comparisons pass. The live companion reads the
fresh log directly. F5 is a separate Hessian defect and is not corrected here.

Wheel SHA-256: `c9bc7c6926b47a29bcf38c834aaf402891d957ce6f39909505e6eac630c431d4`.

Source SHA-256: `f87c5e2eafc732ea1cf40875b6a5520b38ba88badeac96b91be2f00b0c7e1742`.

Native executable SHA-256: `36867a1984c542a59328afa1d4d119754e232e0eab1800e42e49e5a9ae7089f4`.

Independent Luna High review found no blocker: the corrected expression and
adjacent caloric formulas agree with their canonical equations;
`state_at_pressure` inherits the fix through `state`. Numerical execution was
verified by the lead.

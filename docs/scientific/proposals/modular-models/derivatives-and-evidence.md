# Shared derivatives and evidence requirements

[Overview](README.md) · [Model specification](model-specification.md) · [Roadmap](roadmap.md)

## What becomes automatic

Use a common reduced molar basis $a^{res}=A^{res}/(nRT)$, with $n$ the
total amount in moles. A general representation is

$$
a^{res}(\mathbf q;\mathcal M,\boldsymbol\theta)
=\sum_c a_c(\mathbf q,\mathbf z(\mathbf q);\boldsymbol\theta_c),
\qquad
\mathbf R(\mathbf z,\mathbf q;\mathcal M,\boldsymbol\theta)=0.
$$

$\mathcal M$ is the specification; $\boldsymbol\theta$ contains parameters;
$\mathbf q$ represents the independent state and active parameter coordinates.
Temperature, volume and component amounts form one state chart; composition
directions must respect the selected constraint convention. Internal variables
$\mathbf z$ can include site fractions, screening and pairing quantities.

On a smooth branch with nonsingular $\mathbf R_z$,

$$
\frac{d\mathbf z}{d\mathbf q}=-\mathbf R_z^{-1}\mathbf R_q,
\qquad
\frac{da^{res}}{dq}=\sum_c\left[
\frac{\partial a_c}{\partial q}
+\frac{\partial a_c}{\partial\mathbf z}\frac{d\mathbf z}{dq}\right].
$$

The inverse is mathematical notation; evaluation can solve the linear system.
Higher actions differentiate the coupled residual equations again, including
the derivatives of their Jacobians and all shared constitutive relations.
These identities explain requirements; the accepted algorithms remain the owner
of implementation methods, tolerances and failure behavior.

For example,

$$
\mathbf X\longrightarrow\epsilon_r(T,\rho,\mathbf x,\mathbf X)
\longrightarrow a^{DH},a^{Born}
$$

connects association to dielectric response and both ionic terms. Freezing
$\epsilon_r$ while differentiating those terms misses part of the model.
Temperature-dependent diameters and pair parameters need the same treatment.
Silva also shows that introducing concentration dependence at different stages
of the electrostatic derivation produces different models. Complete derivatives
of one chosen expression do not make it equivalent to another derivation.
[Silva2023](references.md#silva2023)

## Pressure closure adds an implicit dependency

At prescribed pressure, density is a root of

$$
F(\rho_m,\mathbf u)=p(T,\rho_m,\mathbf x,\boldsymbol\theta)-p_{set}=0,
$$

where $\mathbf u$ contains the remaining independent coordinates. On a
regular selected branch,

$$
\frac{d\rho_m}{du}=-\frac{F_u}{F_{\rho_m}},
\qquad
\frac{dg}{du}=g_u+g_{\rho_m}\frac{d\rho_m}{du}.
$$

Association and dielectric dependencies are already inside these total state
derivatives. Near a spinodal or singular root, $F_{\rho_m}$ can vanish;
finite actions and unique continuation cannot simply be assumed. Branch
identity, achieved residual and conditioning belong in evidence and diagnostics.
Pressure-closed finite differences must resolve the same branch at perturbed
inputs. Fixed-density evidence does not establish pressure-closed derivatives.

For activity coefficients, derivatives must be taken from the extensive
$A^{res}=nRTa^{res}$, not from $a^{res}$ alone. Fugacity, electrolyte
reference subtraction and concentration-scale conversion must share that same
thermodynamic state and convention. The comparison document explains why this
matters to MIAC interpretation.

## What remains supplied or qualified

Once dependencies are correct, Properties and the one equilibrium driver can
consume the assembled model without separate paper-specific solvers. This
does not supply missing ideal thermal correlations, chemical standard states,
transport properties or reaction chemistry. Missing ideal calorics should leave
valid residual calculations available with explicit total-property limitations.

Automatic differentiation establishes derivatives of implemented equations.
It does not prove source correctness, implicit-root regularity, branch choice,
physical accuracy or validity outside calibration. Numerical and physical
acceptance must remain distinct.

## Numerical qualification matrix

| Claim | Smallest useful falsifier or reference |
|---|---|
| Correct formulation | Independent equation values and justified source/analytic limits; corrected coefficients and unit conversions. |
| Correct association | Mass-action residuals; no-bond limit $X=1$, $a^{assoc}=0$; analytic single-pair case; coupled mixtures with declared sites and edges. |
| Correct internal derivatives | Implicit derivative identities and independent directional differences on regular branches, with step-size/error assessment. |
| Correct shared dielectric chain | Nonzero state/parameter derivatives propagated through dielectric-dependent ionic energies; selected constant-input limits. |
| Correct pressure closure | Pressure residual, branch consistency and derivative actions of pressure-closed observables. |
| Correct ionic limits | Neutral limit, dilute electrolyte behavior and selected Born limiting relations within their mathematical domains. |
| Correct composition and parameter actions | Declared independent directions and active coordinates, including mixed actions; identify exact-zero boundary limitations. |
| Correct higher derivatives | Independently supported actions at each requested order; first/second-order agreement does not qualify third order. |
| Correct equilibrium integration | Applicable balance, stationarity and solved-state action checks through the same driver, with stability/search claims bounded separately. |

Choose tolerances from source precision, numerical conditioning and intended
observables before execution. This proposal introduces no arbitrary universal
tolerance. Persistent tests must falsify numerical claims; prose checks or
API snapshots cannot establish thermodynamic acceptance.

## Physical validation matrix

| Axis | Observations and comparisons | Limit to record |
|---|---|---|
| Temperature | MIAC curvature, density, dielectric response and available caloric quantities | State which observations were fitted and which temperatures were withheld. |
| Density/pressure | Pure-solvent density, salt-induced density increment, compressibility and pressure-dependent activities | Neutral-solvent error and ionic increment error can have different causes. |
| Solvent composition | Pure endpoints, water–methanol/ethanol mixtures, solvation/transfer and MIAC | Mixture refits are not parameter-free transfer predictions. |
| Salt concentration | Dilute limit through concentrated solutions within a stated phase | Near-solubility data are not automatically multiphase or solid-equilibrium validation. |
| Chemical diversity | Independent salts, solvents and, where needed, multivalent/reactive species | Monovalent alkali halides do not establish extraction-complex or amine chemistry. |
| Pairing mechanism | Pairing-sensitive observations alongside bulk thermodynamics | A good MIAC fit alone does not identify pair fractions. |

Compare models using common observations, reference states and residual/error
definitions. Retain property-resolved deviations rather than only one aggregate
objective; an MIAC improvement can conceal worse density or permittivity.
Record parameter uncertainty and identifiability where supported; do not invent
uncertainties absent from a paper.

## Existing association evidence and its exact boundary

At inspected snapshot `0ae7806`, the retained
[association notebook](../../../../analyses/2026-association-native-comparison/notebook.qmd)
compares fixed-depth log-Newton candidates at 8, 16, 24 and 32 corrections with
an independent implicit reference on 110 frozen synthetic neutral states.
The same 107 cases qualify at each depth: 428 of 440 slots, with 12 unknown
because three reference solves exhaust their budget. There are zero false
criteria, 52,429 true criteria and 576 finite-difference rows in the retained
55,309-row comparison. Unknown reference outcomes are not candidate successes.

This supports a bounded fixed-volume numerical comparison, not production
replacement, all-mixture performance, ionic-site networks, pressure closure,
equilibrium or third derivatives. Higher depth adds no qualification gain on
these inputs. The fixed-depth method changes numerical approximation and
differentiation; Yang-style site labels change physical topology. Neither claim
may inherit the other's evidence without testing the affected inputs.

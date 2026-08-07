# PC-SAFT/ePC-SAFT pair-parameter symmetry and meaning

**Evidence date:** 2026-08-07
**Scope:** Standard PC-SAFT and the ePC-SAFT formulations cited below. “Symmetric” means the stored coefficient is invariant under exchanging the two component (or site-endpoint) labels; it does not mean the association topology has no donor/acceptor roles.

## Conclusion

| Family | Physical meaning | Scope and symmetry in the cited formulations | Temperature status | Formulation status |
| --- | --- | --- | --- | --- |
| \(k_{ij}\) | Corrects the Berthelot geometric-mean **cross-dispersion energy**: \(u_{ij}=\sqrt{u_i u_j}(1-k_{ij})\). | One scalar for an unlike **component pair**; not a site-pair or directed quantity. Gross and Sadowski explicitly print \(k_{ij}=k_{ji}\). | Original PC-SAFT uses a scalar \(k_{ij}\). ePC-SAFT may define \(k_{ij}(T)\); this is an empirical formulation choice, not an inherent property. | Core PC-SAFT parameter. |
| \(l_{ij}\) | Corrects the Lorentz arithmetic-mean **cross segment diameter/size**: \(\sigma_{ij}=\tfrac12(\sigma_i+\sigma_j)(1-l_{ij})\). | One scalar for an unlike **component pair**; standard use is unordered and not site-directional. | No \(l_{ij}(T)\) rule was found in the cited primary formulations. Treat it as a fitted constant unless a source defines a correlation. | ePC-SAFT extension; \(l_{ij}=0\) is the uncorrected Lorentz rule. |
| \(k_{ij}^{\mathrm{HB}}\) (also \(k_{ij}^{hb}\), \(k_{ij\_hb}\)) | Corrects the Wolbach–Sandler **cross-association energy**, not the association volume: \(\varepsilon^{A_iB_j}=\tfrac12(\varepsilon^{A_iB_i}+\varepsilon^{A_jB_j})(1-k_{ij}^{\mathrm{HB}})\). | A component-pair scalar applied to the selected association site rule \(A_iB_j\). The cited rule uses one value for the pair and a symmetric pure-energy mean; standard implementations therefore store it reciprocally. Site roles/topology can still distinguish \(A_iB_j\) from another site pairing. It is not an independently directed \(i\!\to\!j\) coefficient. | No \(k_{ij}^{\mathrm{HB}}(T)\) rule was found in the cited primary formulations; the cited values are constants, usually zero except for selected pairs. | A source-specific ePC-SAFT combining-rule correction, absent from original Gross–Sadowski PC-SAFT. |

Therefore, for the standard formulations reviewed, all three are best represented as **unordered component-pair coordinates**. Only association topology/site identity is potentially role-specific. A separate \(ij\) versus \(ji\) value would be a new asymmetric model convention requiring its own source and equations, not an implication of the subscripts.

## Primary-source evidence

1. **Cross dispersion, \(k_{ij}\).** Gross and Sadowski, *Ind. Eng. Chem. Res.* 40 (2001) 1244–1260, p. 1247, Eqs. (23)–(24), define the Lorentz diameter and
   \[
   \epsilon_{ij}=\sqrt{\epsilon_i\epsilon_j}(1-k_{ij}),
   \]
   introducing one binary coefficient to correct unlike-chain segment interactions. Their §3.3 describes one \(k_{ij}\) for any binary component pair. Table 4 (p. 1252) reports one fitted value per binary system over a stated temperature range; Eq. (24) itself has no temperature argument. DOI: [10.1021/ie0003887](https://doi.org/10.1021/ie0003887).

   Gross and Sadowski, *Ind. Eng. Chem. Res.* 41 (2002) 5510–5515, Table 2 footnote, explicitly states \(k_{ij}=k_{ji}\) for \(i\ne j\). The paragraph following their cross-association Eqs. (2)–(3) (p. 5511) also says that no additional binary association correction is introduced: only \(k_{ij}\) corrects dispersion in that formulation. DOI: [10.1021/ie010954d](https://doi.org/10.1021/ie010954d).

2. **Cross size, \(l_{ij}\), and hydrogen-bond correction, \(k_{ij}^{\mathrm{HB}}\).** Held et al., *Chem. Eng. Res. Des.* 92 (2014) 2884–2897, §2.3, p. 2887, Eqs. (15)–(18), gives
   \[
   \sigma_{ij}=\tfrac12(\sigma_i+\sigma_j)(1-l_{ij}),\qquad
   u_{ij}=\sqrt{u_i u_j}(1-k_{ij}),
   \]
   and
   \[
   \varepsilon^{A_iB_j}=\tfrac12(\varepsilon^{A_iB_i}+\varepsilon^{A_jB_j})(1-k_{ij}^{hb}).
   \]
   Their Eq. (18) separately combines the association volume using the Wolbach–Sandler geometric-mean/size factor; it does not identify \(k_{ij}^{hb}\) with that volume. The paper’s nomenclature labels \(k_{ij}\), \(l_{ij}\), and \(k_{ij}^{hb}\) as binary interaction parameters. Table 1 (p. 2888) gives \(l_{ij}=-0.0044\) and \(k_{ij}^{hb}=0.026\) for water/1-butanol. DOI: [10.1016/j.cherd.2014.05.017](https://doi.org/10.1016/j.cherd.2014.05.017).

3. **Temperature dependence and zero defaults.** Ascani, Sadowski, and Held, *J. Chem. Eng. Data* 67 (2022) 1972–1984, Supporting Information, Eqs. (S4)–(S7), writes \(u_{ij}=\sqrt{u_i u_j}(1-k_{ij}(T))\) and explicitly defines
   \[
   k_{ij}(T)=k_{ij,298.15}+k_{ij,T}(T-298.15).
   \]
   The same SI uses \(l_{ij}\) in the cross-diameter rule and \(k_{ij}^{hb}\) in the cross-association-energy rule. It states that \(l_{ij}\) and \(k_{ij}^{hb}\) are usually zero and Table S3 lists one constant value for each selected water/solvent pair. DOI and official SI entry: [10.1021/acs.jced.1c00866](https://doi.org/10.1021/acs.jced.1c00866).

4. **Association-site meaning.** Chapman et al., *Ind. Eng. Chem. Res.* 29 (1990) 1709–1721, §2, defines site-site association as localized, orientation-dependent interaction. It distinguishes association energy \(\epsilon^{AB}\) (square-well depth) from association volume \(\kappa^{AB}\) (bonding volume) and states that \(A_iB_j\) denotes site A on molecule \(i\) interacting with site B on molecule \(j\). DOI: [10.1021/ie00104a021](https://doi.org/10.1021/ie00104a021). Gross and Sadowski 2002 Eqs. (2)–(3) apply this distinction to cross-association: arithmetic-mean energy and a separate volume rule.

## Directionality and notation boundary

The \(A_iB_j\) superscript makes the **association site rule** role-aware, but it does not make the scalar \(k_{ij}^{\mathrm{HB}}\) a directed coefficient. Different site pairings can be separate physical edges in a multi-site model; for one edge, exchanging its two endpoints should preserve the pair value. Thus “site-pair directional” describes topology/role selection, not the usual symmetry of the numeric correction.

Do not confuse \(k_{ij}^{\mathrm{HB}}\) with \(\kappa^{A_iB_i}\) (or \(\kappa^{A_iB_j}\)), the association-volume parameter. Some software/data conventions shorten the latter to a **kHB-like** name. The equation and whether the symbol has component indices are authoritative: \(k_{ij}^{\mathrm{HB}}\) multiplies the cross energy by \(1-k_{ij}^{\mathrm{HB}}\); \(\kappa\) multiplies the association bonding volume.

## Current Engine anchor

The current Engine schema matches the standard symmetry for \(k_{ij}\) and \(l_{ij}\): **PairParameterRecord** canonicalizes component IDs (engine/src/epcsaft/records.py:395-424) and **_pair_matrix** writes the same value to both matrix orientations (engine/src/epcsaft/_resolve.py:168-202). Association energy and volume are separate endpoint-pair families; association endpoints are canonicalized and resolved to reciprocal runtime entries (engine/src/epcsaft/records.py:461-502, engine/src/epcsaft/_resolve.py:205-258). There is currently no first-class \(k_{ij}^{\mathrm{HB}}\) family or native correction array. The native association strength uses separate association energy/volume and the cross size (which includes \(l_{ij}\)) (engine/src/epcsaft/_native/association.hpp:75-83).

## Unresolved facts

- The cited papers do not establish a universal prohibition on asymmetric empirical coefficients outside their own combining rules. Any asymmetric *ij/ji* data must therefore be treated as a separately specified formulation, not silently canonicalized.
- Held et al. and Ascani et al. do not print separate identities \(l_{ij}=l_{ji}\) or \(k_{ij}^{hb}=k_{ji}^{hb}\). Their unordered classification follows the single-value binary-pair equations and the reciprocal pair convention; it is not a universal theorem for arbitrary extensions.
- No cited primary source supplies a universal \(l_{ij}(T)\) or \(k_{ij}^{\mathrm{HB}}(T)\) correlation. The effective cross diameter or association strength can still vary with temperature through other terms, such as pure-component diameters and the association Boltzmann factor.
- The cited \(k_{ij}^{\mathrm{HB}}\) rule is a single component-pair correction for the stated \(A_iB_j\) cross rule. The sources do not determine whether a future multi-site formulation should expose independent \(k_{ij,AB}^{\mathrm{HB}}\) values for every site-type pair.

At 298.15 K and a fixed density of 30000 mol/m³, the manufactured
Li–Mg–Cl–SO₄ state uses 1 kg water and solute amounts (0.2, 0.4, 0.6, 0.2) mol.
The [native reference](../../build/issue-814/reference.csv) gives
$F=\Delta A^{mix}/(RT)=-0.034128$ mol, equal to the hand-derived value,
and $a^{mix}=-0.00059970085135365185$.
The largest scaled exact gradient and Hessian errors are
2.78e-17 and 2.41e-17;
the full residual on-minus-off errors are at most
3.54e-14 and 1.15e-15.
For the zero-sum, charge-neutral direction $(0,1,0,-3,2)$, the first action
is −0.1764 and the second is $-0.0576N$.
Isolated temperature, volume, density and spectator derivatives are zero.
Central differences verify the exact columns; their finest scaled gradient
and Hessian errors are 1.92e-09 and
9.59e-09. Adding sulfate alone leaves $F$
unchanged; setting chloride to zero retains $F=-0.03528$ mol and a nonzero
chloride insertion derivative of 0.00192.

The [chloride comparison](../../build/issue-814/qualification-comparison.json)
retains all 96 predictions and 4,665 numeric values bit-identically,
including signed zeros; both Yao CSV files are byte-identical and each scored
bucket retains 455/455 passing realizations. No coefficients were refitted.
The exact #795 sulfate declaration is admitted with the original Li–Mg–Cl record.
These are numerical verifications of the unchanged equation and admission.
Anion-independent $\theta$ follows Pitzer and Kim (1974), Eqs. 8, 9 and 11;
undeclared $\psi$ values, including sulfate's, are zero as a model assumption.
The $E\theta$ terms remain excluded. Transferring chloride-fitted coefficients
to sulfate is unvalidated; #795 owns its screen and any fit. Independent
delivered review remains with the lead.

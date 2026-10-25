# Retained numerical-allowance corrections

The completed fitted population is preserved at commit `fbed2062`. The pre-cold-scaling allowance records are preserved at `1989de13`. Neither correction below reran a fit or changed the statistical or numerical rules.

The complete-column correction replaces an output-only column norm with the norm of the complete augmented Jacobian, including its exact input-prior rows. Prior-row errors remain zero. The family estimate is consumed only where the inner action diagnostics are absent. This defect was found before a bias summary was inspected. The preceding endpoint and action records are retained beside their replacements with `before-complete-norm-` names.

The cold-column correction applies the ordinary fit's relative-weight square root to the public point evaluator's unweighted optimizer columns before comparing them with the fit's weighted columns. The uncertainty-mode evaluator already standardizes its columns. This defect was found after the first bias summary was inspected. Every earlier endpoint, row error and raw reference diagnostic is retained with a `before-cold-scaling-` name. The native `allowance` entry reconstructs the retained point, exact Jacobian and SVD, evaluates the cold point and recomputes the components without fitting.

An integrity comparison against `fbed2062` found no changed fields in 528 fitted-value records, 264 augmented Jacobians, 264 reconciliation records, 528 work records, 528 diagnostics, 528 failure tables, 528 termination descriptions, two paired-error tables and two signed-bootstrap-mean tables. Only CSV line endings and the additional source-record identifier are ignored in this comparison. Timing values are unchanged too. Its detailed record is `build/issue-510/fit-record-integrity-final.json`.

Source-record identifiers distinguish independent calculations or historical reporting versions whose numerical values happen to agree. Each identifier is the path of that retained table within this leaf. They do not alter numerical fields. Earlier byte representations are preserved in Git.

The final frozen decision passes for MEA and remains inconclusive for ionic. The terminal-reference value term dominates the ionic numerical allowance. No noise, draw, seed, input domain, replicate count, reference definition, fit control, statistical threshold or numerical threshold was tuned after comparison.

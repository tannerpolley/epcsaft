# Stability-return comparison (#750)

Decision **(b)**: stability trial order is not the identity of a thermodynamic minimum. Commit `5048b6d4d` ([PR #654](https://github.com/tannerpolley/ePC-SAFT/pull/654)) replaced the neutral tangent-plane trial search with the native Michelsen search. Published wheel `cc63a135` retains the old ordering; `1cd3ab15` and current main return the new ordering. Both find the same liquid and vapor minima. The central equilibrium result and tangent plane remain within their existing limits.

The comparison now pairs each trial composition with its tangent-plane distance and orders those pairs by composition before applying the unchanged checks. It preserves the number of minima, every composition, and each tangent-plane distance. [The bounded comparison](stability-comparison.json) changes from 0.6814058275270192 by list position to 2.4066920563770346e-11 by composition, within the existing 1e-10 limit.

The complete rerun passes every benchmark case. [rows.csv](rows.csv), [summary.md](summary.md), [environment.json](environment.json), and [the producer log](density-rerun.log) retain its numerical returns, accuracy decisions and timings. Timing values describe this shared-host rerun; no speed claim or physical validation is made. Original numerical references are unchanged.

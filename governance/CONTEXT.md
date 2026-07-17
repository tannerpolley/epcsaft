# Organization Policy Repository Context

`GOVERNANCE.md` revision 2 is the canonical local ecosystem doctrine.

`governance_source_of_truth: true`

`published_to_github: true`

`runtime_source_of_truth: false`

During the active migration, `../ePC-SAFT-migration/MIGRATION.md` owns program
sequencing and its accepted receipts record per-slice authority transfers.
Receipts `promotion-0012-provider-full-candidate` and
`promotion-0015-provider-native-sdk-v1` assign their provider scopes to
`ePC-SAFT/ePC-SAFT`; `promotion-0018-equilibrium-pure-saturation-v1` assigns
its equilibrium scope to `ePC-SAFT/ePC-SAFT-equilibrium`; and
`promotion-0020-regression-methane-saturation-v1` assigns its regression scope
to `ePC-SAFT/ePC-SAFT-regression`. The lab remains the sole transitional
runtime authority for every unpromoted slice.

Consumer Maturity Slice 2 is migration-owned program context, not permanent
organization policy. Exact local and remote commit identities are captured at
each execution gate rather than frozen in this long-lived context file.

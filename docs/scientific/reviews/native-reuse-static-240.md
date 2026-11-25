# Static storage and callbacks in #240

The entries refer to the [base inventory](native-reuse-8078d7d1.md) at
`8078d7d1`. A confined entry is used within its owning call, point or synchronized
context. A false positive stores registration state, a declaration constant or
diagnostic bookkeeping.

| Kind | Entry | Disposition | Reason |
| --- | --- | --- | --- |
| Callback | `state_field` | False positive | The binding captures a pointer to a member, then reads each state supplied by the caller. |
| Callback | `HelmholtzModel::evaluate_locked::complete_hessian` | Confined | Completion takes the live producer explicitly and checks its identity and the destination gradient and association point. |
| Callback | `numbers` | Confined | The returned field reader borrows the current ASL parser or writer for a synchronous record traversal. |
| Callback | `wholes` | Confined | The returned field reader borrows the current ASL parser or writer and owns its fixed size limit. |
| Static | `nanobind__core_module` | False positive | Python module registration metadata. |
| Static | `nanobind__core_slots` | False positive | Python module lifecycle handlers. |
| Static | `callback_error` | False positive | Thread-local text reports the latest callback error. |
| Static | `equilibrium_observable::prepared` | Confined | Exact serialized record lookup and the callback mutex scope model construction and subsequent exact-q rows. |
| Static | `pressure_property::prepared` | Confined | Exact serialized record lookup, callback mutex, branch, q, live mask, order and row failure checks scope the read. |
| Static | `cppad_parallel_enabled` | False positive | CppAD registration flag set before active evaluations. |
| Static | `cppad_slot_exit` | False positive | Thread-exit notification condition variable. |
| Static | `cppad_slots` | False positive | Mutexes reserve CppAD thread slots through thread exit. |
| Static | `cppad_slots_mutex` | False positive | Serializes thread-slot reservations. |
| Static | `cppad_initialized::once` | False positive | One-time CppAD runtime initialization. |
| Static | `PressureSeries::residual::coefficient` | Confined | The current association tangent or factorization overwrites the coefficient buffer before use. |
| Static | `PressureSeries::residual::cross` | Confined | Each entry is assigned from the current composition and site residual before the dot product. |
| Static | `PressureSeries::residual::delta` | Confined | Recomputed from the current density path and geometry before any reads. |
| Static | `PressureSeries::residual::factorization` | Confined | Double evaluations regenerate the factorization before use. A wider evaluation owns its factorization with its refined site point. |
| Static | `PressureSeries::residual::residuals` | Confined | Cleared and assigned from the current mass-action equations at each consumed order. |
| Static | `PressureSeries::residual::rows` | Confined | The current association residual assigns each right-hand-side entry before solving. |
| Static | `PressureSeries::residual::site_series` | Confined | Cleared and rebuilt from the current site solution before the coefficient sweeps. |
| Static | `PressureSeries::residual::strength` | Confined | Cleared and rebuilt from the current density series. |
| Static | `PressureSeries::residual::unit` | Confined | Reassigned from the current factorization before estimating the density floor. |
| Static | `active_trace` | False positive | Selects the current thread's trace sink. |
| Static | `AssociationObservations::current::active` | False positive | RAII associates numerical diagnostics with the current evaluation. |
| Static | `CalculationWork::current::active` | False positive | RAII selects the current work allowance and counters. |
| Static | `property_name::family` | False positive | Immutable property-family names. |
| Static | `property_name::kind` | False positive | Immutable property-kind names. |
| Static | `ipopt_iteration_log::log` | False positive | RAII selects the current iteration log sink. |
| Static | `outside_domains::variables` | False positive | Immutable names for input bounds. |
| Static | `policy::ENTRIES` | False positive | Immutable pointers to declared numerical controls. |
| Static | `interrupted_error` | False positive | Captures a pending Python exception through worker cleanup. |
| Static | `compute_cpu_directory` | False positive | Optional directory for fit CPU counters. |
| Static | `extension_image` | False positive | File descriptor for the native worker's extension image. |
| Static | `fit_host` | False positive | Fit worker launch and interrupt configuration. |

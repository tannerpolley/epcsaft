Native states and references
============================

.. py:currentmodule:: epcsaft._core

.. autoapimodule:: epcsaft._core

Native enum introspection
-------------------------

For every native enumeration, ``name`` is its symbolic string and ``value`` is
its numeric integer value, not a physical quantity.

.. autoapidata:: epcsaft._core.GAS_CONSTANT_J_PER_MOL_K

.. autoapidata:: epcsaft._core.EOS_STANDARD_STATE_ID

.. autoapiclass:: epcsaft._core.ParameterError
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.EosError
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.DensityStructureError
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.ConvergenceError
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.DeclarationError
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.Stop
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.Category
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.Tolerance
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapifunction:: epcsaft._core.tolerances

.. autoapiclass:: epcsaft._core.Convergence
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.Event
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.Span
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.Trace
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.Capture
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.ParameterEntry
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.ParameterDocument
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.Parameters
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapifunction:: epcsaft._core.build_parameters

.. autoapifunction:: epcsaft._core.prepare_thermochemistry

.. autoapifunction:: epcsaft._core.fitted_document

.. autoapiclass:: epcsaft._core.IdealPolynomial
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.IdealShomate
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.IdealNASA7
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.IdealNASA9
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.IdealDIPPR107
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.IdealCorrelation
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.IdealInterval
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.ThermochemistryRecord
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.ReactionLogPolynomial
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.ParameterFamily
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.ActiveParameter
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.DerivativeOrder
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.NeutralProperty
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.NeutralPropertyDerivatives
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.PropertyObservable
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.PropertyActionStatus
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.ParameterActionDirection
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.DensityActionDirection
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.PressureActionDirection
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.PropertyActionResult
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.DensityBranch
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.ThermochemistryStatus
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.PropertyActionDiagnostics
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.PropertyActionBatch
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.State
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.DensityDiagnosticPoint
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.DensityDiagnosticBracket
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.DensityDiagnostics
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.SweepCounters
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.HelmholtzModel
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapifunction:: epcsaft._core.property_asl_record

.. autoapiclass:: epcsaft._core.InfiniteDilutionReference
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.SolventActivity
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapiclass:: epcsaft._core.HenryConstant
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapifunction:: epcsaft._core.henry_constants

.. autoapifunction:: epcsaft._core.solvent_reference

.. autoapifunction:: epcsaft._core.solvent_activity

.. autoapiclass:: epcsaft._core.ElectrolyteActivity
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __int__, __enter__, __exit__

.. autoapifunction:: epcsaft._core.infinite_dilution_reference

.. autoapifunction:: epcsaft._core.electrolyte_activity

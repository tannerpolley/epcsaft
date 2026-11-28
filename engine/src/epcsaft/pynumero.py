"""Opt-in derivative deferral for Engine equilibrium callbacks in PyNumero 6.10.1."""

import os
from typing import Any

import numpy as np
from pyomo.common.env import CtypesEnviron  # type: ignore[import-untyped]
from pyomo.common.tempfiles import TempfileManager  # type: ignore[import-untyped]
from pyomo.contrib.pynumero.asl import AmplInterface  # type: ignore[import-untyped]
from pyomo.contrib.pynumero.interfaces.pyomo_nlp import PyomoNLP  # type: ignore[import-untyped]
from pyomo.solvers.amplfunc_merge import amplfunc_merge  # type: ignore[import-untyped]

_EQUILIBRIUM_FUNCTIONS = {
    b"epcsaft_equilibrium_phase_component_amount",
    b"epcsaft_equilibrium_solid_amount",
    b"epcsaft_equilibrium_system_enthalpy",
}


class _AmplInterface:
    def __init__(self, exact: Any, value: Any, first: Any) -> None:
        self._exact = exact
        self._constraints = np.empty(exact.get_n_constraints(), dtype=np.float64)
        self.eval_f = value.eval_f
        self.eval_g = value.eval_g
        self.eval_deriv_f = first.eval_deriv_f
        self.eval_jac_g = first.eval_jac_g

    def __getattr__(self, name: str) -> Any:
        return getattr(self._exact, name)

    def eval_hes_lag(self, x: Any, lam: Any, hes_lag: Any, obj_factor: float) -> None:
        # PyNumero's sphes uses the exact reader's current point, as set by objval/conval.
        self._exact.eval_f(x)
        self._exact.eval_g(x, self._constraints)
        self._exact.eval_hes_lag(x, lam, hes_lag, obj_factor)


class EquilibriumPyomoNLP(PyomoNLP):  # type: ignore[misc]
    """Use separate ASL readers for equilibrium values, first actions and exact Hessians."""

    _asl: Any

    def __init__(self, pyomo_model: Any, nl_file_options: Any = None) -> None:
        self._original_amplfunc = os.environ.get("AMPLFUNC", "")
        super().__init__(pyomo_model, nl_file_options)

    def _collect_nlp_structure(self) -> None:
        super()._collect_nlp_structure()
        # The Pyomo writer puts F declarations just after the ten-line ASCII NL header.
        with open(self._nl_file, "rb") as nl:
            header = b"".join(nl.readline() for _ in range(10))
            functions = []
            line = nl.readline()
            while line.startswith(b"F"):
                functions.append(line.rstrip(b"\r\n").rsplit(b" ", 1))
                line = nl.readline()
            if not any(name in _EQUILIBRIUM_FUNCTIONS for _, name in functions):
                return
            body = line + nl.read()
        libraries = amplfunc_merge(os.environ | {"AMPLFUNC": self._original_amplfunc})
        readers = []
        for suffix in (b"_value", b"_first_order"):
            path = TempfileManager.create_tempfile(suffix=".nl")
            with open(path, "wb") as nl:
                nl.write(header)
                for prefix, name in functions:
                    alias = name + suffix if name in _EQUILIBRIUM_FUNCTIONS else name
                    nl.write(b" ".join((prefix, alias)) + b"\n")
                nl.write(body)
            with CtypesEnviron(AMPLFUNC=libraries):
                readers.append(AmplInterface(path))
        self._asl = _AmplInterface(self._asl, *readers)

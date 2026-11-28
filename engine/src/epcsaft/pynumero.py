"""Opt-in derivative deferral for Engine equilibrium callbacks in PyNumero 6.10.1."""

import os

import numpy as np
from pyomo.common.env import CtypesEnviron
from pyomo.common.tempfiles import TempfileManager
from pyomo.contrib.pynumero.asl import AmplInterface
from pyomo.contrib.pynumero.interfaces.pyomo_nlp import PyomoNLP
from pyomo.solvers.amplfunc_merge import amplfunc_merge

_EQUILIBRIUM_FUNCTIONS = {
    b"epcsaft_equilibrium_phase_component_amount",
    b"epcsaft_equilibrium_solid_amount",
    b"epcsaft_equilibrium_system_enthalpy",
}


class _AmplInterface:
    def __init__(self, exact, value, first):
        self._exact = exact
        self._constraints = np.empty(exact.get_n_constraints(), dtype=np.float64)
        self.eval_f = value.eval_f
        self.eval_g = value.eval_g
        self.eval_deriv_f = first.eval_deriv_f
        self.eval_jac_g = first.eval_jac_g

    def __getattr__(self, name):
        return getattr(self._exact, name)

    def eval_hes_lag(self, x, lam, hes_lag, obj_factor):
        # PyNumero's sphes uses the exact reader's current point, as set by objval/conval.
        self._exact.eval_f(x)
        self._exact.eval_g(x, self._constraints)
        self._exact.eval_hes_lag(x, lam, hes_lag, obj_factor)


class EquilibriumPyomoNLP(PyomoNLP):
    """Use separate ASL readers for equilibrium values, first actions and exact Hessians."""

    def __init__(self, pyomo_model, nl_file_options=None):
        self._original_amplfunc = os.environ.get("AMPLFUNC", "")
        super().__init__(pyomo_model, nl_file_options)

    def _collect_nlp_structure(self):
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

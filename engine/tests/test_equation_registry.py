from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "docs" / "sync_equation_registry.py"


def _load_sync_module():
    spec = importlib.util.spec_from_file_location("sync_equation_registry", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_equation_parser_requires_unique_complete_metadata(tmp_path: Path) -> None:
    module = _load_sync_module()
    tex = tmp_path / "equations.tex"
    tex.write_text(
        r"""
\section{Residual Helmholtz}
% EqID: residual_helmholtz
% Source: \cite{Figiel2025}, Eq.~(1)
% Status: Source equation
% Description: Sum of active residual contributions.
% Change note: Reconciled to current Provider.
\begin{equation}
a^r=a^{hc}+a^{disp}+a^{assoc}+a^{DH}+a^{Born}
\label{eq:residual_helmholtz}
\end{equation}
""",
        encoding="utf-8",
    )
    entries = module.parse_equations(tex, tmp_path)
    assert [entry.eqid for entry in entries] == ["residual_helmholtz"]
    assert entries[0].status == "Source equation"


def test_equation_parser_rejects_duplicate_eqid(tmp_path: Path) -> None:
    module = _load_sync_module()
    tex = tmp_path / "equations.tex"
    block = r"""
% EqID: residual_helmholtz
% Source: \cite{Figiel2025}, Eq.~(1)
% Status: Source equation
% Description: Sum of active residual contributions.
% Change note: Reconciled to current Provider.
\begin{equation}
a^r=a^{hc}+a^{disp}
\label{eq:residual_helmholtz}
\end{equation}
"""
    tex.write_text(block + block, encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate EqID"):
        module.parse_equations(tex, tmp_path)


def test_source_annotations_link_real_owners_and_require_function_classification(
    tmp_path: Path,
) -> None:
    module = _load_sync_module()
    native = tmp_path / "src" / "epcsaft" / "_native"
    native.mkdir(parents=True)
    owner = native / "owner.hpp"
    owner.write_text(
        """
// EqID: residual_helmholtz
double residual_helmholtz(double a, double b) {
    return a + b;
}

// EqExclude: scalar utility dispatch
double scalar_log(double value) {
    return std::log(value);
}
""",
        encoding="utf-8",
    )
    refs, exclusions = module.parse_source_annotations(
        tmp_path / "src", tmp_path
    )
    assert refs["residual_helmholtz"][0].symbol == "residual_helmholtz"
    assert exclusions[0].reason == "scalar utility dispatch"

    owner.write_text(
        """
double unclassified(double value) {
    return value + 1.0;
}
""",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="unclassified mathematical-owner function"):
        module.parse_source_annotations(tmp_path / "src", tmp_path)


def _entry(module, *, status: str = "Source equation", source: str = r"\cite{Known}"):
    reference = module.SourceReference(
        file="src/owner.hpp",
        line=2,
        comment_line=1,
        symbol="residual_helmholtz",
        context="// EqID: residual_helmholtz",
    )
    return module.EquationEntry(
        eqid="residual_helmholtz",
        section="Residual Helmholtz",
        subsection="",
        label="eq:residual_helmholtz",
        source=source,
        status=status,
        description="Sum of active residual contributions.",
        change_note="Reconciled to current Provider.",
        latex=r"a^r=a^{hc}+a^{disp}",
        tex_file="docs/latex/equations.tex",
        tex_line=10,
        implementation_refs=(reference,),
    )


def test_registry_validation_rejects_stale_links_and_unknown_citations() -> None:
    module = _load_sync_module()
    entry = _entry(module)
    with pytest.raises(ValueError, match="source EqID absent from LaTeX"):
        module.validate_registry([entry], {"unknown": []}, [], {"Known"})

    missing_owner = _entry(module)
    missing_owner = module.EquationEntry(
        **{
            **missing_owner.__dict__,
            "implementation_refs": (),
        }
    )
    with pytest.raises(ValueError, match="has no implementation owner"):
        module.validate_registry([missing_owner], {}, [], {"Known"})

    documentation_only = _entry(module, status="Documentation-only")
    documentation_only = module.EquationEntry(
        **{
            **documentation_only.__dict__,
            "implementation_refs": (),
        }
    )
    module.validate_registry([documentation_only], {}, [], {"Known"})

    with pytest.raises(ValueError, match="unknown citation"):
        module.validate_registry(
            [_entry(module, source=r"\cite{MissingKey}")],
            {"residual_helmholtz": []},
            [],
            {"Known"},
        )


def test_registry_rendering_is_deterministic_and_check_mode_is_read_only(
    tmp_path: Path,
) -> None:
    module = _load_sync_module()
    entry = _entry(module)
    exclusion = module.SourceExclusion(
        file="src/owner.hpp",
        line=8,
        symbol="scalar_log",
        reason="scalar utility dispatch",
    )
    assert module.render_registry([entry], [exclusion]) == module.render_registry(
        [entry], [exclusion]
    )
    markdown = module.render_markdown([entry], [exclusion])
    assert "residual_helmholtz" in markdown
    assert "src/owner.hpp:2" in markdown
    assert "scalar utility dispatch" in markdown

    (tmp_path / "docs" / "latex").mkdir(parents=True)
    (tmp_path / "src" / "epcsaft" / "_native").mkdir(parents=True)
    (tmp_path / "docs" / "latex" / "equations.tex").write_text(
        r"""
\section{Residual Helmholtz}
% EqID: residual_helmholtz
% Source: \cite{Known}
% Status: Source equation
% Description: Sum of active residual contributions.
% Change note: Reconciled to current Provider.
\begin{equation}
a^r=a^{hc}+a^{disp}
\label{eq:residual_helmholtz}
\end{equation}
""",
        encoding="utf-8",
    )
    (tmp_path / "docs" / "latex" / "references.bib").write_text(
        "@article{Known, title={Known source}}\n",
        encoding="utf-8",
    )
    (tmp_path / "src" / "epcsaft" / "_native" / "owner.hpp").write_text(
        """
// EqID: residual_helmholtz
double residual_helmholtz(double a, double b) {
    return a + b;
}
""",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="generated equation documentation is stale"):
        module.synchronize(tmp_path, check=True)
    assert not (tmp_path / "docs" / "equations.md").exists()


def test_repository_equation_registry_is_current_and_complete() -> None:
    module = _load_sync_module()
    module.synchronize(REPO_ROOT, check=True)


def test_required_provider_equation_families_are_present() -> None:
    module = _load_sync_module()
    entries = module.parse_equations(
        REPO_ROOT / "docs" / "latex" / "equations.tex",
        REPO_ROOT,
    )
    eqids = {entry.eqid for entry in entries}
    assert {
        "figiel_bulk_permittivity",
        "ascani_bulk_permittivity",
        "debye_screening_parameter",
        "debye_huckel_helmholtz",
        "original_born_helmholtz",
        "ssm_ds_born_helmholtz",
        "formula_unit_composition",
        "reference_molality_sequence",
        "single_ion_activity_coefficient",
        "mean_ionic_activity_coefficient",
        "molality_activity_conversion",
        "practical_osmotic_coefficient",
        "infinite_dilution_solvation_energy",
        "transfer_free_energy",
        "phase_total_amount",
        "phase_mole_fraction",
        "phase_molar_density",
        "phase_dimensionless_helmholtz",
        "phase_pressure",
        "packing_volume_bounds",
        "neutral_reference_stoichiometric_basis",
    } <= eqids


def test_property_identities_and_cppad_boundary_are_documented() -> None:
    tex = (REPO_ROOT / "docs" / "latex" / "equations.tex").read_text()
    for eqid in (
        "compressibility_factor",
        "pressure",
        "residual_chemical_potential",
        "log_fugacity_coefficient",
        "residual_entropy",
        "residual_enthalpy",
        "residual_heat_capacity_cp",
        "density_residual",
        "mechanical_stability",
    ):
        assert f"% EqID: {eqid}" in tex
    assert "CppAD is the sole nonlinear derivative authority" in tex
    assert "% Status: Analytical derivative" not in tex
    assert "% Status: Manual derivative" not in tex

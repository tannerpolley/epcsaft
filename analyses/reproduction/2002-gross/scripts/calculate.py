"""Declare Gross–Sadowski states; all thermodynamics and root searches run in the Engine."""
import csv
import hashlib
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from time import perf_counter

import epcsaft
from epcsaft import Mixture, Parameters, Pure, _core
from epcsaft.equilibrium import Amounts, Free, Phase, Pinned, Problem, solve_equilibrium, solve_equilibrium_starts, EquilibriumStart, EquilibriumWorkLimits

CASE = Path(__file__).resolve().parents[1]
ROOT = CASE.parents[2]
sys.path.insert(0, str(ROOT / "tools"))
from engine_identity import loaded_identity

RESULTS = CASE / "results/data"
ROW_TOLERANCE = 1e-8


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def write(path, data):
    with path.open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(data[0]) if data else [], lineterminator="\n")
        writer.writeheader()
        writer.writerows(data)


def parameters(names, kij):
    source = {r["component"]: r for r in rows(CASE / "inputs/parameters/pure/any-solvent.csv")}
    source["water"] = source["H2O"]
    # Gross 2002 Figure 2 uses isobutane; the local matrix mistakenly assigns its kij to an alcohol.
    upstream = ROOT / "analyses/reproduction/2001-gross/results/tables/table-002/source/table-002.csv"
    iso = next(r for r in csv.reader(upstream.open()) if r[0] == "isobutane")
    source["isobutane"] = dict(zip(("MW", "m", "s", "e"), (str(float(iso[1])/1000), *iso[2:5])), e_assoc="0", vol_a="0")
    document = _core.ParameterDocument()
    document.component_ids = list(names)
    document.base_choice = "pc-saft"
    document.electrolyte_choice = document.permittivity_choice = document.polar_choice = "none"
    entries, sites, multiplicity, self_edges = [], [], [], []

    def entry(field, index, value, unit, identity=None):
        item = _core.ParameterEntry()
        item.identity = identity or f"{field}/{index}"
        item.field, item.index, item.magnitude, item.unit = field, list(index), float(value), unit
        entries.append(item)

    for i, name in enumerate(names):
        r = source[name]
        for field, key, unit in (("molar_mass", "MW", "kg / mol"), ("segment_count", "m", "dimensionless"),
                                 ("sigma_constant", "s", "angstrom"), ("dispersion_energy_over_k", "e", "kelvin")):
            entry(field, (i,), r[key], unit)
        entry("charge_number", (i,), 0, "elementary-charge")
        if float(r["e_assoc"]):
            a = len(sites)
            sites.extend((i, i))
            multiplicity.extend((1, 1))
            self_edges.append((a, a+1))
            entry("association_energy_over_k", (a, a+1), r["e_assoc"], "kelvin")
            entry("association_volume", (a, a+1), r["vol_a"], "dimensionless")
    if len(names) == 2:
        entry("k_ij", (0, 1), kij, "dimensionless")
    if len(self_edges) == 2:
        a, b = self_edges
        for edge in ((a[0], b[1]), (a[1], b[0])):
            entry("association_combining_rule", (*edge, *a, *b), 0, "dimensionless", "wolbach-sandler")
    document.entries = entries
    document.site_component, document.site_multiplicity = sites, multiplicity
    document.association_choice = "general-site" if sites else "none"
    return Parameters(_core.build_parameters(document), tuple(names),
                      tuple((names[c], "A" if k % 2 == 0 else "B") for k, c in enumerate(sites)))


def figure_tasks():
    tasks = []
    for i, r in enumerate(rows(CASE / "figures/figure-01/source/source-points.csv")):
        tasks.append(dict(figure="figure-01", dataset=r["source_role"], series=r["branch"], names=[r["component"]], kij=0,
                          mode="saturation", T=float(r["temperature_K"]), P=1e5, x=1., quantity="density_kg_m3",
                          reference=float(r["density_kg_m3"]), tolerance=float(r["density_uncertainty_kg_m3"]),
                          input_uncertainty=float(r["temperature_uncertainty_K"]), source=r["source_reference"], index=i))
    for fig in range(2, 11):
        for filename, dataset in (("source-points.csv", "paper_model"), ("literature-points.csv", "experiment")):
            for i, r in enumerate(rows(CASE / f"figures/figure-{fig:02}/source/{filename}")):
                series = r["series"]
                mode = "dew" if "dew" in series or r.get("source_role") == "dew_curve" else "bubble"
                if fig == 2:
                    names, kij, xkey, ykey = ["isobutane", "methanol"], .05, "x_component_1", "y_component_1"
                elif fig == 3:
                    names, kij, xkey, ykey = ["1-propanol", "ethylbenzene"], .023, "x_1_propanol", "y_1_propanol"
                    if dataset == "experiment":
                        xkey, ykey = "x_alcohol", "y_vapor"
                elif fig == 4:
                    names, kij, xkey, ykey = ["1-pentanol", "benzene"], .0135, "x_1_pentanol", "y_1_pentanol"
                    if dataset == "experiment":
                        xkey = "x_alcohol"
                elif fig == 5:
                    alcohol = "2-propanol" if "2-propanol" in series else "1-propanol"
                    names, kij, xkey, ykey = [alcohol, "benzene"], (.021 if alcohol == "2-propanol" else .020), "x_alcohol", "y_alcohol"
                elif fig in (6, 7):
                    names, kij, xkey, ykey = ["butane", "1-butanol" if fig == 6 else "ethanol"], (.015 if fig == 6 else .028), "x_butane", ""
                elif fig in (8, 9):
                    names, kij, xkey, ykey = ["methanol", "cyclohexane" if fig == 8 else "1-octanol"], (.051 if fig == 8 else .020), "x_methanol", "y_methanol"
                else:
                    names, kij, xkey, ykey = ["water", "1-pentanol"], .016, "x_water", "x_water"
                    if series == "upper_vle_curve":
                        mode = "dew"  # PDF Figure 10: upper envelope is the vapor composition.
                if series.startswith("lle"):
                    mode = "lle"
                if series == "vlle_tie_line":
                    mode = "vlle"
                x = float(r[ykey if mode == "dew" else xkey])
                T, P = float(r["T_K"]), float(r.get("P_bar") or r["pressure_bar"]) * 1e5
                quantity = "x" if mode == "lle" else "T_K" if fig in (3, 8, 9, 10) else "P_bar"
                reference = x if quantity == "x" else T if quantity == "T_K" else P / 1e5
                uncertainty_key = "uncertainty_x" if quantity == "x" else "uncertainty_T_C" if quantity == "T_K" else "uncertainty_P_bar"
                record_dataset = "experiment_paper" if r.get("source_role") == "experimental_marker" else dataset
                if dataset == "experiment" and fig in (3, 10):
                    record_dataset = "experiment_supplemental"
                tasks.append(dict(figure=f"figure-{fig:02}", dataset=record_dataset, series=series, names=names, kij=kij, mode=mode,
                                  T=T, P=P, x=x, quantity=quantity, reference=reference,
                                  tolerance=float(r.get(uncertainty_key) or 0),
                                  input_uncertainty=float(r.get("uncertainty_composition") or r.get("uncertainty_x") or 0),
                                  source=r["source_reference"], index=i))
    table = list(csv.reader((CASE / "results/tables/table-001/source/table-001.csv").open()))
    for r in table[1:]:
        if not r[2]:
            continue
        low, high = map(float, r[9].split("-"))
        for i, T in enumerate(((low+high)/2,)):
            tasks.append(dict(figure="table-001", dataset="prediction", series="saturation", names=[r[0]], kij=0,
                              mode="saturation", T=T, P=1e5, x=1., quantity="P_bar", reference=0., tolerance=0.,
                              input_uncertainty=0., source="Gross 2002 Table 1 range; no pointwise reference retained", index=i))
    for fig in (8, 10):
        for i, r in enumerate(rows(CASE / f"figures/figure-{fig:02}/source/experimental-points.csv")):
            names, kij = next((t["names"], t["kij"]) for t in tasks if t["figure"] == f"figure-{fig:02}")
            mode = "lle" if fig == 8 else "dew" if r["series"] == "dew_curve" else "bubble"
            x = float(r["x_methanol"] if fig == 8 else r["x_water"])
            tasks.append(dict(figure=f"figure-{fig:02}", dataset=r["dataset"], series=r["series"], names=names,
                              kij=kij, mode=mode, T=float(r["T_K"]), P=float(r["pressure_bar"])*1e5,
                              x=x, quantity="x" if fig == 8 else "T_K", reference=x if fig == 8 else float(r["T_K"]),
                              tolerance=float(r["uncertainty_x"] if fig == 8 else r["uncertainty_T_C"]),
                              input_uncertainty=float(r["uncertainty_T_C"] if fig == 8 else r["uncertainty_x"]),
                              source=r["source_reference"], index=i))
    return tasks


def calculate(task, declare_only=False, seed=None):
    start = perf_counter()
    out = {**task, "names": "/".join(task["names"]), "value": "", "deviation": "", "absolute_relative_percent": "",
           "within_curve_tolerance": "", "success": False, "status": "", "message": "", "residual": "", "packing_fractions": "", "iterations": "",
           "evaluations": "", "T_solved_K": "", "P_solved_Pa": "", "densities_mol_m3": "", "compositions": "",
           "pressure_aad_percent": "", "liquid_density_aad_percent": "", "liquid_volume_aad_percent": "", "seconds": ""}
    try:
        names, mode, x = task["names"], task["mode"], task["x"]
        endpoint = mode in ("bubble", "dew") and x in (0., 1.)
        if endpoint:
            names = [names[0 if x == 1. else 1]]
        model = (Pure if len(names) == 1 else Mixture)(parameters(names, task["kij"]))
        if len(names) == 1:
            phases = [Phase("L", kind="liquid"), Phase("V", kind="vapor", amount=Pinned(0))]
            feed = Amounts({names[0]: 1.})
        elif mode in ("lle", "vlle"):
            guesses = (.1, .9) if task["figure"] == "figure-08" else (.3, .997)
            phases = [Phase("L1", kind="liquid", composition_guess=(guesses[0], 1-guesses[0])),
                      Phase("L2", kind="liquid", composition_guess=(guesses[1], 1-guesses[1]))]
            if mode == "vlle":
                phases.append(Phase("V", kind="vapor", amount=Pinned(0), composition_guess=(.87, .13)))
            # The water fraction .5 lies outside the source's high-temperature tie lines.
            overall = .8 if task["figure"] == "figure-10" else .5
            feed = Amounts({names[0]: overall, names[1]: 1-overall})
        else:
            phases = [Phase("L", kind="liquid", amount=Pinned(0) if mode == "dew" else Free()),
                      Phase("V", kind="vapor", amount=Pinned(0) if mode == "bubble" else Free())]
            feed = Amounts({names[0]: x, names[1]: 1-x})
        T = Free(task["T"]) if task["quantity"] == "T_K" or mode == "vlle" else task["T"]
        P = Free(task["P"]) if task["quantity"] in ("P_bar", "density_kg_m3") else task["P"]
        problem = Problem(phases=phases, T=T, P=P, feed=feed, row_tolerance=ROW_TOLERANCE)
        if declare_only:
            return model, problem
        if seed is None:
            result = solve_equilibrium(model, problem)
        else:
            attempts = solve_equilibrium_starts(model, problem, [EquilibriumStart(**seed)], EquilibriumWorkLimits())
            result = attempts.results[0]
            if result is None:
                raise RuntimeError(attempts.causes[0])
        out.update(success=result.success and result.requested_tolerance_met is True, status=result.message.partition(":")[0], message=result.message,
                   packing_fractions=json.dumps(result.packing_fractions),
                   residual=result.convergence.residual, iterations=result.iterations, evaluations=result.evaluations,
                   T_solved_K=result.temperature or "", P_solved_Pa=result.pressure or "",
                   densities_mol_m3=json.dumps(result.molar_densities), compositions=json.dumps(result.mole_fractions))
        if not out["success"]:
            return out
        assert result.convergence.residual <= ROW_TOLERANCE
        assert all(rho > 0 for rho in result.molar_densities)
        if mode in ("lle", "vlle"):
            first, second = result.mole_fractions[0], result.mole_fractions[2]
            if abs(first-second) < 1e-6:
                out.update(success=False, message="coalesced liquid phases; no liquid-liquid result")
                return out
            value = min(first, second) if "lean" in task["series"] or "pentanol_rich" in task["series"] else max(first, second)
            if mode == "vlle":
                value = result.temperature
        elif task["quantity"] == "density_kg_m3":
            mass = float(next(r["MW"] for r in rows(CASE / "inputs/parameters/pure/any-solvent.csv") if r["component"] == names[0]))
            value = result.molar_densities[0 if task["series"] == "liquid" else 1] * mass
        elif task["quantity"] == "T_K":
            value = result.temperature
        else:
            value = result.pressure / 1e5
        ref = task["reference"]
        out.update(value=value, deviation=value-ref if task["dataset"] != "prediction" else "", absolute_relative_percent=100*abs(value/ref-1) if ref else "",
                   within_curve_tolerance=abs(value-ref) <= task["tolerance"] if task["dataset"] in ("paper_model", "paper_pc_saft_curve") else "")
        if task["quantity"] == "density_kg_m3" and task["series"] == "liquid" and ref:
            out["liquid_density_aad_percent"] = 100*abs(value/ref-1)
            out["liquid_volume_aad_percent"] = 100*abs(ref/value-1)
        if task["quantity"] == "P_bar":
            out["pressure_aad_percent"] = out["absolute_relative_percent"]
    except Exception as error:
        if declare_only:
            raise
        out["status"] = type(error).__name__
        out["message"] = f"{type(error).__name__}: {error}"
    finally:
        out["seconds"] = perf_counter()-start
    return out


def declarations():
    return [calculate(task, declare_only=True) for task in figure_tasks()]


def main(out=RESULTS):
    out.mkdir(parents=True, exist_ok=True)
    tasks = figure_tasks()
    assert len(tasks) >= 672, len(tasks)
    write(out / "declared-states.csv", tasks)
    started, start = datetime.now(timezone.utc).isoformat(), perf_counter()
    # Preload the editable native extension before starting all CPU workers.
    provenance = loaded_identity()
    with ProcessPoolExecutor(max_workers=os.cpu_count()) as pool:
        results = list(pool.map(calculate, tasks))
    write(out / "default-start-comparisons.csv", results)
    default_seconds = perf_counter()-start
    # The review identified these declarations and demonstrated neighbouring starts.
    neighbours = {21: 20, 22: 20, 23: 20, 41: 40, 42: 39, 43: 39,
                  65: 64, 66: 63, 67: 63, 328: 327, 334: 333, 335: 334,
                  336: 335, 360: 334, 361: 335, 362: 336, 492: 491, 531: 530,
                  558: 557, 722: 724, 723: 724}
    failures = []
    for line, neighbour in neighbours.items():
        task, default, previous = tasks[line-2], results[line-2], results[neighbour-2]
        if default["success"]:  # the default start now returns this state
            continue
        assert previous["success"], (line, neighbour)
        n = len(task["names"])
        fractions = json.loads(previous["compositions"])
        seed = dict(compositions=[fractions[q:q+n] for q in range(0, len(fractions), n)],
                    packing_fractions=json.loads(previous["packing_fractions"]),
                    amounts=[0., 1.] if task["mode"] == "dew" else [1., 0.])
        if task["quantity"] in ("P_bar", "density_kg_m3"):
            seed["pressure"] = previous["P_solved_Pa"]
        else:
            seed["temperature"] = previous["T_solved_K"]
        if line == 558:
            seed["compositions"][0] = [.3, .7]
            seed["packing_fractions"][0] = .4
        recovered = calculate(task, seed=seed)
        assert recovered["success"], (line, recovered["message"])
        model, problem = calculate(task, declare_only=True)
        phases = [dict(name=p.name, kind=p.kind, amount=p.amount.value if isinstance(p.amount, Pinned) else "free",
                       composition_guess=p.composition_guess, packing_guess=p.packing_guess) for p in problem.phases]
        failures.append(dict(csv_line=line, **task, default_status=default["status"], default_message=default["message"],
                             feed_mol=json.dumps(problem.feed.values), phase_spec=json.dumps(phases),
                             temperature_spec=json.dumps(dict(free=isinstance(problem.T, Free),
                                                               value_or_guess=problem.T.guess if isinstance(problem.T, Free) else problem.T)),
                             pressure_spec=json.dumps(dict(free=isinstance(problem.P, Free),
                                                            value_or_guess=problem.P.guess if isinstance(problem.P, Free) else problem.P)),
                             row_tolerance=ROW_TOLERANCE,
                             neighbour_csv_line=neighbour, recovered_start=json.dumps(seed),
                             recovered_status=recovered["status"], recovered_value=recovered["value"],
                             recovered_residual=recovered["residual"]))
        results[line-2] = recovered
    write(out / "default-start-failures.csv", failures)
    write(out / "state-comparisons.csv", results)
    metadata = dict(started_utc=started, calculation_wall_seconds=perf_counter()-start, workers=os.cpu_count(),
                    row_tolerance=ROW_TOLERANCE, overall_water_fraction=.8, default_start_wall_seconds=default_seconds,
                    recovered_default_failures=len(failures), engine=provenance, states=len(tasks),
                    source_pdf=json.loads((RESULTS / 'calculation.json').read_text())['source_pdf'],
                    source_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                                   for p in sorted([*CASE.glob("figures/*/source/*.csv"),
                                                    *CASE.glob("inputs/parameters/**/*.csv"),
                                                    *CASE.glob("results/tables/*/source/*.csv"),
                                                    ROOT / "analyses/reproduction/2001-gross/results/tables/table-002/source/table-002.csv"])} )
    (out / "calculation.json").write_text(json.dumps(metadata, indent=2, default=str)+"\n")
    print(json.dumps({k: metadata[k] for k in ("states", "workers", "calculation_wall_seconds")}))
    print(f"Returned qualified states: {sum(r['success'] for r in results)}")


if __name__ == "__main__":
    main(*map(Path, sys.argv[1:]))

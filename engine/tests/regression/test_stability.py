"""Exact reduced-distance actions on the declared neutral/reacting branches."""
import pytest
from typing import Any
from pathlib import Path
from epcsaft import Parameters,regression as r,equilibrium as eq
from tests.native_inputs import parameters,ETHANE

ROOT=Path(__file__).resolve().parents[3]

def declarations() -> list[Any]:
    neutral=Parameters.from_json(ROOT/'engine/tests/cases/khudaida-idaes-settler/inputs/data/khudaida-2026-figure-2-electrolyte-lle/7/parameters/parameters.json',components=('water','isobutanol'))
    nc=[r.coordinate(neutral,[('k_ij',('water','isobutanol'),1)],origin=-.20,scale=.1,bounds=(-.5,.1))]
    np=eq.Problem(T=293.15,P=1e5,row_tolerance=1e-10,feed=eq.Amounts({'water':.6,'isobutanol':.4}),phases=[eq.Phase('L',kind='liquid',composition_guess=(.6,.4))])
    nr=eq.Problem(T=293.15,P=eq.Free(1e4),row_tolerance=1e-10,feed=np.feed,phases=[eq.Phase('L',kind='liquid',composition_guess=(.6,.4)),eq.Phase('V',kind='vapor',amount=eq.Pinned(0.))])
    reacting=parameters((ETHANE,)*2,pairs=[('k_ij',(0,1),0.)])
    reaction=eq.Reaction({'s0':-1.,'s1':1.},name='isomerisation',correlation=eq.ReactionLogPolynomial(-1.25,90.,0.,0.,250.,150.,350.,eq.EOS_STANDARD_STATE_ID))
    rc=[r.coordinate(reacting,[('k_ij',('s0','s1'),1)],origin=0.,scale=.1,bounds=(-.1,.1)),r.reaction_coordinate('isomerisation','a',origin=-1.25,scale=1.),r.reaction_coordinate('isomerisation','b',origin=90.,scale=250.)]
    rp=eq.Problem(T=200.,P=1e5,row_tolerance=1e-10,feed=eq.Amounts({'s0':.5,'s1':.5}),phases=[eq.Phase('L',kind='liquid')],reactions=[reaction])
    rr=eq.Problem(T=200.,P=eq.Free(2e5),row_tolerance=1e-10,feed=rp.feed,phases=[eq.Phase('L',kind='liquid'),eq.Phase('V',kind='vapor',amount=eq.Pinned(0.))],reactions=[reaction])
    zero=parameters((ETHANE,)*3,pairs=[('k_ij',(0,1),0.),('k_ij',(0,2),0.),('k_ij',(1,2),0.)])
    zreaction=eq.Reaction({'s1':-1.,'s2':1.},name='unfed-isomerisation',correlation=reaction.correlation)
    zc=[r.coordinate(zero,[('dispersion_energy_over_k','s0',1)],origin=191.42,scale=250.,bounds=(150.,250.)),r.reaction_coordinate('unfed-isomerisation','a',origin=-1.25,scale=1.)]
    zp=eq.Problem(T=200.,P=1e5,row_tolerance=1e-10,feed=eq.Amounts({'s0':1.}),phases=[eq.Phase('L',kind='liquid')],reactions=[zreaction])
    return [(neutral,nc,np,nr,2.),(reacting,rc,rp,rr,.5),(zero,zc,zp,None,1.)]


def branch_identity(measure: Any) -> tuple[Any, ...]:
    trial = measure.report.trials[measure.minimizing_trial]
    return (tuple(trial.phase.support), trial.phase.kind, trial.quantities['density_side'],
            tuple(sorted(measure.report.trials[i].seed for i in measure.minima)),
            tuple(p['density_side'] for p in measure.report.present_phases))


@pytest.mark.parametrize("case,relative", [(0, False), (0, True), (1, False), (1, True), (2, False)])
def test_stability_actions_against_resolved_central_differences(case: int, relative: bool) -> None:
    p, coordinates, parent, reference, ratio = declarations()[case]
    row = r.stability(p, parent, scale=.1, pressure_reference=reference if relative else None, pressure_ratio=ratio)
    center = [1.929] if case == 0 else [.5, 0., 0.] if case == 1 else [0., 0.]
    def evaluate(point: list[float]) -> tuple[Any, Any, Any]:
        value = r.differentiate(p, coordinates, optimizer=point, observations=[row])
        prediction = value.stability_predictions[0]
        measure = prediction.stability
        assert measure.available and not measure.event, (measure.cause, measure.event)
        return value, prediction, measure
    value, prediction, measure = evaluate(center)
    for k in range(len(coordinates)):
        for step in (5e-4, 1e-3, 2e-3):
            plus, minus = center.copy(), center.copy()
            plus[k] += step
            minus[k] -= step
            _, _, high = evaluate(plus)
            _, _, low = evaluate(minus)
            assert branch_identity(measure) == branch_identity(high) == branch_identity(low)
            assert (measure.value == 0) == (high.value == 0) == (low.value == 0)
            uncertainty = (high.uncertainty + low.uncertainty) / (2 * step)
            for analytic, difference, error in (
                (prediction.stability_raw_gradient[k], (high.raw_minimum - low.raw_minimum) / (2 * step), uncertainty),
                (100 * value.optimizer_gradient[k], 100 * (high.value - low.value) / (2 * step), 100 * uncertainty if measure.value else 0.),
            ):
                tolerance = 1e-6 + 1e-4 * max(abs(analytic), abs(difference))
                assert error < tolerance, "stencil is inconclusive"
                assert abs(analytic - difference) <= tolerance

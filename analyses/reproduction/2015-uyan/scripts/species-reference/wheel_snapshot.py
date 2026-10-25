import json
import math
from pathlib import Path
import sys
sys.path.insert(0,str(Path.cwd()/'engine'))
from epcsaft import Mixture, ActiveParameter, ParameterFamily
from epcsaft import equilibrium as eq
from tests.native_inputs import parameters, ETHANE, aqueous_electrolyte

def snapshot(mode):
    rows=[]
    for ionic in (False,True):
        p=aqueous_electrolyte() if ionic else parameters((ETHANE,ETHANE),[("k_ij",(0,1),0.)])
        model=Mixture(p,active_parameters=[ActiveParameter(ParameterFamily.DispersionEnergyOverK,[0])])
        for t,pres in ([(298.15,1e5)] if ionic else [(260.,3e6),(280.,5e6)]):
            for common in (False,True):
                for fixed in (False,True):
                    choices={} if mode=='empty' else None
                    if mode=='full':
                        choices={s:eq.SpeciesReference.InfiniteDilution for s in model.component_ids}
                        choices['s0']=eq.SpeciesReference.PureLiquid
                    kwargs={} if mode=='omitted' else {'species_reference':choices}
                    ref=eq.NeutralReference(solvent='s0',solute_molality_ratios={'s1':1.,'s2':1.} if ionic else {'s1':1.},standard_molality_mol_per_kg=1.,**kwargs)
                    basis=eq.ReferenceSourceBasis.CommonMolalityInfiniteDilution if common else eq.ReferenceSourceBasis.RawMoleFractionInfiniteDilution
                    problem=eq.Problem(T=t,P=pres,feed=eq.Amounts({'s0':1.,'s1':1e-6,'s2':1e-6} if ionic else {'s0':2.,'s1':.2}),phases=[eq.Phase('L',kind='liquid')],neutral_reference=ref,
                        reactions=[eq.Reaction({'s0':-1.,'s1':1.,'s2':1.} if ionic else {'s0':-2.,'s1':1.},correlation=eq.ReactionLogPolynomial(math.log(1e-8 if ionic else .3),0.,0.,0.,298.15,250.,350.,'manufactured-reference'),reference=eq.ReactionReference(basis,pres*2 if fixed else None))])
                    result=eq.solve_equilibrium(model,problem)
                    assert result.success,result.message
                    fields={name:getattr(result,name) for name in ('temperature','pressure','mole_fractions','phase_amounts','packing_fractions','residuals','standard_over_rt') if hasattr(result,name)}|{name:[getattr(result,name)] for name in ('temperature','pressure')}
                    d=result.neutral_reference_diagnostics[0]
                    fields['reference']={name:getattr(d,name) for name in ('terminal_contraction','observed_terminal_change','terminal_contraction_pressure_derivative','terminal_contraction_pressure_second_derivative','reference_pressure_pa')}
                    fields['availability']=str(d.status)
                    actions=[]
                    for dt,dp,parameter in [(1.,0.,[0.]),(0.,1e5,[0.]),(0.,0.,[1.])]:
                        direction=eq.SolvedStateActionDirection(dt,dp,[0.]*len(model.component_ids),parameter)
                        item=eq.solved_state_actions(eq.compile_problem(model,problem),result,eq.SolvedStateActionRequest([eq.SolvedStateObservable(eq.SolvedStateObservableKind.PhaseMoleFraction,0,1)],[direction])).results[0]
                        actions.append((str(item.status),item.action))
                    fields['actions']=actions
                    rows.append(fields)
    return rows

def flatten(value, prefix=""):
    if isinstance(value, dict):
        return {k2: x for k, x1 in value.items() for k2, x in flatten(x1, prefix+k+".").items()}
    if isinstance(value, (list, tuple)):
        return {k2: x for i, x1 in enumerate(value) for k2, x in flatten(x1, prefix+str(i)+".").items()}
    return {prefix.rstrip("."): value.hex() if isinstance(value, float) else "unavailable" if value is None else str(value)}

import csv
baseline=Path(__file__).resolve().parents[2]/"results/data/native/species-reference/default-wheel-values.tsv"
with baseline.open() as stream:
    old=list(csv.DictReader(stream,delimiter="\t"))
for mode in ("omitted","empty","full"):
    current=[{"case_index":str(i),**flatten(row)} for i,row in enumerate(snapshot(mode))]
    for actual,expected in zip(current,old,strict=True):
        assert not set(actual)-set(expected),mode
        assert {k:actual.get(k, "not-applicable") for k in expected}==expected,mode
    print(mode,len(current),"rows bitwise identical to base wheel")

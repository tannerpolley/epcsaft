import epcsaft
from epcsaft import equilibrium as eq

record = "data/packets/gross-2001-methane-ethane/3/parameters/parameters.json"
mix = epcsaft.Mixture(epcsaft.Parameters.from_json(record, components=["methane", "ethane"]))
feed = eq.MoleFractions({"methane": 0.3, "ethane": 0.7})
problem = eq.bubble_point(liquid=eq.Phase("L", kind="liquid"), T=200.0, feed=feed)
bubble = eq.solve_equilibrium(mix, problem)          # pressure is free
print(bubble.pressure, bubble.phase_mole_fractions(1))   # Pa; incipient vapor

dT = eq.SolvedStateActionDirection(1.0, 0.0, [0.0, 0.0], [ ])  # unit temperature direction
P = eq.SolvedStateObservable(eq.SolvedStateObservableKind.Pressure)
batch = eq.solved_state_actions(eq.compile_problem(mix, problem), bubble,
                                eq.SolvedStateActionRequest([P], [dT]))
print(batch.results[0].action, batch.results[0].status)       # Pa/K

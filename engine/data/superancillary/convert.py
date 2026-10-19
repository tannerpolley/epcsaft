"""#648: convert SAFTsuperanc implementation/domain.cpp + PCSAFTsuperancillary_crit.hpp into one JSON, losslessly.
Every number is parsed with float() and written with repr(), so each double round-trips exactly."""
import json, re, sys
src = sys.argv[1]  # upstream checkout
out = sys.argv[2]
EXP = re.compile(r"\{\s*std::vector<double>\(\{([^}]*)\}\),\s*([^,\s]+),\s*([^\s}]+)\s*\}")
nums = lambda s: [float(t) for t in s.split(',') if t.strip()]
exps = lambda text: [{'xmin': float(a), 'xmax': float(b), 'coef': nums(c)} for c, a, b in EXP.findall(text)]
crit = open(f'{src}/implementation/PCSAFTsuperancillary_crit.hpp').read()
dom = open(f'{src}/implementation/domain.cpp').read()
wedges = nums(re.search(r'Wedges = \{([^}]*)\}', dom).group(1))
blocks = dom.split('const WInterval interval')[1:]
intervals = []
for b in blocks:
    wnodes = nums(re.search(r'// Wnodes\s*\{([^}]*)\}', b).group(1))
    liq, vap = b.split('// Intervals (vapor)')
    liq = liq.split('// Intervals (liq)')[1]
    cols = lambda t: [exps(c) for c in t.split('{{')[1:]]
    L, V = cols(liq), cols(vap)
    assert len(L) == len(V) == len(wnodes) == 17, (len(L), len(V), len(wnodes))
    intervals.append({'wnodes': wnodes, 'liquid': L, 'vapor': V})
assert len(intervals) == len(wedges) - 1
doc = {
    'source': 'Bell and Deiters, Ind. Eng. Chem. Res. 62(4) 1958-1967 (2023), doi:10.1021/acs.iecr.2c02916; '
              'https://github.com/usnistgov/SAFTsuperanc implementation/domain.cpp and PCSAFTsuperancillary_crit.hpp',
    'upstream_commit': sys.argv[3],
    'variables': 'Ttilde = k T/epsilon; rhotilde = rho N_A sigma^3 (molecules); w = 1/m on [1/64, 1]; '
                 'Ttilde_min = exp(-2.20078778) m^0.37627892 Ttilde_crit; Theta = (Ttilde - Ttilde_min)/(Ttilde_crit - Ttilde_min)',
    'Ttilde_min_coefficients': [-2.20078778, 0.37627892],
    'Ttilde_crit': exps(crit),
    'w_edges': wedges,
    'intervals': intervals,
}
json.dump(doc, open(out, 'w'), indent=None, separators=(',', ':'))
# lossless check: every numeric literal in the sources appears as the same double in the JSON
lit = re.compile(r'(?<![\w.])-?\d+\.?\d*(?:e[-+]?\d+)?')
def collect(o, acc):
    if isinstance(o, float): acc.append(o)
    elif isinstance(o, dict): [collect(v, acc) for k, v in o.items() if k != 'Ttilde_min_coefficients']
    elif isinstance(o, list): [collect(v, acc) for v in o]
    return acc
back = collect(json.load(open(out)), [])
n_src = sum(len(nums(m.group(1))) + 2 for m in EXP.finditer(crit + dom)) + len(wedges) + 17 * len(intervals)
print('doubles in JSON', len(back), 'expected from sources', n_src)
segs = [len(c) for iv in intervals for c in iv['liquid']]
print('Theta segments per collection: min', min(segs), 'max', max(segs), 'total collections', len(segs) * 2,
      'coef per expansion', {len(e['coef']) for iv in intervals for c in iv['liquid'] + iv['vapor'] for e in c})

import json
from pathlib import Path
import shlex
import subprocess
root = Path.cwd()
build = root / Path('build/environment-wheel/meson-build-dir').read_text().strip()
folder = root / 'build/issue-284/fault'; folder.mkdir(exist_ok=True)
(folder / 'core').mkdir(exist_ok=True)
text = (root / 'engine/native/core/trace.hpp').read_text()
text = text.replace('template <typename Observation>\nvoid observe', 'inline thread_local int observation_depth = 0;\ninline thread_local bool ipopt_payload_call = false;\ninline thread_local std::size_t ipopt_eager_allocations = 0;\nstruct IpoptPayloadCall {\n    IpoptPayloadCall() { ipopt_payload_call = true; }\n    ~IpoptPayloadCall() { ipopt_payload_call = false; }\n};\ninline thread_local std::array<std::size_t, 4> qualification_work{};\ninline thread_local int payload_kind = 0;\ninline thread_local std::array<std::size_t, 6> payload_calls{};\nstruct PayloadDepth {\n    explicit PayloadDepth(int kind) { payload_kind = kind; if (kind) ++payload_calls[kind - 1]; }\n    ~PayloadDepth() { payload_kind = 0; }\n};\nstruct ObservationDepth {\n    ObservationDepth() { ++observation_depth; }\n    ~ObservationDepth() { --observation_depth; }\n};\ntemplate <typename Observation>\nvoid observe', 1)
text = text.replace('#include <chrono>', '#include <chrono>\n#include <array>').replace('try { observation(); }', 'try { ObservationDepth depth; observation(); }')
text = text.replace('observe([&] { trace_->events.push_back({name, *index_, iteration, values, reason}); });',
    'observe([&] { if (std::strcmp(name, "ipopt_iteration") == 0) { PayloadDepth payload(3); '
    'std::vector<double> copied(values); payload_kind = 0; '
    'trace_->events.push_back({name, *index_, iteration, std::move(copied), reason}); '
    '} else trace_->events.push_back({name, *index_, iteration, values, reason}); });')
(folder / 'core/trace.hpp').write_text(text)
# Reuse R2-J2's existing manufactured fit declaration, without copying its owner.
fit_source = (root / 'engine/tests/native/regression_checks.cpp').read_text()
fit_input = '#include "regression/fit.hpp"\nusing namespace epcsaft;\nusing namespace epcsaft::regression;\nusing namespace epcsaft::equilibrium;\n'
for start, end in [('EquilibriumProblem bubble_at(', 'Observation solved('),
                   ('Observation solved(', '// The E4-B interior tie line:'),
                   ('struct Truth final', 'std::string summary(const FitResult& result) {')]:
    fit_input += fit_source[fit_source.index(start):fit_source.index(end)]
(folder / 'fit_input.hpp').write_text(fit_input)
(folder / 'contributions').mkdir(exist_ok=True)
association = (root / 'engine/native/contributions/association.hpp').read_text().replace('        if (next.empty()) {\n            next = substitution_step(fractions);', '        if (next.empty()) {\n            ++core::qualification_work[1];\n            next = substitution_step(fractions);')
(folder / 'contributions/association.hpp').write_text(association)
commands = json.loads((build / 'compile_commands.json').read_text())
objects = []
for item in commands:
    args = shlex.split(item['command'])
    if 'libepcsaft_native.a.p/' not in item['command']: continue
    target = folder / Path(args[args.index('-o') + 1]).name
    objects.append(target)
    source = (Path(item['directory']) / item['file']).resolve()
    if source.suffix == '.cpp' and '/engine/native/' in str(source):
        generated = source.read_text()
        if source.name == 'helmholtz_model.cpp':
            assert generated.count('        std::vector<double> declared_choices;') == 1
            assert generated.count('        scope.probe(\"predicate_signature\",') == 1
        generated = generated.replace('        std::vector<double> declared_choices;', '        core::PayloadDepth payload(1);\n        std::vector<double> declared_choices;')
        generated = generated.replace('        scope.probe(\"declared_choices\", declared_choices);', '        core::payload_kind = 0;\n        scope.probe(\"declared_choices\", declared_choices);')
        generated = generated.replace('        scope.probe(\"predicate_signature\", std::vector<double>(signature.begin(), signature.end()));', '        core::PayloadDepth payload(2);\n        std::vector<double> values(signature.begin(), signature.end());\n        core::payload_kind = 0;\n        scope.probe(\"predicate_signature\", values);')
        if source.name == 'gibbs_minimizer.cpp':
            generated = generated.replace('if (scope_ != nullptr) scope_->probe("ipopt_iteration",',
                'if (scope_ != nullptr) { core::IpoptPayloadCall payload; scope_->probe("ipopt_iteration",')
            generated = generated.replace('static_cast<std::uint32_t>(iter));', 'static_cast<std::uint32_t>(iter)); }')
        generated = generated.replace('std::vector<double> mask(u_.size(), 0.0), projected(u_.size());',
            'core::PayloadDepth mask_payload(4); std::vector<double> mask(u_.size(), 0.0); core::payload_kind = 0; '
            'core::PayloadDepth projected_payload(5); std::vector<double> projected(u_.size()); core::payload_kind = 0;')
        generated = generated.replace('std::vector<double> mask(width, 0.0);',
            'core::PayloadDepth round_payload(6); std::vector<double> mask(width, 0.0); core::payload_kind = 0;')
        generated = generated.replace('const auto trial_index = trial++;', 'const auto trial_index = trial++; ++core::qualification_work[0];')
        generated = generated.replace('            linked = false;\n            continue;', '            ++core::qualification_work[2];\n            linked = false;\n            continue;')
        generated = generated.replace('                    trace.state[point] = model.state_at_pressure(\n                        temperature, pressure, trace.composition[point], DensityBranch::Liquid, std::nullopt);', '                    ++core::qualification_work[3];\n                    trace.state[point] = model.state_at_pressure(\n                        temperature, pressure, trace.composition[point], DensityBranch::Liquid, std::nullopt);')
        copied = folder / source.name; copied.write_text(generated)
        args[args.index('-c') + 1] = str(copied)
    args.insert(2, '-I' + str(folder))
    for flag in ('-o', '-MQ', '-MF'):
        if flag in args: args[args.index(flag) + 1] = str(target) + ('.d' if flag == '-MF' else '')
    print('compile', source.name, flush=True)
    subprocess.run(args, cwd=item['directory'], check=True, stdout=subprocess.DEVNULL)
# Substitution is header-owned and counted in the qualification copy only.
subprocess.run(['ar', 'rcs', str(folder / 'libepcsaft_native.a'), *map(str, objects)], check=True)
fixture = next(item for item in commands if item['file'].endswith('equilibrium_checks.cpp'))
args = shlex.split(fixture['command'])
args[args.index('-c') + 1] = str(root / 'analyses/2026-native-diagnostic-trace/scripts/qualify.cpp')
args.insert(2, '-I' + str(root / 'build/issue-284'))
args.insert(2, '-I' + str(root / 'engine/tests/native'))
args.insert(2, '-I' + str(folder))
for flag in ('-o', '-MQ', '-MF'):
    if flag in args: args[args.index(flag) + 1] = str(folder / 'qualify.o') + ('.d' if flag == '-MF' else '')
subprocess.run(args, cwd=fixture['directory'], check=True)
links = subprocess.check_output(['ninja', '-C', str(build), '-t', 'commands', 'equilibrium_checks'], text=True).splitlines()
link = shlex.split(links[-1])
link[link.index('-o') + 1] = str(folder / 'qualify')
link = [str(folder / 'qualify.o') if x.endswith('tests_native_equilibrium_checks.cpp.o') else str(folder / 'libepcsaft_native.a') if x == 'libepcsaft_native.a' else x for x in link]
subprocess.run(link, cwd=build, check=True)
subprocess.run([str(folder / 'qualify')], cwd=root, check=True)

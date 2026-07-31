# Test suite index

- `unit/`: schema normalization, immutable parameter identity, and model resolution.
- `numerical/`: equation contributions, thermodynamic invariants, density closure, and
  independent finite-difference checks.
- `api/`: public constructors, bundle-path I/O, errors, and property-facing behavior.
- `sdk/`: ABI/layout, callback parity, bounded batches, and one installed C++ consumer.

All scientific cases originate in `_raw_cases.py` as compact Python tuples and
dictionaries. Tests neither locate nor import a sibling repository. Bundle files are
created only in pytest temporary directories for path-I/O tests.

# Explicit physical parameter inputs

Kind: computational-architecture

Status: accepted; current interfaces and ownership reconciled 2026-09-06.

The Engine consumes explicitly supplied physical parameters and model choices.
Reusable packets and provenance belong to `data/`; applications explicitly
materialize their selected inputs. The runtime does not discover sibling
repositories or ship a paper-specific parameter catalog.

One native admission owner resolves units, ordered components, pairs, sites,
physical choices and applicability. Python maps those inputs. Missing required
physical information is an identified failure, not an inferred coefficient.
Parameter coordinates have one identity and order through EOS derivatives,
equilibrium actions and fitting. Active families require complete derivative
chains; declaring a coordinate does not establish support.

Physical parameter fingerprints identify normalized scientific values,
component order, options and validity. Source citations and packet membership
retain provenance separately; they do not select runtime physics. Native
behavior follows physical topology and declared formulations, never application,
source, component-name or fingerprint allowlists.

The project map and GitHub issue roadmap own current scope and availability.
Residual polar terms remain excluded by ADR 0006. Historical loader names,
catalog migration records and previous package arrangements are recoverable
from this decision's Git history; they are not current callable requirements.

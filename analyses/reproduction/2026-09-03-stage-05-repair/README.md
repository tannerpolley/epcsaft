# Pinned phase amount and feed sensitivity

This retained study compares the pressure sensitivity to ethane feed against a
central finite difference for methane/ethane at 180 K, with vapor amount pinned
to 0.2 mol and pressure free. The historical directory name preserves its origin.

From the repository root, after building and installing the current Engine wheel:

```bash
.venv/bin/python analyses/reproduction/2026-09-03-stage-05-repair/scripts/reproduce_pinned_feed.py
```

The script reads the explicit methane/ethane packet from `data/packets/` and
prints its comparison. It uses native diagnostic interfaces and must be
requalified when those interfaces or the numerical implementation change.
`results/pinned-feed-sensitivity.json` retains the original calculation; moving
the script does not regenerate or qualify that historical result.

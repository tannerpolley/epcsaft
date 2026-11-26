"""Vector-scale comparison of physical quantities returned by equilibrium solves."""

STATE_TOLERANCE = 1e-10


def scaled_state_differences(a, b):
    scale = max(map(abs, (*a, *b)), default=0.0)
    return [abs(x - y) / scale if scale else 0.0 for x, y in zip(a, b, strict=True)]

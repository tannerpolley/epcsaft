"""Show per-state convergence from retained comparison.csv; no model execution."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

LEAF = Path(__file__).resolve().parents[1]


def main():
    rows = list(csv.DictReader((LEAF / 'results/v1.0.0/comparison.csv').open()))
    fig, ax = plt.subplots(figsize=(10, 2.7), layout='constrained')
    for y, column in enumerate(('previous_converged', 'release_converged')):
        for converged, marker, color in ((True, 'o', '#0072B2'), (False, 'x', '#D55E00')):
            selected = [int(r['index']) for r in rows if (r[column] == 'True') == converged]
            ax.scatter(selected, [y] * len(selected), marker=marker, color=color, s=32,
                       label=('Converged' if converged else 'Unavailable') if y == 0 else None)
    ax.set(yticks=[0, 1], yticklabels=['main-8298c602', 'v1.0.0'], ylim=(-.6, 1.8),
           xlabel='Declared family state index (dimensionless)', xlim=(-1, 73))
    ax.set_title('The release recovers every previously unavailable family state')
    ax.legend(loc='upper center', ncol=2, frameon=False)
    ax.spines[['top', 'right']].set_visible(False)
    fig.savefig(LEAF / 'figures/family-v1.0.0.png', dpi=180)
    plt.close(fig)


if __name__ == '__main__':
    main()

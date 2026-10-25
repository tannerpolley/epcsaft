from pathlib import Path
import sys

CAMPAIGN_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(CAMPAIGN_ROOT / "scripts"))
from shared.render import render_figure  # noqa: E402


if __name__ == "__main__":
    render_figure(2)

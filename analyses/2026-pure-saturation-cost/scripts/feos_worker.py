"""Run the unchanged release FeOs worker."""
from pathlib import Path
import runpy

if __name__ == '__main__':
    runpy.run_path(str(Path(__file__).resolve().parents[3] / 'analyses/2026-release-computational-cost/scripts/feos_worker.py'), run_name='__main__')

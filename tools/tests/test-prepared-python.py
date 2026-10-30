"""Prepared Windows interpreter roots must resolve without creating an environment."""
from pathlib import Path
import sys
import os
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from engine_identity import prepared_python

with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    for relative in ('python.exe', 'Scripts/python.exe', 'venv/Scripts/python.exe'):
        interpreter = root / relative
        interpreter.parent.mkdir(exist_ok=True, parents=True)
        interpreter.touch()
        interpreter.chmod(0o755)
        selected, environment = prepared_python(root, {'PATH': 'selected-only'})
        assert selected == interpreter
        assert environment == {'PATH': 'selected-only', **({'PYTHONUTF8': '1'} if os.name == 'nt' else {})}
        interpreter.unlink()
    try:
        prepared_python(root)
    except ValueError:
        pass
    else:
        raise AssertionError('missing prepared interpreter accepted')
print('Windows prepared interpreter roots and missing-interpreter refusal: pass')

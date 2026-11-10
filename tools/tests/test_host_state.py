from subprocess import CompletedProcess
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from host_state import HostState


class HostStateTest(unittest.TestCase):
    def test_open_pull_pages(self):
        state = HostState.__new__(HostState)
        state.base = "repos/owner/repo"
        state.root = Path.cwd()
        responses = [b'[{"number":1131,"head":{"ref":"issue-1010"}}]',
                     b'[{"number":1133,"head":{"ref":"issue-1124"}}]', b'[]']
        with patch("host_state.subprocess.run", side_effect=[
                CompletedProcess([], 0, stdout=body) for body in responses]) as run:
            items = state.open_items("pulls")
        self.assertEqual([(item["number"], item["head"]["ref"]) for item in items],
                         [(1131, "issue-1010"), (1133, "issue-1124")])
        self.assertEqual([call.args[0] for call in run.call_args_list],
                         [["gh", "api", f"{state.base}/pulls?state=open&per_page=100&page={page}"]
                          for page in (1, 2, 3)])


if __name__ == "__main__":
    unittest.main()

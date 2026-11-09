from io import BytesIO
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from forgejo_state import ForgejoState


class ForgejoStateTest(unittest.TestCase):
    def test_open_pull_pages(self):
        state = ForgejoState.__new__(ForgejoState)
        state.base = "https://forgejo.example/api/v1/repos/owner/repo"
        state.token = "test-token"
        responses = [b'[{"number":1131,"head":{"ref":"issue-1010"}}]',
                     b'[{"number":1133,"head":{"ref":"issue-1124"}}]', b'[]']
        with patch("forgejo_state.urlopen", side_effect=[BytesIO(body) for body in responses]) as urlopen:
            items = state.open_items("pulls")
        self.assertEqual([(item["number"], item["head"]["ref"]) for item in items],
                         [(1131, "issue-1010"), (1133, "issue-1124")])
        self.assertEqual([call.args[0].full_url for call in urlopen.call_args_list],
                         [f"{state.base}/pulls?state=open&limit=50&page={page}" for page in (1, 2, 3)])


if __name__ == "__main__":
    unittest.main()

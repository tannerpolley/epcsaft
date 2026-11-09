"""Read cleanup guards from the Forgejo repository at origin."""
import json
from pathlib import Path
import subprocess
from urllib.parse import urlsplit
from urllib.request import Request, urlopen


class ForgejoState:
    def __init__(self, root):
        remote = subprocess.run(["git", "remote", "get-url", "origin"], cwd=root, check=True,
                                capture_output=True, text=True).stdout.strip()
        url = urlsplit(remote)
        if url.scheme != "https" or not url.netloc or len(url.path.strip("/").removesuffix(".git").split("/")) != 2:
            raise ValueError("invalid origin")
        self.base = f"https://{url.netloc}/api/v1/repos/{url.path.strip('/').removesuffix('.git')}"
        self.token = json.loads((Path.home() / ".local/share/forgejo-cli/keys.json").read_text())["hosts"][url.netloc]["token"]

    def get(self, path):
        request = Request(self.base + path, headers={"Authorization": f"token {self.token}"})
        with urlopen(request, timeout=15) as response:
            return json.load(response)

    def open_items(self, kind):
        items, page = [], 1
        while True:
            batch = self.get(f"/{kind}?state=open&limit=50&page={page}")
            if not isinstance(batch, list):
                raise ValueError("invalid Forgejo list")
            items.extend(batch)
            if not batch:
                return items
            page += 1

    def issue(self, number):
        return self.get(f"/issues/{int(number)}")

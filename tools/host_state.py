"""Read cleanup guards from the host repository at origin."""
import json
import subprocess
from urllib.parse import urlsplit


class HostState:
    def __init__(self, root):
        remote = subprocess.run(["git", "remote", "get-url", "origin"], cwd=root, check=True,
                                capture_output=True, text=True).stdout.strip()
        url = urlsplit(remote)
        if url.scheme != "https" or url.netloc != "github.com" or len(url.path.strip("/").removesuffix(".git").split("/")) != 2:
            raise ValueError("invalid origin")
        self.base = f"repos/{url.path.strip('/').removesuffix('.git')}"
        self.root = root

    def get(self, path):
        return json.loads(subprocess.run(["gh", "api", self.base + path], cwd=self.root, check=True,
                                         capture_output=True, text=True).stdout)

    def open_items(self, kind):
        items, page = [], 1
        while True:
            batch = self.get(f"/{kind}?state=open&per_page=100&page={page}")
            if not isinstance(batch, list):
                raise ValueError("invalid host list")
            items.extend(batch)
            if not batch:
                return items
            page += 1

    def issue(self, number):
        return self.get(f"/issues/{int(number)}")

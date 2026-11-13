"""Read cleanup guards from the host repository at origin and the local worktree list."""
import json
import os
import subprocess
from urllib.parse import urlsplit


class HostState:
    def __init__(self, root):
        remote = subprocess.run(["git", "remote", "get-url", "origin"], cwd=root, check=True,
                                capture_output=True, text=True).stdout.strip()
        url = urlsplit(remote)
        if url.scheme != "https" or url.netloc != "github.com" or len(url.path.strip("/").removesuffix(".git").split("/")) != 2:
            raise ValueError(f"invalid origin {remote!r}")
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


def worktrees(directory):
    """Registered worktrees of the repository holding `directory`, main first: {real path: {field: value}}."""
    listing = subprocess.run(["git", "-C", str(directory), "worktree", "list", "--porcelain", "-z"], check=True,
                             capture_output=True).stdout
    entries = {}
    for record in filter(None, listing.split(b"\0\0")):
        fields = [os.fsdecode(field) for field in record.strip(b"\0").split(b"\0")]
        entries[os.path.realpath(fields[0].removeprefix("worktree "))] = dict(f.partition(" ")[::2] for f in fields[1:])
    return entries


def main_checkout(directory):
    return next(iter(worktrees(directory)))

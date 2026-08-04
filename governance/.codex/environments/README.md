# Codex local environment

The ChatGPT desktop app runs `setup.sh` automatically for a newly created
Codex worktree.

The setup is intentionally narrow:

- create a repository-local Python 3.13 environment;
- generate ignored local Governance state when absent; and
- verify the Governance and access contracts without dirtying the worktree.

It does not mutate Git state, install sibling source trees, publish artifacts,
or change GitHub.

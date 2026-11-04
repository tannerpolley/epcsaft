#!/usr/bin/env bash
# Restore the FeOs packets that source archives omit under their CC BY-NC terms: download the four
# pinned FeOs parameter files, rebuild the feos-*/1 packets (import_feos.py refuses any file whose
# SHA-256 differs from its pin) and regenerate the parameter index. No-op when both are present.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
if [[ -f data/packets/feos-esper-2023/1/packet.toml && -f data/packets/feos-rehner-2023/1/packet.toml ]]; then
    exit 0
fi
commit=$(sed -n 's/^REPO, COMMIT, TREE = "[^"]*", "\([0-9a-f]*\)".*/\1/p' data/tools/import_feos.py)
: "${commit:?no pinned FeOs commit in data/tools/import_feos.py}"
downloads=$(mktemp -d)
trap 'rm -rf "$downloads"' EXIT
for name in gross2001 gross2002 esper2023 rehner2023_binary; do
    curl -fsSL -o "$downloads/$name.json" \
        "https://raw.githubusercontent.com/feos-org/feos/$commit/parameters/pcsaft/$name.json"
done
python3 data/tools/import_feos.py "$downloads"
python3 data/tools/parameter_index.py > data/docs/parameter-index.md

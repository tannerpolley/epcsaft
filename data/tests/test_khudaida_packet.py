from __future__ import annotations

import csv
import tomllib
import unittest
from pathlib import Path

from tools.validate_packets import packet_fingerprint, validate_packet


ROOT = Path(__file__).resolve().parents[1]
PACKET = (
    ROOT
    / "packets"
    / "khudaida-2026-figure-2-electrolyte-lle"
    / "2"
)


class KhudaidaPacketTests(unittest.TestCase):
    def test_packet_is_valid_and_fingerprint_is_stable(self) -> None:
        self.assertEqual(validate_packet(PACKET / "packet.toml"), [])
        self.assertEqual(
            packet_fingerprint(PACKET),
            "ed173d78fcaeca527c89a49af1f34c6880fdc642fc1c08c2d395f131991e53af",
        )

    def test_packet_retains_all_eight_normalized_feeds(self) -> None:
        path = PACKET / "observations" / "figure-02-feed-compositions.csv"
        with path.open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual([int(row["point_id"]) for row in rows], list(range(1, 9)))
        for row in rows:
            total = sum(
                float(row[column])
                for column in (
                    "x_ethanol_salt_free",
                    "x_isobutanol_salt_free",
                    "x_water_salt_free",
                )
            )
            self.assertAlmostEqual(total, 1.0, places=12)

    def test_source_points_cover_reported_tielines_feeds_and_paper_curve(self) -> None:
        path = PACKET / "observations" / "figure-02-source-points.csv"
        with path.open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        counts: dict[str, int] = {}
        for row in rows:
            counts[row["dataset"]] = counts.get(row["dataset"], 0) + 1
        self.assertEqual(counts["experimental_tieline"], 16)
        self.assertEqual(counts["feed"], 8)
        self.assertEqual(counts["paper_epcsaft"], 8)

    def test_provenance_pins_validation_source_without_results(self) -> None:
        provenance = tomllib.loads(
            (PACKET / "provenance" / "validation-origin.toml").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            provenance["origin_commit"],
            "dc1841fc8c7b837d92b7fed9d976a22ba6d5b2e6",
        )
        self.assertIn("Generated results", provenance["scope"])
        self.assertFalse(any(path.name.startswith("results") for path in PACKET.rglob("*")))


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

import csv
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PACKET_ROOT = ROOT / "packets" / "khudaida-2026-figure-2-electrolyte-lle"
PACKET_V2 = PACKET_ROOT / "2"
PACKET = PACKET_ROOT / "3"
PACKET_V4 = PACKET_ROOT / "4"


class KhudaidaPacketTests(unittest.TestCase):
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


if __name__ == "__main__":
    unittest.main()

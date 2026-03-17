import csv
import gzip
import sys
import tempfile
import unittest
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = PROJECT_ROOT / "python"
if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

from select_backbone_by_variation import run  # noqa: E402


class SelectBackboneByVariationTest(unittest.TestCase):
    def test_mock_pipeline_outputs_nested_distance_driven_backbone(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            temp_root = Path(tmpdir)
            meta_path = temp_root / "meta.tsv"
            basic_info_path = temp_root / "基础信息.tsv"
            vcf_path = temp_root / "mock.vcf.gz"
            output_dir = temp_root / "output"

            self._write_meta(meta_path)
            self._write_basic_info(basic_info_path)
            self._write_vcf(vcf_path)

            run(
                meta_path=meta_path,
                vcf_path=vcf_path,
                basic_info_path=basic_info_path,
                output_dir=output_dir,
                tiers=[5, 6],
                distance_mode="alt_hamming",
                seed_mode="deep_lineage_cover",
                deep_lineage_labels_spec="auto",
                min_quality=0.0,
            )

            backbone_5 = self._read_tsv(output_dir / "backbone_5.tsv")
            backbone_6 = self._read_tsv(output_dir / "backbone_6.tsv")
            master = self._read_tsv(output_dir / "05_backbone_selection_master.tsv")
            seed_summary = self._read_tsv(output_dir / "03_deep_lineage_seed_summary.tsv")
            clean_pool = self._read_tsv(output_dir / "02_clean_candidate_pool.tsv")

            self.assertEqual([row["SampleID"] for row in backbone_5], [row["SampleID"] for row in backbone_6[:5]])
            self.assertEqual([row["SampleID"] for row in backbone_6], [row["SampleID"] for row in master[:6]])

            top_five_labels = {row["deep_lineage_label"] for row in backbone_5}
            self.assertTrue({"L0", "L1", "L2", "L3"}.issubset(top_five_labels))

            seed_labels = {row["deep_lineage_label"] for row in seed_summary if row["selected_sample_id"]}
            self.assertTrue({"L0", "L1", "L2", "L3"}.issubset(seed_labels))

            clean_ids = {row["SampleID"] for row in clean_pool}
            self.assertIn("S_A1", clean_ids)
            self.assertNotIn("S_A2", clean_ids)
            self.assertNotIn("S_QLOW", clean_ids)
            self.assertNotIn("S_NOT_FULL", clean_ids)

            selected_rows = {row["SampleID"]: row for row in master}
            self.assertEqual(selected_rows["S_L1"]["selection_stage"], "seed_cover")
            self.assertEqual(selected_rows["S_L2"]["selection_stage"], "seed_cover")
            self.assertEqual(selected_rows["S_L3"]["selection_stage"], "seed_cover")
            self.assertEqual(selected_rows["S_B1"]["selection_stage"], "diversity_expansion")

            self.assertTrue((output_dir / "04_distance_selection_summary.tsv").exists())
            self.assertTrue((output_dir / "README_筛选说明.md").exists())
            self.assertFalse((output_dir / "03_root_haplogroup_summary.tsv").exists())
            self.assertFalse((output_dir / "04_level2_haplogroup_summary.tsv").exists())

    def _read_tsv(self, path: Path):
        with path.open("r", encoding="utf-8", newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def _write_meta(self, path: Path) -> None:
        rows = [
            ["SampleID", "Haplogroup", "Rank", "Quality", "Range"],
            ["S_L0", "L0a1", "1", "1.0000", "1-16569"],
            ["S_L1", "L1b1", "1", "1.0000", "1-16569"],
            ["S_L2", "L2a1", "1", "1.0000", "1-16569"],
            ["S_L3", "L3e1", "1", "1.0000", "1-16569"],
            ["S_A1", "A1", "1", "1.0000", "1-16569"],
            ["S_A1", "A1", "1", "0.9400", "1-16569"],
            ["S_A2", "A2", "1", "1.0000", "1-16569"],
            ["S_B1", "B1", "1", "1.0000", "1-16569"],
            ["S_QLOW", "Q1", "1", "0.9500", "1-16569"],
            ["S_NOT_FULL", "C1", "1", "0.9200", "1-576 16024-16569"],
        ]
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerows(rows)

    def _write_basic_info(self, path: Path) -> None:
        rows = [
            ["ID", "Source"],
            ["S_L0", "PUBLIC_A"],
            ["S_L1", "PUBLIC_A"],
            ["S_L2", "PUBLIC_A"],
            ["S_L3", "PUBLIC_A"],
            ["S_A1", "PUBLIC_B"],
            ["S_A2", "LLT_INTERNAL"],
            ["S_B1", "PUBLIC_C"],
            ["S_QLOW", "PUBLIC_D"],
            ["S_NOT_FULL", "PUBLIC_E"],
        ]
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerows(rows)

    def _write_vcf(self, path: Path) -> None:
        samples = ["S_L0", "S_L1", "S_L2", "S_L3", "S_A1", "S_A2", "S_B1", "S_QLOW"]
        variant_rows = [
            ("1", "A", "G", ["1/1", "0/0", "0/0", "0/0", "0/0", "0/0", "1/1", "0/0"]),
            ("2", "C", "T", ["0/0", "1/1", "0/0", "0/0", "0/0", "./.", "1/1", "0/0"]),
            ("3", "G", "A", ["0/0", "0/0", "1/1", "0/0", "0/0", "0/0", "1/1", "0/0"]),
            ("4", "T", "C", ["0/0", "0/0", "0/0", "1/1", "0/0", "0/0", "1/1", "0/0"]),
            ("5", "A", "C", ["0/0", "0/0", "0/0", "0/0", "1/1", "1/1", "0/0", "0/0"]),
            ("6", "G", "A,T", ["0/0", "0/0", "0/0", "0/0", "0/0", "0/1", "1/2", "0/0"]),
        ]
        with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
            handle.write("##fileformat=VCFv4.2\n")
            handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            header = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", *samples]
            handle.write("\t".join(header) + "\n")
            for pos, ref, alt, genotypes in variant_rows:
                row = ["chrM", pos, ".", ref, alt, ".", "PASS", ".", "GT", *genotypes]
                handle.write("\t".join(row) + "\n")


if __name__ == "__main__":
    unittest.main()

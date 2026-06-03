#!/usr/bin/env python3
"""1-1-tsv_filter.py 的 CLI 默认值测试。"""

import importlib.util
import tempfile
import unittest
from pathlib import Path

import pandas as pd


MODULE_PATH = Path(__file__).resolve().parents[1] / "python" / "1-1-tsv_filter.py"
SPEC = importlib.util.spec_from_file_location("tsv_filter", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class TsvFilterTest(unittest.TestCase):
    def test_cli_default_qc_threshold_matches_pipeline_default(self):
        args = MODULE.build_parser().parse_args(["--tsv", "in.tsv", "--output", "out.tsv"])

        self.assertEqual(args.qc_min, 0.85)

    def test_filter_accepts_prefiltered_input_without_sequence_method(self):
        data = pd.DataFrame(
            {
                "ID": ["keep", "drop"],
                "QC_Other": ["YES", "NO"],
                "QC_Haplogrep": [0.9, 0.9],
                "Haplogroup_17.2": ["A", "B"],
            }
        )
        with tempfile.TemporaryDirectory() as tmp:
            input_path = Path(tmp) / "input.tsv"
            data.to_csv(input_path, sep="\t", index=False)

            observed = MODULE.load_and_filter(str(input_path), qc_min=0.85)

        self.assertEqual(observed["ID"].tolist(), ["keep"])


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3
"""1-2-vcf_fps_select.py 的核心逻辑测试。"""

import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


MODULE_PATH = (
    Path(__file__).resolve().parents[1] / "python" / "1-2-vcf_fps_select.py"
)
SPEC = importlib.util.spec_from_file_location("vcf_fps_select", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class VcfFpsSelectTest(unittest.TestCase):
    def test_segregating_mask_ignores_missing_genotypes(self):
        matrix = np.array(
            [
                [0, 0, -1],
                [0, 1, -1],
                [-1, -1, -1],
                [2, 2, 2],
            ],
            dtype=np.int8,
        )

        observed = MODULE.get_segregating_mask(matrix)

        np.testing.assert_array_equal(observed, [False, True, False, False])

    def test_pairwise_diff_ignores_sites_missing_in_either_sample(self):
        samples_by_sites = np.array(
            [
                [0, 0, -1, 2],
                [1, 0, 1, 2],
                [0, -1, 1, 1],
            ],
            dtype=np.int8,
        )

        observed = MODULE.pairwise_diff(samples_by_sites)

        expected = np.array(
            [
                [0, 1, 1],
                [1, 0, 2],
                [1, 2, 0],
            ],
            dtype=np.int16,
        )
        np.testing.assert_array_equal(observed, expected)

    def test_pairwise_diff_can_run_after_dynamic_import_in_new_process(self):
        code = (
            "import runpy, numpy as np; "
            f"m=runpy.run_path({str(MODULE_PATH)!r}); "
            "d=m['pairwise_diff'](np.array([[0, 1], [1, 1]], dtype=np.int8)); "
            "assert d.tolist() == [[0, 1], [1, 0]]"
        )

        observed = subprocess.run(
            [sys.executable, "-c", code],
            capture_output=True,
            text=True,
        )

        self.assertEqual(observed.returncode, 0, observed.stderr)

    def test_collapse_identical_requires_same_missing_pattern(self):
        samples_by_sites = np.array(
            [
                [0, -1, 1],
                [0, -1, 1],
                [0, 0, 1],
            ],
            dtype=np.int8,
        )

        reps, cluster_sizes = MODULE.collapse_identical(samples_by_sites)

        np.testing.assert_array_equal(reps, [0, 2])
        self.assertEqual(cluster_sizes, {0: 2, 2: 1})

    def test_process_block_keeps_missing_patterns_when_counting_haplotypes(self):
        block_df = pd.DataFrame(
            {
                "ID": ["sample_a", "sample_b"],
                "Haplogroup_17.2": ["A", "A"],
                "QC_Haplogrep": [0.9, 0.9],
            }
        )
        sites_by_samples = np.array([[0, -1], [0, 0]], dtype=np.int8)

        with tempfile.TemporaryDirectory() as tmp:
            _, _, summary = MODULE.process_block(
                "A",
                block_df,
                sites_by_samples,
                np.array([0, 1]),
                min_dist=1,
                max_tips=0,
                temp_dir=tmp,
                numba_threads=1,
            )

        self.assertEqual(summary["n_unique_haplotypes"], 2)

    def test_prepare_resume_rejects_changed_fingerprint(self):
        with tempfile.TemporaryDirectory() as tmp:
            log_dir = Path(tmp) / "log"
            temp_dir = Path(tmp) / "temp"
            first = {"vcf": "a.vcf.gz", "min_dist": 1}
            changed = {"vcf": "a.vcf.gz", "min_dist": 2}

            MODULE.prepare_resume(log_dir, temp_dir, first, overwrite=False)

            with self.assertRaisesRegex(RuntimeError, "运行指纹"):
                MODULE.prepare_resume(log_dir, temp_dir, changed, overwrite=False)


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3
"""校准/UCLD BEAST XML 渲染的最小回归测试。"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "python"))

from beast_common import render_calibrated_beast_xml


class CalibratedBeastXmlTest(unittest.TestCase):
    def test_render_calibrated_xml_replaces_calibration_placeholders(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            template = PROJECT_ROOT / "conf/beast_template_calibrated_backbone.xml"
            destination = tmp / "beast.xml"

            render_calibrated_beast_xml(
                template_path=template,
                taxa_block='\t\t<taxon id="RSRS"/>',
                alignment_block=(
                    '\t\t<sequence>\n'
                    '\t\t\t<taxon idref="RSRS"/>\n'
                    '\t\t\tACGT\n'
                    '\t\t</sequence>'
                ),
                starting_tree_newick="(RSRS:180000,A:180000);",
                calibration_block=(
                    '\t<taxa id="M_clade">\n\t\t<taxon idref="A"/>\n\t</taxa>\n'
                    '\t<taxa id="N_clade">\n\t\t<taxon idref="RSRS"/>\n\t</taxa>\n'
                    '\t<tmrcaStatistic id="tmrca_M" absolute="true">\n'
                    '\t\t<mrca><taxa idref="M_clade"/></mrca>\n'
                    '\t\t<treeModel idref="treeModel"/>\n'
                    '\t</tmrcaStatistic>\n'
                    '\t<tmrcaStatistic id="tmrca_N" absolute="true">\n'
                    '\t\t<mrca><taxa idref="N_clade"/></mrca>\n'
                    '\t\t<treeModel idref="treeModel"/>\n'
                    '\t</tmrcaStatistic>'
                ),
                file_prefix="beast_subtree_0001",
                chain_length=12345,
                log_every=100,
                ucld_mean=2.5e-8,
                popsize=100000.0,
                root_age=180000.0,
                root_stdev=20000.0,
                m_age=50000.0,
                m_stdev=5000.0,
                n_age=57000.0,
                n_stdev=5000.0,
                destination=destination,
            )

            xml = destination.read_text(encoding="utf-8")
            self.assertNotIn("{{", xml)
            self.assertIn('<discretizedBranchRates id="branchRates">', xml)
            self.assertIn('fileName="beast_subtree_0001.trees"', xml)
            self.assertIn('<normalPrior mean="180000" stdev="20000">', xml)
            self.assertIn('<normalPrior mean="50000" stdev="5000">', xml)
            self.assertIn('<normalPrior mean="57000" stdev="5000">', xml)
            self.assertIn('<statistic idref="tmrca_M"/>', xml)
            self.assertIn('<statistic idref="tmrca_N"/>', xml)


if __name__ == "__main__":
    unittest.main()

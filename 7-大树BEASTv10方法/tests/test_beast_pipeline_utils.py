import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "python"))

import beast_pipeline_utils as utils  # noqa: E402


def test_normalize_haplogroup_removes_mutation_suffix():
    assert utils.normalize_haplogroup("C4a1a_T195C!") == "C4a1a"
    assert utils.normalize_haplogroup("M65a_C16311T!!") == "M65a"
    assert utils.normalize_haplogroup(" HV_A73G! ") == "HV"


def test_haplogroup_lineage_membership_uses_phylotree_json():
    haplogroups = {
        "M": {"lineage": ["mt-MRCA(RSRS)", "L3", "M"]},
        "C4a1a": {"lineage": ["mt-MRCA(RSRS)", "L3", "M", "C", "C4", "C4a", "C4a1a"]},
        "H": {"lineage": ["mt-MRCA(RSRS)", "L3", "N", "R", "HV", "H"]},
    }
    assert utils.haplogroup_in_node("C4a1a_T195C!", "M", haplogroups)
    assert utils.haplogroup_in_node("H1", "M", haplogroups) is False


def test_tree_monophyly_reports_extra_taxa(tmp_path):
    tree_path = tmp_path / "tree.nwk"
    tree_path.write_text("((a:1,b:1):1,c:2);", encoding="utf-8")
    membership = {"a", "c"}
    result = utils.check_monophyly(tree_path, membership)
    assert result["is_monophyletic"] is False
    assert result["extra_taxa"] == ["b"]


def test_xml_template_has_required_thorney_blocks():
    xml = utils.render_beast_xml(
        template_path=PROJECT_ROOT / "conf" / "beast_v10_thorney.template.xml",
        taxa=["a", "b"],
        starting_tree="(a:1,b:1):1;",
        data_tree="(a:0.1,b:0.1):0.1;",
        calibration_taxa={"M": ["a"], "N": ["b"]},
        calibration={"root": {"age": 180000, "stdev": 20000}, "M": {"age": 49000, "stdev": 5000}, "N": {"age": 59000, "stdev": 5000}},
        scale=17590,
        clock_rate=1.0e-8,
        chain_length=1000,
        log_every=100,
        tree_log_every=100,
        output_prefix="test",
        skygrid_points=10,
        skygrid_cutoff=198000,
        skygrid_init_logpop=10.0,
    )
    assert "<constrainedTreeModel id=\"treeModel\">" in xml
    assert "<thorneyTreeLikelihood id=\"treeLikelihood\">" in xml
    assert "<bigFastTreeIntervals>" in xml
    assert "<bigFastTreeModel" not in xml
    assert "<treeDataLikelihood" not in xml
    assert "normalPrior id=\"M.prior\"" in xml
    # BEAST v10.5.0 的正确元素名是 simpleMutationBranchMap，而非 mutationBranchMap
    assert "<simpleMutationBranchMap id=\"branchMutationCounts\"" in xml
    assert "<simpleMutationBranchMap idref=\"branchMutationCounts\"" in xml
    assert ">mutationBranchMap" not in xml.replace("simpleMutationBranchMap", "")
    assert "scale=\"17590.0\"" in xml


def test_lsd2_date_file_uses_present_zero_and_negative_bp(tmp_path):
    out = tmp_path / "dates.txt"
    utils.write_lsd2_date_file(
        taxa=["a", "b", "c"],
        calibration_taxa={"M": ["a", "b"], "N": ["c"]},
        node_ages={"M": 49000, "N": 59000},
        output_path=out,
    )
    lines = out.read_text(encoding="utf-8").strip().splitlines()
    assert lines[0] == "a\t0"
    assert "a,b\t-49000" in lines
    assert "c\t-59000" in lines


def test_bayesian_skyline_template_uses_generalized_skyline():
    xml = utils.render_beast_xml(
        template_path=PROJECT_ROOT / "conf" / "beast_v10_thorney_bayesian_skyline.template.xml",
        taxa=["a", "b"],
        starting_tree="(a:1,b:1):1;",
        data_tree="(a:0.1,b:0.1):0.1;",
        calibration_taxa={"M": ["a"], "N": ["b"]},
        calibration={"root": {"age": 180000, "stdev": 20000}, "M": {"age": 49000, "stdev": 5000}, "N": {"age": 59000, "stdev": 5000}},
        scale=17590,
        clock_rate=1.0e-8,
        chain_length=1000,
        log_every=100,
        tree_log_every=100,
        output_prefix="test",
        skygrid_points=10,
        skygrid_cutoff=198000,
        skygrid_init_logpop=10.0,
    )
    assert "<generalizedSkyLineLikelihood id=\"skyline\" linear=\"false\">" in xml
    assert "<exponentialMarkovLikelihood id=\"skyline.popSize.prior\" jeffreys=\"true\">" in xml
    assert "<gmrfSkyGridLikelihood" not in xml
    assert "id=\"skyline.popSize\" dimension=\"10\"" in xml
    assert "upper=\"Infinity\"" in xml


def test_read_calibration_taxa_groups_step1_output(tmp_path):
    table = tmp_path / "calibration_taxa.tsv"
    table.write_text("node\tsample_id\thaplogroup\nM\ta\tM1\nN\tb\tN1\n", encoding="utf-8")
    assert utils.read_calibration_taxa(table) == {"M": ["a"], "N": ["b"]}

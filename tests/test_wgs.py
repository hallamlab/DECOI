"""Scientific linkage and paired-study regression tests."""
import copy
import gzip
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from wgs import allocate_reads, load_genomes, simulate_wgs
from scripts.create_paired_smoke import create
import mock16s_chem as m


def test_cell_and_dna_weighting_and_empty_sample():
    counts = pd.DataFrame({"sample": [20, 20], "blank": [0, 0]}, index=["a", "b"])
    registry = pd.DataFrame({"marker_copies": [2, 1], "genome_length": [2000, 4000]}, index=counts.index)
    cells, fractions, pairs = allocate_reads(counts, registry, 1000, np.random.default_rng(42))
    assert cells["sample"].tolist() == [10, 20]
    assert fractions["sample"].tolist() == [0.2, 0.8]
    assert pairs.sum().tolist() == [1000, 0]
    assert pairs.equals(allocate_reads(counts, registry, 1000, np.random.default_rng(42))[2])
    with pytest.raises(ValueError, match="nonnegative integer"):
        allocate_reads(counts, registry, -1, np.random.default_rng(42))


def test_reference_requires_exact_link_and_copy_number(tmp_path):
    marker = "AACCGGTTACGTACGT"
    (tmp_path / "g.fa").write_text(">contig\n" + "A"*1100 + str(m.Seq(marker).reverse_complement()) + "A"*1100 + "\n")
    manifest = tmp_path / "genomes.tsv"
    manifest.write_text("ASV_ID\tgenome_id\tfasta\tmarker_copies\na\tg\tg.fa\t1\n")
    reg, genomes = load_genomes(tmp_path, {"a": marker})
    assert reg.loc["a", "genome_length"] == 2216
    assert genomes["a"][0].id == "g__contig_1"
    manifest.write_text("ASV_ID\tgenome_id\tfasta\tmarker_copies\na\tg\tg.fa\t2\n")
    with pytest.raises(ValueError, match="exact V4 occurrences"):
        load_genomes(tmp_path, {"a": marker})
    with pytest.raises(ValueError, match="absent"):
        load_genomes(tmp_path, {})


def test_reference_rejects_ambiguous_mapping_and_outside_bundle(tmp_path):
    manifest = tmp_path / "genomes.tsv"
    header = "ASV_ID\tgenome_id\tfasta\tmarker_copies\n"
    manifest.write_text(header + "a\tg\tg.fa\t1\na\th\th.fa\t1\n")
    with pytest.raises(ValueError, match="one-to-one"):
        load_genomes(tmp_path, {"a": "ACGT"})
    manifest.write_text(header + "a\tg\t../outside.fa\t1\n")
    with pytest.raises(ValueError, match="inside"):
        load_genomes(tmp_path, {"a": "ACGT"})


def test_wgs_uses_pre_pcr_state_and_excludes_amplicon_artifacts(tmp_path):
    biological = pd.DataFrame({"s": [10, 30]}, index=["a", "b"])
    final = pd.DataFrame({"s": [900, 1, 5, 200, 300], "blank": [0, 0, 20, 4, 5]},
                         index=["a", "b", "c", "chimera", "mito"])
    reg = pd.DataFrame({"genome_id": ["ga", "gb", "gc"], "marker_copies": [1, 1, 1],
                        "genome_length": [2000, 2000, 2000]}, index=["a", "b", "c"])
    genomes = {a: [m.SeqRecord(m.Seq("ACGT"*500), id=a)] for a in reg.index}
    cfg = {"seed": 42, "wgs": {"read_pairs": 0}}
    def unexpected(*args):
        raise AssertionError("ISS must not be called for empty libraries")
    summary = simulate_wgs(cfg, biological, final, ["c"], reg, genomes, tmp_path,
                           unexpected, m.count_fastq_records, m.stable_uuid)
    shared = pd.read_csv(tmp_path / "wgs/shared_marker_counts.tsv", sep="\t", index_col=0)
    assert shared.s.tolist() == [10, 30, 5]
    assert shared.blank.tolist() == [0, 0, 20]
    assert summary["total_read_pairs"] == 0
    assert summary["samples"] == 2
    for path in (tmp_path / "wgs/fastq").glob("*.gz"):
        with gzip.open(path, "rt") as handle:
            assert handle.read() == ""


@pytest.mark.skipif(os.environ.get("DECOI_RUN_INTEGRATION") != "1", reason="Set DECOI_RUN_INTEGRATION=1 with Cutadapt, R/SparseDOSSA2 and ISS on PATH")
def test_real_paired_study_and_amplicon_regression(tmp_path):
    bundle = tmp_path / "fixture"
    create(bundle)
    cfg = yaml.safe_load((bundle / "config.yaml").read_text())
    m.prepare_reference(cfg)
    paired = bundle / "paired"
    m.simulate(copy.deepcopy(cfg), paired)
    cfg["wgs"]["enabled"] = False
    amplicon = bundle / "amplicon"
    m.simulate(copy.deepcopy(cfg), amplicon)
    for name in ("asv_counts_biological.tsv", "asv_counts_final.tsv", "chemistry.tsv", "sample_metadata.tsv", "ground_truth_feature_registry.tsv"):
        assert (paired / name).read_bytes() == (amplicon / name).read_bytes(), name
    for path in (paired / "fastq").glob("*.gz"):
        with gzip.open(path, "rb") as a, gzip.open(amplicon / "fastq" / path.name, "rb") as b:
            assert a.read() == b.read()
    manifest = pd.read_csv(paired / "assay_manifest.tsv", sep="\t")
    assert manifest.groupby("sample_id").assay.apply(set).tolist() == [{"amplicon", "wgs"}]*3
    assert all((paired / path).exists() for path in manifest.fastq_r1)
    checks = pd.read_csv(paired / "wgs/fastq_validation.tsv", sep="\t")
    assert checks.expected_pairs.tolist() == [100]*3
    assert checks.observed_R1.equals(checks.expected_pairs)
    assert checks.observed_R2.equals(checks.expected_pairs)
    # Extend the existing amplicon study in place, preserving its original data.
    protected={path.relative_to(amplicon): path.read_bytes() for path in amplicon.rglob("*")
               if path.is_file() and path.name not in {"manifest.json", "report.html"}}
    m.add_wgs(cfg, amplicon)
    for relative, original in protected.items():
        assert (amplicon / relative).read_bytes() == original, relative
    assert (amplicon / "wgs/genome_read_pairs.tsv").read_bytes() == (paired / "wgs/genome_read_pairs.tsv").read_bytes()
    with pytest.raises(ValueError, match="already has WGS"):
        m.add_wgs(cfg, amplicon)

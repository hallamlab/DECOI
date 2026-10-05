"""Test genome-derived reference provenance and allele-selection constraints."""
import copy
import os

import pandas as pd
import pytest
import yaml
from Bio.Seq import Seq

from genome_reference import extract_markers, prepare_genome_reference
from scripts.create_paired_smoke import create
import mock16s_chem as m


def test_extract_both_strands_multiple_loci_and_length_bounds():
    forward, reverse, marker = "ACGTRT", "CGTACC", "GATTACAGATTACA"
    molecule = "ACGTAT" + marker + str(Seq(reverse).reverse_complement())
    genome = "T"*100 + molecule + "T"*100 + str(Seq(molecule).reverse_complement()) + "T"*100
    hits = extract_markers(genome, forward, reverse, len(marker), len(marker))
    assert len(hits) == 2
    assert {h[1] for h in hits} == {"+", "-"}
    for sequence, strand, start, end in hits:
        assert sequence == marker
        interval = genome[start:end]
        assert (interval if strand == "+" else str(Seq(interval).reverse_complement())) == marker
    assert extract_markers(genome, forward, reverse, 2, 3) == []


def test_prepare_without_prelinked_asv_list(tmp_path):
    create(tmp_path / "fixture", genome_first=True)
    cfg = yaml.safe_load((tmp_path / "fixture/config.yaml").read_text())
    source = tmp_path / "fixture/genomes"
    assert not (source / "genomes.tsv").exists()
    source_taxonomy=pd.read_csv(source / "taxonomy.tsv",sep="\t")
    source_taxonomy["reference_role"]="biological"
    source_taxonomy.loc[0,"reference_role"]="contaminant"
    source_taxonomy["accession"]=[f"TEST_{i}" for i in range(len(source_taxonomy))]
    source_taxonomy.to_csv(source / "taxonomy.tsv",sep="\t",index=False)
    metadata = prepare_genome_reference(cfg, source, tmp_path / "reference")
    assert metadata["source_type"] == "genomes"
    assert metadata["n_unique_v4"] == 8
    table = pd.read_csv(tmp_path / "reference/genomes/genomes.tsv", sep="\t")
    markers = {r.id: str(r.seq) for r in m.SeqIO.parse(tmp_path / "reference/v4_asv_pool.fasta", "fasta")}
    registry, _ = m.load_genomes(tmp_path / "reference/genomes", markers)
    assert set(registry.index) == set(table.ASV_ID)
    assert registry.marker_copies.eq(1).all()
    assert registry.reference_role.value_counts().to_dict()=={"biological":7,"contaminant":1}
    assert set(registry.accession)==set(source_taxonomy.accession)
    prepared_taxonomy=pd.read_csv(tmp_path / "reference/v4_asv_taxonomy.tsv",sep="\t")
    assert prepared_taxonomy.reference_role.value_counts().to_dict()=={"biological":7,"contaminant":1}
    with pytest.raises(ValueError, match="fresh"):
        prepare_genome_reference(cfg, source, tmp_path / "reference")


def test_exclusions_are_audited(tmp_path):
    create(tmp_path / "fixture", genome_first=True)
    cfg = yaml.safe_load((tmp_path / "fixture/config.yaml").read_text())
    source = tmp_path / "fixture/genomes"
    # Duplicate marker genome, ambiguous genome, and a genome without primers.
    (source / "zz_duplicate.fa").write_bytes((source / "genome_0.fa").read_bytes())
    taxonomy = pd.read_csv(source / "taxonomy.tsv", sep="\t")
    duplicate = taxonomy.iloc[[0]].copy().assign(fasta="zz_duplicate.fa")
    pd.concat([taxonomy, duplicate]).to_csv(source / "taxonomy.tsv", sep="\t", index=False)
    (source / "ambiguous.fa").write_text(">ambiguous\n" + "N"*2000 + "\n")
    (source / "no_marker.fa").write_text(">none\n" + "A"*2000 + "\n")
    # Two different alleles in a single assembly must not silently become one.
    (source / "mixed.fa").write_text((source / "genome_0.fa").read_text() + (source / "genome_1.fa").read_text())
    prepare_genome_reference(cfg, source, tmp_path / "reference")
    audit = pd.read_csv(tmp_path / "reference/genome_selection.tsv", sep="\t").set_index("source_fasta")
    assert audit.loc["zz_duplicate.fa", "status"] == "V4_shared_with_previously_retained_genome"
    assert audit.loc["ambiguous.fa", "status"].startswith("requires_contigs")
    assert audit.loc["no_marker.fa", "status"] == "no_exact_primer_bounded_V4"
    assert audit.loc["mixed.fa", "status"] == "multiple_V4_alleles_not_yet_supported"


@pytest.mark.skipif(os.environ.get("DECOI_RUN_INTEGRATION") != "1", reason="Requires R/SparseDOSSA2 and InSilicoSeq")
def test_genome_first_real_paired_run(tmp_path):
    create(tmp_path / "fixture", genome_first=True)
    cfg = yaml.safe_load((tmp_path / "fixture/config.yaml").read_text())
    m.prepare_reference(cfg)
    output = tmp_path / "output"
    m.simulate(copy.deepcopy(cfg), output)
    links = pd.read_csv(output / "assay_manifest.tsv", sep="\t")
    assert links.groupby("sample_id").assay.apply(set).tolist() == [{"amplicon", "wgs"}]*3
    wgs = pd.read_csv(output / "wgs/fastq_validation.tsv", sep="\t")
    assert wgs.expected_pairs.eq(100).all()
    assert wgs.expected_pairs.equals(wgs.observed_R1)
    assert wgs.expected_pairs.equals(wgs.observed_R2)
    assert pd.read_csv(output / "wgs/genome_registry.tsv", sep="\t").genome_id.str.startswith("GENOME_").all()

import io
import json
from pathlib import Path
import tarfile
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

import mock16s_chem as m
from scripts.build_genome_panel import candidate_plan, digest, replay


def test_reference_roles_keep_contaminants_out_of_community():
    tax = pd.DataFrame({"reference_role": ["biological", "biological", "contaminant"]}, index=["a", "b", "c"])
    seqs = {key: "ACGT" for key in tax.index}
    cfg = {"simulation": {"n_asvs": 2}, "artifacts": {"contaminants": {"enabled": True, "n_asvs": 1, "prevalence": 1, "mean_reads": 10}}}
    chosen = m.choose_asvs(cfg, tax, seqs, np.random.default_rng(42))
    assert set(chosen) == {"a", "b"}
    counts = pd.DataFrame({"s": [20, 30]}, index=chosen)
    result, contaminants = m.add_contaminants(counts, chosen, tax, seqs, cfg, np.random.default_rng(42))
    assert contaminants == ["c"]
    assert result.loc["c", "s"] > 0
    with pytest.raises(ValueError, match="available"):
        m.choose_asvs(cfg, tax, seqs, np.random.default_rng(42), n=3)


def test_candidate_plan_is_deterministic_and_prioritizes_roles(tmp_path):
    header = ["assembly_accession", "group", "assembly_level", "version_status", "excluded_from_refseq", "ftp_path", "genome_size", "species_taxid", "taxid", "organism_name", "refseq_category"]
    rows = []
    lineage = []
    for i, genus in enumerate(["Background", "Streptococcus", "Ralstonia"], 100):
        rows.append([f"GCF_{i}.1", "bacteria", "Complete Genome", "latest", "na", f"https://ftp.ncbi.nlm.nih.gov/genomes/{i}", "2000000", str(i), str(i), genus + " example", "reference genome"])
        lineage.append("\t|\t".join([str(i), genus+" example", "", genus, "Family", "Order", "Class", "Phylum", "", "Bacteria"]) + "\t|\n")
    summary = tmp_path / "summary.tsv"
    summary.write_text("#"+"\t".join(header)+"\n"+"\n".join("\t".join(r) for r in rows)+"\n")
    data = "".join(lineage).encode()
    archive = tmp_path / "taxonomy.tar.gz"
    with tarfile.open(archive, "w:gz") as tar:
        info = tarfile.TarInfo("rankedlineage.dmp"); info.size = len(data)
        tar.addfile(info, io.BytesIO(data))
    policy = {"seed":42, "assembly_levels":["Complete Genome"], "max_genome_bases":12000000,
              "domains":["Bacteria"], "max_assemblies_per_species":2,
              "priority_genera":["Streptococcus"], "contaminant_genera":["Ralstonia"]}
    plan = candidate_plan(policy, summary, archive)
    assert plan == candidate_plan(policy, summary, archive)
    assert [r["Genus"] for r in plan] == ["Ralstonia", "Streptococcus", "Background"]
    assert [r["reference_role"] for r in plan] == ["contaminant", "biological", "biological"]


def test_verify_detects_changed_genome_and_taxonomy(tmp_path):
    genome = tmp_path / "g.fa.gz"; genome.write_bytes(b"fixture bytes")
    taxonomy = tmp_path / "taxonomy.tsv"; taxonomy.write_text("taxonomy fixture\n")
    lock = tmp_path / "panel.lock.json"
    lock.write_text(json.dumps({"genomes":[{"fasta":genome.name,"sha256":digest(genome)}], "taxonomy_table_sha256":digest(taxonomy)}))
    args = SimpleNamespace(lock=lock, output=tmp_path, command="verify")
    replay(args)
    genome.write_bytes(b"changed")
    with pytest.raises(ValueError, match="Checksum mismatch"):
        replay(args)
    genome.write_bytes(b"fixture bytes")
    taxonomy.write_text("changed taxonomy\n")
    with pytest.raises(ValueError, match="Taxonomy table"):
        replay(args)

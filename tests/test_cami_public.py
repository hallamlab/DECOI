import io
import json
import tarfile
from pathlib import Path

import pytest

from scripts import download_cami as cami


def archive(path, name, content=b'test', kind=None):
    with tarfile.open(path, 'w:gz') as handle:
        member = tarfile.TarInfo(name)
        if kind is not None:
            member.type = kind
            member.linkname = '/tmp/outside'
        else:
            member.size = len(content)
        handle.addfile(member, io.BytesIO(content) if kind is None else None)


@pytest.mark.parametrize('name,kind', [('../escape', None), ('/absolute', None), ('link', tarfile.SYMTYPE)])
def test_unsafe_archive_rejected(tmp_path, name, kind):
    path = tmp_path/'setup.tar.gz'
    archive(path, name, kind=kind)
    with pytest.raises(ValueError, match='Unsafe'):
        cami.extract(path, tmp_path/'extracted')
    assert not (tmp_path/'extracted').exists()


def test_cached_archive_integrity(tmp_path, monkeypatch):
    path = tmp_path/'block.setup.tar.gz'
    archive(path, 'source_genomes/g.fa')
    monkeypatch.setattr(cami, 'BUNDLES', {'block': path.stat().st_size})
    monkeypatch.setattr(cami, 'fetch_taxonomy', lambda output: {})
    cami.fetch(tmp_path)
    assert (tmp_path/'block/source_genomes/g.fa').read_bytes() == b'test'
    cami.fetch(tmp_path)
    manifest = json.loads((tmp_path/'download_manifest.json').read_text())
    manifest['archives']['block']['sha256'] = 'wrong'
    (tmp_path/'download_manifest.json').write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match='checksum changed'):
        cami.fetch(tmp_path)


def test_public_sample_assignments():
    for block, count in [('airskinurogenital', 29), ('gastrooral', 20)]:
        numbers = [n for b, ns in cami.SITE_SAMPLES.values() if b == block for n in ns]
        assert sorted(numbers) == list(range(count))


def test_public_preparation(tmp_path, monkeypatch):
    import pandas as pd
    from Bio.Seq import Seq
    from scripts import prepare_cami_study as study
    class SerialPool:
        def __init__(self, **kwargs): pass
        def __enter__(self): return self
        def __exit__(self, *args): pass
        def map(self, fn, args): return map(fn, args)
    monkeypatch.setattr(study, 'ProcessPoolExecutor', SerialPool)
    forward = study.FORWARD.replace('Y', 'C').replace('M', 'A')
    reverse = study.REVERSE.replace('N', 'A').replace('V', 'A').replace('W', 'A')
    locus = forward + 'A'*250 + str(Seq(reverse).reverse_complement())
    for block, count in [('airskinurogenital', 29), ('gastrooral', 20)]:
        bundle = tmp_path/'source'/block
        (bundle/'source_genomes').mkdir(parents=True)
        (bundle/'source_genomes/g.fa').write_text('>genome\n'+locus+'\n')
        (bundle/'genome_to_id.tsv').write_text('g\t/original/server/g.fa\n')
        (bundle/'metadata.tsv').write_text('genome_ID\tNCBI_ID\ng\t9\n')
        (bundle/'taxonomic_profile_0.txt').write_text('@Version:0.9.1\n@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\t_CAMI_GENOMEID\n1\tspecies\t1\tBacteria|P|C|O|F|G|S\t100\tg\n')
        for number in range(count):
            (bundle/f'abundance{number}.tsv').write_text('g\t100\n')
    with tarfile.open(tmp_path/'source'/cami.TAXONOMY_FILE, 'w:gz') as handle:
        parents = ['1', '2', '3', '4', '5', '6', '7']
        ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
        nodes = ''.join(f'{i+2}\t|\t{parent}\t|\t{rank}\t|\n'
                        for i, (parent, rank) in enumerate(zip(parents, ranks)))
        names = ''.join(f'{i+2}\t|\t{name}\t|\t\t|\tscientific name\t|\n'
                        for i, name in enumerate(['Bacteria', 'P', 'C', 'O', 'F', 'G', 'S']))
        for name, content in [('nodes.dmp', nodes), ('names.dmp', names), ('merged.dmp', '9\t|\t8\t|\n')]:
            member = tarfile.TarInfo(name)
            member.size = len(content.encode())
            handle.addfile(member, io.BytesIO(content.encode()))
    study.prepare(tmp_path/'source', tmp_path/'prepared', 1)
    ref = tmp_path/'prepared/reference'
    weights = pd.read_csv(ref/'amplicon_abundances.tsv', sep='\t', index_col=0)
    assert weights.shape == (1, 49)
    assert 'Airways_4' in weights and 'Oral_6' in weights and 'Skin_1' in weights
    assert 'Airways_0' not in weights
    assert pd.read_csv(ref/'v4_asv_taxonomy.tsv', sep='\t').Species.tolist() == ['S']
    assert json.loads((ref/'reference_manifest.json').read_text())['input_layout'] == 'public_cami_setup'


@pytest.mark.parametrize('honor_range', [True, False])
def test_download_resumes_or_restarts_when_range_ignored(tmp_path, monkeypatch, honor_range):
    path = tmp_path/'setup.tar.gz'
    partial = path.with_suffix('.gz.part')
    partial.write_bytes(b'abc')
    class Response(io.BytesIO):
        status = 206 if honor_range else 200
        headers = {'Content-Range': 'bytes 3-5/6'} if honor_range else {}
    def open_url(request, **kwargs):
        assert request.get_header('Range') == 'bytes=3-'
        return Response(b'def' if honor_range else b'abcdef')
    monkeypatch.setattr(cami.urllib.request, 'urlopen', open_url)
    cami.download('https://example.org/setup.tar.gz', path, 6)
    assert path.read_bytes() == b'abcdef'
    assert not partial.exists()


def test_incomplete_download_keeps_partial_for_resume(tmp_path, monkeypatch):
    class Response(io.BytesIO):
        status = 200
    monkeypatch.setattr(cami.urllib.request, 'urlopen', lambda *a, **k: Response(b'abc'))
    path = tmp_path/'setup.tar.gz'
    with pytest.raises(ValueError, match='Incomplete'):
        cami.download('https://example.org/setup.tar.gz', path, 6)
    assert not path.exists()
    assert path.with_suffix('.gz.part').read_bytes() == b'abc'


def test_upstream_checksum_mismatch_rejected_before_extraction(tmp_path, monkeypatch):
    path = tmp_path/'block.setup.tar.gz'
    archive(path, 'genome.fa')
    monkeypatch.setattr(cami, 'BUNDLES', {'block': path.stat().st_size})
    monkeypatch.setattr(cami, 'BUNDLE_MD5', {'block': '0'*32})
    with pytest.raises(ValueError, match='published CAMI MD5'):
        cami.fetch(tmp_path)
    assert not (tmp_path/'block').exists()

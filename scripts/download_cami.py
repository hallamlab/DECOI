#!/usr/bin/env python3
"""Download the official CAMI II toy human microbiome setup bundles."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import tarfile
import urllib.request

BASE_URL = 'https://frl.publisso.de/data/frl:6425518/'
BUNDLES = {'airskinurogenital': 2467407842, 'gastrooral': 1906532721}
# Published at BASE_URL + 'md5sums.tsv'.
BUNDLE_MD5 = {'airskinurogenital': 'e3e9620340d1ef53d751d12fa154267d',
              'gastrooral': 'eebd129ce7f8640683eb744a0e946223'}
TAXONOMY_FILE = 'ncbi-taxonomy_20170222.tar.gz'
TAXONOMY_URL = ('https://raw.githubusercontent.com/CAMI-challenge/CAMISIM/'
                '33421fd0e0de4e0f7d315cc54dd1cfc033fcac69/tools/' + TAXONOMY_FILE)
TAXONOMY_SHA256 = '494c5f06ef7e9b798febd153aff674b602e5c13b2d43aab5db85bedd9949f4a8'
# Original simulation sample numbers, from the public dataset README.md.
SITE_SAMPLES = {
    'Airways': ('airskinurogenital', [4,7,8,9,10,11,12,23,26,27]),
    'Skin': ('airskinurogenital', [1,13,14,15,16,17,18,19,20,28]),
    'Urogenital': ('airskinurogenital', [0,2,3,5,6,21,22,24,25]),
    'Gastrointestinal': ('gastrooral', [0,1,2,3,4,5,9,10,11,12]),
    'Oral': ('gastrooral', [6,7,8,13,14,15,16,17,18,19]),
}


def checksum(path, algorithm='sha256'):
    digest = hashlib.new(algorithm)
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def download(url, path, size):
    """Resume partial downloads; publish only the expected complete byte count."""
    if path.exists():
        if path.stat().st_size != size:
            raise ValueError(f'Unexpected archive size: {path}')
        return
    partial = path.with_suffix(path.suffix + '.part')
    offset = partial.stat().st_size if partial.exists() else 0
    if offset > size:
        raise ValueError(f'Oversized partial download: {partial}')
    if offset < size:
        request = urllib.request.Request(url, headers={'Range': f'bytes={offset}-'} if offset else {})
        with urllib.request.urlopen(request, timeout=120) as response:
            if offset and response.status == 206:
                if not response.headers.get('Content-Range', '').startswith(f'bytes {offset}-'):
                    raise ValueError('Server returned an unexpected download range')
                mode = 'ab'
            elif response.status == 200:
                mode = 'wb'
            else:
                raise ValueError(f'Unexpected HTTP status: {response.status}')
            transferred = offset if mode == 'ab' else 0
            reported = transferred
            with partial.open(mode) as handle:
                for chunk in iter(lambda: response.read(1024 * 1024), b''):
                    handle.write(chunk)
                    transferred += len(chunk)
                    if transferred - reported >= 256 * 1024 * 1024:
                        print(f'  {path.name}: {transferred / size:.0%}', flush=True)
                        reported = transferred
    if partial.stat().st_size != size:
        raise ValueError(f'Incomplete download: {partial}; rerun to resume')
    partial.replace(path)


def extract(archive, destination):
    """Extract regular files/directories only, rejecting paths outside the bundle."""
    if destination.exists():
        raise ValueError(f'Refusing to overwrite {destination}')
    temporary = destination.with_name(destination.name + '.extracting')
    if temporary.exists():
        shutil.rmtree(temporary)
    temporary.mkdir()
    try:
        with tarfile.open(archive, 'r:gz') as bundle:
            for member in bundle:
                target = (temporary / member.name).resolve()
                if not target.is_relative_to(temporary.resolve()) or not (member.isfile() or member.isdir()):
                    raise ValueError(f'Unsafe archive entry: {member.name}')
                if member.isdir():
                    target.mkdir(parents=True, exist_ok=True)
                else:
                    target.parent.mkdir(parents=True, exist_ok=True)
                    with bundle.extractfile(member) as source, target.open('wb') as output:
                        shutil.copyfileobj(source, output)
        temporary.rename(destination)
    except BaseException:
        shutil.rmtree(temporary)
        raise


def fetch_taxonomy(output):
    path = output / TAXONOMY_FILE
    download(TAXONOMY_URL, path, 38530577)
    if checksum(path) != TAXONOMY_SHA256:
        raise ValueError(f'Taxonomy checksum mismatch: {path}')
    return dict(url=TAXONOMY_URL, bytes=38530577, sha256=TAXONOMY_SHA256)


def fetch(output):
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    manifest_path = output / 'download_manifest.json'
    previous = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    records = previous.get('archives', {})
    for name, size in BUNDLES.items():
        archive = output / f'{name}.setup.tar.gz'
        url = f'{BASE_URL}{name}/setup.tar.gz'
        print(f'Downloading/checking {name} ({size / 1e9:.2f} GB)', flush=True)
        download(url, archive, size)
        digest = checksum(archive)
        if name in BUNDLE_MD5 and checksum(archive, 'md5') != BUNDLE_MD5[name]:
            raise ValueError(f'Archive does not match the published CAMI MD5: {archive}')
        if name in records and records[name]['sha256'] != digest:
            raise ValueError(f'Cached archive checksum changed: {archive}')
        destination = output / name
        if destination.exists() and name not in records:
            raise ValueError(f'Untracked extraction directory: {destination}')
        if not destination.exists():
            extract(archive, destination)
        records[name] = dict(url=url, bytes=size, sha256=digest, upstream_md5=BUNDLE_MD5.get(name))
        temporary = manifest_path.with_suffix('.tmp')
        temporary.write_text(json.dumps(dict(dataset_doi='10.4126/FRL01-006425518',
            body_site_source=BASE_URL+'README.md', archives=records,
            checksum_source=BASE_URL+'md5sums.tsv',
            checksum_basis='Published CAMI MD5 and locally computed SHA-256; taxonomy pinned by SHA-256'), indent=2)+'\n')
        temporary.replace(manifest_path)
    manifest = json.loads(manifest_path.read_text())
    manifest['taxonomy'] = fetch_taxonomy(output)
    temporary = manifest_path.with_suffix('.tmp')
    temporary.write_text(json.dumps(manifest, indent=2)+'\n')
    temporary.replace(manifest_path)
    print(f'Ready: {output}', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True, help='Download/cache directory (allow at least 25 GB free)')
    fetch(parser.parse_args().output)

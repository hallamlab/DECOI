"""Resolve CAMI genome taxon IDs using the original CAMISIM taxonomy snapshot."""
import tarfile

RANKS = ('superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')


def lineages(archive, taxids):
    taxids = set(taxids)
    with tarfile.open(archive, 'r:gz') as handle:
        members = {m.name.rsplit('/', 1)[-1]: m for m in handle.getmembers() if m.isfile()}
        def rows(name):
            with handle.extractfile(members[name]) as source:
                for line in source:
                    yield [field.strip() for field in line.decode().split('|')]
        parents = {row[0]: (row[1], row[2]) for row in rows('nodes.dmp')}
        merged = {row[0]: row[1] for row in rows('merged.dmp')} if 'merged.dmp' in members else {}
        paths = {}
        needed = set()
        for original in taxids:
            node = original
            seen = set()
            path = {}
            while node != '1':
                if node in seen:
                    raise ValueError(f'Taxonomy cycle for {original}')
                seen.add(node)
                if node in merged:
                    node = merged[node]
                    continue
                if node not in parents:
                    raise ValueError(f'Taxon {original} missing from CAMISIM taxonomy snapshot')
                parent, rank = parents[node]
                if rank in RANKS:
                    path[rank] = node
                    needed.add(node)
                node = parent
            paths[original] = path
        names = {row[0]: row[1] for row in rows('names.dmp')
                 if row[0] in needed and row[3] == 'scientific name'}
    return {taxid: [names[path[rank]] if rank in path else '' for rank in RANKS]
            for taxid, path in paths.items()}

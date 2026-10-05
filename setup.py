"""Bundle only runtime resources, never downloaded data or manuscript material."""
from pathlib import Path
from setuptools import setup

root_files = ['main.nf', 'nextflow.config', 'mock16s_chem.py', 'genome_reference.py', 'wgs.py']
data = [('share/decoi', root_files)]
for directory in ['scripts', 'config', 'study', 'fixtures', 'references/airway_v1']:
    files = [str(p) for p in sorted(Path(directory).iterdir())
             if p.is_file() and p.suffix in {'.py', '.R', '.yaml', '.json', '.tsv', '.fasta', '.md'}]
    data.append(('share/decoi/' + directory, files))
setup(packages=['decoi'], data_files=data)

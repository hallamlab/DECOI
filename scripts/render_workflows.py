#!/usr/bin/env python3
"""Render DECOI diagrams with the shared MP/ASPIRE module-and-node design."""
from html import escape
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
COLORS = {'data': ('#F5F5F5', '#666666'), 'input': ('#DAE8FC', '#6C8EBF'),
          'output': ('#D5E8D4', '#82B366')}


def draw(name, rows):
    width, height = 1200, 355 + 195 * len(rows)
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-labelledby="title desc">',
             '<title id="title">DECOI workflow overview</title>',
             '<desc id="desc">Conceptual modules connect reference preparation, shared biological truth, chemical and assay effects, read simulation and reproducible outputs. Numbered modules describe concepts, not a serial execution schedule. WGS and DADA2 are optional branches.</desc>',
             f'<rect width="{width}" height="{height}" fill="white"/>',
             '<defs><marker id="arrow" viewBox="0 0 10 10" refX="10" refY="5" markerWidth="7" markerHeight="7" orient="auto"><path d="M0 0 L10 5 L0 10z" fill="#111111"/></marker></defs>',
             '<style>text{font-family:"Times New Roman",Times,serif;fill:#111111}.wire{fill:none;stroke:#111111;stroke-width:1.5;marker-end:url(#arrow)}</style>']

    def rect(x, y, w, h, fill, stroke, radius=16, sw=1.5):
        parts.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{radius}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>')

    def text(x, y, lines, size=22, anchor='middle'):
        if isinstance(lines, str):
            lines = [lines]
        spans = ''.join(f'<tspan x="{x}" dy="{0 if i == 0 else size * 1.18:g}">{escape(line)}</tspan>' for i, line in enumerate(lines))
        parts.append(f'<text x="{x}" y="{y}" font-size="{size}" text-anchor="{anchor}">{spans}</text>')

    def node(x, y, kind):
        fill, stroke = COLORS.get(kind, COLORS['data'])
        if kind == 'compute':
            parts.append(f'<path d="M{x} {y-13} L{x+13} {y} L{x} {y+13} L{x-13} {y}Z" fill="{fill}" stroke="{stroke}" stroke-width="3"/>')
        elif kind == 'module':
            rect(x-13, y-13, 26, 26, '#F5F5F5', '#111111', 0, 3)
        else:
            parts.append(f'<circle cx="{x}" cy="{y}" r="12" fill="{fill}" stroke="{stroke}" stroke-width="3"/>')

    def wire(x1, y1, x2, y2):
        parts.append(f'<path class="wire" d="M{x1} {y1} L{x2} {y2}"/>')

    rect(20, 20, 1160, 80, '#CCCCCC', '#666666')
    text(600, 54, 'DECOI • Simulated microbial studies with Nextflow', 29)
    text(600, 84, 'Mamba environments • Local / Slurm execution • Resources, checkpoints and logs', 21)
    rect(20, 120, 690, 130, *COLORS['input'])
    text(40, 150, 'Inputs', 24, 'start')
    text(40, 180, ['Study design + configuration; SILVA or verified genomes',
                        'Optional: frozen airway panel or linked WGS reference bundle',
                        'Configured assays, biological effects and technical artifacts'], 20, 'start')
    rect(730, 120, 450, 130, '#CCCCCC', '#666666')
    for x, label, kind in [(775, 'Module', 'module'), (865, 'Compute', 'compute'),
                           (955, 'Data', 'data'), (1045, 'Input', 'input'), (1135, 'Output', 'output')]:
        text(x, 157, label, 18)
        node(x, 198, kind)

    for i, (title, source, compute, output, details) in enumerate(rows):
        y = 365 + i * 195
        text(125, y-10, title, 23)
        node(270, y, 'module')
        text(270, y+7, str(i+1), 20)
        node(410, y, 'input' if i == 0 else 'data')
        node(730, y, 'compute')
        node(1060, y, 'output')
        text(410, y-65, source, 21)
        text(730, y-65, compute, 23)
        text(1060, y-65, output, 22)
        text(730, y+43, details, 19)
        wire(284.5, y, 396.5, y)
        wire(423.5, y, 714.879, y)
        wire(745.121, y, 1046.5, y)
        if i < len(rows)-1:
            wire(270, y+14.5, 270, y+195-14.5)
    rect(20, height-85, 1160, 65, *COLORS['output'], radius=12)
    text(600, height-57, 'Matched assays, chemical measurements and recorded ground truth', 23)
    text(600, height-32, 'FASTQs • Metadata • Truth tables • HTML report • Reproducibility manifest', 19)
    parts.append('</svg>')
    (ROOT / 'docs/assets' / name).write_text('\n'.join(parts) + '\n')


def main():
    draw('workflow-brief.svg', [
        (['Reference', 'preparation'], ['SILVA or', 'verified genomes'], ['Extract V4 markers;', 'link taxonomy and genomes'], ['Reference sequences', 'and identity registry'], ['Primer matching / Cutadapt • Deduplication • Checksums and marker-copy links']),
        (['Biological', 'community'], ['Study design and', 'reference identities'], ['Simulate communities;', 'apply biological effects'], ['Shared abundances', 'and chemical truth'], ['SparseDOSSA2 • Cohort / network / batch effects', 'Constructed ASV–compound associations, noise and zeros']),
        (['Assay effects', 'and controls'], ['Shared pre-PCR', 'biological state'], ['Add assay-specific effects', 'and control allocations'], ['Amplicon and', 'optional WGS truth'], ['PCR bias • Contaminants • Mitochondrial fixtures • Chimeras', 'WGS: marker-copy correction and genome-length weighting']),
        (['Read simulation', 'and validation'], ['Assay references', 'and read allocations'], ['Simulate paired reads;', 'check requested counts'], ['Paired FASTQs and', 'validation tables'], ['InSilicoSeq • R1 / R2 count checks', 'Optional amplicon validation: Cutadapt + DADA2']),
        (['Results and', 'reproducibility'], ['Assays, metadata', 'and recorded effects'], ['Assemble the dataset', 'and run summary'], ['Report, manifests', 'and truth tables'], ['Source hashes • Configuration • Random seeds • Software versions'])
    ])
    draw('workflow-main.svg', [
        (['Reference', 'identity'], ['SILVA sequences', 'or genome FASTAs'], ['Extract primer-defined V4;', 'deduplicate and link sources'], ['V4 sequences and', 'reference registry'], ['Cutadapt for SILVA / exact primer matching in genomes', 'Taxonomy • Source hashes • Marker loci and copy numbers']),
        (['Study', 'layout'], ['Study design and', 'reference identities'], ['Expand sample metadata;', 'simulate communities'], ['Sample metadata and', 'initial abundances'], ['Cohorts • Participants • Paired samples', 'SparseDOSSA2 community abundance profiles']),
        (['Biological', 'effects'], ['Initial community', 'abundance profiles'], ['Apply group, network', 'and microbial batch effects'], ['Shared pre-PCR', 'biological state'], ['Differential abundance • Microbial modules • Batch structure', 'Record effects before assay-specific transformations']),
        (['Chemical', 'measurements'], ['Shared biological', 'abundances'], ['Construct signed links;', 'add noise and batch effects'], ['Chemical tables and', 'association truth'], ['ASV–compound coefficients • Zeros • Chemical batch effects', 'Correlations by construction; not mechanistic metabolism']),
        (['Amplicon', 'artifacts'], ['Biological state and', 'artifact settings'], ['Apply PCR bias;', 'inject artifacts and controls'], ['Amplicon truth and', 'read allocations'], ['Contaminants • Mitochondrial marker fixtures • Chimeras', 'Extraction blanks • Recorded artifact identities']),
        (['Matched WGS', '(optional)'], ['Pre-PCR state and', 'linked genomes'], ['Adjust marker copies', 'and genome lengths'], ['WGS truth and', 'read allocations'], ['Use contaminant / control allocations', 'Exclude amplicon-only chimeras and mitochondrial marker fixtures']),
        (['Sequence', 'and validate'], ['References and', 'assay allocations'], ['Simulate paired reads;', 'verify R1 / R2 counts'], ['FASTQs and', 'validation results'], ['InSilicoSeq for both assays • Exact depth checks', 'Optional DADA2: trim, filter, denoise, merge and remove chimeras']),
        (['Dataset and', 'provenance'], ['Assays and all', 'ground-truth tables'], ['Compile metadata,', 'manifests and run report'], ['Reproducible study', 'and HTML report'], ['Reference fixtures • Configuration • Seeds • Hashes • Software versions'])
    ])


if __name__ == '__main__':
    main()

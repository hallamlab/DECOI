"""Build user documentation without importing analysis dependencies."""
import os

project = 'DECOI'
author = 'Hallam Lab and DECOI contributors'
copyright = '2026, DECOI contributors'
extensions = ['myst_parser', 'sphinxcontrib.mermaid']
source_suffix = {'.md': 'markdown'}
root_doc = 'index'
# Explicit publication boundary: local manuscript drafts and audits never enter
# the documentation build, even when they exist in a developer checkout.
include_patterns = ['*.md']
myst_heading_anchors = 4
myst_fence_as_directive = ['mermaid']
html_theme = 'sphinx_rtd_theme'
html_theme_options = {'collapse_navigation': False, 'navigation_depth': 2}
html_title = 'DECOI: Data Emulator for Community Omics tool Benchmarking'
html_baseurl = os.environ.get('READTHEDOCS_CANONICAL_URL', '')
html_static_path = ['assets']
html_css_files = ['docs.css']
mermaid_version = '11.12.1'
mermaid_init_config = {'startOnLoad': False, 'theme': 'neutral', 'flowchart': {'htmlLabels': False}}
mermaid_light_theme = 'neutral'
mermaid_dark_theme = 'neutral'
mermaid_fullscreen = True
mermaid_height = 'auto'

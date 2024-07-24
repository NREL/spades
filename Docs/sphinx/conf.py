import sys

extensions = ["sphinx.ext.mathjax", "breathe", "exhale"]
templates_path = ["_templates"]
source_suffix = ".rst"
master_doc = "index"
project = "SPADES"
copyright = "2023, the Alliance for Sustainable Energy, LLC., through National Renewable Energy Laboratory. All rights reserved"
author = "M. T. Henry de Frahan, E. Young, H. Sitaraman, R. Larsen, D. Vaidhynathan"
version = "0.1"
release = "0.1"
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
pygments_style = "sphinx"
todo_include_todos = False
numfig = True
numfig_format = {"figure": "%s", "table": "%s", "code-block": "%s"}
html_theme = "sphinx_rtd_theme"
# html_static_path = ['_static']
htmlhelp_basename = "spadesdoc"
latex_elements = {}
latex_documents = [
    (master_doc, "spades.tex", "SPADES Documentation", author, "manual"),
]
texinfo_documents = [
    (
        master_doc,
        "spades",
        "SPADES Documentation",
        author,
        "SPADES",
        "Scalable PArallel Discrete Event Simulation.",
        "Miscellaneous",
    ),
]

primary_domain = 'cpp'
highlight_language = 'cpp'

nitpick_ignore_regex = [(r"cpp:identifier", r"amrex.*")]

breathe_default_project = "spades"

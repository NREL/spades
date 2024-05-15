import sys

extensions = [ 'sphinx.ext.mathjax']
templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'
project = u'SPADES'
copyright = u'2023, the Alliance for Sustainable Energy, LLC., through National Renewable Energy Laboratory. All rights reserved'
author = u'E. Young, M. T. Henry de Frahan, H. Sitaraman, R. Larsen'
version = u'0.1'
release = u'0.1'
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
pygments_style = 'sphinx'
todo_include_todos = False
numfig = True
numfig_format = {'figure': '%s', 'table': '%s', 'code-block': '%s'}
html_theme = 'sphinx_rtd_theme'
#html_static_path = ['_static']
htmlhelp_basename = 'spadesdoc'
latex_elements = { }
latex_documents = [
    (master_doc, 'spades.tex', u'SPADES Documentation',
     author, 'manual'),
]
texinfo_documents = [
    (master_doc, 'spades', u'SPADES Documentation',
     author, 'SPADES', 'Scalable PArallel Discrete Events Simulation.',
     'Miscellaneous'),
]

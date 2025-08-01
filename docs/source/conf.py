# -*- coding: utf-8 -*-
#
# NumBAT documentation build configuration file, created by
# sphinx-quickstart on Wed Oct 19 15:23:46 2016.
#
# This file is execfile()d with the current directory set to its
# containing dir.
#
# Note that not all possible configuration values are present in this
# autogenerated file.
#
# All configuration values have a default; values that are commented out
# serve to show the default.

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#

import sys
import re

#MOCK_MODULES = ['scipy', 'scipy.interpolate', 'numpy',
#    'matplotlib', 'matplotlib.pyplot', 'matplotlib.mlab', 'matplotlib.gridspec',
#    'fortran','make_axes_locatable','csv','mpl_toolkits', 'json', 'matplotlib.externals',
#    'matplotlib.cbook', 'matplotlib.axes', 'matplotlib.axes.prop_cycle',
#    'matplotlib.transforms', 'matplotlib.externals.six', 'matplotlib.artist',
#    'matplotlib.axis', 'mpl_toolkits.axes_grid1']


#from unittest.mock import MagicMock

#class Mock(MagicMock):
#    @classmethod
#    def __getattr__(cls, name):
#            return MagicMock()
#
#sys.modules.update((mod_name, Mock()) for mod_name in MOCK_MODULES)




# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#sys.path.insert(0, os.path.abspath('../../backend/'))
#

#sys.path.insert(0, str(Path.home() / 'research/numbat/nb_releases/nb_latest/backend'))

sys.path.insert(0, '../../backend')   # This file is run from docs/source

#sys.path.insert(0, '/home/msteel/research/numbat/latest/backend')
#sys.path.insert(0, r'%HOMEPATH%\numbat\nb_releases\nb_latest\backend')

#sys.path.insert(0, str(Path.home() / 'research/numbat/devel_1/backend'))
#print('path is', sys.path)
#import os
#print(os.getcwd())

from nbversion import *

# -- General configuration ------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.napoleon',   # handle numpy and google docstrings
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'nbsphinx',             # parser for jupyter notebooks
    'sphinxcontrib.bibtex',
    'sphinx_subfigure'

#    'sphinx.ext.coverage',
#    'sphinx.ext.viewcode',
#    'sphinx.ext.githubpages',
#
]

nbsphinx_execute = 'never'

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
#source_suffix = '.rst'

# The encoding of source files.
#
# source_encoding = 'utf-8-sig'

# The master toctree document.
master_doc = 'index'

# General information about the project.
project = u'NumBAT'
copyright = u'2016-2025, Michael Steel, Bjorn Sturmberg, Blair Morrison, Mike Smith, and Christopher Poulton'
author = u'Michael Steel, Bjorn Sturmberg, Blair Morrison, Mike Smith, and Christopher Poulton'

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
#version = u'2.0'
version = NUMBAT_VERSION_STR_MM

# The full version, including alpha/beta/rc tags.
#release = u'2.1.0'
release = NUMBAT_VERSION_STR_MMM

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = 'en'

# There are two options for replacing |today|: either, you set today to some
# non-false value, then it is used:
#
# today = ''
#
# Else, today_fmt is used as the format for a strftime call.
#
today_fmt = '%d %B %Y'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This patterns also effect to html_static_path and html_extra_path
exclude_patterns = []

# The reST default role (used for this markup: `text`) to use for all
# documents.
#
# default_role = None

# If true, '()' will be appended to :func: etc. cross-reference text.
#
# add_function_parentheses = True

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
#
# add_module_names = True

# If true, sectionauthor and moduleauthor directives will be shown in the
# output. They are ignored by default.
#
# show_authors = False

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# A list of ignored prefixes for module index sorting.
# modindex_common_prefix = []

# If true, keep warnings as "system message" paragraphs in the built documents.
# keep_warnings = False

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = False


# -- Options for HTML output ----------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_themes = ['alabaster', 'classic', 'sphinxdox', 'scrolls',
               'agogo', 'nature', 'pyramid', 'haiku', 'traditional',
               'epub', 'bizstyle', 'sphinx_rtd_theme']
html_theme = 'default' #'alabaster'
html_theme = 'pyramid'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
# html_theme_options = {
#         'logo_only': True,
#         'sticky_navigation': False,
#         'display_version': True,
#         }

# Add any paths that contain custom themes here, relative to this directory.
# html_theme_path = []

# The name for this set of Sphinx documents.
# "<project> v<release> documentation" by default.
#
html_title = u'NumBAT - The Numerical Brillouin Analysis Tool'

# A shorter title for the navigation bar.  Default is the same as html_title.
#
#html_short_title = u'NumBAT 2.1'
html_short_title = u'NumBAT ' + NUMBAT_VERSION_STR_MM

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
#
html_logo = 'NumBAT_logo.png'

# The name of an image file (relative to this directory) to use as a favicon of
# the docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
#
html_favicon = 'numbat.ico'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
#html_static_path = ['_static']

# Add any extra paths that contain custom files (such as robots.txt or
# .htaccess) here, relative to this directory. These files are copied
# directly to the root of the documentation.
#
# html_extra_path = []

# If not None, a 'Last updated on:' timestamp is inserted at every page
# bottom, using the given strftime format.
# The empty string is equivalent to '%b %d, %Y'.
#
# html_last_updated_fmt = None

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
#
# html_use_smartypants = True

# Custom sidebar templates, maps document names to template names.
#
# html_sidebars = {}

# Additional templates that should be rendered to pages, maps page names to
# template names.
#
# html_additional_pages = {}

# If false, no module index is generated.
#
# html_domain_indices = True

# If false, no index is generated.
#
# html_use_index = True

# If true, the index is split into individual pages for each letter.
#
# html_split_index = False

# If true, links to the reST sources are added to the pages.
#
# html_show_sourcelink = True

# If true, "Created using Sphinx" is shown in the HTML footer. Default is True.
#
# html_show_sphinx = True

# If true, "(C) Copyright ..." is shown in the HTML footer. Default is True.
#
# html_show_copyright = True

# If true, an OpenSearch description file will be output, and all pages will
# contain a <link> tag referring to it.  The value of this option must be the
# base URL from which the finished HTML is served.
#
# html_use_opensearch = ''

# This is the file name suffix for HTML files (e.g. ".xhtml").
# html_file_suffix = None

# Language to be used for generating the HTML full-text search index.
# Sphinx supports the following languages:
#   'da', 'de', 'en', 'es', 'fi', 'fr', 'hu', 'it', 'ja'
#   'nl', 'no', 'pt', 'ro', 'ru', 'sv', 'tr', 'zh'
#
# html_search_language = 'en'

# A dictionary with options for the search language support, empty by default.
# 'ja' uses this config value.
# 'zh' user can custom change `jieba` dictionary path.
#
# html_search_options = {'type': 'default'}

# The name of a javascript file (relative to the configuration directory) that
# implements a search results scorer. If empty, the default will be used.
#
# html_search_scorer = 'scorer.js'

# Output file base name for HTML help builder.
htmlhelp_basename = 'NumBATdoc'

# -- Options for LaTeX output ---------------------------------------------

latex_maketitle = r'''
\renewcommand{\maketitle}{%
  \noindent\rule{\textwidth}{1pt}\par
    \begingroup % for PDF information dictionary
       \def\endgraf{ }\def\and{\& }%
       \pdfstringdefDisableCommands{\def\\{, }}% overwrite hyperref setup
       \hypersetup{pdfauthor={\@author}, pdftitle={\@title}}%
    \endgroup
  \begin{flushright}
    \sphinxlogo
    \py@HeaderFamily
    {\Huge \@title }\par
    {\itshape\large \py@release \releaseinfo}\par
    \vspace{25pt}
    {\Large
      \begin{tabular}[t]{c}
        \@author
      \end{tabular}}\par
    \vspace{25pt}
    \@date \par
    \py@authoraddress \par
  \end{flushright}
  \@thanks
  \setcounter{footnote}{0}
  \let\thanks\relax\let\maketitle\relax
  %\gdef\@thanks{}\gdef\@author{}\gdef\@title{}
  }
'''

latex_elements = {
        # The paper size ('letterpaper' or 'a4paper').
        #
        'papersize': 'a4paper',

        # The font size ('10pt', '11pt' or '12pt').
        #
        # 'pointsize': '10pt',

        # Additional stuff for the LaTeX preamble.
        #

        # Latex figure (float) alignment
        #
        #'figure_align': 'htbp',
        'figure_align': 'H',

        #     'maketitle': latex_maketitle,

#        'preamble': r'\usepackage{amsmath}' + '\n' +r'\usepackage{amssymb}'
        }
latex_elements['preamble'] = r'\usepackage{amsmath} '
latex_elements['preamble'] += r'\usepackage{amssymb} '
#latex_elements['preamble'] += r'\usepackage{boldsymbol} '

#latex_elements['preamble'] = r'\usepackage{amsmath}' + '\n' +r'\usepackage{amssymb}'

try:
    pngmath_latex_preamble
except NameError:
    pngmath_latex_preamble = ''



# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
        #(master_doc, 'NumBAT.tex', u'NumBAT Documentation',
        (master_doc, 'NumBAT.tex', u'NumBAT - The Numerical Brillouin Analysis Tool',
         u'Michael Steel, Bjorn Sturmberg, Blair Morrison, \\\\ Mike Smith and  Christopher Poulton' , 'manual'),
        ]

# The name of an image file (relative to this directory) to place at the top of
# the title page.
#
latex_logo = "NumBAT_logo.png"

# For "manual" documents, if this is true, then toplevel headings are parts,
# not chapters.
#
# latex_use_parts = False

# If true, show page references after internal links.
#
# latex_show_pagerefs = False

# If true, show URL addresses after external links.
#
# latex_show_urls = False

# Documents to append as an appendix to all manuals.
#
# latex_appendices = []

# It false, will not define \strong, \code, 	itleref, \crossref ... but only
# \sphinxstrong, ..., \sphinxtitleref, ... To help avoid clash with user added
# packages.
#
# latex_keep_old_macro_names = True

# If false, no module index is generated.
#
# latex_domain_indices = True


# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
        (master_doc, 'numbat', u'NumBAT Documentation',
         [author], 1)
        ]

# If true, show URL addresses after external links.
#
# man_show_urls = False


# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
        (master_doc, 'NumBAT', u'NumBAT Documentation',
         author, 'NumBAT', 'One line description of project.',
         'Miscellaneous'),
        ]

# Documents to append as an appendix to all manuals.
#
# texinfo_appendices = []

# If false, no module index is generated.
#
# texinfo_domain_indices = True

# How to display URL addresses: 'footnote', 'no', or 'inline'.
#
# texinfo_show_urls = 'footnote'

# If true, do not generate a @detailmenu in the "Top" node's menu.
#
# texinfo_no_detailmenu = False


# -- Options for Epub output ----------------------------------------------

# Bibliographic Dublin Core info.
epub_title = project
epub_author = author
epub_publisher = author
epub_copyright = copyright

# The basename for the epub file. It defaults to the project name.
# epub_basename = project

# The HTML theme for the epub output. Since the default themes are not
# optimized for small screen space, using the same theme for HTML and epub
# output is usually not wise. This defaults to 'epub', a theme designed to save
# visual space.
#
# epub_theme = 'epub'

# The language of the text. It defaults to the language option
# or 'en' if the language is not set.
#
# epub_language = ''

# The scheme of the identifier. Typical schemes are ISBN or URL.
# epub_scheme = ''

# The unique identifier of the text. This can be a ISBN number
# or the project homepage.
#
# epub_identifier = ''

# A unique identification for the text.
#
# epub_uid = ''

# A tuple containing the cover image and cover page html template filenames.
#
# epub_cover = ()

# A sequence of (type, uri, title) tuples for the guide element of content.opf.
#
# epub_guide = ()

# HTML files that should be inserted before the pages created by sphinx.
# The format is a list of tuples containing the path and title.
#
# epub_pre_files = []

# HTML files that should be inserted after the pages created by sphinx.
# The format is a list of tuples containing the path and title.
#
# epub_post_files = []

# A list of files that should not be packed into the epub file.
epub_exclude_files = ['search.html']

# The depth of the table of contents in toc.ncx.
#
# epub_tocdepth = 3

# Allow duplicate toc entries.
#
# epub_tocdup = True

# Choose between 'default' and 'includehidden'.
#
# epub_tocscope = 'default'

# Fix unsupported image types using the Pillow.
#
# epub_fix_images = False

# Scale large images.
#
# epub_max_image_width = 0

# How to display URL addresses: 'footnote', 'no', or 'inline'.
#
# epub_show_urls = 'inline'

# If false, no index is generated.
#
# epub_use_index = True


# latex_engine = "xelatex"


# Our NumBAT latex macros

with open('latex_macros.sty') as f:
    for macro in f:
        latex_elements['preamble'] += '\n' + macro
        pngmath_latex_preamble += '\n' + macro


# Latex macros in mathjax:
#    See  https://stackoverflow.com/questions/9728292/creating-latex-math-macros-within-sphinx
mathjax3_config = { 'tex': {'macros': {}}}

with open('latex_macros.sty', 'r') as f:
    for line in f:
        macros = re.findall(r'\\(DeclareRobustCommand|newcommand|renewcommand){\\(.*?)}(\[(\d)\])?{(.+)}', line)
        for macro in macros:
            if len(macro[2]) == 0:
                mathjax3_config['tex']['macros'][macro[1]] = "{"+macro[4]+"}"
            else:
                mathjax3_config['tex']['macros'][macro[1]] = ["{"+macro[4]+"}", int(macro[3])]


bibtex_bibfiles = ['numbat_refs.bib']
bibtex_default_style = 'unsrt'
bibtex_reference_style = 'label'


import sphinx.builders.latex.transforms
class DummyTransform(sphinx.builders.latex.transforms.BibliographyTransform):
    def run(self, **kwargs):
        pass
sphinx.builders.latex.transforms.BibliographyTransform = DummyTransform

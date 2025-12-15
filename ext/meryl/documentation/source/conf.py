# -*- coding: utf-8 -*-

import sys
import os

def setup(app):
  app.add_css_file('custom.css')

extensions = [
    'sphinx.ext.todo',
    'sphinx.ext.mathjax',
    'sphinx.ext.ifconfig',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix of source filenames.
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# General information about the project.
project   = u'meryl'
copyright = u'2018-2023, Brian Walenz, Arang Rhie'

# Version info, acts as replacement for |version| and |release|, also used in
# various other places throughout the built documents.
#
version = '1.3'    # The short X.Y version.
release = '1.3'    # The full version, including alpha/beta/rc tags.

# There are two options for replacing |today|: either, you set today to
# some non-false value, then it is used, otherwise, today_fmt is used as
# the format for a strftime call.
#
#today = ''
#today_fmt = '%B %d, %Y'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
#
exclude_patterns = []

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# If true, keep warnings as "system message" paragraphs in the built documents.
keep_warnings = True

# -- Options for HTML output ----------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = 'default'
html_theme_options = {}

#  Build using the RTD theme, if not on RTD.
#    https://read-the-docs.readthedocs.org/en/latest/theme.html
#    https://github.com/snide/sphinx_rtd_theme
#
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

if not on_rtd:
    import sphinx_rtd_theme
    html_theme = 'sphinx_rtd_theme'

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
#html_logo = None

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
#html_favicon = None

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Add any extra paths that contain custom files (such as robots.txt or
# .htaccess) here, relative to this directory. These files are copied
# directly to the root of the documentation.
#html_extra_path = []

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
html_last_updated_fmt = '%b %d, %Y'

# If true, "Created using Sphinx" is shown in the HTML footer. Default is True.
#html_show_sphinx = True

# If true, "(C) Copyright ..." is shown in the HTML footer. Default is True.
#html_show_copyright = True

# Output file base name for HTML help builder.
htmlhelp_basename = 'meryldoc'

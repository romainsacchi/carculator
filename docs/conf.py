# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

# sys.path.insert(0, os.path.abspath('C:\Users\sacchi_r\Documents\GitHub\coarse\coarse'))
sys.path.insert(0, os.path.abspath(".."))


# -- Project information -----------------------------------------------------

project = "Carculator"
copyright = "2019, Paul Scherrer Institut"
author = "Chris Mutel, Brian Cox, Romain Sacchi"

# The full version, including alpha/beta/rc tags
release = "1.6.8"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinxcontrib.bibtex",
    "sphinx_copybutton",
    "sphinx_design",
    "sphinx_immaterial",
]

autoapi_type = "python"
autoapi_dirs = ["../carculator"]

master_doc = "index"

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_logo = "_static/img/mediumsmall.png"
html_favicon = "_static/img/favicon.png"
html_title = "Carculator truck"

html_theme = "sphinx_immaterial"

# material theme options (see theme.conf for more information)
html_theme_options = {
    "icon": {
        "repo": "fontawesome/brands/github",
    },
    "font": {"text": "Fira Sans", "code": "JetBrains Mono"},
    "site_url": "https://carculator.readthedocs.io",
    "repo_url": "https://github.com/romainsacchi/carculator",
    "repo_name": "romainsacchi/carculator_truck",
    "repo_type": "github",
    "edit_uri": "blob/master/docs/",
    "globaltoc_collapse": True,
    "features": ["navigation.top", "search.share", "navigation.tracking", "toc.follow"],
    "palette": [
        {
            "media": "(prefers-color-scheme: light)",
            "scheme": "default",
            "primary": "teal",
            "accent": "amber",
            "toggle": {
                "icon": "material/weather-night",
                "name": "Switch to dark mode",
            },
        },
        {
            "media": "(prefers-color-scheme: dark)",
            "scheme": "slate",
            "primary": "teal",
            "accent": "amber",
            "toggle": {
                "icon": "material/weather-sunny",
                "name": "Switch to light mode",
            },
        },
    ],
    "toc_title_is_page_title": True,
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
html_css_files = ["css/custom.css"]
# html_sidebars = { '**': ['globaltoc.html', 'relations.html', 'sourcelink.html', 'searchbox.html'] }

references = []
bibtex_path = os.path.join(".", "references")

if os.path.exists(bibtex_path):
    for bib_file in os.listdir(bibtex_path):
        if bib_file.endswith(".bib"):
            filename = os.path.basename(bib_file)
            references.append("./references/{}".format(filename))

bibtex_bibfiles = references
bibtex_default_style = "unsrt"
bibtex_reference_style = "author_year"

rst_epilog = """
.. |br| raw:: html 

   <br>

.. |s_caption| raw:: html

    <p style="text-align:center;">

.. |e_caption| raw:: html

    </p>
"""
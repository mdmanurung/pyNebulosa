"""Sphinx configuration for pyNebulosa documentation."""

project = "pyNebulosa"
copyright = "2026, Mikhael Manurung"
author = "Mikhael Manurung"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "myst_nb",
]

# General settings
templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
nitpicky = True
nitpick_ignore = [
    ("py:class", "matplotlib.axes.Axes"),
    ("py:class", "matplotlib.axes._axes.Axes"),
    ("py:class", "matplotlib.figure.Figure"),
    ("py:class", "anndata.AnnData"),
    ("py:class", "anndata._core.anndata.AnnData"),
    ("py:class", "AnnData"),
    ("py:class", "scipy.stats._kde.gaussian_kde"),
    ("py:class", "numpy.ndarray"),
    ("py:class", "np.ndarray"),
]

# Napoleon settings (NumPy-style docstrings)
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_rtype = False

# Autodoc settings
autodoc_member_order = "bysource"
autodoc_typehints = "description"
# Intersphinx mapping
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
    "anndata": ("https://anndata.readthedocs.io/en/stable/", None),
    "scanpy": ("https://scanpy.readthedocs.io/en/stable/", None),
}

# MyST-NB settings (for notebook rendering)
nb_execution_mode = "off"  # Don't re-execute notebooks during build

# HTML output
html_theme = "furo"
html_title = "pyNebulosa"
html_theme_options = {
    "source_repository": "https://github.com/mdmanurung/pyNebulosa",
    "source_branch": "main",
    "source_directory": "docs/",
}
html_static_path = ["_static"]
html_css_files = ["custom.css"]

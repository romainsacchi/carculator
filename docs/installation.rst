.. _install:

Installation
============

Use Python **3.12** (``>=3.12,<3.13``) in a fresh environment.
The shared runtime requires NumPy ``>=1.26.4,<2``.

Install a release
-----------------

After publication, install ``carculator 1.9.7`` from PyPI::

   python3.12 -m venv .venv
   source .venv/bin/activate
   python -m pip install "carculator==1.9.7"

On Windows activate with ``.venv\Scripts\activate``. Alternatively, create a
conda environment with ``conda create -n carculator-release python=3.12 pip``,
activate it, and use the same pip command. Availability of a conda package is
separate from the PyPI release.

Core calculations use bundled resources without a Brightway project, an
ecoinvent installation or network access. The matching ``carculator_utils``
runtime installs ``brightpath>=1.0.0a6,<1.1`` for Brightway Excel, SimaPro CSV
and openLCA JSON-LD export. This uses Brightpath's v1 API, currently an alpha.
Brightpath brings ``bw2io``, XlsxWriter and ``olca-schema``; model and LCIA code
does not import Brightpath or Brightway.

The ``excel`` extra remains a compatibility alias; the Excel writer is already
installed through Brightpath. To select the tested legacy Brightway stack::

   python -m pip install "carculator[brightway]==1.9.7"

The ``brightway`` extra selects ``bw2io<0.9``, ``bw2data<4`` and ``bw2calc<2``.
Brightpath also installs ``bw2io`` without this extra. Exports default to ecoinvent
3.12 cutoff. Older 3.9/3.10 exports require verified supplier mappings and may
stop if those mappings are unavailable. Match external suppliers to the
corresponding background in the destination tool. openLCA exports contain only
foreground processes and require provider and elementary-flow mapping before
calculation. See :doc:`inventory_export` for examples and limitations.

Source checkout and documentation
---------------------------------

This website follows the repository documentation. Changes listed under
``Unreleased`` in the changelog may be newer than the package on PyPI. Use
matching source checkouts to reproduce those changes; record their Git revisions
in addition to package version numbers.

For development, use the matching sibling checkouts and install from this
repository root::

   python -m pip install -e "../carculator_utils[test,excel,brightway]" -e ".[test,docs,excel,brightway]"
   python -m pip check
   python -m pytest
   python -m pip install -r docs/docs_requirements.txt
   python -m sphinx -b html docs docs/_build/html

Use ``docs/docs_requirements.txt`` to install the extensions needed to build
this documentation site.

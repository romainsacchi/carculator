Release notes
=============

1.9.7
-----

When upgrading from an earlier version:

* Create a Python 3.12 environment and install the updated package as described
  in :doc:`installation`. This version requires ``carculator_utils>=1.3.7``.
* Recalculate saved scenarios: changes to energy accounting, battery sizing,
  costs and emissions can affect results.
* Select a supported background scenario: ``SSP2-NPi``,
  ``SSP2-PkBudg1000``, ``SSP2-PkBudg650`` or ``static``. Older ``1150`` and
  ``500`` pathway labels are no longer accepted.

Inventory export now uses Brightpath through ``carculator_utils``. Writers are
runtime dependencies; the ``brightway`` extra selects the tested legacy stack.
Retain one sample before constructing the model and inventory. Brightway
importers remain unlinked; SimaPro consumers must account for its new Latin-1
layout and identifiers. The new openLCA JSON-LD ZIP contains foreground processes
and requires provider and elementary-flow mapping before calculation. See
:doc:`inventory_export` for migration details and return values.

The :download:`changelog <../CHANGELOG.md>` lists the changes in each version.
See :doc:`validity` for calibration evidence and the scope of model validation.

NMC532 chemistry migration
--------------------------

The ecoinvent 3.12 background refresh replaces the former ``NMC-523`` option
with ``NMC-532`` (Ni:Mn:Co = 5:3:2). Use ``NMC-532`` in explicit
``energy_storage`` selections and chemistry-specific custom parameter names.
This changes the inventory chemistry; it is not an alias for the old recipe.
Capacity, cell-mass-share, cycle-life and cost values retain the existing
engineering assumptions under the new name, rather than a new empirical calibration.
Historical validation snapshots retain their original labels. Use matching
updated vehicle and ``carculator_utils`` checkouts for this background refresh.

LCA background update
---------------------

Matching shared utilities now use premise 2.5.4 and ecoinvent 3.12 cutoff, with
a rebuilt A matrix and all 57 background LCIA coefficient matrices. Exports
default to ecoinvent 3.12. Older targets reject suppliers without verified
backward links. Recalculate saved LCA results; the vehicle energy/mass model
is unchanged by this background refresh. Reproduction and validation details
are in carculator_utils' ``docs/background_rebuild.rst``.

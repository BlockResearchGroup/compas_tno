.. _analysis:

********************************************************************************
Analysis
********************************************************************************

.. currentmodule:: compas_tno.optimisers

.. highlight:: python

This tutorial provides a quick tour of the generation of :mod:`Analysis <compas_tno.analysis.Analysis>`.

The :mod:`Analysis <compas_tno.analysis.Analysis>` object stores the main objects of the analysis which are:

* ``Analysis.form``: the :mod:`FormDiagram <compas_tna.diagrams.FormDiagram>` of the analysis,
* ``Analysis.envelope``: the :mod:`Envelope <compas_tna.envelope.Envelope>` of the analysis,
* ``Analysis.optimiser``: the :mod:`Optimiser <compas_tno.optimisers.Optimiser>` of the analysis,

Creating an Analysis object
=============================

The :mod:`Analysis <compas_tno.analysis.Analysis>` object can be created based on one or all of the three elements that it composes.

* :mod:`from_elements <compas_tno.analysis.Analysis.from_elements>`
* :mod:`from_form_and_optimiser <compas_tno.analysis.Analysis.from_form_and_optimiser>`
* :mod:`from_form_and_shape <compas_tno.analysis.Analysis.from_form_and_shape>`

Beyond defining each element separately, a series of functions has been implemented to create the optimiser for a specific problem. Nevertheless, when these are created, the user should check if the `Analysis.optimiser.settings` dictionaty contains the right information for the analysis.

* :mod:`create_lp_analysis <compas_tno.analysis.Analysis.create_lp_analysis>`
* :mod:`create_max_load_analysis <compas_tno.analysis.Analysis.create_max_load_analysis>`
* :mod:`create_maxhrust_analysis <compas_tno.analysis.Analysis.create_maxhrust_analysis>`
* :mod:`create_minthrust_analysis <compas_tno.analysis.Analysis.create_minthrust_analysis>`
* :mod:`create_minthk_analysis <compas_tno.analysis.Analysis.create_minthk_analysis>`
* :mod:`create_compl_energy_analysis <compas_tno.analysis.Analysis.create_compl_energy_analysis>`

Main Analysis Methods
=============================

When the elements are added to the Analysis object, the following methods are useful to assign specific constraints on the problem

* :mod:`apply_selfweight <compas_tno.analysis.Analysis.apply_selfweight>`
* :mod:`apply_envelope <compas_tno.analysis.Analysis.apply_envelope>`
* :mod:`apply_reaction_bounds <compas_tno.analysis.Analysis.apply_reaction_bounds>`
* :mod:`apply_hor_multiplier <compas_tno.analysis.Analysis.apply_hor_multiplier>`

After assigning all modifications necessary, the problem needs to be set up and run with the following commands

* :mod:`set_up_optimiser <compas_tno.analysis.Analysis.set_up_optimiser>`
* :mod:`run <compas_tno.analysis.Analysis.run>`

A full list of relevant methods should be checked at the :mod:`Analysis <compas_tno.analysis.Analysis>` documentation.

Create a minimum thrust analysis
=================================

The code below creates a minimum thrust analysis in a shallow crossvault and using the cross form diagram. This example is discussed detail :ref:`here <example-cross-1>`.

.. code-block:: Python

    from compas_tna.envelope import CrossVaultEnvelope
    from compas_tna.diagrams import FormDiagram
    from compas_masonry.viewers import MasonryViewer
    from compas_tno.analysis import Analysis

    envelope = CrossVaultEnvelope(x_span=(0.0, 10.0),
                                 y_span=(0.0, 10.0),
                                 thickness=0.5,
                                 n=100)
    form = FormDiagram.create_cross_form(discretisation=10)

    analysis = Analysis.create_minthrust_analysis(form, envelope)
    analysis.apply_selfweight()
    analysis.apply_envelope()
    analysis.set_up_optimiser()
    analysis.run()

    view = MasonryViewer(formdiagram=form, envelope=envelope)
    view.setup()
    view.show()


Disclaimer about the solution
=================================

The problem of finding a connected network within the bounds of a masory geometry is nonlinear. Therefore, a solution is not guarantee. For this reason the output of the Analysis should always be checked to see if the exitflag is satisfactory, i.e., if the problem was indeed solved.

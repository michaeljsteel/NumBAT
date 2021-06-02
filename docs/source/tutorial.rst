
.. _chap-tutorial-label:

This chapter provides several resources for learning NumBAT, exploring its applications and validating it against literature results. You should begin by working through the sequence of tutorial exercises which are largely based on literature results. You may then select from examples drawn from a recent tutorial paper by Dr Mike Smith and colleagues, and a range of other literature studies.

Tutorial
--------

In this section we walk through a number of simple simulations that demonstrate the basic use of NumBAT.
:ref:`sec-literature-label` looks at a number of literature examples taken from many of
the well-known groups in this field.
The full Python interface is documented in :ref:`chap-pythonbackend-label`.



Tutorial 1 -- Basic SBS Gain Calculation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This example, contained in ``tutorials/simo-tut_01-first_calc.py`` calculates the backward SBS gain for a rectangular silicon waveguide surrounded by air.

The sequence of operations (annotated in the source code below as Step 1, Step 2, etc) is:

  #. Import NumBAT modules
  #. Define the structure shape and dimensions
  #. Specify the electromagnetic and acoustic modes to be solved for
  #. Construct the waveguide with ``objects.Struct``
  #. Solve the electromagnetic problem. ``mode_calcs.calc_EM_modes()`` returns an object containing electromagnetic mode profiles, propagation constants, and potentially other data which can be accessed through various methods.
  #. Display the propagation constants in units of :math:`\text{m}^{-1}` of the EM modes using ``mode_calcs.kz_EM_all()``
  #. Obtain the effective index of the fundamental mode using ``mode_calcs.neff()``
  #. Identify the desired acoustic wavenumber from the difference of the pump and Stokes propagation constants and solve the acoustic problem.  ``mode_calcs.calc_AC_modes()`` returns an object containing the acoustic mode profiles, frequencies and potentially other data at the propagation constant ``k_AC``.
  #. Display the acoustic frequencies in Hz using ``mode_calcs.nu_AC_all()``.
  #. Calculate the total SBS gain, contributions from photoelasticity and moving boundary effects, and the acoustic loss using ``integration.gain_and_qs()``.


.. literalinclude:: ../../tutorials/simo-tut_01-first_calc.py
    :lines: 0-

.. raw:: latex

    \clearpage


Tutorial 2 -- SBS Gain Spectra
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This example, contained in ``tutorials/simo-tut_02-gain_spectra-npsave.py`` considers the same silicon-in-air structure but adds plotting of fields, gain spectra and techniques for saving and reusing data from earlier calculations. 

Here are some elements to note\:

  #. ``np.savez()`` and ``np.load()`` allow storage of arbitrary data in numpy ``.npz`` files between simulations to accelerate subsequent calculations. Use the flag ``recalc_fields`` to determine whether to recalculate the data from scratch or load data from an existing ``.npz`` file. The data is recovered as a :Simmo: object by calling the ``tolist()`` method. Note that numpy requires the ``allow_pickle=True`` flag for loading array data from file.
  #. This tutorial and many subsequent ones can be made to run faster at the expense of accuracy by appending the argument ``fast=1`` to the command line. This has the effect of specifying a coarser FEM grid. In this case, the output data and fields directory begin with ``ftut_02`` rather than ``tut_02``.
  #. Both electric and magnetic fields can be selected using ``EM_E`` or ``EM_H`` as the value of ``EM_AC`` in ``plotting.mode_fields``. These fields are stored in a folder ``tut_02-fields/`` within the tutorial folder. 
  #. By default, plots are exported as ``png`` format. Pass the option ``pdf_png=pdf`` to plot functions to generate a ``pdf`` output.
  #. Plots of both spectra and modes are generated with a best attempt at font sizes, line widths etc, but the range of potential cases make it impossible to find a selection that works in all cases. Most plot functions therefore support the passing of a ``plotting.Decorator`` object that can vary the settings of these parameters and also pass additional commands to write on the plot axes. See the plotting API for details. This should be regarded as a relatively advanced NumBAT feature.
  #. The ``suppress_imimre`` option suppresses plotting of the :math:`\text{Im}[x]`, :math:`\text{Im}[y]` and :math:`\text{Re}[z]` components of the fields which in a lossless non-leaky problem should normally be zero at all points and therefore not useful to plot.
  #. Vector field plots often require tweaking to get an attractive set of vector arrows.  The ``quiver_points`` option controls the number of arrows drawn along each direction.
  #. The plot functions and the ``Decorator`` class support many options. Consult the API chapter for details on how to fine tune your plots.

.. literalinclude:: ../../tutorials/simo-tut_02-gain_spectra-npsave.py
    :lines: 0-


The following figures show a selection of electromagnetic and acoustic mode profiles produced
in this example.

.. figure:: ../../tutorials/tut_02-fields/EM_E_field_0.png
   :scale: 60 %
   
   Fundamental optical mode fields.


.. figure:: ../../tutorials/tut_02-fields/AC_field_2.png
   :scale: 60 %
   
   Acoustic mode with high gain due to moving boundary effect.


.. figure:: ../../tutorials/tut_02-fields/AC_field_4.png
   :scale: 60 %
   
   Acoustic mode with high gain due to moving boundary effect.

.. raw:: latex

    \clearpage

This example also generates gain spectra.

.. _fig-gainspec1-label:

.. figure:: ../../tutorials/tut_02-gain_spectra-MB_PE_comps.png
   :scale: 35 %
   
   Gain spectra showing gain due to the photoelastic effect, gain due to moving boundary effect, and the total gain.


.. figure:: ../../tutorials/tut_02-gain_spectra-MB_PE_comps_zoom.png
   :scale: 35 %
   
   Zoomed-in gain spectra from :ref:`fig-gainspec1-label`.

.. raw:: latex

    \clearpage


Tutorial 3a -- Investigating Dispersion and np.save/np.load
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This example, contained in ``tutorials/simo-tut_03_1-dispersion-npload.py`` calculates the acoustic dispersion diagram for the problem in the previous tutorial and classifies the modes according to the point group symmetry class.

.. literalinclude:: ../../tutorials/simo-tut_03_1-dispersion-npload.py
    :lines: 0-


.. figure:: ../../tutorials/tut_03_1-dispersion_npload_symmetrised.png
   :scale: 70 %
   
   Acoustic dispersion diagram with modes categorised by symmetry as in Table 1 of "Formal selection rules for Brillouin scattering in integrated waveguides and structured fibers" by C. Wolff, M. J. Steel, and C. G. Poulton ``https://doi.org/10.1364/OE.22.032489``

.. raw:: latex

    \clearpage



Tutorial 3b -- Investigating Dispersion and multiprocessing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This tutorial, contained in ``tutorials/simo-tut_03_2-dispersion-multicore.py`` continues the study of acoustic dispersion and demonstrates the use of Python multiprocessor calls to increase speed.

.. literalinclude:: ../../tutorials/simo-tut_03_2-dispersion-multicore.py
    :lines: 0-


.. figure:: ../../tutorials/tut_03_2-dispersion_multicore.png
   :scale: 70 %
   
   Acoustic dispersion diagram ploted as lines.

.. raw:: latex

    \clearpage


Tutorial 4 -- Parameter Scan of Widths
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This tutorial, contained in ``tutorials/simo-tut_04_scan_widths.py`` demonstrates use of a parameter scan, in this case of the width of the silicon rectangular waveguide, to understand the behaviour of the Brillouin gain.

.. literalinclude:: ../../tutorials/simo-tut_04-scan_widths.py
    :lines: 0-


.. figure:: ../../tutorials/tut_04-gain_spectra-waterfall.png
   :scale: 70 %
   
   Gain spectra as function of waveguide width.

.. raw:: latex

    \clearpage


Tutorial 5 -- Convergence Study
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This tutorial, contained in ``tutorials/simo-tut_05_convergence_study.py`` demonstrates a scan of numerical parameters for our standard silicon-in-air problem to test the convergence of the calculation results.

.. literalinclude:: ../../tutorials/simo-tut_05-convergence_study.py
    :lines: 0-


.. figure:: ../../tutorials/tut_05-convergence-freq_EM.png
   :scale: 50 %
   
   Convergence of optical mode frequencies.


.. figure:: ../../tutorials/tut_05-convergence-freq_AC.png
   :scale: 50 %
   
   Convergence of acoustic mode frequencies.


.. figure:: ../../tutorials/tut_05-convergence-Gain_PE.png
   :scale: 50 %
   
   Convergence of photoelastic gain.


.. figure:: ../../tutorials/tut_05-convergence-Gain_MB.png
   :scale: 50 %
   
   Convergence of moving boundary gain.


.. figure:: ../../tutorials/tut_05-convergence-Gain.png
   :scale: 50 %
   
   Convergence of total gain.

.. raw:: latex

    \clearpage


Tutorial 6 -- Silica Nanowire 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In this tutorial, contained in ``tutorials/simo-tut_06_silica_nanowire.py`` we start to explore some different structures, in this case a silica nanowire surrounded by vacuum.

.. literalinclude:: ../../tutorials/simo-tut_06-silica_nanowire.py
    :lines: 0-


.. figure:: ../../tutorials/tut_06-gain_spectra-MB_PE_comps_SiO2_NW.png
   :scale: 50 %
   
   Gain spectra showing gain due to photoelastic effect, gain due to moving boundary effect, and total gain.

.. raw:: latex

    \clearpage


Tutorial 7 -- Slot Waveguide
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This tutorial, contained in, ``tutorials/simo-tut_07-slot.py`` examines backward SBS in a more complex structure: chalcogenide soft glass (:math:`\text{As}_2\text{S}_3`) embedded in a silicon slot waveguide on a silica slab. 

.. literalinclude:: ../../tutorials/simo-tut_07-slot.py
    :lines: 0-


.. figure:: ../../tutorials/tut_07-gain_spectra-MB_PE_comps_slot.png
   :scale: 50 %
   
   Gain spectra showing gain due to photoelastic effect, gain due to moving boundary effect, and total gain.

.. raw:: latex

    \clearpage


Tutorial 8 -- Slot Waveguide Scan Covering
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This tutorial, contained in, ``tutorials/simo-tut_08-slot_coated-scan.py`` continues the study of the previous slot waveguide, by examining the thickness dependence of a silica capping layer.

.. literalinclude:: ../../tutorials/simo-tut_08-slot_coated-scan.py
    :lines: 0-


.. figure:: ../../tutorials/tut_08-freq_changes.png
   :scale: 50 %
   
   Acoustic frequencies as function of covering layer thickness.

.. raw:: latex

    \clearpage


Tutorial 9 -- Anisotropic Elastic Materials 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This tutorial, contained in, ``tutorials/simo-tut_09-anisotropy.py`` improves the treatment of the silicon rectangular waveguide by accounting for the anisotropic elastic properties of silicon (simply by referencing a different material file for silicon).

.. literalinclude:: ../../tutorials/simo-tut_09-anisotropy.py


.. raw:: latex

    \clearpage


Tutorial 10 -- Multilayered 'Onion'
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This tutorial, contained in, ``tutorials/simo-tut_10-onion.py`` shows how to create multi-layered circular structures.

.. literalinclude:: ../../tutorials/simo-tut_10-onion.py


.. raw:: latex

    \clearpage


Tutorial 11 -- Two-layered 'Onion'
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This tutorial, contained in, ``tutorials/simo-tut_11a-onion2.py`` demonstrates use of the two layered onion structure which generates a more efficient mesh than the full onion template.

.. literalinclude:: ../../tutorials/simo-tut_11a-onion2.py


.. raw:: latex

    \clearpage


Tutorial 12 -- SMF-28 fibre
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This tutorial, contained in, ``tutorials/simo-tut_12-smf28.py`` models backward SBS in the standard SMF28 fibre using the onion2 template.

.. literalinclude:: ../../tutorials/simo-tut_11a-onion2.py


.. raw:: latex

    \clearpage





.. _sec-josab-label:

JOSA B Tutorial
---------------------

Mike Smith et al. have used NumBAT throughout their 2021 SBS tutorial paper,
published in J. Opt. Soc. Am. B.
.. (see
..  M. Smith et al, FIXME
.. `Generation of phonons from electrostriction in small-core optical waveguides 
.. <http://dx.doi.org/10.1063/1.4801936>`_, *JOSA B* **3**, 042109 (2021).
.. )
This tutorial works through backward, forward, and intermodal forward SBS.
The simulation scripts and resultant mode fields are shown below.


BSBS - Circular Waveguide - Silica
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: ../../JOSAB_tutorial/simo-josab-BSBS-1umcylwg-SiO2.py
    :lines: 0-


.. figure:: ../../JOSAB_tutorial/bsbs-josab-1umSiO2fields/EM_E_field_1.png
   :scale: 50 %

   Fundamental optical mode fields.

.. figure:: ../../JOSAB_tutorial/bsbs-josab-1umSiO2fields/EM_E_field_1_Eabs.png
   :scale: 50 %

   Fundamental optical mode fields.

.. figure:: ../../JOSAB_tutorial/bsbs-josab-1umSiO2fields/EM_E_field_1_Et.png
   :scale: 50 %

   Fundamental optical mode fields.

.. figure:: ../../JOSAB_tutorial/bsbs-josab-1umSiO2fields/AC_field_28.png
   :scale: 50 %

   Fundamental acoustic mode fields.

.. figure:: ../../JOSAB_tutorial/bsbs-josab-1umSiO2fields/AC_field_28_uabs.png
   :scale: 50 %

   Fundamental acoustic mode fields.

.. figure:: ../../JOSAB_tutorial/bsbs-josab-1umSiO2fields/AC_field_28_ut.png
   :scale: 50 %

   Fundamental acoustic mode fields.

.. raw:: latex

    \clearpage


BSBS - Rectangular Waveguide - Silicon
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: ../../JOSAB_tutorial/simo-josab-BSBS-450x200nmrectwg-Si.py
    :lines: 0-

.. figure:: ../../JOSAB_tutorial/bsbs-josab-450x200nmSifields/EM_E_field_0.png
   :scale: 50 %

   Fundamental optical mode fields.

.. figure:: ../../JOSAB_tutorial/bsbs-josab-450x200nmSifields/EM_E_field_0_Eabs.png
   :scale: 50 %

   Fundamental optical mode fields.

.. figure:: ../../JOSAB_tutorial/bsbs-josab-450x200nmSifields/EM_E_field_0_Et.png
   :scale: 50 %

   Fundamental optical mode fields.

.. figure:: ../../JOSAB_tutorial/bsbs-josab-450x200nmSifields/AC_field_6.png
   :scale: 50 %

   Fundamental acoustic mode fields.

.. figure:: ../../JOSAB_tutorial/bsbs-josab-450x200nmSifields/AC_field_6_uabs.png
   :scale: 50 %

   Fundamental acoustic mode fields.

.. figure:: ../../JOSAB_tutorial/bsbs-josab-450x200nmSifields/AC_field_6_ut.png
   :scale: 50 %

   Fundamental acoustic mode fields.


Let's also calculate the acoustic dispersion relation for this structure.

.. literalinclude:: ../../JOSAB_tutorial/simo-josab-BSBS-acbands-450x200nmrectwg-Si.py
    :lines: 0-

.. figure:: ../../JOSAB_tutorial/dispersioncurves_classified.png
   :scale: 50 %
   
   Acoustic dispersion diagram with modes categorised by symmetry as in Table 1 of "Formal selection rules for Brillouin scattering in integrated waveguides and structured fibers" by C. Wolff, M. J. Steel, and C. G. Poulton ``https://doi.org/10.1364/OE.22.032489``

.. raw:: latex

    \clearpage



FSBS - Circular Waveguide - Silica
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: ../../JOSAB_tutorial/simo-josab-FSBS-1umcylwg-SiO2.py
    :lines: 0-


.. figure:: ../../JOSAB_tutorial/fsbs-josab-1umSiO2fields/EM_E_field_1.png
   :scale: 50 %

   Fundamental optical mode fields.

.. figure:: ../../JOSAB_tutorial/fsbs-josab-1umSiO2fields/EM_E_field_1_Eabs.png
   :scale: 50 %

   Fundamental optical mode fields.

.. figure:: ../../JOSAB_tutorial/fsbs-josab-1umSiO2fields/EM_E_field_1_Et.png
   :scale: 50 %

   Fundamental optical mode fields.

.. figure:: ../../JOSAB_tutorial/fsbs-josab-1umSiO2fields/AC_field_7.png
   :scale: 50 %

   Fundamental acoustic mode fields.

.. figure:: ../../JOSAB_tutorial/fsbs-josab-1umSiO2fields/AC_field_7_uabs.png
   :scale: 50 %

   Fundamental acoustic mode fields.

.. figure:: ../../JOSAB_tutorial/fsbs-josab-1umSiO2fields/AC_field_7_ut.png
   :scale: 50 %

   Fundamental acoustic mode fields.

.. raw:: latex

    \clearpage



FSBS - Rectangular Waveguide - Silicon
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: ../../JOSAB_tutorial/simo-josab-FSBS-450x200nmrectwg-Si.py
    :lines: 0-


.. figure:: ../../JOSAB_tutorial/fsbs-josab-450x200nmSifields/EM_E_field_0.png
   :scale: 50 %

   Fundamental optical mode fields.

.. figure:: ../../JOSAB_tutorial/fsbs-josab-450x200nmSifields/EM_E_field_0_Eabs.png
   :scale: 50 %

   Fundamental optical mode fields.

.. figure:: ../../JOSAB_tutorial/fsbs-josab-450x200nmSifields/EM_E_field_0_Et.png
   :scale: 50 %

   Fundamental optical mode fields.

.. figure:: ../../JOSAB_tutorial/fsbs-josab-450x200nmSifields/AC_field_6.png
   :scale: 50 %

   Fundamental acoustic mode fields.

.. figure:: ../../JOSAB_tutorial/fsbs-josab-450x200nmSifields/AC_field_6_uabs.png
   :scale: 50 %

   Fundamental acoustic mode fields.

.. figure:: ../../JOSAB_tutorial/fsbs-josab-450x200nmSifields/AC_field_6_ut.png
   :scale: 50 %

   Fundamental acoustic mode fields.

.. raw:: latex

    \clearpage



IFSBS - Circular Waveguide - Silica
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: ../../JOSAB_tutorial/simo-josab-IFSBS-1umcylwg-SiO2.py
    :lines: 0-


.. figure:: ../../JOSAB_tutorial/ifsbs-josab-1umSiO2fields/EM_E_field_0.png
   :scale: 50 %

   Fundamental optical mode fields.

.. figure:: ../../JOSAB_tutorial/ifsbs-josab-1umSiO2fields/EM_E_field_0_Eabs.png
   :scale: 50 %

   Fundamental optical mode fields.

.. figure:: ../../JOSAB_tutorial/ifsbs-josab-1umSiO2fields/EM_E_field_0_Et.png
   :scale: 50 %

   Fundamental optical mode fields.


.. figure:: ../../JOSAB_tutorial/ifsbs-josab-1umSiO2fields/EM_E_field_1.png
   :scale: 50 %

   Second order optical mode fields.

.. figure:: ../../JOSAB_tutorial/ifsbs-josab-1umSiO2fields/EM_E_field_1_Eabs.png
   :scale: 50 %

   Second order optical mode fields.

.. figure:: ../../JOSAB_tutorial/ifsbs-josab-1umSiO2fields/EM_E_field_1_Et.png
   :scale: 50 %

   Second order optical mode fields.


.. figure:: ../../JOSAB_tutorial/ifsbs-josab-1umSiO2fields/AC_field_6.png
   :scale: 50 %

   Fundamental acoustic mode fields.

.. figure:: ../../JOSAB_tutorial/ifsbs-josab-1umSiO2fields/AC_field_6_uabs.png
   :scale: 50 %

   Fundamental acoustic mode fields.

.. figure:: ../../JOSAB_tutorial/ifsbs-josab-1umSiO2fields/AC_field_6_ut.png
   :scale: 50 %

   Fundamental acoustic mode fields.

.. raw:: latex

    \clearpage



IFSBS - Rectangular Waveguide - Silicon
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: ../../JOSAB_tutorial/simo-josab-IFSBS-450x200nmrectwg-Si.py
    :lines: 0-


.. figure:: ../../JOSAB_tutorial/ifsbs-josab-450x200nmSifields/EM_E_field_0.png
   :scale: 50 %

   Fundamental optical mode fields.

.. figure:: ../../JOSAB_tutorial/ifsbs-josab-450x200nmSifields/EM_E_field_0_Eabs.png
   :scale: 50 %

   Fundamental optical mode fields.

.. figure:: ../../JOSAB_tutorial/ifsbs-josab-450x200nmSifields/EM_E_field_0_Et.png
   :scale: 50 %

   Fundamental optical mode fields.

.. figure:: ../../JOSAB_tutorial/ifsbs-josab-450x200nmSifields/EM_E_field_0.png
   :scale: 50 %

   Second order optical mode fields.

.. figure:: ../../JOSAB_tutorial/ifsbs-josab-450x200nmSifields/EM_E_field_0_Eabs.png
   :scale: 50 %

   Second order optical mode fields.

.. figure:: ../../JOSAB_tutorial/ifsbs-josab-450x200nmSifields/EM_E_field_0_Et.png
   :scale: 50 %

   Second order optical mode fields.

.. figure:: ../../JOSAB_tutorial/ifsbs-josab-450x200nmSifields/AC_field_2.png
   :scale: 50 %

   Fundamental acoustic mode fields.

.. figure:: ../../JOSAB_tutorial/ifsbs-josab-450x200nmSifields/AC_field_2_uabs.png
   :scale: 50 %

   Fundamental acoustic mode fields.

.. figure:: ../../JOSAB_tutorial/ifsbs-josab-450x200nmSifields/AC_field_2_ut.png
   :scale: 50 %

   Fundamental acoustic mode fields.

.. raw:: latex

    \clearpage







.. _sec-literature-label:

Literature Examples
---------------------

Having become somewhat familiar with NumBAT, we now set out to replicate a number of examples 
from the recent literature located in the ``lit_examples`` directory.
The examples are presented in chronological order. 
We note the particular importance of examples 5-8 which include experimental and numerical results that are in good agreement.


LitEx 1 -- Laude and Beugnot, *AIP Advances* (2013): BSBS in a silica rectangular waveguide
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This example ``simo-lit_01-Laude-AIPAdv_2013-silica.py``
is based on the calculation of backward SBS
in a small rectangular silica waveguide described in V. Laude and J.-C. Beugnot, 
`Generation of phonons from electrostriction in small-core optical waveguides 
<http://dx.doi.org/10.1063/1.4801936>`_, *AIP Advances* **3**, 042109 (2013).

Observe the use of a material named ``materials.materials_dict["SiO2_2013_Laude"]`` 
specifically modelled on the parameters in this paper.
This technique allows users to easily compare exactly to other authors
without changing their preferred material values for their own samples and experiments.

.. literalinclude:: ../../lit_examples/simo-lit_01-Laude-AIPAdv_2013-silica.py
    :lines: 0-


.. figure:: ../../lit_examples/lit_01-fields/EM_E_field_0.png
   :scale: 50 %
   
   Fundamental optical mode fields.


.. figure:: ../../lit_examples/lit_01-fields/AC_field_4.png
   :scale: 50 %
   
   High gain acoustic mode, marked as C in paper.


.. figure:: ../../lit_examples/lit_01-fields/AC_field_55.png
   :scale: 50 %
   
   High gain acoustic mode, marked as D in paper.


.. figure:: ../../lit_examples/lit_01-gain_spectra-MB_PE_comps-logy.png
   :scale: 50 %
   
   Gain spectra on semilogy axis.
   

.. figure:: ../../lit_examples/lit_01-gain_spectra-MB_PE_comps_zoom.png
   :scale: 50 %
   
   Gain spectra zoomed in on mode D.

.. raw:: latex

    \clearpage

LitEx 2 -- Laude and Beungot, *AIP Advances* (2013): BSBS in a rectangular silicon waveguide 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This example in ``simo-lit_02-Laude-AIPAdv_2013-silicon.py`` again follows the paper of V. Laude and J.-C. Beugnot, 
`Generation of phonons from electrostriction in small-core optical waveguides 
<http://dx.doi.org/10.1063/1.4801936>`_, *AIP Advances* **3**, 042109 (2013),
but this time looks at the *silicon* waveguide case.

.. literalinclude:: ../../lit_examples/simo-lit_02-Laude-AIPAdv_2013-silicon.py
    :lines: 0-


.. figure:: ../../lit_examples/lit_02-fields/AC_field_4.png
   :scale: 50 %
   
   High gain acoustic mode, marked as G in paper.


.. figure:: ../../lit_examples/lit_02-gain_spectra-MB_PE_comps-logy.png
   :scale: 50 %
   
   Gain spectra on semilogy axis.

.. raw:: latex

    \clearpage

LitEx 3 -- Beugnot *et al*, *Nature Communications* (2014): BSBS in a tapered fibre - scanning widths
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This example, in ``simo-lit_03-Beugnot-NatComm_2014.py``,
is based on the calculation of backward SBS
in a micron scale optical fibre described in J.-C. Beugnot *et al.*, 
`Brillouin light scattering from surface acoustic waves in a subwavelength-diameter optical fibre
<http://dx.doi.org/10.1038/ncomms6242>`_, *Nature Communications* **5**, 5242 (2014).

.. literalinclude:: ../../lit_examples/simo-lit_03-Beugnot-NatComm_2014.py
    :lines: 0-


.. figure:: ../../lit_examples/lit_03-gain-width_scan.png
   :scale: 50 %
   
   Full acoustic wave spectrum for silica microwire, as per Fig. 4a in paper.

.. raw:: latex

    \clearpage

LitEx 4 -- Van Laer *et al*, *Nature Photonics* (2015): FSBF in a waveguide on a pedestal
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This example, in ``simo-lit_04-pillar-Van_Laer-NatPhot_2015.py``, 
is based on the calculation of forward SBS
in a pedestal silicon waveguide described in R. Van Laer *et al.*, 
`Interaction between light and highly confined hypersound in a silicon photonic nanowire 
<http://dx.doi.org/10.1038/ncomms6242>`_, *Nature Photonics* **9**, 199 (2015).

Note that the absence of an absorptive boundary in the acoustic model 
causes a problem where the slab layer significantly distorting acoustic modes.
Adding this feature is a priority for a future release of NumBAT.
The following example shows an approximate way to avoid the problem for now.

.. literalinclude:: ../../lit_examples/simo-lit_04-pillar-Van_Laer-NP_2015.py
    :lines: 0-


.. figure:: ../../lit_examples/lit_04-pillar-fields/EM_E_field_0.png
   :scale: 50 %
   
   Fundamental optical mode fields.


.. figure:: ../../lit_examples/lit_04-pillar-fields/AC_field_38.png
   :scale: 50 %
   
   Dominant high gain acoustic mode.
   Note how the absence of an absorptive boundary on the SiO2 slab causes this layer to significantly distorted the acoustic modes.


We may also choose to study the simplified situation where the pedestal is removed.


.. literalinclude:: ../../lit_examples/simo-lit_04-no_pillar-Van_Laer-NP_2015.py
    :lines: 0-

Which gives good agreement for the gain spectrum.

.. figure:: ../../lit_examples/lit_04-no_pillar-gain_spectra-MB_PE_comps.png
   :scale: 50 %
   
   Gain spectrum for the simplified case of a waveguide surrounded by vacuum.


.. raw:: latex

    \clearpage

LitEx 5 -- 2015 - Van Laer *et al*, *New Journal of Physics* (2015): FSBF in a waveguide without a pedestal
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This example, in ``simo-lit_05-Van_Laer-NJP_2015.py``, continues  
the study of forward SBS
in a pedestal silicon waveguide described in R. Van Laer *et al.*, 
`Interaction between light and highly confined hypersound in a silicon photonic nanowire 
<http://dx.doi.org/10.1038/ncomms6242>`_, *Nature Photonics* **9**, 199 (2015).

In this case, we simply remove the pedestal and model the main rectangular waveguide.
This makes the acoustic loss calculation incorrect but avoids the problem of acoustic
energy being excessively concentrated in the substrate.

.. literalinclude:: ../../lit_examples/simo-lit_05-Van_Laer-NJP_2015.py
    :lines: 0-


.. figure:: ../../lit_examples/lit_05-fields/EM_E_field_0.png
   :scale: 50 %
   
   Fundamental optical mode fields.


.. figure:: ../../lit_examples/lit_05-fields/AC_field_6.png
   :scale: 50 %
   
   Dominant high gain acoustic mode.


.. raw:: latex

    \clearpage

LitEx 6 -- Florez *et al*, *Nature Communications* (2016): BSBS self-cancellation in a tapered fibre (:math:`d = 550` nm)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This example, in ``simo-lit_06_1-Florez-NatComm_2016-d550nm.py``,
looks at the phenomenon of Brillouin "self-cancellation" due to 
the electrostrictive and radiation pressure effects acting with opposite sign. 
This was described in O. Florez *et al.*, `Brillouin self-cancellation 
<http://dx.doi.org/10.1038/ncomms11759>`_, *Nature Communications* **7**, 11759 (2016).

.. literalinclude:: ../../lit_examples/simo-lit_06_1-Florez-NatComm_2016-d550nm.py
    :lines: 0-


.. figure:: ../../lit_examples/lit_06_1-fields/AC_field_4.png
   :scale: 50 %
   
   :math:`TR_{21}` acoustic mode fields of a nanowire with diameter 550 nm.


.. figure:: ../../lit_examples/lit_06_1-fields/AC_field_5.png
   :scale: 50 %
   
   :math:`R_{01}` acoustic mode fields of a nanowire with diameter 550 nm.


.. figure:: ../../lit_examples/lit_06_1-gain_spectra-MB_PE_comps.png
   :scale: 50 %
   
   Gain spectra of a nanowire with diameter 550 nm, matching blue curve of Fig. 3b in paper.


.. figure:: ../../lit_examples/lit_06_1-gain_spectra-MB_PE_comps-5.png
   :scale: 50 %

.. figure:: ../../lit_examples/lit_06_1-gain_spectra-MB_PE_comps-6.png
   :scale: 50 %

.. figure:: ../../lit_examples/lit_06_1-gain_spectra-MB_PE_comps-8.png
   :scale: 50 %

.. figure:: ../../lit_examples/lit_06_1-gain_spectra-MB_PE_comps-11.png
   :scale: 50 %
   
   Zoomed in gain spectra around gaint peaks of 550 nm diameter NW.

.. raw:: latex

    \clearpage

LitEx 6b -- Florez *et al*, *Nature Communications* (2016): BSBS self-cancellation in a tapered fibre (:math:`d = 1160` nm)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This example, in ``simo-lit_06_2-Florez-NatComm_2016-1160nm.py``, again looks at the paper 
O. Florez *et al.*, `Brillouin self-cancellation <http://dx.doi.org/10.1038/ncomms11759>`_, *Nature Communications* **7**, 11759 (2016),
but now for a wider core.

.. literalinclude:: ../../lit_examples/simo-lit_06_2-Florez-NatComm_2016-d1160nm.py
    :lines: 0-


.. figure:: ../../lit_examples/lit_06_2-gain_spectra-MB_PE_comps.png
   :scale: 50 %
   
   Gain spectra of a nanowire with diameter 1160 nm, as in Fig. 4 of Florez, showing near perfect cancellation at 5.4 GHz.


.. figure:: ../../lit_examples/lit_06_2-gain_spectra-MB_PE_comps-logy.png
   :scale: 50 %
   
   Gain spectra of a nanowire with diameter 1160 nm, as in Fig. 4 of paper, showing near perfect cancellation at 5.4 GHz.


.. raw:: latex

    \clearpage

LitEx 7 -- Kittlaus *et al*,  *Nature Photonics* (2016), FSBF in a silicon rib waveguide
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This example, in ``../../lit_examples/simo-lit_07-Kittlaus-NatPhot_2016.py``,
explores a first geometry showing large forward SBS in silicon
as described in E. Kittlaus *et al.*, `Large Brillouin amplification in silicon 
<http://dx.doi.org/10.1038/nphoton.2016.112>`_, *Nature Photonics* **10**, 463 (2016).


.. literalinclude:: ../../lit_examples/simo-lit_07-Kittlaus-NatPhot_2016.py
    :lines: 0-


.. figure:: ../../lit_examples/lit_07-fields/EM_E_field_0.png
   :scale: 50 %
   
   Fundamental optical mode fields.


.. figure:: ../../lit_examples/lit_07-fields/AC_field_19.png
   :scale: 50 %
   
   Dominant high gain acoustic mode.


.. figure:: ../../lit_examples/lit_07-gain_spectra-MB_PE_comps.png
   :scale: 50 %
   
   Gain spectra showing gain due to photoelastic effect, gain due to moving boundary effect, and total gain.


.. raw:: latex

    \clearpage

LitEx 8 -- Kittlaus *et al*, *Nature Communications* (2017): Intermodal FSBF in a silicon waveguide 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This example (``simo-lit_08-Kittlaus-NatComm_2017.py``), also from the Yale group,  examines intermode forward Brillouin scattering in silicon.

.. literalinclude:: ../../lit_examples/simo-lit_08-Kittlaus-NatComm_2017.py
    :lines: 0-


.. figure:: ../../lit_examples/lit_08-fields/EM_E_field_0.png
   :scale: 50 %
   
   Fundamental (symmetric TE-like) optical mode fields.


.. figure:: ../../lit_examples/lit_08-fields/EM_E_field_1.png
   :scale: 50 %
   
   2nd lowest order (anti-symmetric TE-like) optical mode fields.


.. figure:: ../../lit_examples/lit_08-fields/AC_field_23.png
   :scale: 50 %
   
   Dominant high gain acoustic mode.


.. figure:: ../../lit_examples/lit_08-gain_spectra-MB_PE_comps.png
   :scale: 50 %
   
   Gain spectra showing gain due to photoelastic effect, gain due to moving boundary effect, and total gain.

.. raw:: latex

    \clearpage


LitEx 9 -- Morrison *et al*, *Optica* (2017):  BSBS in a chalcogenide rib waveguide
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This example, in ``simo-lit_09-Morrison-Optica_2017.py``, from the Sydney group examines backward SBS in a chalcogenide rib waveguide.

.. literalinclude:: ../../lit_examples/simo-lit_09-Morrison-Optica_2017.py
    :lines: 0-


.. figure:: ../../lit_examples/lit_09-gain_spectra-MB_PE_comps.png
   :scale: 50 %
   
   Gain spectra showing gain due to photoelastic effect, gain due to moving boundary effect, and total gain.

.. include:: numbatdefs.txt

.. _chap-intro-label:

****************************
Introduction to |NUMBAT|
****************************


.. role:: raw-math(raw)
    :format: latex html

.. figure:: NumBAT_logo.png
   :scale: 40 %

Introduction
================


|NUMBAT|, the *Numerical Brillouin Analysis Tool*, is a software tool integrating electromagnetic and acoustic mode solvers to calculate the interactions of optical and acoustic waves in waveguides. Most notably, this includes Stimulated Brillouin Scattering (SBS) frequency shifts and optical gains.

This chapter provides some background on the capabilities and techniques used in |NUMBAT|.
If you would like to get straight to computations, jump ahead to the installation and setup instructions in :ref:`chap-install-label`.

Goals
================
|NUMBAT| is designed primarily to calculate the optical gain response from
stimulated Brillouin scattering (SBS) in integrated waveguides. It uses finite element algorithms to solve the electromagnetic and acoustic modes of a wide range of 2D waveguide structures. It can account for photoelastic/electrostriction and moving boundary/radiation pressure effects, as well as arbitrary acoustic anisotropy.

|NUMBAT| also supports user-defined material properties and we hope its creation will drive a community-driven set of standard properties and geometries which will allow all groups to test and validate each other's work.

A full description of |NUMBAT|'s physical and numerical algorithms is available in the article B.C.P Sturmberg at al., "Finite element analysis of stimulated Brillouin scattering in integrated photonic waveguides", *J. Lightwave Technol.*  **37**, 3791-3804 (2019),
available at `<https://dx.doi.org/10.1109/JLT.2019.2920844>`_.


|NUMBAT| is open-source software and the authors welcome additions to the code.  Details for how to contribute are available in :ref:`sec-contribute-label`.

Development team
================

|NUMBAT| was developed by Bjorn Sturmberg, Kokou Dossou, Blair Morrison, Chris
Poulton and Michael Steel in a collaboration between Macquarie University, the
University of Technology Sydney, and the University of Sydney.

We thank Christian Wolff, Mike Smith and Mikolaj Schmidt for contributions.


Citing |NumBAT|
===============
If you use |NumBAT| in published work, we would appreciate a citation
to B.C.P Sturmberg at al.,
"Finite element analysis of stimulated Brillouin scattering in integrated
photonic waveguides", *J. Lightwave Technol.*  **37**, 3791-3804 (2019),
available at `<https://dx.doi.org/10.1109/JLT.2019.2920844>`_
and `<https://arxiv.org/abs/1811.10219>`_,
and a link to the github page at `<https://github.com/michaeljsteel/NumBAT>`_.


.. _sec-contribute-label:

Contributing to NumBAT
================================
NumBAT is open source software licensed under the GPL with all source and documentation available at `github.com <https://github.com/michaeljsteel/NumBAT.git>`_. We welcome additions to NumBAT code, documentation and the materials library. Interested users should fork the standard release from github and make a pull request when ready.  For major changes, we strongly suggest contacting the NumBAT team before starting work at ``michael.steel@mq.edu.au``.


Seeking assistance
================================
We will do our best to support users of |NUMBAT| within the time available.
All requests should be sent to ``michael.steel@mq.edu.au``.

   * For assistance with installing and building |NUMBAT|, please see the instructions
     in :ref:`sec-helpinstall-label` and collect all the required information before writing.

   * For assistance with calculations once |NUMBAT| is working, please send an email
      containing the following information:

      * Your platform (Linux, Windows, MacOS) and specific operating system
      * How recently you have updated your |NUMBAT| code
      * A python script that demonstrates the issue you are trying to solve
        Please make the script as short in both code length and execution time as possible while still exhibiting the issue of concern.



Release notes
=============

Version 2.0
-----------

A number of API changes have been made in |NUMBAT| 2.0 to tidy up the interface and make plotting and analysis simpler and more powerful.
You will need to make some changes to existing files to run in |NUMBAT| 2.0.  Your best guide to new capabilities and API changes is to look through the code in the tutorial examples.

Some key changes you will need to make are as follows:
 * On Linux, the fortran Makefile is now designed to work with a virtual environment python to avoid dependencies on your system python.
 * There is a new core |NUMBAT| module ``numbat`` that should be imported before any other |NUMBAT| modules.
 * It should no longer be necessary to import the ``object`` or ``Numbat`` (note different case) modules.
 * The first call to any |NUMBAT| code should be to create a |NUMBAT| application object by calling ``nbapp = numbat.NumBATApp()``.
 * The default output prefix can now be set as an argument to ``numbat.NumBATApp()``. All output can be directed to a sub-folder of the starting directory with a second argument: ``nbapp = numbat.NumBATApp('tmp', 'tmpdir')``.
 * The waveguide class ``Struct`` has been renamed to ``Structure``.
 * A waveguide is now constructed using ``nbapp.make_waveguide`` rather than ``object.Structure``.
 * The interface for creating materials has changed. You now call the ``materials. make_material(`` *name* ``)`` function. For example ``material_a = materials.make_material('Vacuum')``
 * To access an existing material in an  existing ``Structure`` object (say, in a variable called ``wguide``) use ``wguide.get_material(`` *label* ``)`` For example, ``mat_a = wguide.get_material('b')`` where the allowed labels are ``bkg`` and the letters ``a`` to ``r``.
 * The member name for refractive index in a ``Material`` object has changed from ``n`` to ``refindex_n``.
 * The member name for density in a ``Material`` object has changed from ``n`` to ``rho``.
 * Due to a change in parameters, the function ``plotting.gain_spectra`` is deprecated and replaced by ``plotting.plot_gain_spectra`` with the following changes:
      * The frequency arguments ``freq_min`` and ``freq_max`` should now be passed in units of Hz, not GHz.
      * The argument ``k_AC`` has been removed.
 * In all functions the parameter ``prefix_str`` has been renamed to ``prefix`` for brevity. Using the default output settings in ``NumBATApp()``, these should be rarely needed.
 * All waveguides are now specified as individual plugin classes in the files ``backend/msh/user_waveguides.json`` and ``backend/msh/user_meshes.py``.  These files provide useful examples of how to design and load new waveguide templates. See the following chapter for more details.


About our mascot
================================
The **numbat** (*Myrmecobius fasciatus*) is a delightful insect-eating marsupial
from Western Australia, of which it is the official state animal.
It has two other common Aboriginal names, *noombat* in the Nyungar language,
and *walpurti* in the Pitjantjatjara language.

As a carnivorous marsiupial, they belong to the order Dasyuromorphia, closely related to quolls and the famed thylacines which had similar markings on their lower back.
Once found across southern Australia, numbats are now confined to small local groups
in Western Australia and the species has Endangered status.

.. figure:: ./numbat_face.jpg
   :scale: 40 %

   A numbat at Perth zoo in 2010. `(Creative commons) <https://commons.wikimedia.org/wiki/File:Numbat_Face.jpg>`_.

Apart from the distinctive striped back (which we like to think of as an acoustic wave made flesh), numbats have a number of unique properties. They are the only fully diurnal marsupial.
They are insectivores and eat exclusively termites, perhaps 20000 each day!

To find out how you can support the care and revitalisation of this beautiful animal, check out the work at
`projectnumbat <numbat.org.au>`_ and the `Australian Wildlife Conservancy <https://www.australianwildlife.org/wildlife/numbat>`_.


Acknowledgements
===================

We acknowledge the elders past and present of the following First Nations people,
on whose unceded lands |NUMBAT| has been developed: the Wallumattagal clan of the Dharug
nation, the Cammeraygal people, and the Gadigal people of the Eora nation.


Development of |NumBAT| has been supported in part by the
Australian Research Council under Discovery Projects
DP130100832, DP160101691, DP200101893 and DP220100488.

.. include:: numbatdefs.txt

Introduction
-------------------------------------------

Dr Christian  Wolff and colleagues  have used |NUMBAT| throughout their 2021 SBS tutorial paper
`Brillouin scattering—theory and experiment: tutorial <http://dx.doi.org/10.1364/JOSAB.416747>`_ published in J. Opt. Soc. Am. B.

This set of examples works through their discussions of backward SBS, 
forward Brillouin scattering, and intermodal forward Brillouin scattering.



Example 1 -- Backward SBS in a circular silica waveguide 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Figure 15 in the paper shows the fundamental optical field of a silica 1 micron waveguide, with the gain and other parameters shown in Table 2.
The corresponding results generated with ``simo-josab-01.py`` are as follows:


.. figure:: ../../josab_tutorial/josab_01-fields/EM_E_field_01.png
   :scale: 50 %

   Fundamental optical mode fields.

.. figure:: ../../josab_tutorial/josab_01-fields/AC_field_28.png
   :scale: 50 %

   Elastic mode with largest SBS gain.

.. figure:: ../../josab_tutorial/josab_01-gain_spectra.png
   :scale: 50 %

   Gain spectrum for the silica waveguide.



.. raw:: latex

    \clearpage


Example 2 -- Backward SBS in a rectangular silicon waveguide 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Figure 14 in the paper calculates the backwards SBS properties  of a 
rectangular :math:`450x200` nm silicon waveguide.
The corresponding results generated with ``simo-josab-02.py`` are as follows:


.. figure:: ../../josab_tutorial/josab_02a-gain_spectra.png
.. .. figure:: ../../josab_tutorial/josab_02b-disp-qnu.png
   :scale: 50 %


The fields and gain parameters are as follows:

.. figure:: ../../josab_tutorial/josab_02a-fields/EM_E_field_00.png
   :scale: 50 %

   Fundamental optical mode fields.

.. figure:: ../../josab_tutorial/josab_02a-fields/AC_field_06.png
   :scale: 50 %

   Fundamental elastic mode fields.

We can reproduce Fig. 13 showing the elastic dispersion of this waveguide
silicon waveguide using ``simo-josab-02b-acdisp.py``.

.. figure:: ../../josab_tutorial/josab_02b-disp-qnu.png
   :scale: 50 %
   
   Acoustic dispersion diagram with modes categorised by symmetry as in Table 1 of "Formal selection rules for Brillouin scattering in integrated waveguides and structured fibers" by C. Wolff, M. J. Steel, and C. G. Poulton ``https://doi.org/10.1364/OE.22.032489``

.. raw:: latex

    \clearpage



Example 3 -- Forward Brillouin scattering in a circular silica waveguide 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Figure 16 and Table 3 examine the same waveguides in the case of forward
Brillouin scattering. 

These results can be generated with ``simo-josab-03.py`` and ``simo-josab-04.py``.

Let's see the results for the silica cylinder first:

.. figure:: ../../josab_tutorial/josab_03-gain_spectra.png
   :scale: 50 %

   Gain spectrum for forward SBS of the silica cylinder.

.. figure:: ../../josab_tutorial/josab_03-fields/EM_E_field_01.png
   :scale: 50 %

   Fundamental optical mode field.


.. figure:: ../../josab_tutorial/josab_03-fields/AC_field_06.png
   :scale: 50 %

   Elastic mode of maximum gain.

.. figure:: ../../josab_tutorial/josab_03-fields/AC_field_12.png
   :scale: 50 %

   Elastic mode of second highest gain.




.. raw:: latex

    \clearpage



Example 4 -- Forward Brillouin scattering in a rectangular silicon waveguide 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The corresopnding results for the silicon waveguide can be generated with ``simo-josab-04.py``:

.. figure:: ../../josab_tutorial/josab_04-gain_spectra.png
   :scale: 50 %

   Gain spectrum for forward SBS of the silicon waveguide.

.. figure:: ../../josab_tutorial/josab_04-fields/EM_E_field_00.png
   :scale: 50 %

   Fundamental optical mode field.

.. figure:: ../../josab_tutorial/josab_04-fields/AC_field_06.png
   :scale: 50 %

   Elastic mode of maximum gain.


.. raw:: latex

    \clearpage



Example 5 -- Intermodal Forward Brillouin scattering in a circular silica waveguide 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the problem of intermodal FBS, the paper considers coupling between the two lowest optical modes. The elastic mode of highest gain is actually a degenerate pair:


.. figure:: ../../josab_tutorial/josab_05-gain_spectra.png
   :scale: 50 %

   Gain spectrum for intermodal forward SBS of the silica waveguide.

.. figure:: ../../josab_tutorial/josab_05-fields/EM_E_field_00.png
   :scale: 50 %

   Fundamental optical mode field.

.. figure:: ../../josab_tutorial/josab_05-fields/EM_E_field_01.png
   :scale: 50 %

   Second order optical mode field.


.. figure:: ../../josab_tutorial/josab_05-fields/AC_field_06.png
   :scale: 50 %

   Elastic mode field of maximum gain.


.. raw:: latex

    \clearpage



Example 6 -- Intermodal Forward Brillouin scattering in a rectangular silicon waveguide 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Finally, the silicon waveguide generates extraordinarily high gain when operated
in an intermodal configuration:

.. figure:: ../../josab_tutorial/josab_06-gain_spectra.png
   :scale: 50 %

   Gain spectrum for intermodal forward SBS of the silicon waveguide.


.. figure:: ../../josab_tutorial/josab_06-fields/EM_E_field_00.png
   :scale: 50 %

   Fundamental optical mode field.

.. figure:: ../../josab_tutorial/josab_06-fields/EM_E_field_02.png
   :scale: 50 %

   Second order optical mode field.

.. figure:: ../../josab_tutorial/josab_06-fields/AC_field_02.png
   :scale: 50 %

   Elastic mode field of maximum gain.


.. raw:: latex

    \clearpage






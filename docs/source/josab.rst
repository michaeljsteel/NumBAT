
.. include:: numbatdefs.txt

.. _chap-josab-label:

***********************
JOSA-B Tutorial Paper
***********************

Introduction
-------------------------------------------

Dr Christian  Wolff and colleagues  have used |NUMBAT| throughout their 2021 SBS tutorial paper
`Brillouin scattering—theory and experiment: tutorial <http://dx.doi.org/10.1364/JOSAB.416747>`_ published in J. Opt. Soc. Am. B.
This set of examples works through their discussions of backward SBS,
forward Brillouin scattering, and intermodal forward Brillouin scattering.


As the calculations in this paper used |NUMBAT| with essentially the same code in these
tutorials, we do not bother to include the original figures.



Example 1 -- Backward SBS in a circular silica waveguide
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Figure 15 in the paper shows the fundamental optical field of a silica 1 micron waveguide, with the gain and other parameters shown in Table 2.
The corresponding results generated with ``sim-josab-01.py`` are as follows:


.. figure:: ./images/josab_tutorial/josab_01-fields/EM_E_mode_01.png
   :width: 15cm

   Fundamental optical mode fields.

.. figure:: ./images/josab_tutorial/josab_01-fields/AC_mode_28.png
   :width: 15cm

   Elastic mode with largest SBS gain.

.. figure:: ./images/josab_tutorial/josab_01-gain_spectra.png
   :width: 15cm

   Gain spectrum for the silica waveguide.



.. raw:: latex

    \clearpage


Example 2 -- Backward SBS in a rectangular silicon waveguide
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Figure 14 in the paper calculates the backwards SBS properties  of a
rectangular :math:`450x200` nm silicon waveguide.
The corresponding results generated with ``sim-josab-02.py`` are as follows:


.. figure:: ./images/josab_tutorial/josab_02a-gain_spectra.png
   :width: 15cm


The fields and gain parameters are as follows:

.. figure:: ./images/josab_tutorial/josab_02a-fields/EM_E_mode_00.png
   :width: 15cm

   Fundamental optical mode fields.

.. figure:: ./images/josab_tutorial/josab_02a-fields/AC_mode_02.png

   Fundamental elastic mode fields for mode 2.

.. figure:: ./images/josab_tutorial/josab_02a-fields/AC_mode_06.png
   :width: 15cm

   Fundamental elastic mode fields for mode 6.

We can reproduce Fig. 13 showing the elastic dispersion of this waveguide
silicon waveguide using ``sim-josab-02b-acdisp.py``.

.. figure:: ./images/josab_tutorial/josab_02b-disp-qnu.png
   :width: 15cm

   Acoustic dispersion diagram with modes categorised by symmetry as in Table 1 of "Formal selection rules for Brillouin scattering in integrated waveguides and structured fibers" by C. Wolff, M. J. Steel, and C. G. Poulton ``https://doi.org/10.1364/OE.22.032489``

.. raw:: latex

    \clearpage



Example 3 -- Forward Brillouin scattering in a circular silica waveguide
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Figure 16 and Table 3 examine the same waveguides in the case of forward
Brillouin scattering.

These results can be generated with ``sim-josab-03.py`` and ``sim-josab-04.py``.

Let's see the results for the silica cylinder first:

.. figure:: ./images/josab_tutorial/josab_03-gain_spectra.png
   :width: 15cm

   Gain spectrum for forward SBS of the silica cylinder.

.. figure:: ./images/josab_tutorial/josab_03-fields/EM_E_mode_01.png
   :width: 15cm

   Fundamental optical mode field.


.. figure:: ./images/josab_tutorial/josab_03-fields/AC_mode_06.png
   :width: 15cm

   Elastic mode of maximum gain.

.. figure:: ./images/josab_tutorial/josab_03-fields/AC_mode_12.png
   :width: 15cm

   Elastic mode of second highest gain.




.. raw:: latex

    \clearpage



Example 4 -- Forward Brillouin scattering in a rectangular silicon waveguide
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The corresponding results for the silicon waveguide can be generated with ``sim-josab-04.py``:

.. figure:: ./images/josab_tutorial/josab_04-gain_spectra.png
   :width: 15cm

   Gain spectrum for forward SBS of the silicon waveguide.

.. figure:: ./images/josab_tutorial/josab_04-fields/EM_E_mode_00.png
   :width: 15cm

   Fundamental optical mode field.

.. figure:: ./images/josab_tutorial/josab_04-fields/AC_mode_06.png
   :width: 15cm

   Elastic mode of maximum gain.


.. raw:: latex

    \clearpage



Example 5 -- Intermodal Forward Brillouin scattering in a circular silica waveguide
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the problem of intermodal FBS, the paper considers coupling between the two lowest optical modes. The elastic mode of highest gain is actually a degenerate pair:


.. figure:: ./images/josab_tutorial/josab_05-gain_spectra.png
   :width: 15cm

   Gain spectrum for intermodal forward SBS of the silica waveguide.

.. figure:: ./images/josab_tutorial/josab_05-fields/EM_E_mode_00.png
   :width: 15cm

   Fundamental optical mode field.

.. figure:: ./images/josab_tutorial/josab_05-fields/EM_E_mode_01.png
   :width: 15cm

   Second order optical mode field.


.. figure:: ./images/josab_tutorial/josab_05-fields/AC_mode_06.png
   :width: 15cm

   Elastic mode field of maximum gain.


.. raw:: latex

    \clearpage



Example 6 -- Intermodal Forward Brillouin scattering in a rectangular silicon waveguide
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Finally, the silicon waveguide generates extraordinarily high gain when operated
in an intermodal configuration:

.. figure:: ./images/josab_tutorial/josab_06-gain_spectra.png
   :width: 15cm

   Gain spectrum for intermodal forward SBS of the silicon waveguide.


.. figure:: ./images/josab_tutorial/josab_06-fields/EM_E_mode_00.png
   :width: 15cm

   Fundamental optical mode field.

.. figure:: ./images/josab_tutorial/josab_06-fields/EM_E_mode_02.png
   :width: 15cm

   Second order optical mode field.

.. figure:: ./images/josab_tutorial/josab_06-fields/AC_mode_02.png
   :width: 15cm

   Elastic mode field of maximum gain.


.. raw:: latex

    \clearpage






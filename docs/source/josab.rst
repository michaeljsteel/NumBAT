
.. include:: numbatdefs.txt

Introduction
-------------------------------------------

Dr Mike Smith and colleagues  have used |NUMBAT| throughout their 2021 SBS tutorial paper
`Generation of phonons from electrostriction in small-core optical waveguides <http://dx.doi.org/10.1063/1.4801936>`_ published in J. Opt. Soc. Am. B.

This set of examples works through their discussions of backward SBS, 
forward Brillouin scattering, and intermodal forward Brillouin scattering.


Example 1 -- Backward SBS in a circular silica waveguide 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. .. literalinclude:: ../../JOSAB_tutorial/simo-josab-01-BSBS-1umcylwg-SiO2.py
    :lines: 0-


.. figure:: ../../JOSAB_tutorial/josab_01-fields/EM_E_field_01.png
   :scale: 50 %

   Fundamental optical mode fields.

.. figure:: ../../JOSAB_tutorial/josab_01-fields/AC_field_28.png
   :scale: 50 %

   Fundamental elastic mode fields.

.. raw:: latex

    \clearpage


Example 2 -- Backward SBS in a rectangular silicon waveguide 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. .. literalinclude:: ../../JOSAB_tutorial/simo-josab-02-BSBS-450x200nmrectwg-Si.py
    :lines: 0-

.. figure:: ../../JOSAB_tutorial/josab_02a-fields/EM_E_field_00.png
   :scale: 50 %

   Fundamental optical mode fields.

.. figure:: ../../JOSAB_tutorial/josab_02a-fields/AC_field_06.png
   :scale: 50 %

   Fundamental elastic mode fields.

Let's also calculate the elastic dispersion relation for this structure.

.. .. literalinclude:: ../../JOSAB_tutorial/simo-josab-02b-BSBS-acbands-450x200nmrectwg-Si.py
    :lines: 0-

.. figure:: ../../JOSAB_tutorial/josab_02b-dispersion.png
   :scale: 50 %
   
   Acoustic dispersion diagram with modes categorised by symmetry as in Table 1 of "Formal selection rules for Brillouin scattering in integrated waveguides and structured fibers" by C. Wolff, M. J. Steel, and C. G. Poulton ``https://doi.org/10.1364/OE.22.032489``

.. raw:: latex

    \clearpage



Example 3 -- Forward Brillouin scattering in a circular silica waveguide 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. .. literalinclude:: ../../JOSAB_tutorial/simo-josab-03-FSBS-1umcylwg-SiO2.py
    :lines: 0-


.. figure:: ../../JOSAB_tutorial/josab_03-fields/EM_E_field_01.png
   :scale: 50 %

   Fundamental optical mode fields.


.. figure:: ../../JOSAB_tutorial/josab_03-fields/AC_field_07.png
   :scale: 50 %

   Fundamental elastic mode fields.

.. raw:: latex

    \clearpage



Example 4 -- Forward Brillouin scattering in a rectangular silicon waveguide 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. .. literalinclude:: ../../JOSAB_tutorial/simo-josab-04-FSBS-450x200nmrectwg-Si.py
    :lines: 0-


.. figure:: ../../JOSAB_tutorial/josab_04-fields/EM_E_field_00.png
   :scale: 50 %

   Fundamental optical mode fields.

.. figure:: ../../JOSAB_tutorial/josab_04-fields/AC_field_06.png
   :scale: 50 %

   Fundamental elastic mode fields.

.. raw:: latex

    \clearpage



Example 5 -- Intermodal Forward Brillouin scattering in a circular silica waveguide 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. .. literalinclude:: ../../JOSAB_tutorial/simo-josab-05-IFSBS-1umcylwg-SiO2.py
    :lines: 0-


.. figure:: ../../JOSAB_tutorial/josab_05-fields/EM_E_field_00.png
   :scale: 50 %

   Fundamental optical mode fields.

.. figure:: ../../JOSAB_tutorial/josab_05-fields/EM_E_field_01.png
   :scale: 50 %

   Second order optical mode fields.

.. figure:: ../../JOSAB_tutorial/josab_05-fields/AC_field_06.png
   :scale: 50 %

   Fundamental elastic mode fields.

.. raw:: latex

    \clearpage



Example 5 -- Intermodal Forward Brillouin scattering in a rectangular silicon waveguide 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. .. literalinclude:: ../../JOSAB_tutorial/simo-josab-06-IFSBS-450x200nmrectwg-Si.py
    :lines: 0-


.. figure:: ../../JOSAB_tutorial/josab_06-fields/EM_E_field_00.png
   :scale: 50 %

   Fundamental optical mode fields.

.. figure:: ../../JOSAB_tutorial/josab_06-fields/EM_E_field_01.png
   :scale: 50 %

   Second order optical mode fields.

.. figure:: ../../JOSAB_tutorial/josab_06-fields/AC_field_02.png
   :scale: 50 %

   Fundamental elastic mode fields.

.. raw:: latex

    \clearpage






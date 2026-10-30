
.. include:: numbatdefs.txt

.. _chap-waveguides-label:

******************************************
Materials, waveguides and meshing
******************************************

Now that we have a basic understanding of using |NUMBAT|, this chapter provides detailed information on how to specify a large of range materials and waveguide designs.

We will return to more advanced examples/tutorials in the next chapter.

Materials
----------------------

In order to calculate the modes of a structure we must specify the acoustic and optical properties of all constituent materials.

In |NUMBAT|, this data is read in from human-readable ``.json`` files, which are stored in the directory ``<NumBAT>/backend/material_data``.

These files not only provide the numerical values for optical and acoustic variables, but provide links to the origin of the data. Often they are taken from the literature and the naming convention allows users to select from different parameter values chosen by different authors for the same nominal material.

The intention of this arrangement is to create a library of materials that can serves as
standard reference data within the research community.
They also allow users to check the sensitivity of their results on particular parameters for a given material.

At present, the library contains the following materials:
  * Vacuum (or air)
      * ``Vacuum``
  * The chalcogenide glass Arsenic tri*sulfide
      * ``As2S3_2016_Smith``
      * ``As2S3_2017_Morrison``
      * ``As2S3_2021_Poulton``
  * Fused silica
      * ``SiO2_2013_Laude``
      * ``SiO2_2015_Van_Laer``
      * ``SiO2_2016_Smith``
      * ``SiO2_2021_Smith``
      * ``SiO2_smf28.json``
      * ``SiO2GeO2_smf28.json``
  * Silicon
      * ``Si_2012_Rakich``
      * ``Si_2013_Laude``
      * ``Si_2015_Van_Laer``
      * ``Si_2016_Smith``
      * ``Si_2021_Poulton``
      * ``Si_test_anisotropic``
  * Silicon nitride
      * ``Si3N4_2014_Wolff``
      * ``Si3N4_2021_Steel``
  * Gallium arsenide
      * ``GaAs_2016_Smith``
  * Germanium
      * ``Ge_cubic_2014_Wolff``
  * Lithium niobate
      * ``LiNbO3_2021_Steel``
      * ``LiNbO3aniso_2021_Steel``

Materials can easily be added to this library by copying any of these files as a template and
modifying the properties to suit. The ``Si_test_anisotropic`` file contains all the variables
that |NUMBAT| is setup to read. We ask that stable parameters (particularly those used
for published results) be added to the |NUMBAT| git repository using the same naming convention.


Waveguide Geometries
----------------------

|NUMBAT| encodes different waveguide structures through finite element meshes
constructed using the ``.geo`` language used by the open source tool ``Gmsh``.
Most users will find they can construct all waveguides of interest using the
existing templates. New templates use a ``.geo`` file and a subclass of
``UserGeometryBase``. They can remain in ``<NumBAT>/backend/msh`` using the
legacy file-based mechanism, or live together in an external source directory
that is added to Python's ``sys.path``. This procedure is described in detail
in :ref:`sec-user-defined-waveguides-label`.

All the builtin examples below are
constructed in the same fashion in a parallel ``builtin_meshes.py`` file and
can be used as models for your own designs.

The following figures give some examples of how material types and physical
dimensions are represented in the mesh geometries. In particular, for each
structure template, they identify the interpretation of the dimensional
parameters (``inc_a_x``, ``slab_b_y``, etc), material labels (``material_a``,
``material_b`` etc), and the grid refinement parameters (``lc_bkg``,
``lc_refine_1``, ``lc_refine_2``, etc).  The captions for each structure also
identify the mesh geometry template files in the directory
``<NumBAT>/backend/fortran/msh`` with filenames of the form
``<prefix>_msh_template.geo`` which define the structures and can give ideas
for developing new structure files.


The |NUMBAT| code for creating all these structures can be found  in ``<NumBAT>/docs/source/images/make_meshfigs.py``.

Single inclusion waveguides with surrounding medium
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These structures consist of a single medium inclusion (``mat_a``) with a background material (``mat_bkg``).
The dimensions are set with ``inc_a_x`` and ``inc_a_y``.

.. figure:: images/meshes/wg_rect-mesh-annotated.png
   :width: 16cm

.. figure:: images/meshes/wg_rect-refractive_index.png
   :width: 10cm

   Rectangular waveguide using shape ``rectangular`` (template ``oneincl_msh``).



.. figure:: images/meshes/wg_circ-mesh-annotated.png
   :width: 16cm

.. figure:: images/meshes/wg_circ-refractive_index.png
   :width: 10cm

   Elliptical waveguide using shape ``circular`` (template ``oneincl_msh``).




.. figure:: images/meshes/wg_triangular-mesh-annotated.png
   :width: 16cm

.. figure:: images/meshes/wg_triangular-refractive_index.png
   :width: 10cm

   Triangular waveguide using shape ``triangular``.



Double inclusion waveguides with surrounding medium
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These structures consist of a pair of waveguides with a single common background material.
The dimensions are set by ``inc_a_x/inc_a_y`` and  ``inc_b_x/inc_b_y``. They are separated
horizontally by ``two_inc_sep`` and the right waveguide has a vertical offset of ``y_off``.

.. figure::  images/meshes/wg_twoincl_rect-mesh-annotated.png
   :width: 16cm

.. figure:: images/meshes/wg_twoincl_rect-refractive_index.png
   :width: 10cm

   Coupled rectangular waveguides using shape ``rectangular`` (template ``twoincl_msh``).

.. figure::  images/meshes/wg_twoincl_circ-mesh-annotated.png
   :width: 16cm

.. figure:: images/meshes/wg_twoincl_circ-refractive_index.png
   :width: 10cm


   Coupled circular waveguides using shape ``circular`` (template ``twoincl_msh``).
   There appears to be a bug here!


Rib waveguides
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These structures consist of a rib on one or more substrate layers with zero to two  coating layers.

.. figure::  images/meshes/wg_rib-mesh-annotated.png
   :width: 16cm

.. figure:: images/meshes/wg_rib-refractive_index.png
   :width: 10cm

   A conventional rib waveguide using shape ``rib`` (template ``rib``).



.. figure:: images/meshes/wg_rib_coated-mesh-annotated.png
   :width: 16cm

.. figure:: images/meshes/wg_rib_coated-refractive_index.png
   :width: 10cm

   A coated rib waveguide using shape ``rib_coated`` (template ``rib_coated``).




.. figure:: images/meshes/wg_rib_double_coated-mesh-annotated.png
   :width: 16cm

.. figure:: images/meshes/wg_rib_double_coated-refractive_index.png
   :width: 10cm

   A rib waveguide on two substrates using shape ``rib_double_coated`` (template ``rib_double_coated``).


.. figure:: images/meshes/wg_slot-mesh-annotated.png
   :width: 16cm

Engineered rib  waveguides
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These are examples of more complex rib geometries. These are good examples to study in order to make new designs
using the user-specified waveguide and mesh mechanism.

.. figure:: images/meshes/wg_trapezoidal_rib-mesh-annotated.png
   :width: 16cm

   A trapezoidal rib structure using shape ``trapezoidal_rib``.

.. figure:: images/meshes/wg_pedestal-mesh-annotated.png
   :width: 16cm

   A supported pedestal structure using shape ``pedestal``.


Slot waveguides
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These slot waveguides can be used to enhance the horizontal component of the electric field in the low index region
by the 'slot' effect.

   A slot waveguide using shape ``slot`` (``material_a`` is low index) (template ``slot``).

.. .. figure:: images/meshes/slot_coated_wg-mesh-annotated.png
.. figure:: images/meshes/wg_slot_coated-mesh-annotated.png
   :width: 16cm

   A coated slot waveguide using shape ``slot_coated`` (``material_a`` is low index) (template ``slot_coated``).


Layered circular waveguides
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These waveguides consist of a set of concentric circular rings of a desired
number of layers in either a square or circular outer domain.
Note that ``inc_a_x`` specifies the innermost *diameter*.
The subsequent parameters ``inc_b_x``, ``inc_c_x``, etc specify the annular thickness of each successive layer.



.. figure:: images/meshes/wg_onion2-mesh-annotated.png
   :width: 16cm

.. figure:: images/meshes/wg_onion2-refractive_index.png
   :width: 8cm

   A two-layered concentric structure with background using shape ``onion2`` (template ``onion2``).

.. figure:: images/meshes/wg_onion3-mesh-annotated.png
   :width: 16cm

.. figure:: images/meshes/wg_onion3-refractive_index.png
   :width: 8cm

   A three-layered concentric structure with background using shape ``onion3`` (template ``onion3``).

.. figure:: images/meshes/wg_onion-mesh-annotated.png
   :width: 16cm

.. figure:: images/meshes/wg_onion-refractive_index.png
   :width: 8cm

   A many-layered concentric structure using shape ``onion`` (template ``onion``).


.. figure:: images/meshes/wg_circ_onion2-mesh-annotated.png
   :width: 16cm

.. figure:: images/meshes/wg_circ_onion2-refractive_index.png
   :width: 8cm

   A two-layered concentric structure with a circular outer boundary using shape ``circ_onion2`` (template ``circ_onion2``).


.. figure:: images/meshes/wg_circ_onion3-mesh-annotated.png
   :width: 16cm

.. figure:: images/meshes/wg_circ_onion3-refractive_index.png
   :width: 8cm

   A three-layered concentric structure with a circular outer boundary using shape ``circ_onion3`` (template ``circ_onion3``).


.. figure:: images/meshes/wg_circ_onionN-mesh-annotated.png
   :width: 16cm

.. figure:: images/meshes/wg_circ_onionN-refractive_index.png
   :width: 8cm

   A many-layered concentric structure with a circular outer boundary using shape ``circ_onion`` (template ``circ_onion``).




.. raw:: latex

    \clearpage

.. _sec-user-defined-waveguides-label:

User-defined waveguide geometries
-----------------------------------

Users may incorporate their own waveguide designs in a source-tree directory
outside the NumBAT checkout. The directory contains the Python implementation,
its Gmsh template files, and a manifest that tells NumBAT which waveguides to
load. The complete ``EllipticalRod`` example is in
``examples/user_waveguides`` and its use demonstrated in the tutorial
``sim-tut_21-user-waveguides.py``.

User waveguide directory
^^^^^^^^^^^^^^^^^^^^^^^^

The following is a minimal directory containing one waveguide:

::

    my_waveguides/
    ├── numbat_waveguides.json
    ├── elliptical_rod.py
    └── elliptical_rod_msh_template.geo

A directory may contain any number of Python modules and Gmsh templates. For a
larger collection, it is preferable to use a uniquely named Python package to
avoid collisions with other modules:

::

    my_waveguides/
    ├── numbat_waveguides.json
    └── my_lab_waveguides/
        ├── __init__.py
        ├── elliptical_rod.py
        ├── elliptical_rod_msh_template.geo
        ├── slot_rod.py
        └── slot_rod_msh_template.geo

In your main simulation script, register the directory before creating ``NumBATApp``:

::

    import numbat

    numbat.register_waveguide_directory("/path/to/my_waveguides")
    nbapp = numbat.NumBATApp()

Registration verifies the directory and its manifest, then makes the directory
available to Python's ordinary import system. NumBAT imports the active
waveguide classes when ``NumBATApp`` initializes. The directory can therefore
be kept in a separate Git repository and shared without modifying the NumBAT
checkout.

Waveguide manifest
^^^^^^^^^^^^^^^^^^

Every registered directory must contain a file named
``numbat_waveguides.json``. It contains a ``wguides`` list. Each list item
declares one waveguide:

::

    {
        "wguides": [
            {
                "module": "my_lab_waveguides.elliptical_rod",
                "wg_class": "EllipticalRod",
                "inc_shape": "elliptical_rod",
                "active": true
            },
            {
                "module": "my_lab_waveguides.slot_rod",
                "wg_class": "SlotRod",
                "inc_shape": "lab_slot_rod",
                "active": false
            }
        ]
    }

``module`` is the normal Python import path of the module containing the
class. For the minimal layout above, it would be ``elliptical_rod``. For the
package layout, it is ``my_lab_waveguides.elliptical_rod``.

``wg_class`` is the class name in that module. ``inc_shape`` is the name passed
as ``inc_shape`` to ``nbapp.make_structure`` and must be unique among all
active built-in, legacy, and registered waveguides. Only entries whose
``active`` value is ``true`` are loaded. This allows a single directory to
contain work-in-progress waveguides without making them available to a
simulation.

Python implementation
^^^^^^^^^^^^^^^^^^^^^

Each class named by the manifest must subclass ``UserGeometryBase``. Its module
must import the base class:

::

    from pathlib import Path

    from usermesh import UserGeometryBase

    class MyWaveguide(UserGeometryBase):
        def init_geometry(self):
            self.set_properties("mywaveguide", is_curvi=True)
            self.set_required_parameters(["inc_a_x", "inc_a_y"], 2)
            self.set_allowed_parameters(["lc_refine_1", "lc_refine_2"])

        def gmsh_template_path(self, p_dir_templates):
            return Path(__file__).with_name("mywaveguide_msh_template.geo")

        def apply_parameters(self):
            return [
                ("a1 = 20;", "a1 = %f;", "inc_a_x"),
                ("a1y = 10;", "a1y = %f;", "inc_a_y"),
            ]

The value passed to ``set_properties`` should normally match the manifest's
``inc_shape``. ``init_geometry`` defines the shape name, whether its boundaries
are curvilinear, and the required or optional parameters. Use
``set_parameter_help`` to document any template-specific parameters shown in
NumBAT's help output.

``apply_parameters`` returns string substitutions that replace identifiable
default expressions in the ``.geo`` file with values passed to
``make_structure``. Standard parameters include ``domain_x``, ``domain_y``,
``inc_a_x``, ``inc_a_y``, ``lc_bkg``, ``lc_refine_1``, and ``lc_refine_2``.
The full implementation in ``examples/user_waveguides/elliptical_rod.py`` also
shows substitutions for these common parameters.

For source-tree waveguides, override ``gmsh_template_path`` to return the
template alongside the Python module, as above. An optional
``draw_mpl_frame(self, ax, styles)`` method can add an outline to optical and
elastic mode-profile plots. The built-in ``Circular``, ``Rectangular``, and
``TwoIncl`` classes in ``builtin_meshes.py`` provide useful models.

Gmsh ``.geo`` template
^^^^^^^^^^^^^^^^^^^^^^

Create a Gmsh ``.geo`` file for each geometry. It is usually helpful to start
with ``triangular_msh_template.geo`` or another supplied template. If the
Python implementation above is used, the file must be named
``mywaveguide_msh_template.geo`` and be stored beside the module.

When designing the template, ensure the following:

  - Use appropriate default values for all significant dimensions. This makes
    it easier to inspect the general shape before NumBAT substitutes the exact
    values for a simulation.

  - Make every ``Line`` element unique. Do not create two lines joining the
    same pair of points. Such a design can look correct but produce a malformed
    mesh that fails during a NumBAT calculation.

  - Define the ``Line Loop`` elements for each region with consistent
    handedness. Anticlockwise is the natural convention. Include a minus sign
    when traversing a line in the reverse direction.

  - Group all regions with one material into a single ``Physical Surface``.

  - Group the outer boundary into ``Physical Line`` definitions.

  - Place the coordinate origin at a useful symmetry point near the expected
    fundamental field. This does not affect the calculation, but gives more
    natural axes in output plots.

The supplied mesh templates demonstrate these principles. The
``elliptical_rod_msh_template.geo`` example is a small two-material template
with an elliptical core and surrounding background.

Using a registered waveguide
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After registering the directory and creating ``NumBATApp``, select the
manifest's shape name as usual:

::

    structure = nbapp.make_structure(
        inc_shape="elliptical_rod",
        domain_x=3000,
        domain_y=2000,
        inc_a_x=800,
        inc_a_y=400,
        material_bkg=material_bkg,
        material_a=material_a,
        lc_bkg=0.1,
        lc_refine_1=2,
        lc_refine_2=3,
    )

Tutorial 21, ``examples/tutorials/sim-tut_21-user-waveguides.py``, demonstrates
the complete registration, structure creation, optical-mode, and elastic-mode
workflow.

Migrating a legacy user waveguide
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

NumBAT no longer loads user waveguides from ``backend/msh`` or supports
``wg_impl`` manifest entries. Move each legacy model into a directory outside
the NumBAT checkout and register that directory before creating
``NumBATApp``.

For a legacy model consisting of ``mywaveguides.py``,
``mywaveguide_msh_template.geo``, and an entry in
``backend/msh/user_waveguides.json``:

#. Create a directory such as ``my_waveguides`` outside the NumBAT checkout.
#. Move ``mywaveguides.py`` and ``mywaveguide_msh_template.geo`` into that
   directory.
#. Create ``my_waveguides/numbat_waveguides.json`` with a normal Python module
   entry:

   ::

       {
           "wguides": [
               {
                   "module": "mywaveguides",
                   "wg_class": "MyWaveguide",
                   "inc_shape": "mywaveguide",
                   "active": true
               }
           ]
       }

#. Ensure the geometry class resolves its Gmsh template alongside its module:

   ::

       from pathlib import Path

       def gmsh_template_path(self, p_dir_templates):
           return Path(__file__).with_name("mywaveguide_msh_template.geo")

#. Register the directory before constructing the application:

   ::

       import numbat

       numbat.register_waveguide_directory("/path/to/my_waveguides")
       nbapp = numbat.NumBATApp()

Existing serialized geometries created with the former ``wg_impl`` mechanism
cannot be loaded after this migration. Recreate them from the registered
waveguide model.

Designing and implementing a waveguide should not be daunting, but the first
template can be confusing. Use the supplied examples as models and review the
generated mesh before running a large calculation.


Mesh parameters
----------------------

The parameters ``lc_bkg``, ``lc_refine_1``, ``lc_refine_2`` labelled in the
above figures control the fineness of the FEM mesh and are set when
constructing the waveguide, as discussed in the next chapter.  The first
parameter ``lc_bkg`` sets the reference background mesh size, typically as a
fraction of the length of the outer boundary edge.  A larger ``lc_bkg`` yields
a coarser mesh. Reasonable starting values are ``lc_bkg=0.1`` (10 mesh points
on the outer boundary) to ``lc_bkg=0.05`` (20 mesh points on the outer
boundary).

As well as setting the overall mesh scale with ``lc_bkg``, one can also refine
the mesh near interfaces and near select points in the domain, as may be
observed in the figures in the previous section.  This helps to increase the
mesh resolution in regions where there the electromagnetic and acoustic fields
are likely to be strong and/or rapidly varying.  This is achieved using the
``lc_refine_n`` parameters as follows. At the interface between materials, the
mesh is refined to have characteristic length ``lc_bkg/lc_refine_1``, therefore
a *larger* ``lc_refine_1`` gives a *finer* mesh by a factor of ``lc_refine_1``
at these interfaces.  The meshing program ``Gmsh`` automatically adjusts the
mesh size to smoothly transition from a point that has one mesh parameter to
points that have other meshing parameters. The mesh is typically also refined
in the vicinity of important regions, such as in the center of a waveguide,
which is done with ``lc_refine_2``, which analogously to ``lc_refine_1``,
refines the mesh size at these points as ``lc_bkg/lc_refine_2``.

For more complicated structures, there are additional ``lc_refine_<n>``
parameters.  To see their exact function, look for these expressions in the
particular .geo file.

Choosing appropriate values of ``lc_bkg``, ``lc_refine_1``, ``lc_refine_2`` is
crucial for |NUMBAT| to give accurate results. The appropriate values depend
strongly on the type of structure being studied, and so we strongly recommended
carrying out a convergence test before delving into new structures (see
Tutorial 5 for an example) starting from similar parameters as used in the
tutorial simulations.

As will as giving low accuracy, a structure with too coarse a mesh is often the
cause of the eigensolver failing to converge in which case |NUMBAT| will
terminate with an error.  If you encounter such an error, try the calculation
again with a slightly smaller value for ``lc_bkg``, or slightly higher values
for the ``lc_refine_n`` parameters.

On the other hand, it is wise to begin with relatively coarse meshes. It will
be apparent that the number of  elements scales roughly *quadratically* with
the ``lc_refine`` parameters and so the run-time increases rapidly as the mesh
becomes finer.  For each problem,  some initial experimentation to identify a
mesh resolution that gives reasonable convergence in acceptable simulation is
usually worthwhile.


Viewing the mesh
----------------------

When |NUMBAT| constructs a waveguide, the template ``geo`` file is converted to
a concrete instantiation with the ``lc_refine`` and geometric parameters
adjusted to the requested values. This file is then converted into a ``gmsh``
``.msh`` file. When exploring new structures and their convergence behaviour,
it is a very good idea to view the generated mesh frequently.

You can examine the resolution of your mesh by calling the
``plot_mesh(<prefix>)`` or ``check_mesh()`` methods on a waveguide
``Structure`` object.  The first of these functions saves a pair of images of
the mesh to a ``<prefix>-mesh-annotated.png`` file in the local directory which can be
viewed with your preferred image viewer; the second opens the mesh in a
``gmsh`` window (see Tutorial 1 above).

In addition, the ``.msh`` file generated by |NUMBAT| in any calculation is
stored in ``<NumBAT>/backend/fortran/msh/build`` and can be viewed by running
the command ::

    gmsh <msh_filename>.msh

In some error situations, |NUMBAT| will explicitly suggest viewing the mesh and
will print out the required command to do so.


.. .. raw:: latex

    \clearpage

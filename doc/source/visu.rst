Visualization utilities
=======================

The :mod:`cromosim.visu` module provides two persistent renderers for
microscopic simulations. They share the same simulation layers (people,
contacts, paths, sensors, actual velocities and desired velocities) while using
different graphics backends.

Matplotlib renderer
-------------------

.. currentmodule:: cromosim.visu

.. autoclass:: matplotlibVisualization
   :members:
   :undoc-members:
   :show-inheritance:

A renderer is created once and then updated at every visualization step. This
avoids rebuilding the whole Matplotlib figure during the simulation.

.. code-block:: python

   from cromosim.visu import matplotlibVisualization

   visu = matplotlibVisualization(
       domain,
       cmap="viridis",
       velocity_reference=1.0,
       velocity_arrow_length=1.5,
   )

   visu.update(
       people=people,
       contacts=contacts,
       sensors=sensors,
       time=t,
       plot_people=True,
       plot_contacts=True,
       plot_velocities=True,
       plot_desired_velocities=True,
       plot_paths=True,
       plot_sensors=True,
   )

Fixed velocity-arrow scale
~~~~~~~~~~~~~~~~~~~~~~~~~~

Both renderers use the same fixed conversion from a physical velocity to an
arrow length. If ``velocity_reference`` is :math:`v_{ref}` and
``velocity_arrow_length`` is :math:`L_{ref}`, a velocity vector :math:`v` is
displayed with length

.. math::

   L(v) = \frac{L_{ref}}{v_{ref}}\,\lVert v \rVert.

Consequently, the scale does not depend on the current frame or on the largest
velocity present in that frame. Actual velocities ``U`` and desired velocities
``Vd`` therefore remain directly comparable throughout a simulation.

PyVista renderer
----------------

.. autoclass:: pyvistaVisualization
   :members:
   :undoc-members:
   :show-inheritance:

The PyVista renderer accepts the same velocity scaling parameters. It can also
record frames directly to a movie when ``movie_path`` is supplied.

Legacy plotting functions
-------------------------

The historical plotting helpers remain available for compatibility:

.. autofunction:: plot_people

.. autofunction:: plot_sensors

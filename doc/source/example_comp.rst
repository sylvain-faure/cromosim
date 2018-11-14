Discover the compartment model
===============================

The models presented in this chapter are fundamentally different from all other models in this book.
This class of model does not aim at computing positions and velocities,
but rather at describing the evolution of total amounts of people in specified
compartments. The compartments shall be identified to nodes of a network, those
nodes are connected by edges, which corresponds to paths between nodes.

Reference : [MF2018]_ Chapter 6.

A compartment model example can be find in the directory

.. code-block:: python3

   cromosim/examples/compartments/

and can be run with:

.. code-block:: python3

   python3 compartments.py --json input.json

Results:

.. raw:: html

   <br>
   <br>
   <video style="display:block; margin: 0 auto;" controls poster="_static/compartments.png" width="50%" height="auto">
      <source src="_static/compartments.mp4" />
      <source src="_static/compartments.webm" />
      <source src="_static/compartments.theora.ogv" />
   </video>
   <br>
   <div align=center>
     <i>Compartment model : evacuation</i>
   </div>
   <br>
   <br>

Code:

.. code-block:: python3

   cromosim/examples/compartments/compartments.py


..  literalinclude:: ../../examples/compartments/compartments.py
    :lines: 1-

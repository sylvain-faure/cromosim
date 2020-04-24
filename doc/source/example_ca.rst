Try a cellular automata
========================

Cellular Automata models are microscopic in the sense that agents are
represented individually, but they are also natively discrete in both
space and time. In space, agents are identified to particles located in
cells of a fixed cartesian grid, with an exclusion rule (one particle at
most in each cell). In time, the evolution consists of a succession of steps,
i.e. positions of particles are updated one after the other, in a stochastic
manner.

Reference : [MF2018]_ Chapter 5.

A cellular automata example can be found in the directory

.. code-block:: python3

   cromosim/examples/cellular_automata/

and can be run with:

.. code-block:: python3

   python3 cellular_automata.py --json input.json

Results:

.. raw:: html

   <br>
   <br>
   <video style="display:block; margin: 0 auto;" controls poster="_static/cellular_automata.png" width="50%" height="auto">
      <source src="_static/cellular_automata.mp4" />
   </video>
   <br>
   <div align=center>
     <i>Cellular automata : simple example, an evacuation</i>
   </div>
   <br>
   <br>
   <img src="_static/cellular_automata_paths.png" alt="cellular automata : paths of all individuals" style="margin:0px auto;display:block" width="50%" height="auto"/>
   <br>
   <div align=center>
     <i>The lines represent the paths of all the individuals</i>
   </div>
   <br>
   <br>

Code:

.. code-block:: python3

   cromosim/examples/cellular_automata/cellular_automata.py

..  literalinclude:: ../../examples/cellular_automata/cellular_automata.py
    :lines: 1-

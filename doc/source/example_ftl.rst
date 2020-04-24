Understand one-dimensional microscopic models : the Follow-the-Leader model
============================================================================

The Follow-the-Leader (FTL) model presented here is a one-dimensional
microscopic model. Pedestrians are assumed to walk on a line toward a
common direction and the instantaneous velocity of an individual is a
function of the distance to the next individual.

We consider :math:`N+1` individuals walking on a straight line in a
non-periodic case or in a periodic case which means that the
individual :math:`N+1` is identified to the individual :math:`1`.

We can also consider an inertial version of the Follow-the-Leader
model (an order :math:`2` version) which can be seen as a one-dimensional version of the Social
Force model with asymmetric forcing terms to take into account that
pedestrians are influenced by their immediate neighbor in front of
them.

Reference : [MF2018]_ Chapter 2.

Thanks to several json files given as inputs of the following script,
we can combine periodic or non-periodic domain with inertial or
non-inertial FTL model. Examples can be find in the directory

.. code-block:: python3

   cromosim/examples/follow_the_leader/

To run the **non-inertial FTL model in non-periodic case**:

.. code-block:: python3

   python3 follow_the_leader.py --json input_ftl_order1.json

Results:

.. raw:: html

   <br>
   <br>
   <video style="display:block; margin: 0 auto;" controls poster="_static/ftl_order1.png" width="60%" height="auto">
      <source src="_static/ftl_order1.mp4" />
   </video>
   <br>

To run the **non-inertial FTL model in periodic case**:

   .. code-block:: python3

      python3 follow_the_leader.py --json input_ftl_order1_periodic.json

Results:

.. raw:: html

   <br>
   <br>
   <video style="display:block; margin: 0 auto;" controls poster="_static/ftl_order1_periodic.png" width="60%" height="auto">
      <source src="_static/ftl_order1_periodic.mp4" />
   </video>
   <br>


To run the **inertial FTL model in non-periodic case**:

.. code-block:: python3

   python3 follow_the_leader.py --json input_ftl_order2.json

Results:

.. raw:: html

      <br>
      <br>
      <video style="display:block; margin: 0 auto;" controls poster="_static/ftl_order2.png" width="60%" height="auto">
         <source src="_static/ftl_order2.mp4" />
      </video>
      <br>

To run the **inertial FTL model in periodic case**:

.. code-block:: python3

         python3 follow_the_leader.py --json input_ftl_order2_periodic.json

Results:

.. raw:: html

      <br>
      <br>
      <video style="display:block; margin: 0 auto;" controls poster="_static/ftl_order2_periodic.png" width="60%" height="auto">
         <source src="_static/ftl_order2_periodic.mp4" />
      </video>
      <br>

Code:

.. code-block:: python3

   cromosim/examples/follow_the_leader/follow_the_leader.py


..  literalinclude:: ../../examples/follow_the_leader/follow_the_leader.py
    :lines: 1-

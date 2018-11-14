Make a microscopic simulation
===============================

In the two following models, individuals are seen as spherical particles. They
move at a desired velocity to reach a goal (typically a door). Doing nothing
more would lead to many overlaps between individuals hence the two families of
models below can prevent these overlaps through two different approaches. For
the Social Force models some forces are added to act against overlaps, and for
the Granular model the velocities are projected into a permissible velocity
space which ensures the absence of overlaps.

The domain used for both models is obtained in using the following background
image:

.. raw:: html

   <br>
   <img src="_static/event.png" alt="Computational domain" style="margin:0px auto;display:block" width="50%" height="auto"/>
   <br>



Social force model
-------------------

The Social Force model introduced in the 90’s. Pedestrians are identified
with inertial particles, submitted to a forcing term which implements the
individual tendencies, and extra forces which account for interactions with
other pedestrians (typically the tendency to preserve a certain distance
with neighbors).

Reference : [MF2018]_ Chapter 3.

An example can be find in the directory

.. code-block:: python3

   cromosim/examples/micro/social

and can be launched with

.. code-block:: python3

   python3 micro_social.py --json input.json

.. raw:: html

  <br>
  <br>
  <video style="display:block; margin: 0 auto;" controls poster="_static/micro_social.png" width="50%" height="auto">
    <source src="_static/micro_social.mp4" />
    <source src="_static/micro_social.webm" />
    <source src="_static/micro_social.theora.ogv" />
  </video>
  <br>
  <div align=center>
    <i>Social force model : an evacuation in complex geometry</i>
  </div>
  <br>
  <br>
  <img src="_static/micro_social_sensor_0_5929.png" alt="Results of sensor 1" style="margin:0px auto;display:block" width="50%" height="auto"/>
  <br>
  <div align=center>
    <i>Sensor 1</i>
  </div>
  <br>
  <img src="_static/micro_social_sensor_1_5929.png" alt="Results of sensor 2" style="margin:0px auto;display:block" width="50%" height="auto"/>
  <br>
  <div align=center>
    <i>Sensor 2</i>
  </div>
  <br>
  <br>

..  literalinclude:: ../../examples/micro/social/micro_social.py
    :lines: 1-


Granular model
---------------

The Granular model comes from crowd motion models of the granular type : each
individual is identified to a hard disk of a prescribed size, subject to a
non-overlapping constraint with their neighbors. The approach relies on a
desired velocity for each individual (the velocity they would take if they were
alone), and the global velocity field shall be defined as the closest to the
desired one among all those feasible fields (fields which do not lead to
overlapping of disks).

Reference : [MF2018]_ Chapter 4.

An example can be find in the directory

.. code-block:: python3

   cromosim/examples/micro/granular

and can be launched with

.. code-block:: python3

   python3 micro_granular.py --json input.json


.. raw:: html

   <br>
   <br>
   <video style="display:block; margin: 0 auto;" controls poster="_static/micro_granular.png" width="50%" height="auto">
     <source src="_static/micro_granular.mp4" />
     <source src="_static/micro_granular.webm" />
     <source src="_static/micro_granular.theora.ogv" />
   </video>
   <br>
   <div align=center>
     <i>Granular model : an evacuation in complex geometry</i>
   </div>
   <br>
   <br>
   <img src="_static/micro_granular_sensor_0_4097.png" alt="Results of sensor 1" style="margin:0px auto;display:block" width="50%" height="auto"/>
   <br>
   <div align=center>
     <i>Sensor 1</i>
   </div>
   <br>
   <img src="_static/micro_granular_sensor_1_4097.png" alt="Results of sensor 2" style="margin:0px auto;display:block" width="50%" height="auto"/>
   <br>
   <div align=center>
     <i>Sensor 2</i>
   </div>
   <br>
   <br>


..  literalinclude:: ../../examples/micro/granular/micro_granular.py
    :lines: 1-

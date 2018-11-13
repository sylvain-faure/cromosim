Make a microscopic simulation
===============================

.. raw:: html

   <img src="_static/event.png" alt="Computational domain" style="margin:0px auto;display:block" width="50%" height="auto"/>
   <br>



Social force model
-------------------

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

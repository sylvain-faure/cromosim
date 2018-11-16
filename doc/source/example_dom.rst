Build a domain
===============

A domain object contains:

* a geometry defined by a background image (completely white by
  default) and Matplotlib shapes that can be added (lines, circles, ellipses, rectangles or polygons). These elements are walls if they are black. To represent a door serving as an exit to the individuals that will be injected into the domain, a red line will be used.
* a 2D array of the size of the image representing the distance to the nearest wall (black color). This means that each pixel in the image has its distance to the nearest wall.
* a 2D array of the size of the image representing the distance to the nearest exit (red color). The opposite of the gradient of this distance represents the desired direction for an individual desiring to reach this outcome. This desired direction is also stored in an array.

Reference : [MF2018]_ Chapter 8.

Manually
---------

An example can be found in the directory

.. code-block:: python3

   cromosim/examples/domain/

and can be launched with

.. code-block:: python3

   python3 domain_manually_computed.py

Manually means that, from a completely white background image, the
domain will be defined and drawn using Matplotlib shapes (Line2D,...). Hence, the
first lines of the script build the boundary of the room (wall in black color) and
the doors (in red color):

..  literalinclude:: ../../examples/domain/domain_manually_computed.py
    :lines: 1-46

which gives the following domain:

.. raw:: html

   <img src="_static/room_400x300_domain.png" alt="A domain manually defined using lines" style="margin:0px auto;display:block" width="50%" height="auto"/>
   <br>


Then, the desired velocity and the wall distance can also be computed manually
and added to the domain object:

..  literalinclude:: ../../examples/domain/domain_manually_computed.py
    :lines: 47-

.. raw:: html

   <img src="_static/domain_manually_wall_distance.png" alt="Domain : wall distance" style="margin:0px auto;display:block" width="60%" height="auto"/>
   <div align=center>
     <i>Wall distance</i>
   </div>
   <br>
   <img src="_static/domain_manually_desired_velocity.png" alt="Domain : desired velocity" style="margin:0px auto;display:block" width="60%" height="auto"/>
   <div align=center>
     <i>Desired velocity (to reach the door)</i>
   </div>

Automatically
--------------

Starting from a background image drawn using three colors (red for the
doors, black for the walls and white):

.. raw:: html

   <img src="_static/room_400x300.png" alt="A domain manually defined using lines" style="margin:0px auto;display:block" width="50%" height="auto"/>
   <br>


It is still possible to add Matplotlib shapes:

..  literalinclude:: ../../examples/domain/domain_auto_computed.py
    :lines: 1-44

and then compute automatically the desired velocity and the wall distance thanks to a fast-marching method.

..  literalinclude:: ../../examples/domain/domain_auto_computed.py
    :lines: 45-

.. raw:: html

   <img src="_static/domain_with_shapes.png" alt="Domain with shapes" style="margin:0px auto;display:block" width="60%" height="auto"/>
   <div align=center>
    <i>Shapes have been add to the background image</i>
   </div>
   <br>
   <img src="_static/domain_wall_distance.png" alt="Domain : wall distance" style="margin:0px auto;display:block" width="60%" height="auto"/>
   <div align=center>
    <i>Wall distance</i>
   </div>
   <br>
   <img src="_static/domain_desired_velocity.png" alt="Domain : desired velocity" style="margin:0px auto;display:block" width="60%" height="auto"/>
   <div align=center>
    <i>Desired velocity (to reach the door)</i>
   </div>

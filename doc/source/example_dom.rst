Build a domain
===============

A `Domain` object contains:

* a geometry defined by a background image (completely white by
  default) and `Matplotlib` shapes that can be added (lines, circles, ellipses,
  rectangles or polygons). Depending on the color chosen for a shape, its
  meaning changes. For example, you can decide that a red line corresponds to
  one door and a green line represents another.
  To build a computational domain, we therefore use a background image
  (drawing, map of the site, satellite photo, etc.) or we leave an entirely
  white background, then we optionally add colored shapes, and finally we
  define a color code. In general the walls are coded in black. All colors
  are taken in rgb format: `[r, g, b]` with `r`, `g` and `b` integers between
  `0` and `255`.
* a 2D array of the size of the image representing the distance to the nearest
  wall (coded by all the colors defined in `wall_colors`). This means that each
  pixel in the image has its distance to the nearest wall.
* A set of possible destinations attributable to people. A destination is
  defined by a `Destination` object which contains all the elements useful for
  its construction: the color codes representing the destination on the domain
  (red for a door represented by a red line for example), exclusion zones also
  defined via colors (black for an impassable wall for example), the names of
  the next destination and the next domain (we can chain several destinations
  to go up one floor: one `Domain` object per floor, one `Destination` object
  leading to the stairs, another to climb the stairs, etc ...).
  Each destination will contain a 2D array of the size of the image containing
  the distance to the objective (the aimed door for example). The opposite of
  the gradient of this distance represents the desired direction (or named
  desired velocity) for an individual desiring to reach this outcome. This
  desired direction is also stored in an array.

Reference : [MF2018]_ Chapter 8.

All the following examples can be found in the directory

.. code-block:: python

   cromosim/examples/domain/


The simplest example
--------------------

.. code-block:: python

   python domain_room.py

This domain uses an image as background and two forms: a red line used to
define a door and a black disc to add an obstacle a few meters before the door.

.. raw:: html

  <br>
  <table style="width:100%">
    <tr>
      <td>
        <center>
          <em><small>The image that serves as the background</small></em>
          <img src="_static/room.png" alt="The image that serves as the background" style="margin:0px" width="80%" height="auto"/>
        </center>
      </td>
      <td>
        <center>
          <em><small>The domain: background, black disk and red line added</small></em>
          <img src="_static/room_domain.png" alt="The domain : background + a shape (black circle)" style="margin:0px auto;display:block" width="80%" height="auto"/>
        </center>
      </td>
    </tr>
    <tr>
      <td>
        <center>
          <em><small>Distance to walls and its gradient</small></em>
          <img src="_static/room_wall_distance.png" alt="Distance to walls and its gradient" style="margin:0px auto;display:block" width="80%" height="auto"/>
        </center>
      </td>
      <td>
        <center>
          <em><small>Distance to destination (red door) and desired velocity</small></em>
          <img src="_static/room_desired_velocity.png" alt="Distance to destination (red door) and the desired velocity" style="margin:0px auto;display:block" width="80%" height="auto"/>
        </center>
      </td>
    </tr>
  </table>
  <br>

..  literalinclude:: ../../examples/domain/domain_room.py
  :caption: domain_room.py
  :linenos:
  :language: python
  :lines: 18-


A circular domain
-----------------

.. code-block:: python

  python domain_stadium.py

This domain contains black walls `rgb=[0,0,0]`, a red line `rgb=[255,0,0]` and
a green line `rgb=[0,255,0]`, all included in the given background.

If we calculate the desired velocity leading to the red line by solving the
Eikonal equation using a Fast-Marching method, then we get arrows directed
towards this line at any point on the track. For example, a person just to the
right of the line will go to the left and a person just to the left will go to
the right, so people will not turn around the track.

To correct this problem, we add a fictitious wall (green line) used only during
the Fast-Marching method: thanks to it a person located to the left of the
red line will go to the left and go around the track.

Finally if we want people to do several laps, we impose a desired velocity
equal to `(-1.0)` on the green and red pixels using the
`desired_velocity_from_color` option of the `Destination` object.

.. raw:: html

  <br>
  <table style="width:100%">
    <tr>
      <td>
        <center>
          <em><small>The image that serves as the background</small></em>
          <img src="_static/stadium.png" alt="The image that serves as the background" style="margin:0px" width="80%" height="auto"/>
        </center>
      </td>
      <td>
        <center>
          <em><small>The domain, identical to the background image</small></em>
          <img src="_static/stadium_domain.png" alt="The domain, identical to the background image" style="margin:0px auto;display:block" width="92%" height="auto"/>
        </center>
      </td>
    </tr>
    <tr>
      <td>
        <center>
          <em><small>Distance to walls and its gradient</small></em>
          <img src="_static/stadium_wall_distance.png" alt="Distance to walls and its gradient" style="margin:0px auto;display:block" width="80%" height="auto"/>
        </center>
      </td>
      <td>
        <center>
          <em><small>Distance to destination (red door) and desired velocity</small></em>
          <img src="_static/stadium_desired_velocity.png" alt="Distance to destination (red door) and the desired velocity" style="margin:0px auto;display:block" width="92%" height="auto"/>
        </center>
      </td>
    </tr>
  </table>
  <br>

..  literalinclude:: ../../examples/domain/domain_stadium.py
  :caption: domain_stadium.py
  :linenos:
  :language: python
  :lines: 16-


A domain with several destinations
----------------------------------

.. code-block:: python

  python domain_shibuya_crossing.py


Shibuya Crossing is a popular scramble crossing in Shibuya, Tokyo, Japan.
Here some walls have been added (magenta walls) to channel pedestrians on the pedestrian
crossings and we have defined five destinations leading to five disks of
different colors (blue, brown, cyan, green and pink).

.. raw:: html

  <br>
  <table style="width:100%">
    <tr>
      <td>
        <center>
          <em><small>The domain, identical to the background image</small></em>
          <img src="_static/shibuya_crossing_domain.png" alt="The domain, identical to the background image" style="margin:0px auto;display:block" width="80%" height="auto"/>
        </center>
      </td>
      <td>
        <center>
          <em><small>Distance to walls and its gradient</small></em>
          <img src="_static/shibuya_crossing_wall_distance.png" alt="Distance to walls and its gradient" style="margin:0px auto;display:block" width="80%" height="auto"/>
        </center>
      </td>
    </tr>
    <tr>
      <td>
        <center>
          <em><small>Distance to destination (blue disk) and desired velocity</small></em>
          <img src="_static/shibuya_crossing_blue.png" alt="Distance to destination (blue disk) and the desired velocity" style="margin:0px auto;display:block" width="80%" height="auto"/>
        </center>
      </td>
      <td>
        <center>
          <em><small>Distance to destination (brown disk) and desired velocity</small></em>
          <img src="_static/shibuya_crossing_brown.png" alt="Distance to destination (brown disk) and the desired velocity" style="margin:0px auto;display:block" width="80%" height="auto"/>
        </center>
      </td>
    </tr>
    <tr>
      <td>
        <center>
          <em><small>Distance to destination (cyan disk) and desired velocity</small></em>
          <img src="_static/shibuya_crossing_cyan.png" alt="Distance to destination (cyan disk) and the desired velocity" style="margin:0px auto;display:block" width="80%" height="auto"/>
        </center>
      </td>
      <td>
        <center>
          <em><small>Distance to destination (green disk) and desired velocity</small></em>
          <img src="_static/shibuya_crossing_green.png" alt="Distance to destination (green disk) and the desired velocity" style="margin:0px auto;display:block" width="80%" height="auto"/>
        </center>
      </td>
    </tr>
    <tr>
      <td>
        <center>
          <em><small>Distance to destination (pink disk) and desired velocity</small></em>
          <img src="_static/shibuya_crossing_pink.png" alt="Distance to destination (pink disk) and the desired velocity" style="margin:0px auto;display:block" width="80%" height="auto"/>
        </center>
      </td>
      <td>
      </td>
    </tr>
  </table>
  <br>


..  literalinclude:: ../../examples/domain/domain_shibuya_crossing.py
  :caption: domain_shibuya_crossing.py
  :linenos:
  :language: python
  :lines: 15-


A domain build from a json file
-------------------------------

.. code-block:: python

  python domain_from_json.py --json input_room.json
  python domain_from_json.py --json input_stadium.json
  python domain_from_json.py --json input_shibuya_crossing.json

JSON (JavaScript Object Notation) is a lightweight data-interchange format
which can be used to describe one or more domains. This means in our case
that with a single python script taking as argument a json file we will be
able to build several domains.

..  literalinclude:: ../../examples/domain/domain_from_json.py
  :caption: domain_from_json.py
  :linenos:
  :language: python
  :lines: 20-

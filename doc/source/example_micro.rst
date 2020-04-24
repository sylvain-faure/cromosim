Make a microscopic simulation
===============================

In the two following models, individuals are seen as spherical particles. They
move at a desired velocity to reach a goal (typically a door). Doing nothing more
would lead to many overlaps between individuals. Hence the two families of
models below can prevent these overlaps through two different approaches. For
the `Social Force models`, and, for some forces are added to act against overlaps, and for
the `Granular model` the velocities are projected into a permissible velocity
space which ensures the absence of overlaps.


Social force model
-------------------

The `Social Force model` has been introduced in the 90’s. Pedestrians are identified
with inertial particles submitted to a forcing term which implements the
individual tendencies and extra forces which account for interactions with
other pedestrians (typically the tendency to preserve a certain distance
with neighbors).

Reference : [MF2018]_ Chapter 3.

Some examples can be found in the directory

.. code-block:: python

   cromosim/examples/micro/social

and can be launched with

.. code-block:: python

   python micro_social.py --json input_room.json
   python micro_social.py --json input_event.json
   python micro_social.py --json input_stadium.json
   python micro_social.py --json input_shibuya_crossing.json
   python micro_social.py --json input_stairs.json

.. raw:: html

  <table style="width:100%">
    <tr>
      <th>
        <center>
          <em><small>Evacuation of a room, visualization of trajectories</small></em>
        </center>
      </th>
    </tr>
    <tr>
      <td>
        <center>
          <video style="display:block; margin: 0 auto;" controls poster="_static/micro_social_room_paths.png" width="40%" height="auto">
            <source src="_static/micro_social_room_paths.mp4" />
          </video>
        </center>
      </td>
    </tr>
  </table>
  <br>
  <table style="width:100%">
    <tr>
      <th>
        <center>
          <em><small>Circular track around a stadium</small></em>
        </center>
      </th>
    </tr>
    <tr>
      <td>
        <center>
          <video style="display:block; margin: 0 auto;" controls poster="_static/micro_social_stadium.png" width="40%" height="auto">
            <source src="_static/micro_social_stadium.mp4" />
          </video>
        </center>
      </td>
    </tr>
  </table>
  <br>
  <table style="width:100%">
    <tr>
      <th>
        <center>
          <em><small>Evacuation of an exhibition hall: two groups of people with the same destination</small></em>
        </center>
      </th>
      <th>
        <center>
          <em><small>Evacuation of an exhibition hall: sensors (green lines) results</small></em>
        </center>
      </th>
    </tr>
    <tr>
      <td>
        <center>
          <video style="display:block; margin: 0 auto;" controls poster="_static/micro_social_event.png" width="40%" height="auto">
            <source src="_static/micro_social_event.mp4" />
          </video>
        </center>
      </td>
      <td>
        <center>
          <img src="_static/micro_social_event_sensor.png" alt="Evacuation of an exhibition hall : Sensors (green lines) results" style="margin:0px auto;display:block" width="100%" height="auto"/>
        </center>
      </td>
    </tr>
  </table>
  <br>
  <table style="width:100%">
    <tr>
      <th>
        <center>
          <em><small>Shibuya crossing (Japan) : five groups of people and five different destinations</small></em>
        </center>
      </th>
    </tr>
    <tr>
      <td>
        <center>
          <video style="display:block; margin: 0 auto;" controls poster="_static/micro_social_shibuya_crossing.png" width="40%" height="auto">
            <source src="_static/micro_social_shibuya_crossing.mp4" />
          </video>
        </center>
      </td>
    </tr>
  </table>
  <br>
  <table style="width:100%">
    <tr>
      <th>
        <center>
          <em><small>Two floors of a building: a group of people goes up and another goes down
          <br>Floor 1 on the left and 0 on the right</small></em>
        </center>
      </th>
    </tr>
    <tr>
      <td>
        <center>
          <video style="display:block; margin: 0 auto;" controls poster="_static/micro_social_stairs_0_1.png" width="100%" height="auto">
            <source src="_static/micro_social_stairs_0_1.mp4" />
          </video>
        </center>
      </td>
    </tr>
  </table>

..  literalinclude:: ../../examples/micro/social/micro_social.py
  :caption: micro_social.py
  :linenos:
  :language: python
  :lines: 1-


Granular model
---------------

The `Granular model` comes from crowd motion models of the granular type : each
individual is identified to a hard disk of a prescribed size, subject to a
non-overlapping constraint with their neighbors. The approach relies on a
desired velocity for each individual (the velocity they would take if they were
alone), and the global velocity field shall be defined as the closest to the
desired one among all those feasible fields (fields which do not lead to
overlapping of disks).

Reference : [MF2018]_ Chapter 4.

An example can be find in the directory

.. code-block:: python

   cromosim/examples/micro/granular

and can be launched with


.. code-block:: python

   python micro_granular.py --json input_room.json
   python micro_granular.py --json input_event.json
   python micro_granular.py --json input_stadium.json
   python micro_granular.py --json input_shibuya_crossing.json
   python micro_granular.py --json input_stairs.json

.. raw:: html

  <table style="width:100%">
    <tr>
      <th>
        <center>
          <em><small>Evacuation of a room, visualization of trajectories</small></em>
        </center>
      </th>
    </tr>
    <tr>
      <td>
        <center>
          <video style="display:block; margin: 0 auto;" controls poster="_static/micro_granular_room_paths.png" width="40%" height="auto">
            <source src="_static/micro_granular_room_paths.mp4" />
          </video>
        </center>
      </td>
    </tr>
  </table>
  <br>
  <table style="width:100%">
    <tr>
      <th>
        <center>
          <em><small>Circular track around a stadium</small></em>
        </center>
      </th>
    </tr>
    <tr>
      <td>
        <center>
          <video style="display:block; margin: 0 auto;" controls poster="_static/micro_granular_stadium.png" width="40%" height="auto">
            <source src="_static/micro_granular_stadium.mp4" />
          </video>
        </center>
      </td>
    </tr>
  </table>
  <br>
  <table style="width:100%">
    <tr>
      <th>
        <center>
          <em><small>Evacuation of an exhibition hall: two groups of people with the same destination</small></em>
        </center>
      </th>
      <th>
        <center>
          <em><small>Evacuation of an exhibition hall: sensors (green lines) results</small></em>
        </center>
      </th>
    </tr>
    <tr>
      <td>
        <center>
          <video style="display:block; margin: 0 auto;" controls poster="_static/micro_granular_event.png" width="40%" height="auto">
            <source src="_static/micro_granular_event.mp4" />
          </video>
        </center>
      </td>
      <td>
        <center>
          <img src="_static/micro_granular_event_sensor.png" alt="Evacuation of an exhibition hall : Sensors (green lines) results" style="margin:0px auto;display:block" width="100%" height="auto"/>
        </center>
      </td>
    </tr>
  </table>
  <br>
  <table style="width:100%">
    <tr>
      <th>
        <center>
          <em><small>Shibuya crossing (Japan) : five groups of people and five different destinations</small></em>
        </center>
      </th>
    </tr>
    <tr>
      <td>
        <center>
          <video style="display:block; margin: 0 auto;" controls poster="_static/micro_granular_shibuya_crossing.png" width="40%" height="auto">
            <source src="_static/micro_granular_shibuya_crossing.mp4" />
          </video>
        </center>
      </td>
    </tr>
  </table>
  <br>
  <table style="width:100%">
    <tr>
      <th>
        <center>
          <em><small>Two floors of a building: a group of people goes up and another goes down
          <br>Floor 1 on the left and 0 on the right</small></em>
        </center>
      </th>
    </tr>
    <tr>
      <td>
        <center>
          <video style="display:block; margin: 0 auto;" controls poster="_static/micro_granular_stairs_0_1.png" width="100%" height="auto">
            <source src="_static/micro_granular_stairs_0_1.mp4" />
          </video>
        </center>
      </td>
    </tr>
  </table>

..  literalinclude:: ../../examples/micro/granular/micro_granular.py
  :caption: micro_granular.py
  :linenos:
  :language: python
  :lines: 1-

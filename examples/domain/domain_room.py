# Authors:
#     Sylvain Faure <sylvain.faure@universite-paris-saclay.fr>
#     Bertrand Maury <bertrand.maury@universite-paris-saclay.fr>
#
#     cromosim/examples/domain/domain_room.py
#
#     python domain_room.py
#
# License: GPL

"""
This script allows you to create a very simple domain: a room with a single door.
The layout of the room is given by the image "room.png". Then we add a
door using a red line (Line2D object from Matplotlib) and an obstacle using a
black circle (Circle object from Matplotlib).
"""

import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.lines import Line2D

from cromosim.domain import Domain
from cromosim.domain import Destination

# To create a Domain object from a background image
dom = Domain(name='room', background='room.png', pixel_size=0.1)

# To define the color for the walls
wall_color = [0, 0, 0]

# To add an obstacle using a matplotlib shape colored with wall_color :
#     Circle( (center_x,center_y), radius )
circle = Circle((20.0, 7.0), 1.0)
dom.add_shape(circle, outline_color=wall_color, fill_color=wall_color)

# To define the color for the issue of the room
door_color = [255, 0, 0]

# To add a door using a matplotlib shape :
#     Line2D(xdata, ydata, linewidth)
line = Line2D([17.0, 23.0], [3.1, 3.1], linewidth=2)
dom.add_shape(line, outline_color=door_color, fill_color=door_color)

# To build the domain :
dom.build_domain()

# To plot the domain : backgroud + added shapes
dom.plot(id=1, title="Domain")

# To create a Destination object towards the door
dest = Destination(name='door', colors=[door_color],
                   excluded_colors=[wall_color])
dom.add_destination(dest)

# To plot the wall distance and its gradient
dom.plot_wall_dist(id=2, step=20,
                   title="Distance to walls and its gradient",
                   savefig=False, filename="room_wall_distance.png")

# To plot the distance to the red door and the correspondant
# desired velocity
dom.plot_desired_velocity('door', id=3, step=20,
                          title="Distance to the destination and desired velocity",
                          savefig=False, filename="room_desired_velocity.png")

print("===> Domain: ", dom)

plt.show()

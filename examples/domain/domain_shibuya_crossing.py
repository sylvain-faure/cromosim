# Authors:
#     Sylvain Faure <sylvain.faure@universite-paris-saclay.fr>
#     Bertrand Maury <bertrand.maury@universite-paris-saclay.fr>
#
#     cromosim/examples/domain/domain_shibuya_crossing.py
#
#     python domain_shibuya_crossing.py
#
# License: GPL

"""
This script allows you to make a domain with multiple destinations.
"""

import matplotlib.pyplot as plt

from cromosim.domain import Domain
from cromosim.domain import Destination


# The walls are represented here by the black and magenta colors
black = [0, 0, 0]
magenta = [148, 33, 146]
wall_colors = [black, magenta]

# To create a Domain object from a background image
dom = Domain(name='shibuya_crossing',
             pixel_size=0.033,
             background="shibuya_crossing.png",
             wall_colors=wall_colors)

# To build the domain :
dom.build_domain()

# To plot the domain
dom.plot(id=1, title="Domain")

# To define all the destinations

# --> towards the blue disk
blue = [4, 51, 255]
dest_blue = Destination(name='blue', colors=[blue], excluded_colors=wall_colors)
dom.add_destination(dest_blue)

# --> towards the brown disk
brown = [170, 121, 66]
dest_brown = Destination(name='brown', colors=[brown],
                         excluded_colors=wall_colors)
dom.add_destination(dest_brown)

# --> towards the cyan disk
cyan = [0, 253, 255]
dest_cyan = Destination(name='cyan', colors=[cyan],
                        excluded_colors=wall_colors)
dom.add_destination(dest_cyan)

# --> towards the green disk
green = [0, 249, 0]
dest_green = Destination(name='green', colors=[green],
                         excluded_colors=wall_colors)
dom.add_destination(dest_green)

# --> towards the pink disk
pink = [255, 64, 255]
dest_pink = Destination(name='pink', colors=[pink],
                        excluded_colors=wall_colors)
dom.add_destination(dest_pink)

print("===> Domain: ", dom)

# To plot the wall distance and its gradient
dom.plot_wall_dist(id=3, step=40,
                   title="Distance to walls and its gradient",
                   savefig=False, filename="shibuya_crossing_wall_distance.png")

# To plot all the desired velocity fields
dom.plot_desired_velocity('blue', id=4, step=40,
                          title='Distance to the destination "blue" and desired velocity',
                          savefig=False, filename="shibuya_crossing_blue.png")

dom.plot_desired_velocity('brown', id=5, step=40,
                          title='Distance to the destination "brown" and desired velocity',
                          savefig=False, filename="shibuya_crossing_brown.png")

dom.plot_desired_velocity('cyan', id=6, step=40,
                          title='Distance to the destination "cyan" and desired velocity',
                          savefig=False, filename="shibuya_crossing_cyan.png")

dom.plot_desired_velocity('green', id=7, step=40,
                          title='Distance to the destination "green" and desired velocity',
                          savefig=False, filename="shibuya_crossing_green.png")

dom.plot_desired_velocity('pink', id=8, step=40,
                          title='Distance to the destination "pink" and desired velocity',
                          savefig=False, filename="shibuya_crossing_pink.png")

plt.show()

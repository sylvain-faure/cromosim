# Authors:
#     Sylvain Faure <sylvain.faure@universite-paris-saclay.fr>
#     Bertrand Maury <bertrand.maury@universite-paris-saclay.fr>
#
#     cromosim/examples/domain/domain_stadium.py
#
#     python domain_stadium.py
#
# License: GPL

"""
This script allows to make a circular domain, here an athletics track
around a stadium.
"""

from cromosim import *

## To create a Domain object from a background image
## By default the black color is for the walls
dom = Domain(name = 'stadium', pixel_size = 0.1, background = "stadium.png")

## To build the domain :
dom.build_domain()

## To plot the domain : backgroud + added shapes (none here, they were added
## directly to the drawing using drawing software)
dom.plot(id=1,title="Domain")

## To create a Destination object towards the red line
## Since the domain is circular here, we want people to go to the red line
## and then cross it and continue to turn. For this, the color green is used
## to block one side of the track when calculating the desired velocity (as a
## temporary wall will do). Once this calculation is done, the desired speed
## is imposed at [-1,0] on the red and green pixels..
dest = Destination(name='running_track',
                   colors=[[255,0,0]],
                   excluded_colors = [[0,0,0],[0,255,0]],
                   desired_velocity_from_color =
                        [[255,0,0, -1,0],[0,255,0, -1,0]])
dom.add_destination(dest)

## To plot the wall distance and its gradient
dom.plot_wall_dist(id=2,step=30,
    title="Distance to walls and its gradient",
    savefig=False, filename="stadium_wall_distance.png")

## To plot the distance to the red line and the correspondant desired velocity
dom.plot_desired_velocity('running_track',id=3,step=30,
    title="Distance to the destination and desired velocity",
    savefig=False, filename="stadium_desired_velocity.png")

print("===> Domain: ",dom)

plt.show()

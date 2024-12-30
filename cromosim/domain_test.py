import pytest
from matplotlib.patches import Circle
from matplotlib.lines import Line2D
import numpy as np

from cromosim.domain import Domain
from cromosim.domain import Destination

def test_create_domain():
    # To create a Domain object without background image

    width = 600
    height = 400

    dom = Domain(
       name='room', 
       pixel_size=0.1, 
       width=width, height=height,
       wall_colors=[[0, 0, 0]])

    # To add an obstacle using a matplotlib shape colored with wall_color :
    #     Circle( (center_x,center_y), radius )
    circle = Circle((30.0, 20.0), 2.0)
    dom.add_shape(circle, outline_color=[0, 0, 0], fill_color=[0, 0, 0])

    # To add a door using a matplotlib shape :
    #     Line2D(xdata, ydata, linewidth)
    line = Line2D([17.0, 23.0], [3.1, 3.1], linewidth=2)
    dom.add_shape(line, outline_color=[255, 0, 0], fill_color=[255, 0, 0])

    # To build the domain :
    dom.build_domain()

    # To create a Destination object towards the door
    dest = Destination(
        name='door', colors=[[255, 0, 0]],
        excluded_colors=[[0, 0, 0]])
    dom.add_destination(dest)

    assert dom.image.shape == (height, width, 3)
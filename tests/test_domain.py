import numpy as np
from matplotlib.patches import Circle
from matplotlib.lines import Line2D

from cromosim.domain import Domain
from cromosim.domain import Destination


def test_create_domain1():
    # To create a Domain object without background image

    width = 600
    height = 400

    dom = Domain(
        name="room",
        pixel_size=0.1,
        width=width,
        height=height,
        wall_colors=[[0, 0, 0]],
    )

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
        name="door",
        colors=[[255, 0, 0]],
        excluded_colors=[[0, 0, 0]],
    )
    dom.add_destination(dest)

    assert dom.image.shape == (height, width, 3)



def test_create_domain(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    width = 600
    height = 400

    dom = Domain(
        name="room",
        pixel_size=0.1,
        width=width,
        height=height,
        wall_colors=[[0, 0, 0]],
    )

    circle = Circle((30.0, 20.0), 2.0)
    dom.add_shape(circle, outline_color=[0, 0, 0], fill_color=[0, 0, 0])

    line = Line2D([17.0, 23.0], [3.1, 3.1], linewidth=2)
    dom.add_shape(line, outline_color=[255, 0, 0], fill_color=[255, 0, 0])

    dom.build_domain()

    dest = Destination(
        name="door",
        colors=[[255, 0, 0]],
        excluded_colors=[[0, 0, 0]],
    )
    dom.add_destination(dest)

    assert dom.image.shape == (height, width, 3)
    assert dom.wall_distance.shape == (height, width)
    assert dest.distance.shape == (height, width)
    assert dest.desired_velocity_X.shape == (height, width)
    assert dest.desired_velocity_Y.shape == (height, width)


def test_people_fields_use_domain_origin_velocity_scale_and_radius():
    dom = Domain(
        name="offset-domain",
        pixel_size=1.0,
        xmin=10.0,
        ymin=20.0,
        width=3,
        height=2,
    )
    dom.wall_distance = np.array([[5.0, 6.0, 7.0], [8.0, 9.0, 10.0]])

    dest = Destination(name="door", colors=[[255, 0, 0]], velocity_scale=2.0)
    dest.distance = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    dest.desired_velocity_X = np.array([[1.0, 0.0, -1.0], [0.5, 0.0, -0.5]])
    dest.desired_velocity_Y = np.array([[0.0, 1.0, 0.0], [0.0, -1.0, 0.0]])
    dom.destinations = {"door": dest}

    xyrv = np.array(
        [
            [10.5, 20.5, 0.1, 1.5],
            [12.5, 21.5, 0.2, 2.0],
        ]
    )
    destinations = np.array(["door", "door"])

    rows, cols, desired = dom.people_desired_velocity(xyrv, destinations)
    np.testing.assert_array_equal(rows, [0, 1])
    np.testing.assert_array_equal(cols, [0, 2])
    np.testing.assert_allclose(desired, [[3.0, 0.0], [-2.0, 0.0]])

    rows2, cols2, target_distance = dom.people_target_distance(
        xyrv, destinations, rows, cols
    )
    np.testing.assert_array_equal(rows2, rows)
    np.testing.assert_array_equal(cols2, cols)
    np.testing.assert_allclose(target_distance, [0.9, 5.8])

    rows3, cols3, wall_distance = dom.people_wall_distance(xyrv, rows, cols)
    np.testing.assert_array_equal(rows3, rows)
    np.testing.assert_array_equal(cols3, cols)
    np.testing.assert_allclose(wall_distance, [4.9, 9.8])

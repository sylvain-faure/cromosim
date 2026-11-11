# Authors:
#     Sylvain Faure <sylvain.faure@universite-paris-saclay.fr>
#     Bertrand Maury <bertrand.maury@universite-paris-saclay.fr>
#
#      cromosim/examples/cellular_automata/cellular_automata.py
#      python cellular_automata.py --json input.json
#
# License: GPL

import sys
import os
import json
import numpy as np
from optparse import OptionParser
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, Rectangle, Polygon
from matplotlib.lines import Line2D

import logging
logging.basicConfig(
    level=logging.INFO,
    # format="%(asctime)s | %(levelname)-8s | %(name)s | %(message)s",
    format="%(levelname)-8s | %(name)s | %(message)s",
)
LOGGER = logging.getLogger(__name__)
# LOG_FORMAT = "%(asctime)s | %(levelname)-8s | %(name)s | %(message)s"

from cromosim.ca import plot_people_according_to_initial_door_distance
from cromosim.ca import exit, plot_people_paths
from cromosim.ca import sequential_update, parallel_update
from cromosim.domain import Domain
from cromosim.domain import Destination

"""
    python cellular_automata.py --json input.json
"""

parser = OptionParser(usage="usage: %prog [options] filename",
                      version="%prog 1.0")
parser.add_option('--json', dest="jsonfilename", default="input.json",
                  type="string", action="store", help="Input json filename")
opt, remainder = parser.parse_args()

LOGGER.info(f"Loading JSON configuration: {opt.jsonfilename}")

with open(opt.jsonfilename) as json_file:
    try:
        config = json.load(json_file)
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid JSON in {opt.jsonfilename}: {exc}") from exc

os.makedirs(config["prefix"], exist_ok=True)

"""
    Get parameters from json file :
    For the domain :
    |    name: string
    |        Domain name
    |    background: string
    |        Image file used as background
    |    px: float
    |        Pixel size in meters (also called space step)
    |    width: integer
    |        Domain width (equal to the width of the background image)
    |    height: integer
    |        Domain height (equal to the height of the background image)
    |    wall_colors: list
    |        rgb colors for walls
    |        [ [r,g,b],[r,g,b],... ]
    |    shape_lines: list
    |        Used to define the Matplotlib Polyline shapes,
    |        [
    |          {
    |             "xx": [x0,x1,x2,...],
    |             "yy": [y0,y1,y2,...],
    |             "linewidth": float,
    |             "outline_color": [r,g,b],
    |             "fill_color": [r,g,b]
    |          },...
    |        ]
    |    shape_circles: list
    |        Used to define the Matplotlib Circle shapes,
    |        [
    |           {
    |             "center_x": float,
    |             "center_y": float,
    |             "radius": float,
    |             "outline_color": [r,g,b],
    |             "fill_color": [r,g,b]
    |            },...
    |        ]
    |    shape_ellipses: list
    |        Used to define the Matplotlib Ellipse shapes,
    |        [
    |           {
    |             "center_x": float,
    |             "center_y": float,
    |             "width": float,
    |             "height": float,
    |             "angle_in_degrees_anti-clockwise": float (degre),
    |             "outline_color": [r,g,b],
    |             "fill_color": [r,g,b]
    |            },...
    |        ]
    |    shape_rectangles: list
    |        Used to define the Matplotlib Rectangle shapes,
    |        [
    |           {
    |             "bottom_left_x": float,
    |             "bottom_left_y": float,
    |             "width": float,
    |             "height": float,
    |             "angle_in_degrees_anti-clockwise": float (degre),
    |             "outline_color": [r,g,b],
    |             "fill_color": [r,g,b]
    |            },...
    |        ]
    |    shape_polygons: list
    |        Used to define the Matplotlib Polygon shapes,
    |        [
    |           {
    |             "xy": float,
    |             "outline_color": [r,g,b],
    |             "fill_color": [r,g,b]
    |            },...
    |        ]
    |    destinations: list
    |        Used to define the Destination objects,
    |        [
    |           {
    |             "name": string,
    |             "colors": [[r,g,b],...],
    |             "excluded_colors": [[r,g,b],...],
    |             "desired_velocity_from_color": [] or
    |             [
    |                {
    |                   "color": [r,g,b],
    |                   "gradient": [ex,ey]
    |                },...
    |             ],
    |             "velocity_scale": float,
    |             "next_destination": [null] or [string],
    |             "next_domain": [null] or [string],
    |             "next_transit_box": [null] or [ [[x0,y0],...,[x3,y3]] ],
    |             "next_transit_box_xshift": [0.0],
    |             "next_transit_box_yshift": [0.0],
    |            },...
    |        ]
    |--------------------
    prefix: string
        Folder name to store the results
    seed: integer
        Random seed which can be used to reproduce a random selection if >0
    Np: integer
        Number of persons
    kappa: float
        Parameter for Static Floor Field
    update_strategy: string 
        "sequential" or "parallel" strategy to compute the new positions of people
    Tf: float
        Final time
    dt: float
        Time step
    plot_period: integer
        The results will be displayed every "plot_period" iterations
"""


"""
    Build the Domain objects
"""

domain_config = config["domain"]

name = domain_config["name"]
LOGGER.info(f"Start building domain: '{name}'")
if domain_config["background"] == "":
    domain = Domain(
        name=name,
        pixel_size=domain_config["px"],
        width=domain_config["width"],
        height=domain_config["height"],
        wall_colors=domain_config["wall_colors"],
    )
else:
    domain = Domain(
        name=name,
        background=domain_config["background"],
        pixel_size=domain_config["px"],
        wall_colors=domain_config["wall_colors"],
    )
    
LOGGER.info(
    f"Building domain: '{name}' => "
    f"Add geometric objects described in the JSON "
    f"configuration")

for shape_config in domain_config["shape_lines"]:
    shape = Line2D(
        shape_config["xx"], shape_config["yy"],
        linewidth=shape_config["linewidth"],
    )
    domain.add_shape(
        shape,
        outline_color=shape_config["outline_color"],
        fill_color=shape_config["fill_color"],
    )

for shape_config in domain_config["shape_circles"]:
    shape = Circle(
        (shape_config["center_x"], shape_config["center_y"]),
        shape_config["radius"],
    )
    domain.add_shape(
        shape,
        outline_color=shape_config["outline_color"],
        fill_color=shape_config["fill_color"],
    )

for shape_config in domain_config["shape_ellipses"]:
    shape = Ellipse(
        (shape_config["center_x"], shape_config["center_y"]),
        shape_config["width"], shape_config["height"],
        shape_config["angle_in_degrees_anti-clockwise"],
    )
    domain.add_shape(
        shape,
        outline_color=shape_config["outline_color"],
        fill_color=shape_config["fill_color"],
    )

for shape_config in domain_config["shape_rectangles"]:
    shape = Rectangle(
        (shape_config["bottom_left_x"], shape_config["bottom_left_y"]),
        shape_config["width"], shape_config["height"],
        shape_config["angle_in_degrees_anti-clockwise"],
    )
    domain.add_shape(
        shape,
        outline_color=shape_config["outline_color"],
        fill_color=shape_config["fill_color"],
    )

for shape_config in domain_config["shape_polygons"]:
    shape = Polygon(shape_config["xy"])
    domain.add_shape(
        shape,
        outline_color=shape_config["outline_color"],
        fill_color=shape_config["fill_color"],
    )
    
LOGGER.info(f"Building domain: {name} => build_domain()")    
domain.build_domain()

LOGGER.info(
    f"Building domain: {name} => "
    f"Add the destinations described in the JSON "
    f"configuration")

for destination_index, destination_config in enumerate(
        domain_config["destinations"]
    ):
    desired_velocity_from_color = [
        np.concatenate((cdv["color"], cdv["desired_velocity"]))
        for cdv in destination_config["desired_velocity_from_color"]
    ]
    destination = Destination(
        name=destination_config["name"], 
        colors=destination_config["colors"],
        excluded_colors=destination_config["excluded_colors"],
        desired_velocity_from_color=desired_velocity_from_color,
        velocity_scale=destination_config["velocity_scale"],
        next_destination=destination_config["next_destination"],
        next_domain=destination_config["next_domain"],
        next_destination_weights=destination_config["next_destination_weights"],
        next_transit_box=destination_config["next_transit_box"],
        next_transit_box_xshift=destination_config["next_transit_box_xshift"],
        next_transit_box_yshift=destination_config["next_transit_box_yshift"],
    )
    LOGGER.info(
            f"Building domain: {name} => "
            f"Add destination {destination_index}: {destination_config['name']}")
    domain.add_destination(destination)
    
    domain.plot_desired_velocity(
        destination_config["name"], 
        id=10+destination_index, 
        step=20,
        title=f"Destination '{destination_config["name"]}' "
                f"in '{name}' domain")
    
    domain.plot(
        id=0, 
        title=f"Geometry of the '{name}' domain")
    domain.plot_wall_dist(
        id=1, 
        step=20,
        title=f"Distance from walls for the '{name}' domain")


"""
    Initialization
"""

# Current time
t = 0

people = np.ma.MaskedArray(np.zeros((domain.height, domain.width), dtype=int),
                        mask=domain.wall_mask)

# Number of cells
Nc = domain.height*domain.width - domain.wall_id[0].shape[0]
LOGGER.info(f"===> Number of cells = {Nc}")

# People initialisation taking to account masked positions
rng = np.random.RandomState()
if (config["seed"] > 0):
    rng = np.random.RandomState(config["seed"])
LOGGER.info(f"===> seed  = {rng.get_state()[1][0]}")

s = 0
if (Nc < domain.height*domain.width):
    imin = domain.wall_id[0].min()+1
    imax = domain.wall_id[0].max()-1
    jmin = domain.wall_id[1].min()+1
    jmax = domain.wall_id[1].max()-1
else:
    imin = 0
    imax = domain.height
    jmin = 0
    jmax = domain.width
while (s != config["Np"]):
    # people.data[rng.randint(imin,imax+1,Np-s),
    #             rng.randint(jmin+int(0.25*(jmax-jmin)),
    #                         jmax-int(0.25*(jmax-jmin))+1,Np-s)] = 1
    people.data[
        rng.randint(imin, imax+1, config["Np"]-s), 
        rng.randint(jmin, jmax+1, config["Np"]-s)
    ] = 1
    s = np.sum(people)
ind = np.where(people == 1)
people_ij = -np.ones((config["Np"], 2), dtype=int)  # i(0<=i<=height-1) j(0<=j<=width-1)
people_ij[:, 0] = ind[0]
people_ij[:, 1] = ind[1]
# print("===> Init : people_ij = ",people_ij)

# Array to store the results
results = np.copy(people_ij).reshape((config["Np"], 2, 1))

# Static Floor Field
weight = np.exp(-config["kappa"]*domain.destinations["door"].distance)

cc = 0
iter = 0

# plot_people_according_to_current_door_distance(1, people, dom,
#     savefig=True, filename=prefix+'/cellular_automata_'+str(iter).zfill(6)+'.png')
plot_people_according_to_initial_door_distance(1, people, domain,
                                            results, savefig=True,
                                            filename=config["prefix"]+'/cellular_automata_'+str(iter).zfill(6)+'.png')

while (np.sum(people) > 0):
    if (config["update_strategy"] == "parallel"):
        people, people_ij = parallel_update(people, people_ij, weight,
                                            friction=0)
    elif (config["update_strategy"] == "sequential"):
        people, people_ij = sequential_update(people, people_ij, weight,
                                            shuffle='random', randomstate=rng)
    else:
        print("Bad value for update strategy... EXIT")
        sys.exit()
    people, people_ij = exit(domain, people, people_ij)
    results = np.concatenate((results, people_ij.reshape((config["Np"], 2, 1))), axis=2)
    t += config["dt"]
    cc += 1
    iter += 1
    LOGGER.info(f"========> time = {t}, number of persons = {np.sum(people)}")
    if (cc >= config["plot_period"]):
        # plot_people_according_to_current_door_distance(1, people, dom,
        #     savefig=True,
        #     filename=prefix+'/cellular_automata_'+str(iter).zfill(6)+'.png')
        plot_people_according_to_initial_door_distance(1, people, domain,
                                                    results, savefig=True,
                                                    filename=config["prefix"]+'/cellular_automata_'+str(iter).zfill(6)+'.png')
        plt.pause(0.05)
        cc = 0

# plot_people_according_to_current_door_distance(1, people, dom)
# plot_people_according_to_initial_door_distance(2, people, dom, results)
# plot_people_according_to_exit_times(3, config["dt"], people, dom, results)
plot_people_paths(2, config["dt"], domain.pixel_size, people, domain, results,
                savefig=True, filename=config["prefix"]+'/cellular_automata_paths.png')

plt.show()

# Authors:
#     Sylvain Faure <sylvain.faure@universite-paris-saclay.fr>
#     Bertrand Maury <bertrand.maury@universite-paris-saclay.fr>
#
#      cromosim/examples/compartments/compartments.py
#      python compartments.py --json input.json
#
# License: GPL

import sys
import os
import json
from optparse import OptionParser
import numpy as np
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

from cromosim.domain import Domain
from cromosim.domain import Destination
from cromosim.comp import plot_compt, iteration

"""
    python compartments.py --json input.json
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
    |                   "desired_velocity": [ex,ey]
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
    Np_rooms: integer
        Number of persons per rooms
    area: float
        for visualization...
    CircAngles: list
        for visualization (orientation): [0,0,-90,...]
    RoomNames: list of strings
        all the room names
    Np_rooms: list
        number of persons in each room
    RoomInlets:
        all the room links: [[1,2], [2,3,4],...]
    "TravelTimeWeights": list 
        weights for each way: [[1,1], [1,2,1],...]
    DoorCenters:
        coordinates of the door centers
    RoomCenters:
        coordinates of the room centers 
    DoorRoomCapacity: list
        all the door capacities (number of persons per minute)
    Nseconds: int
        final time in seconds
"""

jdom = config["domain"]
LOGGER.info(f"===> JSON data used to build the domain: {jdom}" )
prefix = config["prefix"]
if not os.path.exists(prefix):
    os.makedirs(prefix)
seed = config["seed"]
Np_rooms = config["Np_rooms"]
area = config["area"]
RoomNames = config["RoomNames"]
Np_rooms = np.array(config["Np_rooms"])
RoomInlets = config["RoomInlets"]
TravelTimeWeights = config["TravelTimeWeights"]
DoorCenters = np.array(config["DoorCenters"])
RoomCenters = np.array(config["RoomCenters"])
CircAngles = np.array(config["CircAngles"])
DoorRoomCapacity = config["DoorRoomCapacity"]
Nseconds = config["Nseconds"]
LOGGER.info(f"===> Number of persons per rooms: {Np_rooms}")
Np = 0
for n in Np_rooms:
    Np += n
LOGGER.info(f"===> Number of persons: {Np}")
Nrooms = len(RoomNames)
LOGGER.info(f"===> Number of rooms: {Nrooms}")
LOGGER.info(f"===> Capacity of each exit door: {DoorRoomCapacity}")

"""
    Build the Domain objects
"""
jname = jdom["name"]
LOGGER.info(f"===> Build domain: {jname}")
jbg = jdom["background"]
jpx = jdom["px"]
jwidth = jdom["width"]
jheight = jdom["height"]
jwall_colors = jdom["wall_colors"]
if (jbg == ""):
    dom = Domain(name=jname, pixel_size=jpx, width=jwidth, height=jheight,
                 wall_colors=jwall_colors)
else:
    dom = Domain(name=jname, background=jbg, pixel_size=jpx,
                 wall_colors=jwall_colors)
# To add lines : Line2D(xdata, ydata, linewidth)
for sl in jdom["shape_lines"]:
    line = Line2D(sl["xx"], sl["yy"], linewidth=sl["linewidth"])
    dom.add_shape(line, outline_color=sl["outline_color"],
                  fill_color=sl["fill_color"])
# To add circles : Circle( (center_x,center_y), radius )
for sc in jdom["shape_circles"]:
    circle = Circle((sc["center_x"], sc["center_y"]), sc["radius"])
    dom.add_shape(circle, outline_color=sc["outline_color"],
                  fill_color=sc["fill_color"])
# To add ellipses : Ellipse( (center_x,center_y), width, height,
#                            angle_in_degrees_anti-clockwise )
for se in jdom["shape_ellipses"]:
    ellipse = Ellipse((se["center_x"], se["center_y"]),
                      se["width"], se["height"],
                      se["angle_in_degrees_anti-clockwise"])
    dom.add_shape(ellipse, outline_color=se["outline_color"],
                  fill_color=se["fill_color"])
# To add rectangles : Rectangle( (bottom_left_x,bottom_left_y), width, height,
#                                 angle_in_degrees_anti-clockwise )
for sr in jdom["shape_rectangles"]:
    rectangle = Rectangle((sr["bottom_left_x"], sr["bottom_left_y"]),
                          sr["width"], sr["height"],
                          sr["angle_in_degrees_anti-clockwise"])
    dom.add_shape(rectangle, outline_color=sr["outline_color"],
                  fill_color=sr["fill_color"])
# To add polygons : Polygon( [[x0,y0],[x1,y1],...] )
for spo in jdom["shape_polygons"]:
    polygon = Polygon(spo["xy"])
    dom.add_shape(polygon, outline_color=spo["outline_color"],
                  fill_color=spo["fill_color"])
# To build the domain : background + shapes
dom.build_domain()
# To add all the available destinations
for j, dd in enumerate(jdom["destinations"]):
    desired_velocity_from_color = []
    for gg in dd["desired_velocity_from_color"]:
        desired_velocity_from_color.append(
            np.concatenate((gg["color"], gg["desired_velocity"])))
    dest = Destination(name=dd["name"], colors=dd["colors"],
                       excluded_colors=dd["excluded_colors"],
                       desired_velocity_from_color=desired_velocity_from_color,
                       velocity_scale=dd["velocity_scale"],
                       next_destination=dd["next_destination"],
                       next_domain=dd["next_domain"],
                       next_transit_box=dd["next_transit_box"])
    LOGGER.info(f"===> Destination: {dest}")
    dom.add_destination(dest)

    dom.plot_desired_velocity(dd["name"], id=10+j, step=20)

LOGGER.info(f"===> Domain: {dom}")

dom.plot(id=0)
dom.plot_wall_dist(id=1, step=20)


# Maximal number of inlets for a room
NiorMax = 0
# Number of inlets for each room
Nior = []
for ri in RoomInlets:
    NiorMax = max(NiorMax, len(ri))
    Nior.append(len(ri))
LOGGER.info(f"===> Maximal number of inlets for a room: {NiorMax}")
LOGGER.info(f"===> Number of inlets for each room: {Nior}")

# Room numbers associate to the inlets for each room
List_iOr = np.zeros([Nrooms, NiorMax], dtype=int)
for i, ri in enumerate(RoomInlets):
    List_iOr[i, :len(ri)] = ri
# LOGGER.info(f"===> Room numbers associate to the inlets for each room: {List_iOr}")

# Compute travel times :
T_iOr = np.zeros([Nrooms, NiorMax], dtype=int)
for ir in np.arange(Nrooms):
    for i in np.arange(Nior[ir]):
        jr = List_iOr[ir, i]
        ij_exit = np.rint(DoorCenters[ir, :2]/dom.pixel_size).astype(int)
        ij_entrance = np.rint(DoorCenters[jr, :2]/dom.pixel_size).astype(int)
        T_iOr[ir, i] = np.ceil(
            dom.destinations["door"].distance.data[ij_entrance[1], ij_entrance[0]] -
            dom.destinations["door"].distance.data[ij_exit[1], ij_exit[0]])
        T_iOr[ir, i] = TravelTimeWeights[ir][i]*T_iOr[ir, i]
LOGGER.info(f"Travel times: {T_iOr}")

NPrir = np.zeros([Nrooms, NiorMax])

# Number of persons for each room
NPir = np.zeros([Nseconds, Nrooms])
NPir[0, :] = Np_rooms

# Number of persons who's waiting to leave each room
NPwir = np.zeros([Nseconds, Nrooms])
NPwir[0, :] = NPir[0, :]

# Nrooms of persons upstream the exit
Flux = np.zeros([Nseconds, Nrooms])

LOGGER.info(f"--> NPir at t=0: {NPir[0, :].sum()}")

it = 0
plot_compt(5, RoomNames, RoomCenters, DoorCenters, CircAngles, NPir[it, :],
           NPrir, NPwir[it, :], Nior, List_iOr, dom, area, savefig=True,
           filename=prefix+"fig_"+str(it).zfill(6)+".png",
           title='time = '+str(it)+' s')

for it in np.arange(1, Nseconds):

    Flux, NPir, NPwir, NPrir = iteration(it, Nrooms, DoorRoomCapacity, NPir, NPrir,
                                         NPwir, Nior, List_iOr, T_iOr, Flux)
    LOGGER.info(f"---------> it = {it}")
    # LOGGER.info(f"--> Flux : ", Flux)
    LOGGER.info(f"--> NPir : {NPir[it, :]}")
    # LOGGER.info(f"--> NPwir : {NPwir[it, :]}")
    # LOGGER.info(f"--> NPrir : {NPrir}")
    plot_compt(5, RoomNames, RoomCenters, DoorCenters, CircAngles, NPir[it, :],
               NPrir, NPwir[it, :], Nior, List_iOr, dom, area, savefig=True,
               filename=prefix+"fig_"+str(it).zfill(6)+".png",
               title='time = '+str(it)+' s')
    plt.pause(0.1)

plt.show()

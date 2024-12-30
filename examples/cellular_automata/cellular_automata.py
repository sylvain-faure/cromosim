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

from cromosim.ca import plot_people_according_to_initial_door_distance
from cromosim.ca import exit, plot_people_paths
from cromosim.ca import sequential_update, parallel_update
from cromosim.domain import Domain
from cromosim.domain import Destination

"""
    python cellular_automata.py --json input.json
"""

parser = OptionParser(usage="usage: %prog [options] filename", version="%prog 1.0")
parser.add_option('--json', dest="jsonfilename", default="input.json", type="string",
                  action="store", help="Input json filename")
opt, remainder = parser.parse_args()
print("===> JSON filename = ", opt.jsonfilename)
with open(opt.jsonfilename) as json_file:
    try:
        input = json.load(json_file)
    except json.JSONDecodeError as msg:
        print(msg)
        print("Failed to load json file ", opt.jsonfilename)
        print("Check its content : ")
        print("(https://fr.wikipedia.org/wiki/JavaScript_Object_Notation)")
        sys.exit()

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
    |             "next_destination": null or string,
    |             "next_domain": null or string,
    |             "next_transit_box": null or [[x0,y0],...,[x3,y3]]
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
    Tf: float
        Final time
    dt: float
        Time step
    drawper: integer
        The results will be displayed every "drawper" iterations
"""

jdom = input["domain"]
print("===> JSON data used to build the domain : ", jdom)
prefix = input["prefix"]
if not os.path.exists(prefix):
    os.makedirs(prefix)
seed = input["seed"]
Np = input["Np"]
update_strategy = input["update_strategy"]
kappa = input["kappa"]
Tf = input["Tf"]
dt = input["dt"]
drawper = input["drawper"]
print("===> Number of persons = ", Np)
print("===> Final time, Tf = ", Tf)
print("===> Time step, dt = ", dt)
print("===> To draw the results each drawper iterations, drawper = ", drawper)

"""
    Build the Domain objects
"""
jname = jdom["name"]
print("===> Build domain ", jname)
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
    if (j > 0) or (dd["name"] != "door"):
        print("===> EXIT : For the moment the cellular automata implemented does")
        print("     not make it possible to direct individuals to several")
        print("     destinations, nor to use several Domain object. The json")
        print("     file must therefore only contain one domain and one destination")
        print("     necessarily called \"door\".")
        sys.exit()
    desired_velocity_from_color = []
    for gg in dd["desired_velocity_from_color"]:
        desired_velocity_from_color.append(
            np.concatenate((gg["color"], gg["gradient"])))
    dest = Destination(name=dd["name"], colors=dd["colors"],
                       excluded_colors=dd["excluded_colors"],
                       desired_velocity_from_color=desired_velocity_from_color,
                       velocity_scale=dd["velocity_scale"],
                       next_destination=dd["next_destination"],
                       next_domain=dd["next_domain"],
                       next_transit_box=dd["next_transit_box"])
    print("===> Destination : ", dest)
    dom.add_destination(dest)

    dom.plot_desired_velocity(dd["name"], id=10+j, step=20)

print("===> Domain : ", dom)

dom.plot(id=0)
dom.plot_wall_dist(id=1, step=20)


"""
    Initialization
"""

# Current time
t = 0

people = np.ma.MaskedArray(np.zeros((dom.height, dom.width), dtype=int),
                           mask=dom.wall_mask)

# Number of cells
Nc = dom.height*dom.width - dom.wall_id[0].shape[0]
print("===> Number of cells = ", Nc)

# People initialisation taking to account masked positions
rng = np.random.RandomState()
if (seed > 0):
    rng = np.random.RandomState(seed)
print("===> seed  = ", rng.get_state()[1][0])

s = 0
if (Nc < dom.height*dom.width):
    imin = dom.wall_id[0].min()+1
    imax = dom.wall_id[0].max()-1
    jmin = dom.wall_id[1].min()+1
    jmax = dom.wall_id[1].max()-1
else:
    imin = 0
    imax = dom.height
    jmin = 0
    jmax = dom.width
while (s != Np):
    # people.data[rng.randint(imin,imax+1,Np-s),
    #             rng.randint(jmin+int(0.25*(jmax-jmin)),
    #                         jmax-int(0.25*(jmax-jmin))+1,Np-s)] = 1
    people.data[rng.randint(imin, imax+1, Np-s), rng.randint(jmin, jmax+1, Np-s)] = 1
    s = np.sum(people)
ind = np.where(people == 1)
people_ij = -np.ones((Np, 2), dtype=int)  # i(0<=i<=height-1) j(0<=j<=width-1)
people_ij[:, 0] = ind[0]
people_ij[:, 1] = ind[1]
# print("===> Init : people_ij = ",people_ij)

# Array to store the results
results = np.copy(people_ij).reshape((Np, 2, 1))

# Static Floor Field
weight = np.exp(-kappa*dom.destinations["door"].distance)

cc = 0
iter = 0

# plot_people_according_to_current_door_distance(1, people, dom,
#     savefig=True, filename=prefix+'/cellular_automata_'+str(iter).zfill(6)+'.png')
plot_people_according_to_initial_door_distance(1, people, dom,
                                               results, savefig=True,
                                               filename=prefix+'/cellular_automata_'+str(iter).zfill(6)+'.png')

while (np.sum(people) > 0):
    if (update_strategy == "parallel"):
        people, people_ij = parallel_update(people, people_ij, weight,
                                            friction=0)
    elif (update_strategy == "sequential"):
        people, people_ij = sequential_update(people, people_ij, weight,
                                              shuffle='random', randomstate=rng)
    else:
        print("Bad value for update strategy... EXIT")
        sys.exit()
    people, people_ij = exit(dom, people, people_ij)
    results = np.concatenate((results, people_ij.reshape((Np, 2, 1))), axis=2)
    t += dt
    cc += 1
    iter += 1
    print("========> time = ", t, " number of persons = ", np.sum(people))
    if (cc >= drawper):
        # plot_people_according_to_current_door_distance(1, people, dom,
        #     savefig=True,
        #     filename=prefix+'/cellular_automata_'+str(iter).zfill(6)+'.png')
        plot_people_according_to_initial_door_distance(1, people, dom,
                                                       results, savefig=True,
                                                       filename=prefix+'/cellular_automata_'+str(iter).zfill(6)+'.png')
        plt.pause(0.05)
        cc = 0

# plot_people_according_to_current_door_distance(1, people, dom)
# plot_people_according_to_initial_door_distance(2, people, dom, results)
# plot_people_according_to_exit_times(3, dt, people, dom, results)
plot_people_paths(2, dt, dom.pixel_size, people, dom, results,
                  savefig=True, filename=prefix+'/cellular_automata_paths.png')

plt.show()

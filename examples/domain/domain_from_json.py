# Authors:
#     Sylvain Faure <sylvain.faure@universite-paris-saclay.fr>
#     Bertrand Maury <bertrand.maury@universite-paris-saclay.fr>
#
#     cromosim/examples/domain/domain_from_json.py
#
#     python domain_from_json.py --json input.json
#
# License: GPL

"""
It is possible to describe the domain to be built as well as these
destinations in a JSON file. This file therefore contains all the elements
for drawing the domain and calculating the desired velocity fields.
Below there are explanations of its content.
To create a domain, you just have to write a json file and supply it as an
argument to this script.
"""

import sys, os
from cromosim import *
from matplotlib.patches import Circle
from matplotlib.lines import Line2D
from optparse import OptionParser
import json

## Parse the arguments of the python command. Then the json file is loaded.
parser = OptionParser(usage="usage: %prog [options] filename",
                      version="%prog 1.0")
parser.add_option('--json',dest="jsonfilename",
                  default="input.json",type="string",
                  action="store",help="Input json filename")
opt, remainder = parser.parse_args()
print("===> JSON filename = ",opt.jsonfilename)
with open(opt.jsonfilename) as json_file:
    try:
        input = json.load(json_file)
    except json.JSONDecodeError as msg:
        print(msg)
        print("Failed to load json file ",opt.jsonfilename)
        print("Check its content :")
        print("https://fr.wikipedia.org/wiki/JavaScript_Object_Notation)")
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
    |                   "desired_velocity": [ex,ey]
    |                },...
    |             ],
    |             "velocity_scale": float,
    |             "next_destination": null or string,
    |             "next_domain": null or string,
    |             "next_transit_box": null or [[x0,y0],...,[x3,y3]]
    |            },...
    |        ]
    |--------------------
"""

jdom = input["domain"]
print("===> JSON data used to build the domain : ",jdom)

"""
    Build the Domain object
"""

jname = jdom["name"]
print("===> Build domain ",jname)
jbg = jdom["background"]
jpx = jdom["px"]
jwidth = jdom["width"]
jheight = jdom["height"]
jwall_colors = jdom["wall_colors"]
if (jbg==""):
    dom = Domain(name=jname, pixel_size=jpx, width=jwidth, height=jheight,
                 wall_colors=jwall_colors)
else:
    dom = Domain(name=jname, background=jbg, pixel_size=jpx,
                 wall_colors=jwall_colors)
## To add lines : Line2D(xdata, ydata, linewidth)
for sl in jdom["shape_lines"]:
    line = Line2D(sl["xx"],sl["yy"],linewidth=sl["linewidth"])
    dom.add_shape(line,outline_color=sl["outline_color"],
                  fill_color=sl["fill_color"])
## To add circles : Circle( (center_x,center_y), radius )
for sc in jdom["shape_circles"]:
    circle = Circle( (sc["center_x"], sc["center_y"]), sc["radius"] )
    dom.add_shape(circle,outline_color=sc["outline_color"],
                  fill_color=sc["fill_color"])
## To add ellipses : Ellipse( (center_x,center_y), width, height,
##                            angle_in_degrees_anti-clockwise )
for se in jdom["shape_ellipses"]:
    ellipse = Ellipse( (se["center_x"], se["center_y"]),
                        se["width"], se["height"],
                        se["angle_in_degrees_anti-clockwise"])
    dom.add_shape(ellipse,outline_color=se["outline_color"],
                  fill_color=se["fill_color"])
## To add rectangles : Rectangle( (bottom_left_x,bottom_left_y), width,
##                                 height, angle_in_degrees_anti-clockwise )
for sr in jdom["shape_rectangles"]:
    rectangle = Rectangle( (sr["bottom_left_x"],sr["bottom_left_y"]),
                           sr["width"], sr["height"],
                           sr["angle_in_degrees_anti-clockwise"])
    dom.add_shape(rectangle,outline_color=sr["outline_color"],
                  fill_color=sr["fill_color"])
## To add polygons : Polygon( [[x0,y0],[x1,y1],...] )
for spo in jdom["shape_polygons"]:
    polygon = Polygon(spo["xy"])
    dom.add_shape(polygon,outline_color=spo["outline_color"],
                  fill_color=spo["fill_color"])

## To build the domain : background + shapes
dom.build_domain()

## To add all the available destinations
for j,dd in enumerate(jdom["destinations"]):
    desired_velocity_from_color=[]
    for gg in dd["desired_velocity_from_color"]:
        desired_velocity_from_color.append(
            np.concatenate((gg["color"],gg["desired_velocity"])))
    dest = Destination(name=dd["name"],colors=dd["colors"],
    excluded_colors=dd["excluded_colors"],
    desired_velocity_from_color=desired_velocity_from_color,
    velocity_scale=dd["velocity_scale"],
    next_destination=dd["next_destination"],
    next_domain=dd["next_domain"],
    next_transit_box=dd["next_transit_box"])
    print("===> Destination : ",dest)
    dom.add_destination(dest)

    ## To plot all the desired velocity field
    dom.plot_desired_velocity(dd["name"],id=10+j,step=20)

print("===> Domain: ",dom)

## To plot the domain
dom.plot(id=0)

## To plot the wall distance and its gradient
dom.plot_wall_dist(id=1,step=20)

plt.show()

# Authors:
#     Sylvain Faure <sylvain.faure@math.u-psud.fr>
#     Bertrand Maury <bertrand.maury@math.u-psud.fr>
#
#      cromosim/examples/compartments/compartments.py
#      python compartments.py --json input.json
#
# License: GPL

import sys, os
from cromosim import *
from cromosim.comp import *
from optparse import OptionParser
import json

"""
    python3 compartments.py --json input.json
"""
parser = OptionParser(usage="usage: %prog [options] filename",version="%prog 1.0")
parser.add_option('--json',dest="jsonfilename",default="input.json",type="string",
                  action="store",help="Input json filename")
opt, remainder = parser.parse_args()
print("===> JSON filename = ",opt.jsonfilename)
with open(opt.jsonfilename) as json_file:
    input = json.load(json_file)


"""
    Get parameters from json file :

    name: string
        Domain name
    prefix: string
        Folder name to store the results
    background: string
        Image file used as background
    px: float
        Pixel size in meters (also called space step)
    width: integer
        Domain width (equal to the width of the background image)
    height: integer
        Domain height (equal to the height of the background image)
    wall_lines : list of numpy arrays
        Polylines used to build walls, [ [[x0,x1,x2,...],[y0,y1,y2,...]],... ]
    wall_ellipses : list of numpy arrays
        Ellipses used to build walls, [ [x_center,y_center, width, height, \
        angle_in_degrees_anti-clockwise],... ]
    wall_polygons : list of numpy arrays
        Polygons used to build walls, [ [[x0,x1,x2,...],[y0,y1,y2,...]],... ]
    door_lines: list of numpy arrays
        Polylines used to build doors, [ [[x0,x1,x2,...],[y0,y1,y2,...]],... ]
    seed: integer
        Random seed which can be used to reproduce a random selection if >0
    Np_rooms: integer
        Number of persons per rooms
    area : float
    Nrooms
    RoomNames
    RoomInlets
    DoorCenters
    RoomCenters
    DoorRoomCapacity
    Nsecondes
"""
name = input["name"]
prefix = input["prefix"]
if not os.path.exists(prefix):
    os.makedirs(prefix)
background = input["background"]
px = input["px"]
width = input["width"]
height = input["height"]
wall_lines = input["wall_lines"]
wall_ellipses = input["wall_ellipses"]
wall_polygons = input["wall_polygons"]
door_lines = input["door_lines"]
seed = input["seed"]
Np_rooms = input["Np_rooms"]
area = input["area"]
RoomNames = input["RoomNames"]
Np_rooms = sp.array(input["Np_rooms"])
RoomInlets = input["RoomInlets"]
TravelTimeWeights = input["TravelTimeWeights"]
DoorCenters = sp.array(input["DoorCenters"])
RoomCenters = sp.array(input["RoomCenters"])
CircAngles = sp.array(input["CircAngles"])
DoorRoomCapacity = input["DoorRoomCapacity"]
Nsecondes = input["Nsecondes"]
print("===> Number of persons per rooms : ",Np_rooms)
Np = 0
for n in Np_rooms:
    Np += n
print("===> Number of persons : ",Np)
Nrooms = len(RoomNames)
print("===> Number of rooms : ",Nrooms)
print("===> Capacity of each exit door : ",DoorRoomCapacity)


"""
    Build the Domain
"""

## To create an Domain object
if (background==""):
    dom = Domain(name=name, pixel_size=px, width=width, height=height)
else:
    dom = Domain(name=name, background=background, pixel_size=px)
## To add lines : Line2D(xdata, ydata, linewidth)
for xy in wall_lines:
    line = Line2D( xy[0],xy[1], linewidth=1)
    dom.add_wall(line)
## To add ellipses : Ellipse( (x_center,y_center), width, height,
##                             angle_in_degrees_anti-clockwise )
for e in wall_ellipses:
    ellipse = Ellipse( (e[0], e[1]), e[2], e[3], e[4])
    dom.add_wall(ellipse)
## To add polygons : Polygon( xy )
for p in wall_polygons:
    polygon = Polygon(p)
    dom.add_wall(polygon)
## To add doors :
for xy in door_lines:
    line = Line2D( xy[0],xy[1], linewidth=1)
    dom.add_door(line)
## To build the domain : background + shapes
dom.build_domain()
## To compute the distance to the walls
dom.compute_wall_distance()
## To compute the desired velocity
dom.compute_desired_velocity()
## To show the domain dimensions
print("===> Domain : ",dom)
print("===> Wall lines : ",wall_lines)
print("===> Door lines : ",door_lines)


## Maximal number of inlets for a room
NiorMax = 0
## Number of inlets for each room
Nior = []
for ri in RoomInlets:
    NiorMax = max(NiorMax,len(ri))
    Nior.append(len(ri))
print("===> Maximal number of inlets for a room : ",NiorMax)
print("===> Number of inlets for each room : ",Nior)

## Room numbers associate to the inlets for each room
List_iOr = sp.zeros([Nrooms,NiorMax],dtype=int)
for i,ri in enumerate(RoomInlets):
    List_iOr[i,:len(ri)] = ri
#print("===> Room numbers associate to the inlets for each room",List_iOr)

## Compute travel times :
T_iOr = sp.zeros([Nrooms,NiorMax],dtype=int)
for ir in sp.arange(Nrooms):
    for i in sp.arange(Nior[ir]):
        jr = List_iOr[ir,i]
        ij_exit = sp.rint(DoorCenters[ir,:2]/px).astype(int)
        ij_entrance = sp.rint(DoorCenters[jr,:2]/px).astype(int)
        T_iOr[ir,i] = sp.ceil( \
                      dom.door_distance.data[ij_entrance[1],ij_entrance[0]] - \
                      dom.door_distance.data[ij_exit[1],ij_exit[0]] )
        T_iOr[ir,i] = TravelTimeWeights[ir][i]*T_iOr[ir,i]
print("Travel times : ",T_iOr)

NPrir = sp.zeros([Nrooms,NiorMax])

## Number of persons for each room
NPir = sp.zeros([Nsecondes , Nrooms])
NPir[0,:] = Np_rooms

## Number of persons who's waiting to leave each room
NPwir = sp.zeros([Nsecondes , Nrooms])
NPwir[0,:]=NPir[0,:]

## Nrooms of persons upstream the exit
Flux = sp.zeros([Nsecondes , Nrooms])

print("--> NPir at t=0 : ",NPir[0,:].sum())

it=0
plot_compt(5, RoomNames, RoomCenters, DoorCenters, CircAngles, NPir[it,:],
           NPrir, NPwir[it,:], Nior, List_iOr, dom, area, savefig=True,
           filename=prefix+"fig_"+str(it).zfill(6)+".png",
           title='time = '+str(it)+' s')

for it in sp.arange(1,Nsecondes):

    Flux,NPir,NPwir,NPrir = iteration(it, Nrooms, DoorRoomCapacity, NPir, NPrir,
                                      NPwir, Nior, List_iOr, T_iOr, Flux)
    print("---------> it = ",it)
    #print("--> Flux : ",Flux)
    print("--> NPir : ",NPir[it,:])
    #print("--> NPwir : ",NPwir[it,:])
    #print("--> NPrir : ",NPrir)
    plot_compt(5, RoomNames, RoomCenters, DoorCenters, CircAngles, NPir[it,:],
               NPrir, NPwir[it,:], Nior, List_iOr, dom, area, savefig=True,
               filename=prefix+"fig_"+str(it).zfill(6)+".png",
               title='time = '+str(it)+' s')

plt.show()

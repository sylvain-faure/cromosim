# Authors:
#     Sylvain Faure <sylvain.faure@math.u-psud.fr>
#     Bertrand Maury <bertrand.maury@math.u-psud.fr>
#
#      cromosim/examples/cellular_automata/cellular_automata.py
#      python cellular_automata.py --json input.json
#
# License: GPL

import sys, os
from cromosim import *
from cromosim.ca import *
from optparse import OptionParser
import json

"""
    python3 cellular_automata.py --json input.json
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
    wall_lines : list of numpy arrays
        Polylines used to build walls, [ [[x0,x1,x2,...],[y0,y1,y2,...]],... ]
    door_lines: list of numpy arrays
        Polylines used to build doors, [ [[x0,x1,x2,...],[y0,y1,y2,...]],... ]
    update_strategy: string
            Rules used to move the individuals : 'sequential' or 'parallel'
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
Np = input["Np"]
update_strategy = input["update_strategy"]
kappa = input["kappa"]
Tf = input["Tf"]
dt = input["dt"]
drawper = input["drawper"]
print("===> Number of persons = ",Np)
print("===> Final time, Tf = ",Tf)
print("===> Time step, dt = ",dt)
print("===> To draw the results each drawper iterations, drawper = ",drawper)

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


"""
    Initialization
"""

## Current time
t = 0

people = sp.ma.MaskedArray(sp.zeros((height,width),dtype=int), mask=dom.mask)

## Number of cells
Nc = height*width - dom.mask_id[0].shape[0]
print("===> Number of cells = ",Nc)

## People initialisation taking to account masked positions
rng = sp.random.RandomState()
if (seed>0):
    rng = sp.random.RandomState(seed)
print("===> seed  = ",rng.get_state()[1][0])

s = 0
if (Nc<height*width):
    imin = dom.mask_id[0].min()+1
    imax = dom.mask_id[0].max()-1
    jmin = dom.mask_id[1].min()+1
    jmax = dom.mask_id[1].max()-1
else:
    imin = 0; imax = height
    jmin = 0; jmax = width
while (s != Np):
    #people.data[rng.randint(imin,imax+1,Np-s),
    #            rng.randint(jmin+int(0.25*(jmax-jmin)),
    #                        jmax-int(0.25*(jmax-jmin))+1,Np-s)] = 1
    people.data[rng.randint(imin,imax+1,Np-s),rng.randint(jmin,jmax+1,Np-s)] = 1
    s = sp.sum(people)
ind = sp.where(people==1)
people_ij = -sp.ones((Np,2),dtype=int)# i(0<=i<=height-1) j(0<=j<=width-1)
people_ij[:,0] = ind[0]
people_ij[:,1] = ind[1]
#print("===> Init : people_ij = ",people_ij)

## Array to store the results
results = sp.copy(people_ij).reshape((Np,2,1))

## Static Floor Field
weight = sp.exp(-kappa*dom.door_distance)

cc = 0
iter = 0

plot_people_according_to_current_door_distance(1, people, dom,
    savefig=True, filename=prefix+'/cellular_automata_'+str(iter).zfill(6)+'.png')
#plot_people_according_to_initial_door_distance(2, people, dom, results)

while (sp.sum(people)>0):
    if (update_strategy == "parallel"):
        people, people_ij = parallel_update(people, people_ij, weight,
                                            friction = 0)
    elif (update_strategy == "sequential"):
        people, people_ij = sequential_update(people, people_ij, weight,
                                              shuffle='random', randomstate = rng)
    else:
        print("Bad value for update strategy... EXIT")
        sys.exit()
    people, people_ij = exit(dom, people, people_ij)
    results = sp.concatenate((results,people_ij.reshape((Np,2,1))), axis=2)
    t += dt
    cc += 1
    iter += 1
    print("========> time = ",t," number of persons = ",sp.sum(people))
    if (cc>=drawper):
        plot_people_according_to_current_door_distance(1, people, dom,
            savefig=True, filename=prefix+'/cellular_automata_'+str(iter).zfill(6)+'.png')
        #plot_people_according_to_initial_door_distance(2, people, dom, results)
        cc = 0

#plot_people_according_to_current_door_distance(1, people, dom)
#plot_people_according_to_initial_door_distance(2, people, dom, results)
#plot_people_according_to_exit_times(3, dt, people, dom, results)
plot_people_paths(2, dt, px, people, dom, results,
                  savefig=True, filename=prefix+'/cellular_automata_paths.png')

plt.show()

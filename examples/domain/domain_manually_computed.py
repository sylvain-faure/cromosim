# Authors:
#     Sylvain Faure <sylvain.faure@math.u-psud.fr>
#     Bertrand Maury <bertrand.maury@math.u-psud.fr>
#
#      cromosim/examples/domain/domain_manually_computed.py
#      python domain_manually_computed.py
#
# License: GPL

from cromosim import *
import scipy as sp
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Circle, Rectangle, Polygon, Arrow
from matplotlib.patches import Circle
from matplotlib.lines import Line2D

## Domain without a background image
pixel_size = 0.1
dom = Domain(name = 'room_400x300', pixel_size = pixel_size, width = 400, height = 300)
print(dom)

##  (x0,y1)------------------------------(x3,y1)
##    I                                     I
##    I                                     I
##  (x0,y0)-----(x1,y0)-door-(x2,y0)-----(x3,y0)
##
##                     (xd,yd)
x0 = 5.0
x1 = 19.0
x2 = 21.0
x3 = 35.0
y0 = 5.0
y1 = 25.0
xd = 20.0
yd = 2.5

## Line2D(xdata, ydata, linewidth)
line = Line2D( [x2,x3,x3,x0,x0,x1],[y0,y0,y1,y1,y0,y0], 2)
dom.add_wall(line)

## Line2D(xdata, ydata, linewidth)
line = Line2D( [x1,x2],[y0,y0], linewidth=2)
dom.add_door(line)

dom.build_domain()

## Compute manually a mask
mask = (dom.X<x0)+(dom.X>x3)+(dom.Y<y0)+(dom.Y>y1)

## Compute manually a desired velocity
desired_velocity_X = xd - dom.X
desired_velocity_Y = yd - dom.Y
norm = sp.sqrt( desired_velocity_X**2+desired_velocity_Y**2 )
desired_velocity_X = sp.ma.MaskedArray(desired_velocity_X/norm,mask=mask)
desired_velocity_Y = sp.ma.MaskedArray(desired_velocity_Y/norm,mask=mask)
dom.door_distance = norm
dom.desired_velocity_X = desired_velocity_X
dom.desired_velocity_Y = desired_velocity_Y

## Compute manually the wall distance
##  I-------------wall2----------------I
## wall3                             wall1
##  I---wall0---pL--door--pR---wall0---I
dist_wall2 = y1-dom.Y
dist_wall1 = x3-dom.X
dist_wall3 = dom.X-x0
dist_wall0 = (dom.X<=x1)*(dom.Y-y0) + \
             (dom.X>=x2)*(dom.Y-y0) + \
             (dom.X>x1)*(dom.X<x2)*sp.minimum( \
                sp.sqrt((dom.X-x1)**2+(dom.Y-y0)**2), \
                sp.sqrt((dom.X-x2)**2+(dom.Y-y0)**2)  \
             )
wall_distance = sp.ma.MaskedArray(sp.minimum(sp.minimum(sp.minimum( \
                    dist_wall0,dist_wall1),dist_wall2),dist_wall3), mask=mask \
                )
grad = sp.gradient(wall_distance,edge_order=2)
grad_X = grad[1]/pixel_size
grad_Y = grad[0]/pixel_size
norm = sp.sqrt(grad_X**2+grad_Y**2)
wall_grad_X = grad_X/norm
wall_grad_Y = grad_Y/norm
dom.wall_distance = wall_distance
dom.wall_grad_X = wall_grad_X
dom.wall_grad_Y = wall_grad_Y


dom.plot()
dom.plot_wall_dist(id=2)
dom.plot_desired_velocity(id=3)

## Custom plot function
fig = plt.figure(10)
ax1 = fig.add_subplot(111)
ax1.imshow(dom.image,interpolation='nearest',
           extent=[dom.xmin,dom.xmax,dom.ymin,dom.ymax],
           origin='lower')
ax1.contour(dom.X, dom.Y, dom.wall_distance, 50, alpha=1.0)
#ax1.imshow(wall_distance,interpolation='nearest',
#           extent=[dom.xmin,dom.xmax,dom.ymin,dom.ymax],
#           alpha=0.7, origin='lower')
#step = 10
#ax1.quiver(dom.X[::step, ::step],
#           dom.Y[::step, ::step],
#           wall_grad_X[::step, ::step],
#           wall_grad_Y[::step, ::step])
#ax1.quiver(dom.X[::step, ::step],
#           dom.Y[::step, ::step],
#           desired_velocity_X[::step, ::step],
#           desired_velocity_Y[::step, ::step])
plt.draw()
plt.show()

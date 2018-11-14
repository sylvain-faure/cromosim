# Authors:
#     Sylvain Faure <sylvain.faure@math.u-psud.fr>
#     Bertrand Maury <bertrand.maury@math.u-psud.fr>
#
#     cromosim/examples/domain/domain_auto_computed.py
#     python domain_auto_computed.py
#
# License: GPL

from cromosim import *
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Circle, Rectangle, Polygon, Arrow
from matplotlib.patches import Circle
from matplotlib.lines import Line2D

## Domain from a background image
dom = Domain(name = 'room_400x300', background = 'room_400x300.png', pixel_size = 0.1)

## Circle( (x_center,y_center), radius )
circle = Circle((20.0,10.0), 1.0)
dom.add_wall(circle)

## Ellipse( (x_center,y_center), width, height, angle_in_degrees_anti-clockwise )
ell = Ellipse( (20.0, 20.0), 3.0, 1.5, 10.0)
dom.add_wall(ell)

## Rectangle( (x_lower_left,y_lower_left), width, height, angle_in_degrees_anti-clockwise)
rect = Rectangle( (10.0, 15.0), 1.0, 4.0, 45)
dom.add_wall(rect)

## Polygon( xy )
poly = Polygon([[30.0, 5.0],  [32.0, 6.0],  [25.0, 8.0]])
dom.add_wall(poly)

## Line2D(xdata, ydata, linewidth)
line = Line2D( [6.0,8.0,7.0,5.0],[15.0,10.0,6.0,8.0], 2)
dom.add_wall(line)

## Line2D(xdata, ydata, linewidth)
line = Line2D( [25.0,29.0],[10.0,10.0], linewidth=2)
dom.add_door(line)
dom.build_domain()
dom.compute_wall_distance()
dom.compute_desired_velocity()
print(dom)

## Plot functions include in vap.py
dom.plot()
dom.plot_wall_dist(id=2)
dom.plot_desired_velocity(id=3)

## Custom plot function
fig = plt.figure(10)
ax1 = fig.add_subplot(111)
ax1.imshow(dom.image,interpolation='nearest',
           extent=[dom.xmin,dom.xmax,dom.ymin,dom.ymax],
           origin='lower')
ax1.contour(dom.X, dom.Y, dom.wall_distance,50,alpha=1.0)
#ax1.imshow(dom.wall_distance,interpolation='nearest',
#           extent=[dom.xmin,dom.xmax,dom.ymin,dom.ymax],
#           alpha=0.7, origin='lower')
#step = 10
#ax1.quiver(dom.X[::step, ::step],dom.Y[::step, ::step],
#           dom.wall_grad_X[::step, ::step],
#           dom.wall_grad_Y[::step, ::step])
plt.draw()

plt.show()

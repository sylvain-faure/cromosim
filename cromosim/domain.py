# Authors:
#     Sylvain Faure <sylvain.faure@math.u-psud.fr>
#     Bertrand Maury <bertrand.maury@math.u-psud.fr>
# License: GPL

import scipy as sp
<<<<<<< HEAD
=======
import numpy as np
>>>>>>> 680862a0447594077e44b73b1b1f0908dfdc3a94
try:
    from scipy.misc import imread
except:
    from imageio import imread
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Circle, Rectangle, Polygon
from matplotlib.lines import Line2D
import PIL
from PIL import Image
from PIL import ImageDraw
import skfmm
import sys

class Destination():

    def __init__(self,
                 name = 'destination',
                 colors = [],
                 excluded_colors = [],
                 gradient_from_color=[],
                 speed=None,
                 next_destination='',
                 next_transit_box=None,
                 next_domain_id = None):
        self.name = name
        self.colors = colors
        self.excluded_colors = excluded_colors
        self.gradient_from_color = gradient_from_color
        self.speed = speed
        self.next_destination = next_destination
        self.next_transit_box = next_transit_box
        self.next_domain_id = next_domain_id
        self.virtual_people_indexes = []
        self.X = None
        self.Y = None
        self.distance = None
        self.desired_velocity_X = None
        self.desired_velocity_Y = None

    def __str__(self):
        return "--> Destination : "+  \
        "\n    name : "+str(self.name) + \
        "\n    colors : "+str(self.colors) + \
        "\n    excluded_colors : "+str(self.excluded_colors) + \
        "\n    gradient_from_color : "+str(self.gradient_from_color) + \
        "\n    next_destination : "+str(self.next_destination)



class Domain():
    """
    To define the computational domain :
        - a background : empty (white) or a PNG image which only \
        contains the colors white, red (for the doors) and black \
        (for the walls)
        - supplementary doors represented by matplotlib shapes : \
        line2D
        - supplementary walls represented by matplotlib shapes : \
        line2D, circle, ellipse, rectangle or polygon

    To compute the obstacle distances and the desired velocities

    Attributes
    ----------

    pixel_size : float
        size of a pixel in meters
    width : int
        width of the background image (number of pixels)
    height : int
        height of the background image (number of pixels)
    xmin : float
        x coordinate of the origin (bottom left corner)
    xmax : float
        xmax = xmin + width*pixel_size
    ymin : float
        y coordinate of the origin (bottom left corner)
    ymax : float
        ymax = ymin + height*pixel_size
    X : numpy array
        x coordinates (meshgrid)
    Y : numpy array
        y coordinates (meshgrid)
    image : numpy array
        pixel array (r,g,b,a)
        The Pillow image is converted to a numpy arrays, then \
        using flipud
        the origin of the image is put it down left instead the \
        top left
    image_red : numpy array
        red values of the image (r,g,b,a)
    image_green : numpy array
        green values of the image (r,g,b,a)
    image_blue : numpy array
        blue values of the image (r,g,b,a)
    wall_mask : numpy array
        boolean array : true for wall pixels
    wall_mask_id : numpy array
        wall pixel indices
    wall_distance : numpy array
        distance (m) to the wall
    wall_grad_X : numpy array
        gradient of the distance to the wall (first component)
    wall_grad_Y : numpy array
        gradient of the distance to the wall (second component)
    door_distance : numpy array
        distance (m) to the door
    desired_velocity_X : numpy array
        opposite of the gradient of the distance to the door : desired velocity \
         (first component)
    desired_velocity_Y : numpy array
        opposite of the gradient of the distance to the door : desired velocity \
        (second component)

    Examples
    --------

    - An example with a domain built using only shape elements :
      cromosim/examples/domain_manually_computed.py

    - An example with a domain built using an image and shape elements :
      cromosim/examples/domain_auto_computed.py
    """

    def __init__(self, name = 'Domain', background = 'White', pixel_size = 1.0,
                 xmin = 0.0, width = 100,
                 ymin = 0.0, height = 100,
                 npzfile = None):
        """
        Constructor of a Domain object

        Parameters
        ----------

        name : string
            domain name (default : 'Domain')
        background : string
            name of the background image (default : 'White', no image)
        pixel_size : float
            size of a pixel in meters (default : 1.0)
        xmin : float
            x coordinate of the origin, bottom left corner (default : 0.0)
        ymin : float
            y coordinate of the origin, bottom left corner (default : 0.0)
        width : int
            width of the background image (default : 100 pixels)
        height : int
            height of the background image (default : 100 pixels)
        npzfile : string
            to build domain from a npz file which contains all variables
        """
<<<<<<< HEAD
        if (npzfile is None):
            self.__shapes = []
            self.__outline_color_shapes = []
            self.__fill_color_shapes = []
            self.__image_filename = ''
            self.name = name
            self.__background = background
            self.destinations = None
            self.pixel_size = pixel_size
            self.xmin, self.ymin = [xmin, ymin]
            if (self.__background != 'White'):
                image = imread(self.__background)
                self.width = image.shape[1]  ## width = number of columns
                self.height = image.shape[0] ## height = number of rows
            else:
                self.width, self.height = [width, height]
            self.xmax = self.xmin + self.width*pixel_size
            self.ymax = self.ymin + self.height*pixel_size
            self.X, self.Y = sp.meshgrid(sp.arange(self.width), sp.arange(self.height))
            self.X = 0.5*self.pixel_size + self.xmin + self.X*self.pixel_size
            self.Y = 0.5*self.pixel_size + self.ymin + self.Y*self.pixel_size
            self.wall_colors = [[0,0,0]] # walls are black by default
            self.wall_mask = None
            self.wall_id = None
            self.wall_distance = None
            self.wall_grad_X = None
            self.wall_grad_Y = None
            if (self.__background != 'White'):
                self.image = Image.open(self.__background)
            else:
                self.image =  Image.new("RGB", (self.width, self.height), "white")
            self.draw = ImageDraw.Draw(self.image)
=======
        self.__walls = []
        self.__doors = []
        self.__image_filename = ''
        self.__name = name
        self.__background = background
        self.pixel_size = pixel_size
        self.xmin, self.ymin = [xmin, ymin]
        if (self.__background != 'White'):
            image = imread(self.__background)
            self.width = image.shape[1]  ## width = number of columns
            self.height = image.shape[0] ## height = number of rows
        else:
            self.width, self.height = [width, height]
        self.xmax = self.xmin + self.width*pixel_size
        self.ymax = self.ymin + self.height*pixel_size
        self.X, self.Y = np.meshgrid(np.arange(self.width), np.arange(self.height))
        self.X = 0.5*self.pixel_size + self.xmin + self.X*self.pixel_size
        self.Y = 0.5*self.pixel_size + self.ymin + self.Y*self.pixel_size
        self.wall_distance = None
        self.wall_grad_X = None
        self.wall_grad_Y = None
        self.desired_velocity_X = None
        self.desired_velocity_Y = None
        self.door_distance = None
>>>>>>> 680862a0447594077e44b73b1b1f0908dfdc3a94

        else:
            data = sp.load(npzfile,allow_pickle=True)
            #print("read shapes : ",data["__shapes"])
            shapes=data["__shapes"].tolist()
            #print("read color shapes : ",data["__color_shapes"])
            color_shapes=data["__color_shapes"].tolist()
            #print("read name : ",data["name"])
            self.name=str(data["name"])
            #print("read __background : ",data["__background"])
            self.__background=str(data["__background"])
            #print("read destinations : ",data["destinations"])
            self.destinations=dict(data["destinations"].tolist())
            #print("read pixel_size : ",data["pixel_size"])
            self.pixel_size=data["pixel_size"]
            #print("read xmin : ",data["xmin"])
            self.xmin=data["xmin"]
            #print("read ymin : ",data["ymin"])
            self.ymin=data["ymin"]
            #print("read xmax : ",data["xmax"])
            self.xmax=data["xmax"]
            #print("read ymax : ",data["ymax"])
            self.ymax=data["ymax"]
            #print("read width : ",data["width"])
            self.width=data["width"]
            #print("read height : ",data["height"])
            self.height=data["height"]
            #print("read X : ",data["X"])
            self.X=data["X"]
            #print("read Y : ",data["Y"])
            self.Y=data["Y"]
            #print("read wall_colors : ",data["wall_colors"])
            self.wall_colors=data["wall_colors"]
            #print("read wall_mask : ",data["wall_mask"])
            self.wall_mask=data["wall_mask"]
            #print("read wall_id : ",data["wall_id"])
            self.wall_id=data["wall_id"]
            #print("read wall_distance : ",data["wall_distance"])
            self.wall_distance=data["wall_distance"]
            #print("read wall_grad_X : ",data["wall_grad_X"])
            self.wall_grad_X=data["wall_grad_X"]
            #print("read wall_grad_Y : ",data["wall_grad_Y"])
            self.wall_grad_Y=data["wall_grad_Y"]
            ## Build again the image with shapes... (which can not be saved in a npz file...)
            if (self.__background != 'White'):
                self.image = Image.open(self.__background)
            else:
                self.image =  Image.new("RGB", (self.width, self.height), "white")
            self.draw = ImageDraw.Draw(self.image)
            # Add shapes
            self.__shapes = []
            self.__outline_color_shapes = []
            self.__fill_color_shapes = []
            for i,shape in enumerate(shapes):
                 self.add_shape(shape,
                    outline_color=outline_color_shapes[i],
                    fill_color=fill_color_shapes[i])
            self.build_domain()

    def save(self, outfile):
        """
        To save the content of the domain in a file

        Parameters
        ----------

        outfile : string
            output filename
        """
        sp.savez(outfile,
        __wall_shapes=self.__wall_shapes,
        __wall_color_shapes=self.__wall_color_shapes,
        name=self.name,
        __background=self.__background,
        pixel_size=self.pixel_size,
        xmin=self.xmin,
        ymin=self.ymin,
        xmax=self.xmax,
        ymax=self.ymax,
        width=self.width,
        height=self.height,
        X=self.X,
        Y=self.Y,
        wall_colors=self.wall_colors,
        wall_mask=self.wall_mask,
        wall_id=self.wall_id,
        wall_distance=self.wall_distance,
        wall_grad_X=self.wall_grad_X,
        wall_grad_Y=self.wall_grad_Y,
        destinations=self.destinations)

    def add_shape(self, shape, outline_color = [0,0,0], fill_color = [255,255,255]):
        """
        To add a matplotlib shape : \
        line2D, circle, ellipse, rectangle or polygon

        Parameters
        ----------

        shape : matplotlib shape
            line2D, circle, ellipse, rectangle or polygon
        """
        self.__shapes.append(shape)
        self.__outline_color_shapes.append(outline_color)
        self.__fill_color_shapes.append(fill_color)
        if ( isinstance(shape, Circle) or isinstance(shape, Ellipse) or
             isinstance(shape, Rectangle) or isinstance(shape, Polygon)):
             xy = shape.get_verts()/self.pixel_size
             xy[:,1] = self.height - xy[:,1]
             self.draw.polygon(sp.around(xy.flatten()).tolist(),
                outline="rgb("+str(outline_color[0])+","+
                               str(outline_color[1])+","+
                               str(outline_color[2])+")",
                fill="rgb("+str(fill_color[0])+","+
                            str(fill_color[1])+","+
                            str(fill_color[2])+")")
             linewidth = shape.get_linewidth()
             self.draw.line(sp.around(xy.flatten()).tolist(),
                width=int(linewidth),
                fill="rgb("+str(outline_color[0])+","+
                           str(outline_color[1])+","+
                           str(outline_color[2])+")")
        elif  isinstance(shape, Line2D):
             linewidth = shape.get_linewidth()
             xy = shape.get_xydata()/self.pixel_size
             xy[:,1] = self.height - xy[:,1]
             self.draw.line(sp.around(xy.flatten()).tolist(),
                width=int(linewidth),
                fill="rgb("+str(outline_color[0])+","+
                            str(outline_color[1])+","+
                            str(outline_color[2])+")")

    # def add_door_shape(self, shape, color = [255,0,0]):
    #     """
    #     To add a door represented by matplotlib shapes : \
    #     line2D (only)
    #
    #     Parameters
    #     ----------
    #
    #     shape : matplotlib shape
    #         line2D
    #     """
    #     self.__door_shapes.append(shape)
    #     self.__door_color_shapes.append(color)
    #     linewidth = shape.get_linewidth()
    #     xy = shape.get_xydata()/self.pixel_size
    #     xy[:,1] = self.height - xy[:,1]
    #     self.draw.line(sp.around(xy.flatten()).tolist(), width=int(linewidth),
    #         fill="rgb("+str(color[0])+","+str(color[1])+","+str(color[2])+")")

    def build_domain(self):
        """
        To build the domain : reads the background image (if supplied) \
        and initializes all the color arrrays
        """
<<<<<<< HEAD
        self.__image_filename = self.name+'_domain.png'
        self.image.save(self.__image_filename)
=======
        if (self.__background != 'White'):
            image = Image.open(self.__background)
        else:
            image =  Image.new("RGB", (self.width, self.height), "white")
        draw = ImageDraw.Draw(image)
        for iw in self.__walls:
            if ( isinstance(iw, Circle) or isinstance(iw, Ellipse) or
                 isinstance(iw, Rectangle) or isinstance(iw, Polygon)):
                 xy = iw.get_verts()/self.pixel_size
                 xy[:,1] = self.height - xy[:,1]
                 draw.polygon(np.around(xy.flatten()).tolist(),
                              outline="rgb(0, 0, 0)", fill="rgb(0, 0, 0)")
            elif  isinstance(iw, Line2D):
                 linewidth = iw.get_linewidth()
                 xy = iw.get_xydata()/self.pixel_size
                 xy[:,1] = self.height - xy[:,1]
                 draw.line(np.around(xy.flatten()).tolist(), width=int(linewidth),
                           fill="rgb(0, 0, 0)")
        for id in self.__doors:
            linewidth = id.get_linewidth()
            xy = id.get_xydata()/self.pixel_size
            xy[:,1] = self.height - xy[:,1]
            draw.line(np.around(xy.flatten()).tolist(), width=int(linewidth),
                      fill="rgb(255, 0, 0)")
        self.__image_filename = self.__name+'_domain.png'
        image.save(self.__image_filename)
>>>>>>> 680862a0447594077e44b73b1b1f0908dfdc3a94
        ## Easy way to convert a Pillow image to numpy arrays...
        ## The origin of img is at the top (left) and flipud allows to put it down...
        self.image = np.flipud(imread(self.__image_filename))
        self.image_red   = self.image[:,:,0]
        self.image_green = self.image[:,:,1]
        self.image_blue  = self.image[:,:,2]
<<<<<<< HEAD
        self.wall_mask = sp.zeros_like(self.image_red)
        for c in self.wall_colors:
            self.wall_mask += (self.image_red == c[0]) \
                        *(self.image_green == c[1]) \
                        *(self.image_blue == c[2])
        self.wall_id = sp.where( self.wall_mask>0 )
        ## Compute wall distances : wall = "wall_colors" pixels
        if (self.wall_id[0].size>0):
            self.compute_wall_distance()
        else:
            print("WARNING : Failed to compute wall distance!")
            print("WARNING : Wall colors are ",self.wall_colors)
            print("WARNING : Check that there are pixels with these colors!")
            sys.exit()
=======
        self.mask =  (self.image_red == 0) \
                    *(self.image_green == 0) \
                    *(self.image_blue == 0)
        self.mask_id = np.where( self.mask )
>>>>>>> 680862a0447594077e44b73b1b1f0908dfdc3a94

    def compute_wall_distance(self):
        """
        To compute the geodesic distance to the walls in using \
        a fast-marching method
        """
<<<<<<< HEAD
        phi = sp.ones(self.image_red.shape)
        if (len(self.wall_id[0])>0):
            phi[self.wall_id] = 0
=======
        phi = np.ones(self.image_red.shape)
        if (len(self.mask_id[0])>0):
            phi[self.mask_id] = 0
>>>>>>> 680862a0447594077e44b73b1b1f0908dfdc3a94
            self.wall_distance = skfmm.distance(phi, dx=self.pixel_size)
            grad = np.gradient(self.wall_distance,edge_order=2)
            grad_X = grad[1]/self.pixel_size
            grad_Y = grad[0]/self.pixel_size
            norm = np.sqrt(grad_X**2+grad_Y**2)
            norm = (norm>0)*norm+(norm==0)*0.001
            self.wall_grad_X = grad_X/norm
            self.wall_grad_Y = grad_Y/norm
        else:
            self.wall_distance = 1.0e99*np.ones(self.image_red.shape)

    def add_destination(self, dest):
        """
        To compute the geodesic distance to the doors in using \
        a fast-marching method. The opposite of the gradient of this distance
        corresponds to the desired velocity which permits to reach the closest
        door

        Returns
        -------

        door_distance : numpy array
            distance to the closest door
        desired_velocity_X : numpy array
            opposite of the gradient of the door distance, x component
        desired_velocity_Y : numpy array
            opposite of the gradient of the door distance, y component
        """
<<<<<<< HEAD

        ## Compute desired velocities to this destination, in this domain...

        # Compute mask for dest.excluded_colors (which can be for example the
        # wall colors...)
        excluded_color_mask = sp.zeros_like(self.image_red)
        for c in dest.excluded_colors:
            excluded_color_mask += (self.image_red == c[0]) \
                                  *(self.image_green == c[1]) \
                                  *(self.image_blue == c[2])
        excluded_color_id = sp.where( excluded_color_mask>0 )

        # Compute mask for dest.colors (which correspond which correspond to
        # the destination of the people)
        dest_mask = sp.zeros_like(self.image_red)
        for ic,rgb in enumerate(dest.colors):
            dest_mask = (self.image_red == rgb[0]) \
                       *(self.image_green == rgb[1]) \
                       *(self.image_blue == rgb[2])
        mask_id = sp.where( (dest_mask>=1) )

        # Define rhs of the Eikonal equation (i.e. the speed)
        if (dest.speed is not None):
            if (dest.speed.shape != self.image_red.shape):
                print("Bad speed shape ! Failed to compute the destination distance...")
                sys.exit()
        else:
            dest.speed = sp.ones_like(self.image_red)


        if (mask_id[0].size>0):
            dest_mask[mask_id] = True
            phi = sp.ones(self.image_red.shape)
            phi[mask_id] = 0

            phi = sp.ma.MaskedArray(phi, mask=excluded_color_mask)
            dest.distance = skfmm.travel_time(phi, dest.speed, dx=self.pixel_size)
            #dest.distance = skfmm.distance(phi, dx=self.pixel_size)
            if (excluded_color_id[0].size>0):
                tmp_dist = dest.distance.filled(9999)
            else:
                tmp_dist = dest.distance
            grad = sp.gradient(tmp_dist,edge_order=2)
        else:
            dest.distance = -sp.ones_like(self.image_red)
            grad = sp.gradient(sp.zeros_like(self.image_red),edge_order=2)

        for l,rgbgrad in enumerate(dest.gradient_from_color):
            test = (self.image_red == int(rgbgrad[0])) \
                  *(self.image_green == int(rgbgrad[1])) \
                  *(self.image_blue == int(rgbgrad[2]))
            indices = sp.where( test == True )
            grad[1][indices] = -rgbgrad[3]
            grad[0][indices] = -rgbgrad[4]

=======
        mask_red = (self.image_red == 255) \
                  *(self.image_green == 0) \
                  *(self.image_blue == 0)
        ind_red = np.where( mask_red )
        phi = np.ones(self.image_red.shape)
        phi[ind_red] = 0
        phi = np.ma.MaskedArray(phi, mask=self.mask)
        self.door_distance = skfmm.distance(phi, dx=self.pixel_size)
        tmp_dist = self.door_distance.filled(9999)
        grad = np.gradient(tmp_dist,edge_order=2)
>>>>>>> 680862a0447594077e44b73b1b1f0908dfdc3a94
        grad_X = -grad[1]/self.pixel_size
        grad_Y = -grad[0]/self.pixel_size
        norm = np.sqrt(grad_X**2+grad_Y**2)
        norm = (norm>0)*norm+(norm==0)*0.001
        dest.desired_velocity_X = grad_X/norm
        dest.desired_velocity_Y = grad_Y/norm
        try:
            self.destinations[dest.name] = dest
        except:
            self.destinations = { dest.name: dest }

    def people_desired_velocity(self, people, dest_name='default', I=None, J=None):
        """
        This function determines people desired velocities from the desired \
        velocity array computed by Domain thanks to a fast-marching method.

        Parameters
        ----------
        people: numpy array
            people coordinates and radius : x,y,r
        dest_name: string
            name of people destination
        I : numpy array (None by default)
            people index i
        J : numpy array (None by default)
            people index j

        Returns
        -------
        I : numpy array
            people index i
        J : numpy array
            people index j
        Vd : numpy array
            people desired velocity
        """
        if ((I is None) or (J is None)):
            I = sp.floor((people[:,1]-self.ymin-0.5*self.pixel_size)/self.pixel_size).astype(int)
            J = sp.floor((people[:,0]-self.xmin-0.5*self.pixel_size)/self.pixel_size).astype(int)
        Vd = sp.zeros( (people.shape[0],2) )
        Vd[:,0] = self.destinations[dest_name].desired_velocity_X[I,J]
        Vd[:,1] = self.destinations[dest_name].desired_velocity_Y[I,J]
        return I,J,Vd

    def people_target_distance(self, people, dest_name='default', I = None, J = None):
        """
        This function determines distances to the current target for all people

        Parameters
        ----------
        people: numpy array
            people coordinates and radius : x,y,r
        dest_name: string
            name of people destination
        I : numpy array (None by default)
            people index i
        J : numpy array (None by default)
            people index j
        Returns
        -------
        I : numpy array
            people index i
        J : numpy array
            people index j
        D : numpy array
            distances to the current target
        """

        if ((I is None) or (J is None)):
            I = sp.floor((people[:,1]-self.ymin-0.5*self.pixel_size)/self.pixel_size).astype(int)
            J = sp.floor((people[:,0]-self.xmin-0.5*self.pixel_size)/self.pixel_size).astype(int)
        D = self.destinations[dest_name].distance[I,J]-people[:,2]
        return I,J,D


    def people_wall_distance(self, people, I = None, J = None):
        """
        This function determines distances to the nearest wall for all people

        Parameters
        ----------
        people: numpy array
            people coordinates and radius : x,y,r
        I : numpy array (None by default)
            people index i
        J : numpy array (None by default)
            people index j

        Returns
        -------
        I : numpy array
            people index i
        J : numpy array
            people index j
        D : numpy array
            distances to the nearest wall
        """
        if ((I is None) or (J is None)):
            I = sp.floor((people[:,1]-self.ymin-0.5*self.pixel_size)/self.pixel_size).astype(int)
            J = sp.floor((people[:,0]-self.xmin-0.5*self.pixel_size)/self.pixel_size).astype(int)
        D = self.wall_distance[I,J]-people[:,2]
        return I,J,D

    def people_update_destination(self, people, people_destination, people_domain_id, dest_name='default'):

        I = sp.floor((people[:,1]-self.ymin-0.5*self.pixel_size)/self.pixel_size).astype(int)
        J = sp.floor((people[:,0]-self.xmin-0.5*self.pixel_size)/self.pixel_size).astype(int)
        D = self.destinations[dest_name].distance[I,J]-people[:,2]
        ind = sp.where(D<=-0.01)[0]
        #print(self.name," : dest = ",dest_name," change destination for ind = ",ind)
        if (ind.shape[0]>0) and (self.destinations[dest_name].next_domain_id is not None):
            people_destination[ind] = self.destinations[dest_name].next_destination
            people_domain_id[ind] = self.destinations[dest_name].next_domain_id
        return people,people_destination,people_domain_id


    def plot(self,id=1,dpi=150):
        """
        To plot the computational domain

        Parameters
        ----------

        id : integer
            Figure id (number)
        dpi : integer
            Figure resolution
        """
        fig = plt.figure(id)
        ax1 = fig.add_subplot(111)
        ax1.imshow(self.image,interpolation='nearest',extent=[self.xmin,self.xmax,
                   self.ymin,self.ymax], origin='lower')
        #plt.savefig('.png',dpi=dpi)
        # ax1.axes.get_xaxis().set_visible(False)
        # ax1.axes.get_xaxis().set_visible(False)
        ax1.set_axis_off()
        #fig.add_axes(ax1)
        plt.draw()

    def plot_wall_dist(self,
                       step=10,
                       scale=10,
                       scale_units='inches',
                       id=1,
                       dpi=150):
        """
        To plot the wall distances

        Parameters
        ----------

        id : integer
            Figure id (number)
        dpi : integer
            Figure resolution
        """
        fig = plt.figure(id)
        ax1 = fig.add_subplot(111)
        ax1.imshow(self.image,interpolation='nearest',
                   extent=[self.xmin,self.xmax,self.ymin,self.ymax], origin='lower')
        ax1.imshow(self.wall_distance,interpolation='nearest',
                   extent=[self.xmin,self.xmax,self.ymin,self.ymax],alpha=0.7,
                   origin='lower')
        ax1.quiver(self.X[::step, ::step],self.Y[::step, ::step],
                   self.wall_grad_X[::step, ::step],
                   self.wall_grad_Y[::step, ::step],
                   scale=scale, scale_units=scale_units)
        ax1.set_axis_off()
        #plt.savefig('.png',dpi=dpi)
        plt.draw()

    def plot_desired_velocity(self,
                              destination_name,
                              step=10,
                              scale=10,
                              scale_units='inches',
                              id=1,
                              dpi=150):
        """
        To plot the desired velocity

        Parameters
        ----------
        destination_name : string
            name of the considered destination
        step : integer
            draw an arrow every step pixels
        scale : integer
            scaling for the quiver arrows
        scale_units : string
            unit name for quiver arrow scaling
        id : integer
            Figure id (number)
        dpi : integer
            Figure resolution
        """
        fig = plt.figure(id)
        ax1 = fig.add_subplot(111)
        ax1.imshow(self.image,interpolation='nearest',
                   extent=[self.xmin,self.xmax,self.ymin,self.ymax], origin='lower')
        ax1.imshow(self.destinations[destination_name].distance,interpolation='nearest',
                   extent=[self.xmin,self.xmax,self.ymin,self.ymax],alpha=0.7,
                   origin='lower')
        ax1.quiver(self.X[::step, ::step],self.Y[::step, ::step],
                   self.destinations[destination_name].desired_velocity_X[::step, ::step],
                   self.destinations[destination_name].desired_velocity_Y[::step, ::step],
                   scale=scale, scale_units=scale_units)
        #plt.savefig('.png',dpi=dpi)
        ax1.set_axis_off()
        plt.draw()

    def __str__(self):
        """
        To print the main caracteristics of a Domain object
        """
        return "--> "+self.name+  \
        " :\n    dimensions : ["+str(self.xmin)+","+str(self.xmax)+"]x["+ \
                               str(self.ymin)+","+str(self.ymax)+"]"+ \
        "\n    width : "+str(self.width)+" height : "+str(self.height)+ \
        "\n    background image : "+str(self.__background)+ \
        "\n    image of the domain : "+str(self.__image_filename)+ \
        "\n    walls : "+str(self.__wall_shapes)+ \
        "\n    doors : "+str(self.__door_shapes)+ \
        "\n    destinations : "+str(self.destinations)

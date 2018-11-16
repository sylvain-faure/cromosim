# Authors:
#     Sylvain Faure <sylvain.faure@math.u-psud.fr>
#     Bertrand Maury <bertrand.maury@math.u-psud.fr>
# License: GPL

import scipy as sp
from scipy.misc import imread
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Circle, Rectangle, Polygon
from matplotlib.lines import Line2D
import PIL
from PIL import Image
from PIL import ImageDraw
import skfmm

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
    mask : numpy array
        boolean array : true for black pixels
    mask_id : numpy array
        black pixel indices
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
                 ymin = 0.0, height = 100):
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
        """
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
        self.X, self.Y = sp.meshgrid(sp.arange(self.width), sp.arange(self.height))
        self.X = 0.5*self.pixel_size + self.xmin + self.X*self.pixel_size
        self.Y = 0.5*self.pixel_size + self.ymin + self.Y*self.pixel_size
        self.wall_distance = None
        self.wall_grad_X = None
        self.wall_grad_Y = None
        self.desired_velocity_X = None
        self.desired_velocity_Y = None
        self.door_distance = None

    def add_wall(self, shape):
        """
        To add a wall represented by matplotlib shapes : \
        line2D, circle, ellipse, rectangle or polygon

        Parameters
        ----------

        shape : matplotlib shape
            line2D, circle, ellipse, rectangle or polygon
        """
        self.__walls.append(shape)

    def add_door(self, shape):
        """
        To add a door represented by matplotlib shapes : \
        line2D (only)

        Parameters
        ----------

        shape : matplotlib shape
            line2D
        """
        self.__doors.append(shape)

    def build_domain(self):
        """
        To build the domain : reads the background image (if supplied) \
        and initializes all the color arrrays
        """
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
                 draw.polygon(sp.around(xy.flatten()).tolist(),
                              outline="rgb(0, 0, 0)", fill="rgb(0, 0, 0)")
            elif  isinstance(iw, Line2D):
                 linewidth = iw.get_linewidth()
                 xy = iw.get_xydata()/self.pixel_size
                 xy[:,1] = self.height - xy[:,1]
                 draw.line(sp.around(xy.flatten()).tolist(), width=int(linewidth),
                           fill="rgb(0, 0, 0)")
        for id in self.__doors:
            linewidth = id.get_linewidth()
            xy = id.get_xydata()/self.pixel_size
            xy[:,1] = self.height - xy[:,1]
            draw.line(sp.around(xy.flatten()).tolist(), width=int(linewidth),
                      fill="rgb(255, 0, 0)")
        self.__image_filename = self.__name+'_domain.png'
        image.save(self.__image_filename)
        ## Easy way to convert a Pillow image to numpy arrays...
        ## The origin of img is at the top (left) and flipud allows to put it down...
        self.image = sp.flipud(imread(self.__image_filename))
        self.image_red   = self.image[:,:,0]
        self.image_green = self.image[:,:,1]
        self.image_blue  = self.image[:,:,2]
        self.mask =  (self.image_red == 0) \
                    *(self.image_green == 0) \
                    *(self.image_blue == 0)
        self.mask_id = sp.where( self.mask )

    def compute_wall_distance(self):
        """
        To compute the geodesic distance to the walls in using \
        a fast-marching method
        """
        phi = sp.ones(self.image_red.shape)
        if (len(self.mask_id[0])>0):
            phi[self.mask_id] = 0
            self.wall_distance = skfmm.distance(phi, dx=self.pixel_size)
            grad = sp.gradient(self.wall_distance,edge_order=2)
            grad_X = grad[1]/self.pixel_size
            grad_Y = grad[0]/self.pixel_size
            norm = sp.sqrt(grad_X**2+grad_Y**2)
            norm = (norm>0)*norm+(norm==0)*0.001
            self.wall_grad_X = grad_X/norm
            self.wall_grad_Y = grad_Y/norm
        else:
            self.wall_distance = 1.0e99*sp.ones(self.image_red.shape)

    def compute_desired_velocity(self):
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
        mask_red = (self.image_red == 255) \
                  *(self.image_green == 0) \
                  *(self.image_blue == 0)
        ind_red = sp.where( mask_red )
        phi = sp.ones(self.image_red.shape)
        phi[ind_red] = 0
        phi = sp.ma.MaskedArray(phi, mask=self.mask)
        self.door_distance = skfmm.distance(phi, dx=self.pixel_size)
        tmp_dist = self.door_distance.filled(9999)
        grad = sp.gradient(tmp_dist,edge_order=2)
        grad_X = -grad[1]/self.pixel_size
        grad_Y = -grad[0]/self.pixel_size
        norm = sp.sqrt(grad_X**2+grad_Y**2)
        norm = (norm>0)*norm+(norm==0)*0.001
        self.desired_velocity_X = grad_X/norm
        self.desired_velocity_Y = grad_Y/norm
        return self.door_distance, self.desired_velocity_X, self.desired_velocity_Y

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
        plt.draw()

    def plot_wall_dist(self,id=1,dpi=150):
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
        step = 10
        ax1.quiver(self.X[::step, ::step],self.Y[::step, ::step],
                   self.wall_grad_X[::step, ::step],
                   self.wall_grad_Y[::step, ::step])
        #plt.savefig('.png',dpi=dpi)
        plt.draw()

    def plot_desired_velocity(self,id=1,dpi=150):
        """
        To plot the desired velocity

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
        ax1.imshow(self.door_distance,interpolation='nearest',
                   extent=[self.xmin,self.xmax,self.ymin,self.ymax],alpha=0.7,
                   origin='lower')
        step = 10
        ax1.quiver(self.X[::step, ::step],self.Y[::step, ::step],
                   self.desired_velocity_X[::step, ::step],
                   self.desired_velocity_Y[::step, ::step])
        #plt.savefig('.png',dpi=dpi)
        plt.draw()

    def __str__(self):
        """
        To print the main caracteristics of a Domain object
        """
        return "--> "+self.__name+  \
        " :\n    dimensions : ["+str(self.xmin)+","+str(self.xmax)+"]x["+ \
                               str(self.ymin)+","+str(self.ymax)+"]"+ \
        "\n    width : "+str(self.width)+" height : "+str(self.height)+ \
        "\n    background image : "+str(self.__background)+ \
        "\n    image of the domain : "+str(self.__image_filename)+ \
        "\n    walls : "+str(self.__walls)+ \
        "\n    doors : "+str(self.__doors)

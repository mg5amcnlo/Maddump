from __future__ import absolute_import
from __future__ import print_function
from copy import copy, deepcopy
import numpy as np
from scipy.interpolate import splrep,splint,interp1d
from scipy.optimize import bisect
from random import randint, random
import multiprocessing
import os 
from six.moves import range
from six.moves import zip

pjoin = os.path.join

from . import fit2D_card as fit2D
import madgraph.various.banner as banner_mod
import madgraph.various.lhe_parser as lhe_parser

def plot (infile,outfile=None):
    """ Plot the canvas with the 2D mesh."""
    import matplotlib
#    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Rectangle

    x,y,width,height = np.loadtxt(infile, unpack = True)
    
    fig, ax = plt.subplots()

    a = 0.
    b = 0.

    min_x = min(x)
    min_y = min(y)
    rectangles = []
    for i in range(len(x)):
        rectangles.append(Rectangle((x[i], y[i]),
                                    width[i],height[i],fill=False))
        # if y[i]==min_y:
        #     a += width[i]
        # if x[i]==min_x:
        #     b += height[i]
        tmp_x=x[i]+width[i]
        if tmp_x > a:
            a = tmp_x
        tmp_y=y[i]+height[i]
        if tmp_y > b:
            b = tmp_y
        
    pc = PatchCollection(rectangles, facecolor='w', edgecolor='b' )
        
    ax.add_collection(pc)

    plt.xlim((min_x,a))
    plt.ylim((min_y,b))

    if not outfile:
        plt.show()
    else:
        plt.savefig(outfile)

    plt.close('all')
#===============================================================================
# Basic Geometrical Point in 2-dimensions
#===============================================================================
class Point (object):
    
    def __init__(self,x=0,y=0):
        self.x = x
        self.y = y
    
    def __str__ (self):
        return "({0}, {1})".format(self.x, self.y)

    def __add__(self,other):
        """ Overwrite the sum operator by applying the
        sum component-by-component. The result is a new point."""
        return Point(self.x+other.x, self.y+other.y)

    def exchange_xy (self):
        """ Exchange components: x<->y."""
        self.x, self.y = self.y, self.x

    def midpoint(self,other):
        """ Calculate the middle point of the point with another one.
        The result is a new point."""
        sum = self+other
        return Point(0.5*sum.x,0.5*sum.y)

# overwrite the inequality relations
# ordering priority:  first in x variable and in the y one

    def cmp (self,other):
        if self.x > other.x: return 1
        if self.x < other.x: return -1
        if self.y > other.y: return 1
        if self.y < other.y: return -1

    def __eq__(self, other):
        return self.cmp(other) == 0
    def __le__(self, other):
        return self.cmp(other) <= 0
    def __ge__(self, other):
        return self.cmp(other) >= 0
    def __gt__(self, other):
        return self.cmp(other) > 0
    def __lt__(self, other):
        return self.cmp(other) < 0
    def __ne__(self, other):
        return self.cmp(other) != 0
        

#===============================================================================
# Point with weight
#===============================================================================
class WeightedPoint (Point):

    def __init__(self,x=0,y=0,weight=1):
        self.weight = weight
        super(WeightedPoint,self).__init__(x,y)

    def isincell(self,cell):
        """ Return true if the coordinates are within the given cell."""
        if self.x < cell.corner.x or \
           self.x > cell.corner.x + cell.width : return False
        if self.y < cell.corner.y or \
           self.y > cell.corner.y + cell.height : return False
        return True


#===============================================================================
# Cell class
#===============================================================================
class Cell (object):
    """Class to define the basic element of the 2Dmesh."""
    
    def __init__(self,pos,width,height,x=0):
        self.corner = pos
        self.width = width
        self.height = height
        
        self.weight = 0

        if x == 0:
            self.pts = []
            self.npts = 0
        else:
            self.pts = []
            for pt in x:
                if pt.isincell(self):
                    self.pts.append(pt)
                    self.weight += pt.weight
            self.npts = len(self.pts)
            self.weight

    def addpts(self,x):
        """ Fill the cell with points given as list of WeightedPoints.
        Only points whose coordinates are within the cell are added."""
        for pt in x:
            if pt.isincell(self):
                self.pts.append(pt)
                self.weight += pt.weight 
        self.npts += len(self.pts)
        

    def error(self):
        """ Compute the total error from the weight of each point."""
        tmp = 0 
        for pt in self.pts:
            tmp += pt.weight**2
        return np.sqrt(tmp)
        
    def split_vertically(self,xsplit,x1=0,x2=0):
        """ Divide a cell into two new cells with a vertical line 
        starting at a given x coordinate. If the mother cell contains
        some points, they will be assigned to the daughter cells."""
        if x1 == 0 and x2 == 0:
            lcell = Cell(self.corner,xsplit-self.corner.x,self.height)
            rcell = Cell(Point(xsplit,self.corner.y),
                         self.width-xsplit+self.corner.x,self.height)
            x1=[]
            x2=[]
            for pt in self.pts:
                if pt.isincell(lcell): 
                    x1.append(pt)
                else:
                    x2.append(pt)
            lcell.addpts(x1)
            rcell.addpts(x2)
        else:
            lcell = Cell(self.corner,xsplit-self.corner.x,self.height,x1)
            rcell = Cell(Point(xsplit,self.corner.y),
                         self.width-xsplit+self.corner.x,self.height,x2)
        return lcell, rcell

    
    def split_horizontally(self,ysplit,x1=0,x2=0):
        """ Divide a cell into two new cells with a horizontal line 
        starting at a given y coordinate. If the mother cell contains
        some points, they will be assigned to the daughter cells."""
        if x1 == 0 and x2 == 0:
            dcell = Cell(self.corner,self.width,ysplit-self.corner.y)
            ucell = Cell(Point(self.corner.x,ysplit),self.width,
                         self.height-ysplit+self.corner.y)
            x1=[]
            x2=[]
            for pt in self.pts:
                if pt.isincell(dcell): 
                    x1.append(pt)
                else:
                    x2.append(pt)
            dcell.addpts(x1)
            ucell.addpts(x2)
        else:
            dcell = Cell(self.corner,self.width,ysplit-self.corner.y,x1)
            ucell = Cell(Point(self.corner.x,ysplit),self.width,
                         self.height-ysplit+self.corner.y,x2)
        return dcell, ucell


    def multi_split_vertically(self,isplit,xsplit):
        """ Divide a cell into  new cells with a vertical line 
        starting at a given x coordinate. If the mother cell contains
        some points, they will be assigned to the daughter cells."""
        self.pts.sort()
        subcell = []
        
        corner = self.corner
        x=self.pts[:isplit[1]-1]

        for i in range(len(isplit)):
            subcell.append(Cell(corner,xsplit[i]-corner.x,self.height,x))
            corner = Point(xsplit[i],corner.y)
            if i != len(isplit)-1:
                x=self.pts[isplit[i]-1:isplit[i+1]-1]
                
        x=self.pts[isplit[-1]-1:]
        subcell.append(Cell(corner,self.corner.x+self.width-xsplit[-1],self.height,x))
            
        return subcell
        
    def multi_split_horizontally(self,isplit,ysplit):
        """ Divide a cell into new cells with a vertically line 
        starting at a given x coordinate. If the mother cell contains
        some points, they will be assigned to the daughter cells."""

        for pt in self.pts:
            pt.exchange_xy()
        self.pts.sort()
        for pt in self.pts:
            pt.exchange_xy()

        subcell = []
        
        corner = self.corner
        x=self.pts[:isplit[1]-1]

        for i in range(len(isplit)):
            subcell.append(Cell(corner,self.width,ysplit[i]-corner.y,x))
#            subcell.append(Cell(corner,xsplit[i]-corner.x,self.height,x))
            corner = Point(corner.x,ysplit[i])
#            corner = Point(xsplit[i],corner.y)
            if i != len(isplit)-1:
                x=self.pts[isplit[i]-1:isplit[i+1]-1]
                
        x=self.pts[isplit[-1]-1:]
        subcell.append(Cell(corner,self.width,self.corner.y+self.height-ysplit[-1],x))
#        subcell.append(Cell(corner,self.corner.x+self.width-xsplit[-1],self.height,x))
            
        return subcell


    
    def __str__(self):
        s = "Cell : corner={0}, width={1}, height={2}".format(self.corner,
                    self.width, self.height)
        if self.npts==0:
            s = s +"\n It is empty"
        else:
            s = s + "\n It contains {0} points: ".format(self.npts)
            for pt in self.pts:
                s = s + "\n {0}".format(pt)
        return s

#===============================================================================
# 2D Histogram as collection of Cells 
#===============================================================================
class CellHistogram(object):
    """Class to define the 2Dmesh as a collection of cells fitted starting 
    from the data sample."""
    
    def __init__(self, pos, width, height, npts_exit=25 ,suffix=None):
        self.cells = []
        self.intervals_x = []        
        self.intervals_y = []        
        self.canvas = Cell(pos,width,height)
        self.npts = 0
        self.weight = 0
        self.width=width
        self.height=height
        self.npts_exit=npts_exit
        if not suffix:
            self.suffix=''
        else:
            self.suffix='-'+str(suffix)
            
        
    def add_pts(self,x):
        """ Fill the canvas with points given as list of WeightedPoints.
        Only points whose coordinates are within the canvas region 
        are added."""
        if self.npts > 0:
            print("Error: Histogram is not empty!")
        elif False in [isinstance(x[i],WeightedPoint) for i in range(len(x))]:
            print("Input Error: Points must be WeightedPoint objects!")
        else:
            self.canvas.addpts(x)
            self.npts = self.canvas.npts
            self.weight = self.canvas.weight

    def write_cell(self,cell,file):
        
        # check if the cell is almost empty
        corner_x = cell.corner.x
        corner_y = cell.corner.y
        width = cell.width
        height = cell.height


        peripheral = False
        # check whether the cell is on the boundary or not
        
        fac=0.98
        epstol = 0.001
        wgt = cell.weight
        cell.pts.sort()
        if self.canvas.corner.x != 0.:
            if( abs(corner_x /self.canvas.corner.x -1.) < epstol ):
                peripheral = True
                weight=0
                for pt in cell.pts[::-1]:
                    weight += pt.weight
                    if (weight > fac*wgt):
                        isplit = cell.pts.index(pt)
                        break
                width = width -(cell.pts[isplit+1].x-corner_x)
                corner_x = cell.pts[isplit+1].x

        if(abs( (corner_x+width)/(self.canvas.corner.x+self.canvas.width) -1.) < epstol):
            peripheral = True
            weight=0
            for pt in cell.pts:
                weight += pt.weight
                if (weight > fac*wgt):
                    isplit = cell.pts.index(pt)
                    break
            width = cell.pts[isplit-1].x-corner_x
            
        for pt in cell.pts:
            pt.exchange_xy()
        cell.pts.sort()
        for pt in cell.pts:
            pt.exchange_xy()
            
        if self.canvas.corner.y != 0.:
            if( abs( corner_y/self.canvas.corner.y -1. ) < epstol ):
                peripheral = True
                weight=0
                for pt in cell.pts[::-1]:
                    weight += pt.weight
                    if (weight > fac*wgt):
                        isplit = cell.pts.index(pt)
                        break
                height = height -(cell.pts[isplit+1].y-corner_y)
                corner_y = cell.pts[isplit+1].y
            
        if( abs( (corner_y+height)/(self.canvas.corner.y+self.canvas.height) -1. ) < epstol) :
            peripheral = True
            weight=0
            for pt in cell.pts:
                weight += pt.weight
                if (weight > fac*wgt):
                    isplit = cell.pts.index(pt)
                    break
            height = cell.pts[isplit-1].y-corner_y

        if not peripheral:
            fac=0.98
            cell.pts.sort()
            wgt = cell.weight
            weight=0
            for pt in cell.pts:
                weight += pt.weight
                if (weight > fac*wgt):
                    isplit = cell.pts.index(pt)
                    break
            while (True):
                if( cell.pts[0].x-corner_x < width/2. ):
                    if ( cell.pts[isplit-1].x-corner_x < width/2. ):
                        width=width/2.
                        continue
                if( cell.pts[1].x-corner_x > width/2. ):
                    corner_x = corner_x + width/2.
                    width=width/2.
                    continue
                break

            for pt in cell.pts:
                pt.exchange_xy()
            cell.pts.sort()
            for pt in cell.pts:
                pt.exchange_xy()
            #wgt = cell.weight
            weight=0
            for pt in cell.pts:
                weight += pt.weight
                if (weight > fac*wgt):
                    isplit = cell.pts.index(pt)
                    break
            while (True):
                if( cell.pts[0].y-corner_y < height/2. ):
                    if ( cell.pts[isplit-1].y-corner_y < height/2. ):
                        height = height/2.
                        continue
                if( cell.pts[1].y-corner_y > height/2. ):
                    corner_y = corner_y + height/2.
                    height = height/2.
                    continue
                break
        

        s = str(corner_x) + "\t" + str(corner_y) + "\t" + \
            str(width) + "\t" + str(height) + "\n" 
        file.write(s)

        
    def fit(self,fit_type='cell',start_direction='horizontally',ncores=1):
        """ Fit 2Dcell/1Dinterval mesh according to data points distribution. """
        method_name = 'fit_' + str(fit_type)
        fit_method = getattr(self, method_name)

        local_canvas = deepcopy(self.canvas)
        canvas=[]
        tmp_canvas = []

        nstart=1
        self.flag_inv_order=False
        
        if (fit_type=='cell'):
            if start_direction == 'horizontally':
                split = [self.equalweight_split_horizontally,self.equalweight_split_vertically]
            elif start_direction == 'vertically':
                split = [self.equalweight_split_vertically,self.equalweight_split_horizontally]
            
            if ncores in [2**j for j in range(1,10)]:
                canvas.append(local_canvas)
                n = int(np.log(ncores)/np.log(2))
                if (n%2==0):
                    for i in range(n/2):
                        for cell in canvas:
                            cell1,cell2 = split[0](cell)
                            tmp_canvas.append(cell1)
                            tmp_canvas.append(cell2)
                            
                        canvas = deepcopy(tmp_canvas)
                        tmp_canvas= []
                            
                        for cell in canvas:
                            cell1,cell2 = split[1](cell)
                            tmp_canvas.append(cell1)
                            tmp_canvas.append(cell2)
                            
                        canvas = deepcopy(tmp_canvas)
                        tmp_canvas= []

                else:
                    cell1,cell2 = split[0](local_canvas)
                    tmp_canvas.append(cell1)
                    tmp_canvas.append(cell2)

                    canvas = deepcopy(tmp_canvas)
                    tmp_canvas= []

                    for i in range((n-1)/2):
                        for cell in canvas:
                            cell1,cell2 = split[1](cell)
                            tmp_canvas.append(cell1)
                            tmp_canvas.append(cell2)
                            
                        canvas = deepcopy(tmp_canvas)
                        tmp_canvas= []
                            
                        for cell in canvas:
                            cell1,cell2 = split[0](cell)
                            tmp_canvas.append(cell1)
                            tmp_canvas.append(cell2)
                            
                        canvas = deepcopy(tmp_canvas)
                        tmp_canvas= []

                    self.flag_inv_order = True
                    
            elif (ncores%2==0):
                print('Warning: multicores has been tested only for number of cores = 2^j, j integer')
                tmp_canvas = self.multi_equalweight_split_vertically(local_canvas,ncores/2)
                canvas = tmp_canvas
                tmp_canvas= []
                for cell in canvas:
                    cell1,cell2 = self.equalweight_split_horizontally(cell)
                    tmp_canvas.append(cell1)
                    tmp_canvas.append(cell2)
                canvas = tmp_canvas
                tmp_canvas= []

                flag_inv_order = True
                        
            else:
                print('Warning: multicores has been tested only for number of cores = 2^j, j integer')
                canvas = self.multi_equalweight_split_horizontally(local_canvas,ncores)

        else:
            canvas = self.multi_equalweight_split_vertically(local_canvas,ncores)
                
        if (not canvas): exit(-1)
        
        jobs = []
        k=0
        for i in range(ncores):
            p = multiprocessing.Process(target=fit_method, args=(canvas[i],k,nstart,start_direction))
            jobs.append(p)
            p.start()
            k +=1
                
        for p in jobs:
            p.join()


        method_name = 'combine_result_' + str(fit_type)
        combine_method = getattr(self, method_name)

        return combine_method(ncores)


    def combine_result_cell(self,ncores):
        filenames = ['cell_fortran_'+str(i)+'.dat' for i in range(ncores)]
        with open('cell_fortran'+str(self.suffix)+'.dat', 'w') as outfile:
            for fname in filenames:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)
                os.remove(fname)

    
    def combine_result_1D_x(self,ncores):          
        filenames = ['ehist_'+str(i)+'.dat' for i in range(ncores)]
        filetmp = 'ehist_tmp.dat'
        with open(filetmp, 'w') as outfile:
            for fname in filenames:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)
                os.remove(fname)

        out = open('ehist'+str(self.suffix)+'.dat', 'w')
        tmp = np.loadtxt(filetmp, unpack = True)
        tmp_T = np.transpose(tmp)
        sort = tmp_T[tmp_T[:,0].argsort()]
        
        # x,y,width,height = np.loadtxt('cell_fortran'+str(self.suffix)+'.dat', unpack = True)

        # xmin=min(x)
        # if xmin > sort[0,0]:
        #     sort[0,0] = xmin

        # c = [x[i]+width[i] for i in range(len(x))]
        # xmax = max(c)
        # if xmax < sort[-1,1]:
        #     sort[-1,1]= xmax 

        for line in sort:
            out.write(str(line)[1:-1]+"\n")

        os.remove(filetmp)

    # def combine_result_1D_y(self,ncores):          
    #     filenames = ['thetahist_'+str(i)+'.dat' for i in range(ncores)]
    #     filetmp = 'thetahist_tmp.dat'
    #     with open(filetmp, 'w') as outfile:
    #         for fname in filenames:
    #             with open(fname) as infile:
    #                 for line in infile:
    #                     outfile.write(line)
    #             os.remove(fname)

    #     out = open('thetahist.dat', 'w')
    #     tmp = np.loadtxt(filetmp, unpack = True)
    #     tmp_T = np.transpose(tmp)
    #     sort = tmp_T[tmp_T[:,0].argsort()]
    #     for line in sort:
    #         out.write(str(line)[1:-1]+"\n")

    #     os.remove(filetmp)
        

    def fit_cell(self,canvas,i,ncores,start_direction):

        npts=self.npts
        weight_exit = self.npts_exit*self.weight/npts

        file = 'cell_fortran_'+str(i)+'.dat'
        out = open(file,'w')

#        if ncores!=1 :
#            cells = self.multi_equalweight_split_vertically(canvas,ncores)
#        else:

        cells = [canvas]

        if start_direction == 'horizontally':
            split = [self.equalweight_split_horizontally,self.equalweight_split_vertically]
        elif start_direction == 'vertically':
            split = [self.equalweight_split_vertically,self.equalweight_split_horizontally]
        
        if self.flag_inv_order==False:
            while ( (not cells) == False):
                for cell in cells:
                    if(cell.weight > weight_exit):
                        cell1, cell2= split[0](cell)
                    
                        if (cell1.weight > weight_exit):
                            subcell1, subcell2 = split[1](cell1)
                            cells.append(subcell1)
                            cells.append(subcell2)

                        else:
                            self.write_cell(cell1,out)

                        if (cell2.weight > weight_exit):
                            subcell3, subcell4 = split[1](cell2)
                            cells.append(subcell3)
                            cells.append(subcell4)
                        
                        else:
                            self.write_cell(cell2,out)

                    else:
                        self.write_cell(cell,out)

                    cells.remove(cell)

        else:
            while ( (not cells) == False):
                for cell in cells:
                    if(cell.weight > weight_exit):
                        cell1, cell2= split[1](cell)
                    
                        if (cell1.weight > weight_exit):
                            subcell1, subcell2 = split[0](cell1)
                            cells.append(subcell1)
                            cells.append(subcell2)

                        else:
                            self.write_cell(cell1,out)


                        if (cell2.weight > weight_exit):
                            subcell3, subcell4 = split[0](cell2)
                            cells.append(subcell3)
                            cells.append(subcell4)
                        
                        else:
                            self.write_cell(cell2,out)

                    else:
                        self.write_cell(cell,out)
                            
                    cells.remove(cell)
                
        print('end fit')

        out.close()

        
    def fit_1D_x(self,canvas,i,nstart,start_direction):

        file = 'ehist_'+str(i)+'.dat'
        out = open(file,'w')
        
        cells = [canvas]

        npts=self.npts
#        weight_exit = self.weight/50.
        weight_exit = self.weight/50.

        split = [self.equalweight_split_vertically]

        while ( (not cells) == False):
            for cell in cells:
                if(cell.weight > weight_exit):
                    cell1, cell2= split[0](cell)
                    cells.append(cell1)
                    cells.append(cell2)
                else:
                    s = str(cell.corner.x) + "\t" + str(cell.corner.x + cell.width) + "\t" + \
                        str(cell.weight/cell.width) + "\t" + str(cell.error()/cell.width) + "\n" 
                    out.write(s)

                    
                cells.remove(cell)

        print('end fit')

        out.close()

        
    # def fit_1D_y(self):

    #     file = 'thetahist_'+str(i)+'.dat'
    #     out = open(file,'w')

    #     cells = [canvas]

    #     npts=self.npts
    #     weight_exit = self.weight/50
        
    #     while ( (not cells) == False):
    #         for cell in cells:
    #             if(cell.weight > weight_exit):
    #                 cell1, cell2= self.equalweight_split_horizontally(cell)
    #                 cells.append(cell1)
    #                 cells.append(cell2)
    #             else:
    #                 s = str(cell.corner.y) + "\t" + str(cell.corner.y + cell.height) + "\t" + \
    #                     str(cell.weight/cell.height) + "\t" + str(cell.error()/cell.height) + "\n" 

    #             cells.remove(cell)
    #     print('end fit')


    def equalweight_split_vertically(self,cell):
        """ Split a cell vertically with the criterion of equal weight."""
        cell.pts.sort()

        isplit = len(cell.pts)/2
        weight=0
        for pt in cell.pts:
            weight += pt.weight
            if (weight > cell.weight/2):
                isplit = cell.pts.index(pt)
                break
                
        xsplit = (cell.pts[isplit-1].midpoint(cell.pts[isplit])).x

        lcell, rcell = cell.split_vertically(xsplit,
                                            cell.pts[:isplit-1],cell.pts[isplit-1:])
        return lcell, rcell


    def equalweight_split_horizontally(self,cell):
        """ Split a cell horizontally with the criterion of equal weight."""
        for pt in cell.pts:
            pt.exchange_xy()
        cell.pts.sort()
        for pt in cell.pts:
            pt.exchange_xy()

        isplit = len(cell.pts)/2
        weight=0
        for pt in cell.pts:
            weight += pt.weight
            if (weight > cell.weight/2):
                isplit = cell.pts.index(pt)
                break
                
        ysplit = (cell.pts[isplit-1].midpoint(cell.pts[isplit])).y
 
        dcell, ucell = cell.split_horizontally(ysplit,
                                            cell.pts[:isplit-1],cell.pts[isplit-1:])

        return dcell,ucell

    def multi_equalweight_split_vertically(self,cell,npart):
        """ Multi-Split (npart) a cell vertically with the criterion 
        of equal weight."""

        if(npart==1):
            return [cell]
        
        cell.pts.sort()
        weight=0
        isplit=[]
        i=1
        for pt in cell.pts:
            weight += pt.weight
            if (weight > i*cell.weight/npart):
                isplit.append(cell.pts.index(pt))
                i+=1
            if (i==npart):
                break

        xsplit=[]
        for i in isplit:
            xsplit.append( (cell.pts[i-1].midpoint(cell.pts[i])).x )

        return  cell.multi_split_vertically(isplit,xsplit)

    def multi_equalweight_split_horizontally(self,cell,npart):
        """ Multi-Split (npart) a cell horizontally with the criterion 
        of equal weight."""

        if(npart==1):
            return [cell]
        
        for pt in cell.pts:
            pt.exchange_xy()
        cell.pts.sort()
        for pt in cell.pts:
            pt.exchange_xy()

        weight=0
        isplit=[]
        i=1
        for pt in cell.pts:
            weight += pt.weight
            if (weight > i*cell.weight/npart):
                isplit.append(cell.pts.index(pt))
                i+=1
            if (i==npart):
                break

        ysplit=[]
        for i in isplit:
            ysplit.append( (cell.pts[i-1].midpoint(cell.pts[i])).y )

        return  cell.multi_split_horizontally(isplit,ysplit)

    
    def equalpts_split_vertically(self,i):
        """ Split a cell vertically with the criterion of equal number of points."""
        y = copy.deepcopy(self.cells[i].pts)
        y.sort()
        if len(y)%2 == 0:
            isplit = len(y)/2
        else:
            isplit = (len(y)-1)/2+randint(0,1)
        xsplit = (y[isplit-1].midpoint(y[isplit])).x

        lcell, rcell = self.cells[i].split_vertically(xsplit,
                                            y[:isplit],y[isplit:])
        
        self.cells[i:i+1] = lcell, rcell        

    def equalpts_split_horizontally(self,i):
        """ Split a cell horizontally with the criterion of equal number of points."""
        y = copy.deepcopy(self.cells[i].pts)
        for pt in y:
            pt.exchange_xy()
        y.sort()
        for pt in y:
            pt.exchange_xy()
        if len(y)%2 == 0:
            isplit = len(y)/2
        else:
            isplit = (len(y)-1)/2+randint(0,1)

        ysplit = (y[isplit-1].midpoint(y[isplit])).y
        dcell, ucell = self.cells[i].split_horizontally(ysplit,
                                            y[:isplit],y[isplit:])
        self.cells[i:i+1] = dcell, ucell        
        

    def plot (self,outfile=None):
        """ Plot the canvas with the 2D mesh."""
        import matplotlib
#        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.collections import PatchCollection
        from matplotlib.patches import Rectangle

        fig, ax = plt.subplots()
        rectangles = []
        if not self.cells:
            rectangles.append(Rectangle((self.canvas.corner.x, self.canvas.corner.y),
                                   self.canvas.width,self.canvas.height,fill=False))
            ax.add_patch(rectangles[0])
        else:
            for cell in self.cells:
                rectangles.append(Rectangle((cell.corner.x, cell.corner.y),
                                            cell.width,cell.height,fill=False))
            pc = PatchCollection(rectangles, facecolor='w', edgecolor='b' )

            ax.add_collection(pc)
        plt.xlim((self.canvas.corner.x,self.width))
        plt.ylim((self.canvas.corner.y,self.height))
        
        if not outfile:
            plt.show()
        else:
            plt.savefig(outfile)


class fit2D_energy_theta(CellHistogram):
    
    lpp2 = {'electron' : 0, 'DIS' : 1}
    ebeam2 = {'electron' : 0.000511, 'DIS' : 0.938}
    
    def __init__(self,proc_characteristics,input_lhe_evts,interaction_channel):
        # load proc info and fit params 
        self.proc_characteristics = proc_characteristics
        self.fit2D_card = fit2D.Fit2DCard(pjoin('fit2D_card.dat'))
        self.interaction_channel = interaction_channel
        #check for corrupted events
        #self.fixLHE(input_lhe_evts)        
        #store and reweight evts
        self.npass,self.E_min,self.E_max,self.theta_min,self.theta_max,self.data = \
                    self.store_reweight(input_lhe_evts)
        
        super(fit2D_energy_theta,self).__init__(Point(self.E_min,self.theta_min), \
                                                self.E_max-self.E_min,self.theta_max-self.theta_min,self.fit2D_card['npoints_cell'])

        self.add_pts(self.data)
        self.resol_fac = 3.
        
    def fixLHE(self,inputLHE):
        evtsLHEin = lhe_parser.EventFile(inputLHE)
        evtsLHEout = lhe_parser.EventFile(path=inputLHE+'-tmp.gz',mode='w')
        banner = evtsLHEin.get_banner()
        banner.write(evtsLHEout, close_tag=False)
        for event in evtsLHEin:
            corrupted = False
            for particle in event:
                p = lhe_parser.FourMomentum(particle)
                if(np.isnan(p.E)):
                    corrupted = True
            if not corrupted:
                evtsLHEout.write_events(event)
        os.system('mv '+inputLHE+'-tmp.gz '+inputLHE)

    def heaviside(self,x):
        if x>0:
            return 1
        else:
            return 0
        
    def max_travel_distance(self,ctheta,stheta,cphi,sphi):
        depth = self.fit2D_card['depth']
        z1 = self.fit2D_card['d_target_detector']        
        z2 = z1 + depth

        if (ctheta<0.):
            return 0.
        else:
            theta = np.arccos(ctheta)

        # off-axis
        # if (self.fit2D_card['off_axis']):
        #     thetac = self.fit2D_card['thetac']
        #     delta_theta = self.fit2D_card['theta_aperture']
        #     yc = z1*np.sin(thetac)
        #     rcone_proj = z1*(np.sin(thetac+delta_theta)-np.sin(thetac))

        #     if abs(theta -thetac) > delta_theta:
        #         return 0.
            
        #     r = z1*stheta            
        #     sphi_star = (r**2-rcone_proj**2+yc**2)/(2.*r*yc)

        #     if sphi < sphi_star:
        #         return 0.
        #     else:            
        #         return self.fit2D_card['off_axis_depth']

        if (self.fit2D_card['off_axis']):
            
            yc = self.fit2D_card['yc']
            r_proj = self.fit2D_card['radius']

            thetac = np.arctan(yc/z1)
            thetah = np.arctan((yc+r_proj)/z1)
            thetal = np.arctan((yc-r_proj)/z1)

            if theta>thetah or theta<thetal:
                return 0.

            r = z1*np.tan(theta)
            cphi_star = (r**2-r_proj**2+yc**2)/(2.*r*yc)

            #print sphi,cphi_star
            
            if sphi > 0 and sphi > cphi_star:
                return depth/ctheta
            else:            
                return 0.

            
        # cylinder detector
        if (self.fit2D_card['cylinder']):
            # theta_min = self.fit2D_card['theta_min']
            # if theta_min!=0.:
            #     print('Error: theta_min != 0. do not supported!')
            #     exit(-1)
            theta_max = self.fit2D_card['theta_max']
            if theta > theta_max:
                return 0.
            radius = z1*np.tan(theta_max)
            theta_star = np.arctan(radius/z2)
            if (theta<theta_star):
                return depth/ctheta
            else:
                return z1*(np.tan(theta_max)-np.tan(theta))/stheta

        # parallelepiped detector
        if (self.fit2D_card['parallelepiped']): 
            xmin = -self.fit2D_card['x_side']/2.
            xmax = self.fit2D_card['x_side']/2.
            ymin = -self.fit2D_card['y_side']/2.
            ymax = self.fit2D_card['y_side']/2.
            tgth = np.tan(theta)
            x1 = z1*tgth*cphi
            y1 = z1*tgth*sphi
            in_z1 = self.heaviside(x1-xmin)*self.heaviside(xmax-x1) \
                    *self.heaviside(y1-ymin)*self.heaviside(ymax-y1)
            if in_z1 == 0.:
                return 0.
            x2 = z2*tgth*cphi
            y2 = z2*tgth*sphi
            in_z2 = self.heaviside(x2-xmin)*self.heaviside(xmax-x2) \
                    *self.heaviside(y2-ymin)*self.heaviside(ymax-y2)
            if in_z2 != 0.:
                return np.sqrt((x2-x1)**2+(y2-y1)**2+depth**2)
            if x2 > xmax:
                x3 =  xmax
                y3 =  xmax*sphi/cphi
                z3 =  xmax/tgth/cphi 
                return np.sqrt((x3-x1)**2+(y3-y1)**2+(z3-z1)**2)
            if x2 < xmin:
                x3 =  xmin          
                y3 =  xmin*sphi/cphi
                z3 =  xmin/tgth/cphi 
                return np.sqrt((x3-x1)**2+(y3-y1)**2+(z3-z1)**2)
            if y2 > ymax:
                x3 =  ymax*cphi/sphi
                y3 =  ymax
                z3 =  ymax/tgth/sphi
                return np.sqrt((x3-x1)**2+(y3-y1)**2+(z3-z1)**2)
            if y2 < ymin:
                x3 =  ymin*cphi/sphi
                y3 =  ymin
                z3 =  ymin/tgth/sphi
                return np.sqrt((x3-x1)**2+(y3-y1)**2+(z3-z1)**2)

        
    def eff_function(self,E,ctheta,stheta,cphi,sphi):
        depth = self.fit2D_card['depth']
        tmp =  self.max_travel_distance(ctheta,stheta,cphi,sphi)*ctheta/depth
        if tmp > 0.:
            return tmp
        else:
            return 0.
        
    def store_reweight(self,input_lhe_evts):
        E_min = 1e20
        E_max = 0.
        theta_min = np.pi/2.
        theta_max = 0.
        npass = 0
        data = []
        lhe_evts = lhe_parser.EventFile(input_lhe_evts) 
        self.testplot = self.fit2D_card['testplot']
        if self.testplot:
            scatterdata = open('in_DM.dat','w')
            scatterdata.write('# E \t theta \t cos(phi) \t sin(phi) \n')
        # flux normalization 
        norm = self.fit2D_card['flux_norm']
        prod_xsec_flag = self.fit2D_card['prod_xsec_in_norm']
        if prod_xsec_flag:
            norm = norm*lhe_evts.cross
        detector_density = self.fit2D_card['detector_density']
        depth = self.fit2D_card['depth']
        fac_cm2_pb = 10**36
        norm = norm*(detector_density*6.022e23)*depth/fac_cm2_pb
        if self.interaction_channel == 'electron':
            norm = norm * self.fit2D_card['Z_average']/self.fit2D_card['A_average']
        nevt=self.fit2D_card['nevts_norm']
        if nevt < 0:
            nevt = len(lhe_evts)
        for event in lhe_evts:
            for particle in event:
                if len(event) == 1:
                    if particle.status == 1:
                        continue
                else:
                    if particle.status != 1: # stable final state
                        continue
                    
                if abs(particle.pid) == int(self.proc_characteristics['DM']):
                    p = lhe_parser.FourMomentum(particle)
                    ctheta = (p.pz/p.norm)
                    stheta = np.sqrt(1.-ctheta**2)
                    if ((p.px==0) and (p.py==0)):
                        cphi=1.
                        sphi=0.
                    else:
                        cphi = p.px/np.sqrt(p.px**2+p.py**2)
                        sphi = p.py/np.sqrt(p.px**2+p.py**2)
                    weight = self.eff_function(p.E,ctheta,stheta,cphi,sphi)
                    if weight != 0.:
                        npass+=1
                        theta = np.arccos(ctheta)
                        data.append(WeightedPoint(p.E,theta, \
                                                  norm/nevt*weight))
                        if self.testplot:
                            s = str(p.E) + '\t' + str(theta) + '\t' + str(cphi) + \
                                '\t' + str(sphi) + '\t' + str(norm/nevt*weight) +'\n'  
                            scatterdata.write(s)
                        if p.E < E_min:
                            E_min = p.E
                        if p.E > E_max:
                            E_max = p.E
                        if theta < theta_min:
                            theta_min = theta
                        if theta > theta_max:
                            theta_max = theta
        return npass,E_min,E_max,theta_min,theta_max,data
                            #print(E_min,theta_min,E_max,theta_max,len(data))

                            
    def do_fit(self):
        self.ncores = self.fit2D_card['ncores']
        # Generate the 2DMesh, output file: cell_fortran.dat
        self.fit(ncores=self.ncores)
        if self.testplot:
            plot('cell_fortran.dat','mesh2D.png')
        # 3 step: Generate 1D distribution (integrated in angles)
        # and of the output ehist.dat
        self.fit('1D_x',ncores=self.ncores)
        # Update parameters in the run card
        nevts = self.fit2D_card['nevts_interaction']
        if nevts < 0:
            nevts = int (3 * self.npass)
        run_card = banner_mod.RunCard(pjoin('run_card.dat'))
        run_card['nevents'] = nevts
        if self.interaction_channel:
            run_card['lpp2'] = self.lpp2[self.interaction_channel]
        run_card['ebeam1'] = self.E_max
        if self.interaction_channel:
            run_card['ebeam2'] = self.ebeam2[self.interaction_channel]

        else:
            run_card['ebeam2'] = 0.938

        run_card['use_syst'] = False
        run_card.write(pjoin('run_card.dat'))


    def check_consistency(self):
        meshL= CellHistogram(Point(self.E_min,self.theta_min), \
                             self.E_max-self.E_min,self.theta_max-self.theta_min,self.fit2D_card['npoints_cell']/2.,'L')
        meshL.add_pts(self.data)
        meshL.fit(ncores=self.ncores)
        plot('cell_fortran-L.dat','mesh2D-L.png')
        meshH= CellHistogram(Point(self.E_min,self.theta_min), \
                            self.E_max-self.E_min,self.theta_max-self.theta_min,self.fit2D_card['npoints_cell']*2.,'H')
        meshH.add_pts(self.data)
        meshH.fit(ncores=self.ncores)
        plot('cell_fortran-H.dat','mesh2D-H.png')

        meshfiles=['cell_fortran-L.dat','cell_fortran.dat','cell_fortran-H.dat']
        genfiles=['gen-L.dat','gen.dat','gen-H.dat']
        cmpfiles=['cmp-L.dat','cmp.dat','cmp-H.dat']

        outstat=open('cmp_stat.dat','w')
                
        a1,a2,a3,a4 = np.loadtxt('cell_fortran.dat',unpack=True)
        cells = np.transpose(np.array([a1,a2,a3,a4]))
        norm = len(a1)
        
        x=[]
        y=[]
        pdata=[]
        i=0
        npoints=100000
        while(i < npoints):
            # note: points are generated in the whole rectangular domain
            # in which the original data points lie. We rejected here the case in
            # which there are no data points
            tmp_x = self.E_min + random() * (self.E_max - self.E_min)
            tmp_y = self.theta_min + random() * (self.theta_max - self.theta_min)
            tmp_prob = self.get_prob(cells,tmp_x,tmp_y)/float(norm)
            if(tmp_prob > 0.):
                x.append(tmp_x)
                y.append(tmp_y)
                pdata.append(tmp_prob)
                i+=1
        
        for j in range(len(genfiles)):
            a1,a2,a3,a4 = np.loadtxt(meshfiles[j],unpack=True)
            cells = np.transpose(np.array([a1,a2,a3,a4]))
            norm = len(a1)

            out = open(cmpfiles[j],'w')
            pgen=[]
            for i in range(npoints):
                pgen.append(self.get_prob(cells,x[i],y[i])/float(norm))
                out.write(str(x[i])+'\t'+str(y[i])+'\t'+str(pdata[i])+'\t'+str(pgen[i])+'\t'+str(abs(pdata[i]-pgen[i]))+'\n')
            out.close()

            mean1,error1,mean2,error2 = self.do_mean_error(cmpfiles[j])
                
            outstat.write(cmpfiles[j].replace('.dat','')+'\t'+str(mean1)+'\t'+str(error1)+'\t\t'+str(mean2)+'\t'+str(error2)+'\n')
            
        outstat.close()

        self.do_theta_check()
        
                
    # def find_index(self,x,y,xarr,yarr):
    #     # find the 2d bin in the set (xarr,yarr)
    #     # in which the given point (x,y) lies
    #     # (binary search algorithm)

    #     a=0
    #     b=len(xarr)-1
    #     while( b-a > 1):
    #         i=int((b+a)/2.)
    #         if xarr[i] > x:
    #             b=i
    #         else:
    #             a=i
                
    #     if xarr[a]< x:
    #         if xarr[b]>x :
    #             i=a
    #         else:
    #             print('Something wrong 1',xarr[a],xarr[b],x)
    #             print(xarr)
    #     elif xarr[b] < x:
    #         if xarr[b+1]>x :
    #             i=b
    #         else:
    #             print('Something wrong 2',xarr[b],xarr[b+1],x)
                
    #     a=0
    #     b=len(yarr)-1
    #     while( b-a > 1):
    #         j=int((b+a)/2.)
    #         if yarr[j] > y:
    #             b=j
    #         else:
    #             a=j

    #     if yarr[a]<y:
    #         if yarr[b]>y :
    #             j=a
    #         else:
    #             print('Something wrong 3',yarr[a],yarr[b],y)
    #     elif yarr[b]<y:
    #         if yarr[b+1]>y :
    #             j=b
    #         else:
    #             print('Something wrong',yarr[b],yarr[b+1],y)

    #     return i,j

    
    def get_prob(self,cells,x,y):
        cells=cells[cells[:,0].argsort()]

        if x <= min(cells[:,0]):
            return 0.
        if x >= max(np.array(cells[:,0]+cells[:,2])):
            return 0.

        if x >= max(cells[:,0]):
            tmp=np.transpose(np.array([cells[:,0],cells[:,1],cells[:,0]+cells[:,2],cells[:,3]]))
        else:
            a=0
            b=len(cells[:,0])-1
            while( b-a > 1):            
                i=int((b+a)/2.)
                if cells[i,0] > x:
                    b=i
                else:
                    a=i

            if cells[a,0] < x and cells[b,0] > x:
                i=a
            else:
                print(('get_prob error 1',x,cells[a,0],cells[b,0]))
                return 0.

            tmp=np.transpose(np.array([cells[:i+1,0],cells[:i+1,1],cells[:i+1,0]+cells[:i+1,2],cells[:i+1,3]]))

        tmp=tmp[tmp[:,2].argsort()]

        if y <= min(tmp[:,2]):
            selected_cells = tmp
        else:
            a=0
            b=len( tmp[:,2] )-1
            while( b-a > 1):
                i=int((b+a)/2.)
                if tmp[i,2] > y:
                    b=i
                else:
                    a=i

            if tmp[a,2] < y and tmp[b,2] > y:
                i=b
            else:
                print(('get_prob error 2',y,tmp[a,2],tmp[b,2]))
                return 0.

            selected_cells = np.transpose(np.array([tmp[i:,0],tmp[i:,1],tmp[i:,2],tmp[i:,3]]))

        for cell in selected_cells:
            if y>cell[1] and y<cell[1]+cell[3]:
                return 1./(cell[2]-cell[0])/cell[3]

        return 0.


    def do_regenerate(self,energyFile=None,outfile=None,
                      mode='2d',E=None,npoints=30000):

        if mode == '2d':
            xmin,xmax,y,yerr= np.loadtxt('ehist.dat',unpack=True)
            x = (xmin+xmax)/2.
            x=np.insert(x,0,(xmin[0]+x[0])/2.)
            x=np.insert(x,0,xmin[0])
            x=np.append(x,(xmax[-1]+x[-1])/2.)
            x=np.append(x,xmax[-1])
            y=np.insert(y,0,y[0]/2.)
            y=np.insert(y,0,y[0]/2.)
            y=np.append(y,y[-1]/2.)
            y=np.append(y,y[-1]/2.)
            yerr=np.insert(yerr,0,yerr[0]/2.)
            yerr=np.insert(yerr,0,yerr[0]/2.)
            yerr=np.append(yerr,yerr[-1]/2.)
            yerr=np.append(yerr,yerr[-1]/2.)

            resfac=max(y)
            y=y/resfac
            yerr=1./yerr*resfac
            xb=min(xmin)
            xe=max(xmax)
            kx=3
            s=len(x)
            tck = splrep(x, y, yerr, xb, xe, kx, 0,  s)

            E=[]
            theta=[]
            Emin = xb
            Emax = xe
            norm = splint(xb,xe,tck)

            out = open(outfile,'w')
            eps=self.resol_fac* min(self.gen_cells[:,2])
            for i in range(npoints):
                E.append( bisect(self.Esolve,Emin,Emax,args=(tck,Emin,random(),norm)) )
                r=[random() for j in range(2)]
                theta.append(self.pick_theta(E[i],self.gen_cells,eps,r))
                out.write(str(E[i])+'\t'+str(theta[i])+'\n')

        elif mode == '1d':
            eps=self.resol_fac* min(self.gen_cells[:,2])
            
            theta=[]
            if E:
                for i in range(npoints):
                    r=[random() for j in range(2)]
                    theta.append(self.pick_theta(E,self.gen_cells,eps,r))
            else:
                print('Error')

            if outfile:
                out=open(outfile,'w')
                for x in theta:
                    out.write(str(x)+'\n')
            else:
                return np.array(theta)

    def Esolve(self,x, t,x0,r,norm):
        return float(-splint(x0,x,t) + norm*(1.-r) )


    def pick_cells(self,cells,Em,Ep):
        selected_cells = []
        for cell in cells:
            if Em < cell[0]:
                if Ep < cell[0]:
                    continue
                else :
                    selected_cells.append(cell)
            elif Em > cell[0] and Em < cell[0]+cell[2]:
                selected_cells.append(cell)
            else:
                continue
        return selected_cells
    
   
    def pick_theta(self,E,cells,eps,r):
        weights=[]
        wtot=0.
        Em = E-eps
        Ep = E+eps

        selected_cells = np.array(self.pick_cells(cells,Em,Ep))
        
        for i in range(len(selected_cells)):
            Emin =  selected_cells[i,0]
            Emax =  selected_cells[i,0] + selected_cells[i,2]
            d = Emax-Emin
            if Em < Emin:
                if Ep < Emax:
                    w=(Ep-Emin)/d
                    weights.append(w)
                    wtot += w
                elif Ep > Emax:
                    w=1.
                    weights.append(w)
                    wtot += w
            elif Em > Emin:
                if Ep < Emax:
                    w=(Ep-Em)/d
                    weights.append(w)
                    wtot += w
                elif Ep > Emax:
                    w=(Emax-Em)/d
                    weights.append(w)
                    wtot += w

        #pick a theta value according to the hit cells and their weights;
        # once a cell is selected, a value of theta is taken uniformly inside it
        
        s=0.
        for i in range(len(selected_cells)):
            s += weights[i]/wtot
            if r[0] < s:
                theta = selected_cells[i,1] + r[1]*selected_cells[i,3]
                return theta


    def do_mean_error(self,infile):
        
        x,y,pdata,pgen,diff = np.loadtxt(infile,unpack=True)
        norm = np.sum(pdata)

        N=len(x)
        d1=0.
        d2=0.
        for i in range(N):
            if pdata[i] > 0. or pgen[i] > 0.: 
                d1 += pdata[i]/(pdata[i]+pgen[i]) / float(N)
                d2 += pdata[i]/(pdata[i]+pgen[i]) * pdata[i]/norm

        dvar1=0.
        dvar2=0.
        for i in range(N):
            if pdata[i] > 0. or pgen[i] > 0.: 
                dvar1 += (d1-pdata[i]/(pdata[i]+pgen[i]))**2/float(N-1)
                dvar2 += (d2-pdata[i]/(pdata[i]+pgen[i]))**2 * pdata[i]/norm

        derr1=np.sqrt(dvar1)
        derr2=np.sqrt(dvar2)

        return d1,derr1,d2,derr2

    
    def do_theta_check(self,meshfile='cell_fortran.dat'):
        E,theta,cpsi,spsi,weight= np.loadtxt('in_DM.dat',unpack=True)        
        a1,a2,a3,a4 = np.loadtxt(meshfile,unpack=True)
        self.gen_cells = np.transpose(np.array([a1,a2,a3,a4]))

        Emin=min(a1)
        Emax=max(a1+a3)
        Estr = self.fit2D_card['E_arr']
        Earr = [ float(Estr.split(',')[i]) for i in range(len(Estr.split(',')))]

        npoints= 5 * len(E)
        nevtscore=npoints/self.ncores 

        # generate points using the mesh
        jobs=[]
        for i in range(self.ncores):
            outfile='regen-'+str(i)+'.dat'
            p = multiprocessing.Process(target=getattr(self,'do_regenerate'),
                                        args=('ehist.dat',outfile,'2d',None,nevtscore))
            jobs.append(p)
            p.start()

        for p in jobs:
            p.join()

        # combine results
        filenames = ['regen-'+str(i)+'.dat' for i in range(self.ncores)]
        with open('regen.dat', 'w') as outfile:
            for fname in filenames:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)
                os.remove(fname)

        # compare the angular distribution 
        jobs = []
        Em=[]
        Ep=[]
        k=0
        for E0 in Earr:
            resol_par = self.resol_fac*min(self.gen_cells[:,2])
            selected_cells = np.array(self.pick_cells(self.gen_cells,E0-resol_par,E0+resol_par))

            # re-adjust the energy bin size
            resol_par = 2.*min(selected_cells[:,2])
            Em.append(E0-resol_par)
            Ep.append(E0+resol_par)

            if k < self.ncores:
                p = multiprocessing.Process(target=getattr(self,'theta_check'), args=(E0,Em[Earr.index(E0)],Ep[Earr.index(E0)]))
                jobs.append(p)
                p.start()
                k += 1
            else:
                for p in jobs:
                    p.join()
                jobs=[]
                k=0
                p = multiprocessing.Process(target=getattr(self,'theta_check'), args=(E0,Em[Earr.index(E0)],Ep[Earr.index(E0)]))
                jobs.append(p)
                p.start()
                k += 1

        for p in jobs:
            p.join()

        for i in range(len(Earr)):
            infile='theta_check'+str(round(Earr[i],1))+'.dat'
            self.plot_theta_check(infile,Em[i],Ep[i])

            
    def theta_check(self,E0,Em,Ep):

        E,theta,cpsi,spsi,weight= np.loadtxt('in_DM.dat',unpack=True)        
        theta_stripe = []
        weight_stripe = []
        for i in range(len(E)):
            if E[i]>Em and E[i]<Ep:
                theta_stripe.append(theta[i])
                weight_stripe.append(weight[i])

        nbins = self.fit2D_card['nbins']
        tmp = int(len(theta_stripe)/nbins)
        if(tmp < 15):
            nbins= int(len(theta_stripe)/15)+3

        #outfile = 'checkTheta-'+str(E0)+'.dat'
        #theta_gen=self.do_regenerate(mode='1d',E=E0)
            
        Hst,bin_edges = np.histogram(theta_stripe,bins=nbins,weights=weight_stripe)
        x = [(bin_edges[i]+bin_edges[i+1])/2. for i in range(len(Hst)) ]

        E_gen,theta_gen = np.loadtxt('regen.dat',unpack=True)
        theta_stripe_gen = []
        for i in range(len(E_gen)):
            if E_gen[i]>Em and E_gen[i]<Ep:
                theta_stripe_gen.append(theta_gen[i])

        #E_gen,theta_gen = np.loadtxt('checkSpline.dat',unpack=True)
        Hst_gen,bin_edges_gen = np.histogram(theta_stripe_gen,bins=bin_edges,density=True)

        # compute the Poisson error associated to the data Histogram
        sumw2 = []
        for left, right in zip(bin_edges, bin_edges[1:]): 
            ix = np.where((theta_stripe >= left) & (theta_stripe <= right))[0]
            sumw2.append(np.sum([weight_stripe[i] ** 2 for i in ix]))

        sumw2 = np.array(sumw2)
        d=[]
        for left, right in zip(bin_edges, bin_edges[1:]):
            d.append(right-left)

        d=np.array(d)
        norm = np.sum(Hst*d)
        Hst = Hst/norm
        err = np.sqrt(sumw2)/norm

        y = np.interp(x,x,Hst,left=0.,right=0.)
        yp = np.interp(x,x,Hst+err,left=0.,right=0.)
        ym = np.interp(x,x,Hst-err,left=0.,right=0.)

        outfile='theta_check'+str(round(E0,1))+'.dat'
        out=open(outfile,'w')
        for i in range(len(x)):
            s=str(x[i])+'\t'+str(y[i])+'\t'+str(ym[i])+'\t'+str(yp[i])+'\t'+str(Hst_gen[i])+'\n'
            out.write(s)
        out.close()
        
    def plot_theta_check(self,infile,Em,Ep):
        import matplotlib
#        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.ticker import MultipleLocator
    
        x,y,ym,yp,Hst_gen = np.loadtxt(infile,unpack=True)

        ### Begin Plot settings ###

        font = {'family': 'serif',
                'color':  'black',
                'weight': 'normal',
                'size': 14,
        }

        font_title = {'family': 'serif',
                      'color':  'black',
                      'weight': 'normal',
                      'size': 20,
        }

        plt.rcParams['xtick.labelsize']=10
        plt.rcParams['ytick.labelsize']=10

        fig=plt.figure()
        ax1 = fig.add_axes([0.1,0.4,0.8,0.5],
                            xticklabels=[], ylim=(0.,max(yp)+max(yp)/10.) )
        ax2 = fig.add_axes([0.1,0.1,0.8,0.3],
                            ylim=(0.,2.4 ) )

        #ax1.set_title(r'$\theta$ distribution for E= '+str(round(E0,1))+' GeV',fontdict=font_title)
        ax2.set_xlabel(r'$\theta$',fontdict=font)
        ax1.set_ylabel('pdf',fontdict=font)
        ax2.set_ylabel('ratio',fontdict=font)

        # minor ticks
        ax2.yaxis.set_major_locator(MultipleLocator(0.5))
        ax2.yaxis.set_minor_locator(MultipleLocator(0.1))
        ax2.tick_params(axis='y',which='minor',direction='in',right='on')
        ax2.grid(True)

        ### End Plot settings ###
        
        ax1.plot(x,y)
        ax1.fill_between(x,ym,yp)
        ax1.plot(x,Hst_gen)

        par_min=min(y)/10.
        ax2.plot(x,y/(y+par_min))
        ax2.fill_between(x,ym/(y+par_min),yp/(y+par_min))
        ax2.plot(x,Hst_gen/(y+par_min))

        # adjust aspect-ratio 
        ratio=0.7
        ax1.set_aspect(1.0/ax1.get_data_ratio()*ratio)
        ratio=0.42
        ax2.set_aspect(1.0/ax2.get_data_ratio()*ratio)

        for label in ax2.xaxis.get_ticklabels()[::2]:
            label.set_visible(False)

        textstr=str(round(Em,1))+'GeV<E<'+str(round(Ep,1))+'GeV'
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax1.text(0.05, 0.95, textstr, transform=ax1.transAxes, fontsize=12,
                 verticalalignment='top', bbox=props)

        plt.savefig(infile.replace('.dat','.pdf'))
        plt.close(fig)

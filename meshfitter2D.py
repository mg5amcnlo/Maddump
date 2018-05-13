from copy import copy, deepcopy
import numpy as np
from random import randint 
import multiprocessing
import os 

pjoin = os.path.join

import fit2D_card as fit2D
import madgraph.various.banner as banner_mod
import madgraph.various.lhe_parser as lhe_parser

def plot (infile,outfile=None):
    """ Plot the canvas with the 2D mesh."""
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
    
    def __init__(self, pos, width, height,npts_exit=25):
        self.cells = []
        self.intervals_x = []        
        self.intervals_y = []        
        self.canvas = Cell(pos,width,height)
        self.npts = 0
        self.weight = 0
        self.width=width
        self.height=height
        self.npts_exit=npts_exit
        
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
        wgt = cell.weight
        cell.pts.sort()
        if(corner_x == self.canvas.corner.x):
            peripheral = True
            weight=0
            for pt in cell.pts[::-1]:
                weight += pt.weight
                if (weight > fac*wgt):
                    isplit = cell.pts.index(pt)
                    break

            width = width -(cell.pts[isplit+1].x-corner_x)
            corner_x = cell.pts[isplit+1].x

        if(corner_x+width == self.canvas.corner.x+self.canvas.width):
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

        if(corner_y == self.canvas.corner.y):
            peripheral = True
            weight=0
            for pt in cell.pts[::-1]:
                weight += pt.weight
                if (weight > fac*wgt):
                    isplit = cell.pts.index(pt)
                    break
            height = height -(cell.pts[isplit+1].y-corner_y)
            corner_y = cell.pts[isplit+1].y
            
        if(corner_y+height == self.canvas.corner.y+self.canvas.height):
            peripheral = True
            weight=0
            for pt in cell.pts:
                weight += pt.weight
                if (weight > fac*wgt):
                    isplit = cell.pts.index(pt)
                    break
            height = cell.pts[isplit-1].y-corner_y
            
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
                    cell1,cell2 = self.split[0](local_canvas)
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
        with open('cell_fortran.dat', 'w') as outfile:
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

        out = open('ehist.dat', 'w')
        tmp = np.loadtxt(filetmp, unpack = True)
        tmp_T = np.transpose(tmp)
        sort = tmp_T[tmp_T[:,0].argsort()]
        
        x,y,width,height = np.loadtxt('cell_fortran.dat', unpack = True)

        xmin=min(x)
        if xmin > sort[0,0]:
             sort[0,0] = xmin

        c = [x[i]+width[i] for i in range(len(x))]
        xmax = max(c)
        if xmax < sort[-1,1]:
            sort[-1,1]= xmax 

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
        weight_exit = self.weight/25

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
        

    def plot (self):
        """ Plot the canvas with the 2D mesh."""
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
        plt.show()


class fit2D_energy_theta(CellHistogram):
    
    lpp2 = {'electron' : 0, 'DIS' : 1}
    ebeam2 = {'electron' : 0.000511, 'DIS' : 0.938}
    
    def __init__(self,proc_characteristics,input_lhe_evts,interaction_channel):
        # load proc info and fit params 
        self.proc_characteristics = proc_characteristics
        self.fit2D_card = fit2D.Fit2DCard(pjoin('fit2D_card.dat'))
        self.interaction_channel = interaction_channel
        #store and reweight evts
        self.npass,self.E_min,self.E_max,self.theta_min,self.theta_max,self.data = \
                    self.store_reweight(input_lhe_evts)
        
        super(fit2D_energy_theta,self).__init__(Point(self.E_min,self.theta_min), \
                                                self.E_max-self.E_min,self.theta_max-self.theta_min,25)

        self.add_pts(self.data)

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
        # cylinder detector
        if (self.fit2D_card['cylinder']):
            theta_min = self.fit2D_card['theta_min']
            if theta_min!=0.:
                print('Error: theta_min != 0. do not supported!')
                exit(-1)
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
        return self.max_travel_distance(ctheta,stheta,cphi,sphi)*ctheta/depth

    
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
        target_density = self.fit2D_card['target_density']
        depth = self.fit2D_card['depth']
        fac_cm2_pb = 10**36
        norm = norm*(target_density*6.022e23)*depth/fac_cm2_pb
        nevt = len(lhe_evts)
        for event in lhe_evts:
            for particle in event:
                if particle.status == 1: # stable final state 
                    if abs(particle.pid) == int(self.proc_characteristics['DM']):
                        p = lhe_parser.FourMomentum(particle)
                        ctheta = (p.pz/p.norm)
                        stheta = np.sqrt(1.-ctheta**2)
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
                                    '\t' + str(sphi) + '\n'  
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
        ncores = self.fit2D_card['ncores']
        # Generate the 2DMesh, output file: cell_fortran.dat
        self.fit(ncores=ncores)
        if self.testplot:
            plot('cell_fortran.dat','mesh2D.png')
        # 3 step: Generate 1D distribution (integrated in angles)
        # and of the output ehist.dat
        self.fit('1D_x',ncores=ncores)
        # Update parameters in the run card
        run_card = banner_mod.RunCard(pjoin('run_card_default.dat'))
        run_card['lpp2'] = self.lpp2[self.interaction_channel]
        run_card['ebeam1'] = self.E_max
        run_card['ebeam2'] = self.ebeam2[self.interaction_channel]
        run_card['use_syst'] = False
        run_card.write(pjoin('run_card.dat'))

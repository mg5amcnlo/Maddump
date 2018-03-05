from copy import copy, deepcopy
from numpy import sqrt,loadtxt,log,transpose
from random import randint 
import multiprocessing
from os import remove


def plot (file,nrun=-1):
    """ Plot the canvas with the 2D mesh."""
    import matplotlib.pyplot as plt
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Rectangle

    x,y,width,height = loadtxt(file, unpack = True)
    
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

    if nrun<0:
        plt.show()
    else:
        s='plot_'+str(nrun)+'.png'
        plt.savefig(s)

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
        return sqrt(tmp)
        
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

        # #        twocells= False
        cell.pts.sort()
        wgt = cell.weight
        loop=True
        while(loop):
            weight=0
            for pt in cell.pts:
                weight += pt.weight
                if (weight > 0.95*wgt):
                    isplit = cell.pts.index(pt)
                    break
            if ( (cell.pts[isplit-1].x-cell.corner.x)  < .6*width ):
                #print((cell.pts[isplit].x-cell.corner.x),cell.width)
                width= cell.pts[isplit-1].x-cell.corner.x
                wgt=0.95*wgt
            else:
                loop=False

        loop=True
#        wgt = cell.weight
        while(loop):
            weight=0
            for pt in cell.pts[::-1]:
                weight += pt.weight
                if (weight > 0.95*wgt):
                    isplit = cell.pts.index(pt)
                    break
            if ( (width+corner_x-cell.pts[isplit-1].x)  < .6*width ):
                #print((cell.pts[isplit].x-cell.corner.x),cell.width)
                #width= cell.pts[isplit-1].x-cell.corner.x
                width = width+corner_x-cell.pts[isplit-1].x
                corner_x = cell.pts[isplit-1].x
                wgt=0.95*wgt
            else:
                loop=False

                
        # #d1 = cell.pts[1].x-cell.corner.x
        # #d2 = cell.width - cell.pts[-1].x
        # d12 = cell.pts[-1].x - cell.pts[-2].x
        # d12ref = cell.pts[-2].x - cell.corner.x
        
        # # if d1 > cell.width/2.:
        # #     corner_x = cell.pts[1].x
        # #     corner_y = cell.corner.y
        # #     width = cell.width - cell.pts[1].x
        # #     height = cell.height
        # # if d2 > cell.width/2.:
        # #     corner_x = cell.corner.x
        # #     corner_y = cell.corner.y
        # #     width = cell.pts[-1].x
        # #     height = cell.height

        # if d12>2*d12ref:
        #     width = cell.pts[-2].x - cell.corner.x 
                
        for pt in cell.pts:
            pt.exchange_xy()
        cell.pts.sort()
        for pt in cell.pts:
            pt.exchange_xy()

        loop=True
        wgt = cell.weight
        while(loop):
            weight=0
            for pt in cell.pts:
                weight += pt.weight
                if (weight > 0.95*wgt):
                    isplit = cell.pts.index(pt)
                    break
            if ( (cell.pts[isplit-1].y-cell.corner.y) < .6*height ):
                height= cell.pts[isplit-1].y-cell.corner.y
                wgt=0.95*wgt
            else:
                loop=False

                               
        loop=True
#        wgt = cell.weight
        while(loop):
            weight=0
            for pt in cell.pts[::-1]:
                weight += pt.weight
                if (weight > 0.95*wgt):
                    isplit = cell.pts.index(pt)
                    break
            if ( (height+corner_y - cell.pts[isplit-1].y) < .6*height ):
                height= height+corner_y - cell.pts[isplit-1].y
                corner_y = cell.pts[isplit-1].y
                wgt=0.95*wgt
            else:
                loop=False

        # #d1 = cell.pts[1].y-cell.corner.y
        # #d2 = cell.height - cell.pts[-1].y

        # d12 = cell.pts[-1].y - cell.pts[-2].y
        # d12ref = cell.pts[-2].y - cell.corner.y

        # # if d1 > cell.height/2.:
        # #     if d1==dmax:
        # #         corner_x = cell.corner.x
        # #         corner_y = cell.pts[1].y
        # #         width = cell.width
        # #         height = cell.height - cell.pts[1].y
        # #     elif d2==dmax:
        # #         corner_x = cell.corner.x
        # #         corner_y = cell.corner.y
        # #         width = cell.width
        # #         height = cell.pts[-1].y

        # if d12>2*d12ref:
        #     height = cell.pts[-2].y - cell.corner.y 
                    
        s = str(corner_x) + "\t" + str(corner_y) + "\t" + \
            str(width) + "\t" + str(height) + "\n" 
        file.write(s)
            
    def fit(self,fit_type='cell',ncores=1):
        """ Fit 2Dcell/1Dinterval mesh according to data points distribution. """
        method_name = 'fit_' + str(fit_type)
        fit_method = getattr(self, method_name)

        local_canvas = deepcopy(self.canvas)
        canvas=[]
        tmp_canvas = []

        nstart=1
        flag_inv_order=False
        
        if (fit_type=='cell'):
            
            if ncores in [2**j for j in range(1,10)]:
                canvas.append(local_canvas)
                n = int(log(ncores)/log(2))
                if (n%2==0):
                    for i in range(n/2):
                        for cell in canvas:
                            cell1,cell2 = self.equalweight_split_horizontally(cell)
                            tmp_canvas.append(cell1)
                            tmp_canvas.append(cell2)
                            
                        canvas = tmp_canvas
                        tmp_canvas= []
                            
                        for cell in canvas:
                            cell1,cell2 = self.equalweight_split_vertically(cell)
                            tmp_canvas.append(cell1)
                            tmp_canvas.append(cell2)
                            
                        canvas = tmp_canvas
                        tmp_canvas= []

                else:
                    cell1,cell2 = self.equalweight_split_horizontally(local_canvas)
                    tmp_canvas.append(cell1)
                    tmp_canvas.append(cell2)

                    canvas = tmp_canvas
                    tmp_canvas= []

                    for i in range((n-1)/2):
                        for cell in canvas:
                            cell1,cell2 = self.equalweight_split_vertically(cell)
                            tmp_canvas.append(cell1)
                            tmp_canvas.append(cell2)
                            
                        canvas = tmp_canvas
                        tmp_canvas= []
                            
                        for cell in canvas:
                            cell1,cell2 = self.equalweight_split_horizontally(cell)
                            tmp_canvas.append(cell1)
                            tmp_canvas.append(cell2)
                            
                        canvas = tmp_canvas
                        tmp_canvas= []

                    flag_inv_order = True
                    
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
            if (not flag_inv_order):
                p = multiprocessing.Process(target=fit_method, args=(canvas[i],k,nstart))
            else:
                p = multiprocessing.Process(target=fit_method, args=(canvas[i],k,nstart,True))
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
                remove(fname)

    
    def combine_result_1D_x(self,ncores):          
        filenames = ['ehist_'+str(i)+'.dat' for i in range(ncores)]
        filetmp = 'ehist_tmp.dat'
        with open(filetmp, 'w') as outfile:
            for fname in filenames:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)
                remove(fname)

        out = open('ehist.dat', 'w')
        tmp = loadtxt(filetmp, unpack = True)
        tmp_T = transpose(tmp)
        sort = tmp_T[tmp_T[:,0].argsort()]
        for line in sort:
            out.write(str(line)[1:-1]+"\n")

        remove(filetmp)

    # def combine_result_1D_y(self,ncores):          
    #     filenames = ['thetahist_'+str(i)+'.dat' for i in range(ncores)]
    #     filetmp = 'thetahist_tmp.dat'
    #     with open(filetmp, 'w') as outfile:
    #         for fname in filenames:
    #             with open(fname) as infile:
    #                 for line in infile:
    #                     outfile.write(line)
    #             remove(fname)

    #     out = open('thetahist.dat', 'w')
    #     tmp = loadtxt(filetmp, unpack = True)
    #     tmp_T = transpose(tmp)
    #     sort = tmp_T[tmp_T[:,0].argsort()]
    #     for line in sort:
    #         out.write(str(line)[1:-1]+"\n")

    #     remove(filetmp)
        

    def fit_cell(self,canvas,i,ncores,flag_inv_order=False):

        npts=self.npts
        weight_exit = self.npts_exit*self.weight/npts

        file = 'cell_fortran_'+str(i)+'.dat'
        out = open(file,'w')

#        if ncores!=1 :
#            cells = self.multi_equalweight_split_vertically(canvas,ncores)
#        else:
        cells = [canvas]
        
        if flag_inv_order==False:
            while ( (not cells) == False):
                for cell in cells:
                    if(cell.weight > weight_exit):
                        cell1, cell2= self.equalweight_split_horizontally(cell)
                    
                        if (cell1.weight > weight_exit):
                            subcell1, subcell2 = self.equalweight_split_vertically(cell1)
                            cells.append(subcell1)
                            cells.append(subcell2)

                        else:
                            # s = str(cell1.corner.x) + "\t" + str(cell1.corner.y) + "\t" + \
                            #     str(cell1.width) + "\t" + str(cell1.height) + "\n" 
                            # out.write(s)
                            self.write_cell(cell1,out)

                        if (cell2.weight > weight_exit):
                            subcell3, subcell4 = self.equalweight_split_vertically(cell2)
                            cells.append(subcell3)
                            cells.append(subcell4)
                        
                        else:
                            # s = str(cell2.corner.x) + "\t" + str(cell2.corner.y) + "\t" + \
                            #     str(cell2.width) + "\t" + str(cell2.height) + "\n" 
                            # out.write(s)
                            self.write_cell(cell2,out)

                    else:
                        # s = str(cell.corner.x) + "\t" + str(cell.corner.y) + "\t" + \
                        #     str(cell.width) + "\t" + str(cell.height) + "\n" 
                        # out.write(s)
                        self.write_cell(cell,out)


                    cells.remove(cell)

        else:
            while ( (not cells) == False):
                for cell in cells:
                    if(cell.weight > weight_exit):
                        cell1, cell2= self.equalweight_split_vertically(cell)
                    
                        if (cell1.weight > weight_exit):
                            subcell1, subcell2 = self.equalweight_split_horizontally(cell1)
                            cells.append(subcell1)
                            cells.append(subcell2)

                        else:
                            # s = str(cell1.corner.x) + "\t" + str(cell1.corner.y) + "\t" + \
                            #     str(cell1.width) + "\t" + str(cell1.height) + "\n" 
                            # out.write(s)
                            self.write_cell(cell1,out)


                        if (cell2.weight > weight_exit):
                            subcell3, subcell4 = self.equalweight_split_horizontally(cell2)
                            cells.append(subcell3)
                            cells.append(subcell4)
                        
                        else:
                            # s = str(cell2.corner.x) + "\t" + str(cell2.corner.y) + "\t" + \
                            #     str(cell2.width) + "\t" + str(cell2.height) + "\n" 
                            # out.write(s)
                            self.write_cell(cell2,out)

                    else:
                        # s = str(cell.corner.x) + "\t" + str(cell.corner.y) + "\t" + \
                        #     str(cell.width) + "\t" + str(cell.height) + "\n" 
                        # out.write(s)
                        self.write_cell(cell,out)
                            

                    cells.remove(cell)
                
        print('end fit')

        out.close()
        
    def fit_1D_x(self,canvas,i,nstart):

        file = 'ehist_'+str(i)+'.dat'
        out = open(file,'w')

        local_canvas = deepcopy(self.canvas)
        
        x,y,width,height = loadtxt('cell_fortran.dat', unpack = True)

        xmin=min(x)
        if xmin > local_canvas.corner.x:
            local_canvas.corner.x = xmin

        c = [x[i]+width[i] for i in range(len(x))]
        xmax = max(c)
        if xmax < local_canvas.corner.x+local_canvas.width:
            local_canvas.width = xmax - local_canvas.corner.x 
        
        cells = [local_canvas]

        npts=self.npts
        weight_exit = self.weight/50

        while ( (not cells) == False):
            for cell in cells:
                if(cell.weight > weight_exit):
                    cell1, cell2= self.equalweight_split_vertically(cell)
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

        maxGapFac = 5
        for subcell in [lcell,rcell]:
            subcell.pts.sort()
            d1L = subcell.pts[1].x-subcell.corner.x
            d12L = subcell.pts[2].x -subcell.pts[1].x

            d1R = subcell.width - (subcell.pts[-1].x-subcell.corner.x)
            d12R = subcell.pts[-1].x - subcell.pts[-2].x

            if d1L > maxGapFac*d12L:
                d12L /= 10
                subcell.width = subcell.width - (subcell.pts[1].x - d12L - subcell.corner.x)
                subcell.corner.x = subcell.pts[1].x - d12L

            if d1R > maxGapFac*d12R:
                d12R /= 10
                subcell.width = subcell.pts[-1].x + d12R - subcell.corner.x

            for pt in subcell.pts:
                pt.exchange_xy()

            subcell.pts.sort()
            for pt in subcell.pts:
                pt.exchange_xy()

            d1D = subcell.pts[1].y-subcell.corner.y
            d12D = subcell.pts[2].y -subcell.pts[1].y

            d1U = subcell.height - (subcell.pts[-1].y-subcell.corner.y)
            d12U = subcell.pts[-1].y -subcell.pts[-2].y

            if d1D > maxGapFac*d12D:
                d12D /= 10
                subcell.height = subcell.height - (subcell.pts[1].y-d12D -subcell.corner.y)
                subcell.corner.y = subcell.pts[1].y - d12D
                #                corner_y = cell.corner.y
                #                height = cell.height
            if d1U > maxGapFac*d12U:
                d12U /= 10
                subcell.height = subcell.pts[-1].y + d12U - subcell.corner.y

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

        maxGapFac = 5
        for subcell in [dcell,ucell]:

            subcell.pts.sort()
            d1L = subcell.pts[1].x-subcell.corner.x
            d12L = subcell.pts[2].x -subcell.pts[1].x

            d1R = subcell.width - (subcell.pts[-1].x-subcell.corner.x)
            d12R = subcell.pts[-1].x - subcell.pts[-2].x

            if d1L > maxGapFac*d12L:
                d12L /= 10
                subcell.width = subcell.width - (subcell.pts[1].x-d12L -subcell.corner.x)
                subcell.corner.x = subcell.pts[1].x-d12L
                #                corner_y = cell.corner.y
                #                height = cell.height
            if d1R > maxGapFac*d12R:
                d12R /= 10
                subcell.width = subcell.pts[-1].x + d12R - subcell.corner.x

            for pt in subcell.pts:
                pt.exchange_xy()
            subcell.pts.sort()
            for pt in subcell.pts:
                pt.exchange_xy()

            d1D = subcell.pts[1].y-subcell.corner.y
            d12D = subcell.pts[2].y -subcell.pts[1].y

            d1U = subcell.height - (subcell.pts[-1].y-subcell.corner.y)
            d12U = subcell.pts[-1].y -subcell.pts[-2].y

            if d1D > maxGapFac*d12D:
                d12D /= 10
                subcell.height = subcell.height - (subcell.pts[1].y-d12D -subcell.corner.y)
                subcell.corner.y = subcell.pts[1].y-d12D
                #                corner_y = cell.corner.y
                #                height = cell.height
            if d1U > maxGapFac*d12U:
                d12U /= 10
                subcell.height = subcell.pts[-1].y + d12U - subcell.corner.y

                
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

    # def export_histogram_cell(self,file,cells):
    #     out = open(file,'w')
    #     for cell in cells:
    #         s = str(cell.corner.x) + "d0" + "\t" + str(cell.corner.y) + "d0" + "\t" + \
    #             str(cell.width) + "d0" + "\t" + str(cell.height) + "d0" + "\n" 
    #         out.write(s)
    #     out.close()


    # def export_histogram(self,file,fit_type='cell'):
    #     """ Export the 2d/1D mesh in raw text output file."""

    #     method_name = 'export_histogram_' + str(fit_type)
    #     export_histogram_method = getattr(self, method_name) 
    #     return export_histogram_method(file)

    # def export_histogram_cell(self,file):
    #     if self.npts == 0:
    #         if len(self.cells)<2:
    #             print "The histogram is empty, nothing to export"
    #             return 
    #     out = open(file,'w')
    #     for cell in self.cells:
    #         s = str(cell.corner.x) + "\t" + str(cell.corner.y) + "\t" + \
    #             str(cell.width) + "\t" + str(cell.height) + "\n" 
    #         out.write(s)
    #     out.close()

    # def export_histogram_1D_x(self,file):
    #     if self.npts == 0:
    #         if len(self.intervals_x)<2:
    #             print "The histogram is empty, nothing to export"
    #             return     
    #     out = open(file,'w')
    #     nintervals_x = len(self.intervals_x)
    #     tmp = [self.intervals_x[i].corner for i in range(nintervals_x)]
    #     tmp_sort = [self.intervals_x[i].corner for i in range(nintervals_x)]
    #     tmp_sort.sort()
    #     for i in range(nintervals_x-1):
    #         bin_width = tmp_sort[i+1].x - tmp_sort[i].x
    #         s = str(tmp_sort[i].x) + "\t" + str(tmp_sort[i+1].x) + "\t" \
    #             + str(self.intervals_x[tmp.index(tmp_sort[i])].weight/bin_width) + "\t" \
    #             + str(self.intervals_x[tmp.index(tmp_sort[i])].error()/bin_width) + "\n" 
    #         out.write(s)
    #     last_index = tmp.index(tmp_sort[nintervals_x-1])
    #     bin_width = self.intervals_x[last_index].width
    #     s = str(tmp_sort[nintervals_x-1].x) + "\t" + str(tmp_sort[nintervals_x-1].x \
    #         + self.intervals_x[last_index].width) \
    #         + "\t" + str(self.intervals_x[last_index].weight/bin_width) \
    #         + "\t" + str(self.intervals_x[last_index].error()/bin_width) + "\n" 
    #     out.write(s)
    #     out.close()

    # def export_histogram_1D_y(self,file):
    #     if self.npts == 0:
    #         if len(self.intervals_y)<2:
    #             print "The histogram is empty, nothing to export"
    #             return     
    #     out = open(file,'w')
    #     nintervals_y = len(self.intervals_y)
    #     tmp = [self.intervals_y[i].corner for i in range(nintervals_y)]
    #     tmp_sort = [self.intervals_y[i].corner for i in range(nintervals_y)]
    #     tmp_sort.sort()
    #     for i in range(nintervals_y-1):
    #         bin_width = tmp_sort[i+1].y - tmp_sort[i].y
    #         s = str(tmp_sort[i].y) + "\t" + str(tmp_sort[i+1].y) + "\t" \
    #             + str(self.intervals_y[tmp.index(tmp_sort[i])].weight/bin_width) + "\t" \
    #             + str(self.intervals_y[tmp.index(tmp_sort[i])].error()/bin_width) + "\n" 
    #         out.write(s)
    #     last_index = tmp.index(tmp_sort[nintervals_y-1])
    #     bin_width = self.intervals_x[last_index].height
    #     s = str(tmp_sort[nintervals_y-1].y) + "\t" + str(tmp_sort[nintervals_y-1].y \
    #         + self.intervals_x[last_index].width) \
    #         + "\t" + str(self.intervals_y[last_index].weight/bin_width) \
    #         + "\t" + str(self.intervals_y[last_index].error()/bin_width) + "\n" 
    #     out.write(s)
    #     out.close()
    

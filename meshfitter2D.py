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
        from numpy import sqrt
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
            
    def fit(self,fit_type='cell'):
        """ Fit 2Dcell/1Dinterval mesh according to data points distribution. """
        method_name = 'fit_' + str(fit_type)
        fit_method = getattr(self, method_name) 
        return fit_method()
        
    def fit_cell(self):
        import copy
        canvas = copy.copy(self.canvas)
        
        npts=self.npts
        weight_exit = self.npts_exit*self.weight/npts

        cells = [canvas]
        
        while ( (not cells) == False):
            for cell in cells:
                if(cell.weight > weight_exit):
                    cell1, cell2= self.equalweight_split_horizontally(cell)
                    
                    if (cell1.weight > weight_exit):
                        subcell1, subcell2 = self.equalweight_split_vertically(cell1)
                        cells.append(subcell1)
                        cells.append(subcell2)

                    else:
                        self.cells.append(cell1)
                        
                    if (cell2.weight > weight_exit):
                        subcell3, subcell4 = self.equalweight_split_vertically(cell2)
                        cells.append(subcell3)
                        cells.append(subcell4)
                        
                    else:
                        self.cells.append(cell2)
                    
                else:
                    self.cells.append(cell)

                cells.remove(cell)
        print('end fit')

    def fit_1D_x(self):
        import copy
        canvas = copy.copy(self.canvas)

        npts=self.npts
        weight_exit = self.weight/50

        cells = [canvas]
        
        while ( (not cells) == False):
            for cell in cells:
                if(cell.weight > weight_exit):
                    cell1, cell2= self.equalweight_split_horizontally(cell)
                    cells.append(cell1)
                    cells.append(cell2)
                else:
                    self.intervals_x.append(cell)

                cells.remove(cell)
        print('end fit')

    def fit_1D_y(self):
        import copy
        canvas = copy.copy(self.canvas)

        npts=self.npts
        weight_exit = self.weight/50

        cells = [canvas]
        
        while ( (not cells) == False):
            for cell in cells:
                if(cell.weight > weight_exit):
                    cell1, cell2= self.equalweight_split_vertically(cell)
                    cells.append(cell1)
                    cells.append(cell2)
                else:
                    self.intervals_y.append(cell)

                cells.remove(cell)
        print('end fit')


    def equalweight_split_horizontally(self,cell):
        """ Split a cell horizontally with the criterion of equal weight."""

        import copy
        from random import randint 
        y = copy.deepcopy(cell.pts)
        y.sort()
        
        weight=0
        for pt in y:
            weight += pt.weight
            if (weight > cell.weight/2):
                isplit = y.index(pt)-1
                break
                
        xsplit = (y[isplit-1].midpoint(y[isplit])).x

        lcell, rcell = cell.split_vertically(xsplit,
                                            y[:isplit],y[isplit:])

        return lcell, rcell


    def equalweight_split_vertically(self,cell):
        """ Split a cell vertically with the criterion of equal weight."""

        import copy
        from random import randint
        y = copy.deepcopy(cell.pts)
        for pt in y:
            pt.exchange_xy()
        y.sort()
        for pt in y:
            pt.exchange_xy()
        
        weight=0
        for pt in y:
            weight += pt.weight
            if (weight > cell.weight/2):
                isplit = y.index(pt)-1
                break
                
        ysplit = (y[isplit-1].midpoint(y[isplit])).y
 
        dcell, ucell = cell.split_horizontally(ysplit,
                                            y[:isplit],y[isplit:])
        
        return dcell,ucell


    def equalpts_split_horinzotally(self,i):
        """ Split a cell horizontally with the criterion of equal number of points."""
        import copy
        from random import randint 
#        i = self.cells.index(cell) 
#        y = copy.deepcopy(cell.pts)
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

    def equalpts_split_vertically(self,i):
        """ Split a cell vertically with the criterion of equal number of points."""
        import copy
        from random import randint
#        i = self.cells.index(cell) 
#        y = copy.deepcopy(cell.pts)
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

        
    def export_histogram(self,file,fit_type='cell'):
        """ Export the 2d/1D mesh in raw text output file."""

        method_name = 'export_histogram_' + str(fit_type)
        export_histogram_method = getattr(self, method_name) 
        return export_histogram_method(file)

    def export_histogram_cell(self,file):
        if self.npts == 0:
            if len(self.cells)<2:
                print "The histogram is empty, nothing to export"
                return 
        out = open(file,'w')
        for cell in self.cells:
            s = str(cell.corner.x) + "\t" + str(cell.corner.y) + "\t" + \
                str(cell.width) + "\t" + str(cell.height) + "\n" 
            out.write(s)
        out.close()

    def export_histogram_1D_x(self,file):
        if self.npts == 0:
            if len(self.intervals_x)<2:
                print "The histogram is empty, nothing to export"
                return     
        out = open(file,'w')
        nintervals_x = len(self.intervals_x)
        tmp = [self.intervals_x[i].corner for i in range(nintervals_x)]
        tmp_sort = [self.intervals_x[i].corner for i in range(nintervals_x)]
        tmp_sort.sort()
        for i in range(nintervals_x-1):
            bin_width = tmp_sort[i+1].x - tmp_sort[i].x
            s = str(tmp_sort[i].x) + "\t" + str(tmp_sort[i+1].x) + "\t" \
                + str(self.intervals_x[tmp.index(tmp_sort[i])].weight/bin_width) + "\t" \
                + str(self.intervals_x[tmp.index(tmp_sort[i])].error()/bin_width) + "\n" 
            out.write(s)
        last_index = tmp.index(tmp_sort[nintervals_x-1])
        bin_width = self.intervals_x[last_index].width
        s = str(tmp_sort[nintervals_x-1].x) + "\t" + str(tmp_sort[nintervals_x-1].x \
            + self.intervals_x[last_index].width) \
            + "\t" + str(self.intervals_x[last_index].weight/bin_width) \
            + "\t" + str(self.intervals_x[last_index].error()/bin_width) + "\n" 
        out.write(s)
        out.close()

    def export_histogram_1D_y(self,file):
        if self.npts == 0:
            if len(self.intervals_y)<2:
                print "The histogram is empty, nothing to export"
                return     
        out = open(file,'w')
        nintervals_y = len(self.intervals_y)
        tmp = [self.intervals_y[i].corner for i in range(nintervals_y)]
        tmp_sort = [self.intervals_y[i].corner for i in range(nintervals_y)]
        tmp_sort.sort()
        for i in range(nintervals_y-1):
            bin_width = tmp_sort[i+1].y - tmp_sort[i].y
            s = str(tmp_sort[i].y) + "\t" + str(tmp_sort[i+1].y) + "\t" \
                + str(self.intervals_y[tmp.index(tmp_sort[i])].weight/bin_width) + "\t" \
                + str(self.intervals_y[tmp.index(tmp_sort[i])].error()/bin_width) + "\n" 
            out.write(s)
        last_index = tmp.index(tmp_sort[nintervals_y-1])
        bin_width = self.intervals_x[last_index].height
        s = str(tmp_sort[nintervals_y-1].y) + "\t" + str(tmp_sort[nintervals_y-1].y \
            + self.intervals_x[last_index].width) \
            + "\t" + str(self.intervals_y[last_index].weight/bin_width) \
            + "\t" + str(self.intervals_y[last_index].error()/bin_width) + "\n" 
        out.write(s)
        out.close()

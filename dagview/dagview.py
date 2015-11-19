# -*- coding: utf-8 -*-
import sys
import math
import time
import wx
import wx.lib.scrolledpanel

class GFIntermediate:
    """Class defining the intermediates in a GeoFold DAG"""
        
    def __init__(self,number=0, center=(0,0), radius=0.,dagfile=''):
        """Initialize given information from imagemap"""
        self.number = number
        self.radius = radius
        self.center = center
        self.dagfile = dagfile
        if dagfile != '':
            success = self.read_dagfile(dagfile)
            assert success == True, 'Could not read dagfile: %s'%(self.dagfile)
        else:
            self.iflag = ""
            #self.state = 0
            self.sym = 0
            self.sas = 0.
            self.entropy = 0.
            self.voids = 0
            self.hbonds = 0
            self.concentration = 0.        
            self.barrels = []
            self.barrelflags = []
                
    
    def read_dagfile(self,dagfile):
        """Given the GFIntermediate initialized with data from imagemap
        open its parent dagfile and read in remaining information.  Returns
        False if something went wrong"""
        try:
            readDAG = open(dagfile,'r')
        except IOError:
            sys.stderr.write('\n\nError: DAG file %s could not be opened\n'%(dagfile))
            sys.stderr.flush()
            return False
        while 1:
            line = readDAG.readline()
            if line == '':
                sys.stderr.write("Error: End of file reached")
                sys.stderr.flush()
                return False
            #find ISEGMT lines            
            if line[0:6] == 'ISEGMT':
                #find matching ISEGMT number
                try:
                    iseg_num = int(line[6:14])
                    if iseg_num == self.number:
                        #Read in remaining information
                        line
                        self.iflag = readDAG.readline().split()[0]
                        line = line.split()
                        self.sym = int(line[3])
                        self.sas = float(line[4])
                        self.entropy = float(line[5])
                        self.voids = int(line[6])
                        self.hbonds = int(line[7])
                        self.concentration = float(line[8])
                        self.barrels = []
                        self.barrelflags = []
                        for i in range(9,17):
                            self.barrels.append(int(line[i]))
                        for barrel in self.barrels:
                            if barrel != 0:
                                self.barrelflags.append(self.setbarrelflags(readDAG,barrel))
                            if barrel == 0:
                                self.barrelflags.append(('',''))
                        return True
                except Exception as e:
                    sys.stderr.write("Error: "+e.message)
                    sys.stderr.flush()
                    return False
        readDAG.close()
        
    def setbarrelflags(self,readDAG,barrel):
        '''Given a non-zero barrel.  return it's u1flags and u2flags for this intermediate'''
        #initialize flags
        u1flag = ''
        u2flag = ''
        #find the barrel number
        barrel_num = len(self.barrelflags)+1
        foundFlags = False
        while not foundFlags:
            line = readDAG.readline()
            if line == '':
                sys.stderr.write("setbarrelflags::Error: End of file reached")
                sys.stderr.flush()
                raise IOError
            line = line.split()
            if line[0] == 'BARREL' and int(line[1]) == barrel_num:
                    #read until we find the right seam
                    foundSeam = False
                    while not foundFlags:
                        line = readDAG.readline()
                        if line == '':
                            sys.stderr.write("setbarrelflags::Error: End of file reached")
                            sys.stderr.flush()
                            raise IOError
                        line = line.split()
                        if line[0] == 'SEAM' and int(line[1]) == barrel:
                            #read until you find u1flag
                            while not foundFlags:
                                line = readDAG.readline()
                                if line == '':
                                    sys.stderr.write("setbarrelflags::Error: End of file reached")
                                    sys.stderr.flush()
                                    raise IOError
                                line = line.split()
                                if line[0] == 'U1FLAG':
                                    u1flag = line[1]
                                    #u2flag is on the next line
                                    u2flag = readDAG.readline().split()[1]
                                    foundFlags = True
        for i in range(0,len(self.iflag)):
            if self.iflag[i] == '.':
                u1flag = u1flag[:i]+'.'+u1flag[i+1:]
                u2flag = u2flag[:i]+'.'+u2flag[i+1:]
        return (u1flag,u2flag)
                            
    def contains_point(self,(x,y)):
        """Returns True if the given coordinate is within the space on the map
        defined by this intermediate (e.g. (x,y) lies within self.radius of 
        self.center)"""
        center_x, center_y = self.center
        if math.sqrt((center_x-x)**2+(center_y-y)**2) <= self.radius:
            return True
        else:
            return False
            
    def show(self):
        ''' Will display the intermediate in the pymol window'''
        import pymol 
        #residues = self.get_residues()
        residues = get_flag_residues(self.iflag)
        u1res = []
        u2res = []
        for (u1,u2) in self.barrelflags:
            if u1 != '':
                u1res.append(get_flag_residues(u1))
                u2res.append(get_flag_residues(u2))
        intermediate  = ' intermediate_%s'%(str(self.number))
        pymol.cmd.hide()
        pymol.cmd.select(intermediate,residues)
        pymol.cmd.show_as('cartoon',intermediate)
        pymol.cmd.color('purple',intermediate)
        if len(u1res)!=0:
            for i in range(0,len(u1res)):
                u1label = 'i_%d_barrel_%d_u1'%(self.number,i)
                u2label = 'i_%d_barrel_%d_u2'%(self.number,i)
                if u1res[i].split() != ['resi']:
                    pymol.cmd.select(u1label,u1res[i])
                    pymol.cmd.color('yellow',u1label)
                if u2res[i].split() != ['resi']:
                    pymol.cmd.select(u2label,u2res[i])
                    pymol.cmd.color('green',u2label)
        pymol.cmd.deselect()

def get_flag_residues(flag):
    '''gets the residue labeling for given flag (iflag,u1flag,u2flag)'''
    residues = []        
    for i in range(0,len(flag)):
        if flag[i] != '.':
            residues.append(str(i+1))
    residues = 'resi %s' %(','.join(residues))
    return residues
    
class GFTransition:
    """Class defining the transition states in a GeoFold DAG"""
    def __init__(self,number = 0, coords = ((0,0),(0,0),(0,0),(0,0)), dagfile = ''):
        #info from imgmap
        self.number = number
        self.coords = coords
        self.dagfile = dagfile
        #info from dag
        if dagfile != '':
            success = self.read_dagfile(dagfile)
            assert success == True, 'Could not read dagfile: %s'%(dagfile)            
        else:
            self.f = 0
            self.u1 = 0
            self.u2 = 0
            self.entropy = 0.
            self.cuttype = ''
            self.iseam = 0
            self.traffic = 0.
            
    def read_dagfile(self,dagfile):
        try:
            readDAG = open(dagfile,'r')
        except IOError:
            sys.stderr.write('\n\nGFTransition::Error: Could not open file: %s\n'%(dagfile))
            sys.stderr.flush()
            return False
        while 1:
            line = readDAG.readline()
            if line == '':
                sys.stderr.write("\n\nError: End of file reached\n")
                sys.stderr.flush()
                return False
            #Find TSTATE line
            if line[0:6] == 'TSTATE':
                line = line.split()
                if int(line[1]) == self.number:
                    try:
                        self.f = int(line[2])
                        self.u1 = int(line[3])
                        self.u2 = int(line[4])
                        self.entropy = float(line[5])
                        self.cuttype = line[6]
                        self.iseam = int(line[7])
                        self.traffic = float(line[8])
                    except Exception as e:
                        sys.stderr.write("\n\nError: %s\n"%(e.message))
                        sys.stderr.flush()
                        return False
                    else:
                        return True
        
        
            
        

    def contains_point(self,(x,y)):
        """Returns True if the given coordinate is within the space on the map
        defined by this intermediate (e.g. (x,y) lies within the box bounded by
        self.coords).  This uses the left-hand test"""
        ((x1,y1),(x2,y2),(x3,y3),(x4,y4)) = self.coords
        #1,2
        if self.isLeft((x1,y1),(x2,y2),(x,y)):
            return False
        #2,3
        if self.isLeft((x2,y2),(x3,y3),(x,y)):
            return False
        #3,4
        if self.isLeft((x3,y3),(x4,y4),(x,y)):
            return False
        #4,5
        if self.isLeft((x4,y4),(x1,y1),(x,y)):
            return False
        return True
        
    def isLeft(self,(x1,y1),(x2,y2),(x,y)):
        A = -(y2-y1)
        B = x2-x1
        C = -(A*x1 + B*y1)
        D = A*x + B*y + C
        return D > 0
        
    def show(self,intermediates):
        '''Displays the transition state on the pymol viewer'''
        import pymol
        if self.u2 == 0:
            u2 = 0
        for intermediate in intermediates:
            if intermediate.number == self.f:
                f = intermediate
            if intermediate.number == self.u1:
                u1 = intermediate
            if intermediate.number == self.u2:
                u2 = intermediate
        u1res = get_flag_residues(u1.iflag)
        if u2 != 0:
            u2res = get_flag_residues(u2.iflag)
        #pymol.cmd.select('f',fres)
        f.show()
        pymol.cmd.select('u1',u1res)
        if u2 != 0:
            pymol.cmd.select('u2',u2res)
        #pymol.cmd.hide()        
        #pymol.cmd.show_as('cartoon','f')
        if u2 != 0:
            pymol.cmd.color('red','u1')
            pymol.cmd.color('blue','u2')
            pymol.cmd.deselect()
        #is seam
        else:
            pymol.cmd.hide()
            u1.show()
        

def parseImgMap(mapFile,dag=''):
    """This function takes the html imagemap file generated by GeoFold
    and uses it to create a list of Intermediate and transition states"""
    
    transitions = []
    intermediates = []
    readMap = open(mapFile,'r')
    for line in readMap:
        if "<area shape" in line:
            line = line.split('"')            
            querystring = line[5].split('=')
            #intermediate
            if line[1] == 'circle':
                print 'Found intermediate'
                number = int(querystring[1].strip('n&amp;dag'))
                if dag == '':
                    dagfile = querystring[2].strip('&amp;')
                else:
                    dagfile = dag
                for i in range(6,len(line)):
                    if line[i].strip() == 'coords=':
                        coords = line[i+1].split(',')
                        center = (int(coords[0]),int(coords[1]))
                        radius = int(coords[2])
                        intermediates.append(GFIntermediate(number,center,radius,dagfile))
            #Transition
            else:
                number = int(querystring[1].strip('t&amp;dagbphsm'))
                if dag == '':
                    dagfile = querystring[2].strip('&amp;')
                else:
                    dagfile = dag
                for i in range(6,len(line)):
                    if line[i].strip() == 'coords=':
                        coord = line[i+1].split(',')
                        coord = [int(j) for j in coord]
                        coords = ((coord[0],coord[1]),(coord[2],coord[3]),(coord[4],coord[5]),(coord[6],coord[7]))
                        transitions.append(GFTransition(number,coords,dagfile))
    return (intermediates,transitions)

def findFiles(folderPath):
    None    
    
class dagPanel(wx.lib.scrolledpanel.ScrolledPanel):
    def __init__(self,parent, dagImg,dagFile,dagMap):        
        wx.lib.scrolledpanel.ScrolledPanel.__init__(self,parent,-1)
        self.intermediates,self.transitions = parseImgMap(dagMap,dagFile)
        vbox = wx.BoxSizer(wx.VERTICAL)
        img = wx.StaticBitmap(self, -1, wx.Bitmap(dagImg, wx.BITMAP_TYPE_ANY))
        vbox.Add(img)
        
        img.Bind(wx.EVT_MOUSE_EVENTS,self.onClick)
        self.SetSizer(vbox)
        self.SetupScrolling()

    def onClick(self,event):
        if event.GetClickCount() == 1 and event.ButtonUp():
            print event.GetPosition()
            (x,y) = event.GetPosition()
            x = int(x)
            y = int(y)
            print x
            print y
            print len(self.intermediates)
            print len(self.transitions)
            print self.intermediates[0].center
            print self.intermediates[0].radius
            notFound = True
            for intermediate in self.intermediates:
                if intermediate.contains_point((x,y)):
                    print "intermediate %d"%(intermediate.number)
                    notFound = False
                    break
            if notFound:
                for transition in self.transitions:
                    if transition.contains_point((x,y)):
                        print "transition %d"%(transition.number)
                        notFound = False
                        break
            if notFound:
                print "notFound"
        
def startPyMOL(pdb):
    '''starts PyMOL for us.  Only for testing.  PyMOL should already be opened
    by InteractiveROSETTA'''
    import __main__
    __main__.pymol_argv = ["pymol", "-qhxi"]
    import pymol
    pymol.finish_launching()
    pymol.cmd.load(pdb)
    pymol.cmd.show_as('cartoon')
    pymol.cmd.hide()
    pymol.cmd.show_as('cartoon','resi 1,2,3,10,11,12,13,14,15')
    
    
if __name__ == '__main__':
    #intermediates,transitions = parseImgMap('2b3p_florynewtest.21846_1.dag.html','2b3p_florynewtest.21846_1.dag.out')
    #print 'parsed map!'
    #startPyMOL('2b3p_florynewtest.21846.pdb')
    '''for intermediate in intermediates:
        intermediate.show()
        time.sleep(0.01)
    transitions[0].show(intermediates)
    for transition in transitions:
        transition.show(intermediates)
        time.sleep(0.01)
    for intermediate in intermediates:
        if intermediate.number == 302:
            intermediate.show()'''
    app = wx.App(0)
    frame = wx.Frame(None,-1)
    testPanel = dagPanel(frame,'2b3p_florynewtest.21846_1.dag.png','2b3p_florynewtest.21846_1.dag.html','2b3p_florynewtest.21846_1.dag.out')
    frame.Show()
    app.MainLoop()
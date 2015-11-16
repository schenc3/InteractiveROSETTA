# -*- coding: utf-8 -*-
import sys
import math

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
                        for i in range(9,17):
                            self.barrels.append(int(line[i]))
                        return True
                except Exception as e:
                    sys.stderr.write("Error: "+e.message)
                    sys.stderr.flush()
                    return False
        readDAG.close()
    
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
            pass

    
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

def parseImgMap(mapFile):
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
                dagfile = querystring[2].strip('&amp;')
                for i in range(6,len(line)):
                    if line[i].strip() == 'coords=':
                        coords = line[i+1].split(',')
                        center = (int(coords[0]),int(coords[1]))
                        radius = int(coords[2])
                        intermediates.append(GFIntermediate(number,center,radius,dagfile))
            #Transition
            else:
                number = int(querystring[1].strip('t&amp;dagbphsm'))
                dagfile = querystring[2].strip('&amp;')
                for i in range(6,len(line)):
                    if line[i].strip() == 'coords=':
                        coord = line[i+1].split(',')
                        coord = [int(j) for j in coord]
                        coords = ((coord[0],coord[1]),(coord[2],coord[3]),(coord[4],coord[5]),(coord[6],coord[7]))
                        transitions.append(GFTransition(number,coords,dagfile))
    return (intermediates,transitions)

def findFiles(folderPath):
    None    
    
def startPyMOL(pdb):
    '''starts PyMOL for us.  Only for testing.  PyMOL should already be opened
    by InteractiveROSETTA'''
    import __main__
    __main__.pymol_argv = ["pymol", "-qhxi"]
    import pymol
    pymol.finish_launching()
    pymol.cmd.load(pdb)
    pymol.cmd.show_as('cartoon')
    
    
if __name__ == '__main__':
    startPyMOL('1LMB1.pdb')
    
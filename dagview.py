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
            self.read_dagfile(dagfile)
        else:
            self.iflag = ""
            self.state = 0
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
            sys.stderr.write('Error: DAG file %s could not be opened'%(dagfile))
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
                        for i in range(9,17):
                            self.barrels.append(int(line[i]))
                        return True
                except Exception as e:
                    sys.stderr.write("Error: "+e)
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

    
class GFTransition:
    """Class defining the transition states in a GeoFold DAG"""
    def __init__(self):
        #info from imgmap
        self.number = 0
        self.coords = ((0,0),(0,0),(0,0),(0,0))
        self.dagfile = ""
        #info from dag
        

    def contains_point(self,(x,y)):
        """Returns True if the given coordinate is within the space on the map
        defined by this intermediate (e.g. (x,y) lies within the box bounded by
        self.coords)."""

def parseImgMap(mapFile):
    """This function takes the html imagemap file generated by GeoFold
    and uses it to create a list of Intermediate and transition states"""

def findFiles(folderPath):
    None    
# -*- coding: utf-8 -*-
import sys
import math
import time
import wx
import shutil
import wx.lib.scrolledpanel
import platform
import webbrowser
from tools import *

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


    def show(self,pymol,boundaries):
        ''' Will display the intermediate in the pymol window'''
        #import pymol
        #residues = self.get_residues()
        #residues = '(%s) AND %s'%(get_flag_residues(self.iflag,boundaries),self.IDs)
        residues = 'model %s and (%s)'%(self.IDs[0][0],get_flag_residues(self.iflag,boundaries))
        logInfo('residues: %s'%(residues))
        u1res = []
        u2res = []
        for (u1,u2) in self.barrelflags:
            if u1 != '':
                #u1res.append('(%s) and %s'%(get_flag_residues(u1),self.ID))
                #u2res.append('(%s) and %s'%(get_flag_residues(u2),self.ID))
                u1res.append('model %s and (%s)'%(self.IDs[0][0],get_flag_residues(u1,boundaries)))
                u2res.append('model %s and (%s)'%(self.IDs[0][0],get_flag_residues(u2,boundaries)))
        intermediate  = ' intermediate_%s'%(str(self.number))
        pymol.cmd.hide('ribbon','Native')
        pymol.cmd.hide('cartoon','Native')
        pymol.cmd.show_as('ribbon','Native')
        #pymol.cmd.color('white',self.ID)
        pymol.cmd.set('ribbon_color','white','Native')
        pymol.cmd.select(intermediate,residues)
        pymol.cmd.show_as('cartoon',intermediate)
        #pymol.cmd.color('purple',intermediate)
        pymol.cmd.set("cartoon_color",'purple',intermediate)
        if len(u1res)!=0:
            for i in range(0,len(u1res)):
                u1label = 'i_%d_barrel_%d_u1'%(self.number,i)
                u2label = 'i_%d_barrel_%d_u2'%(self.number,i)
                #if u1res[i] != '(resi ) and %s'%(self.ID):
                if u1res[i] != 'model %s and ()'%(self.IDs[0][0]):
                    logInfo('u1res[%i]: %s'%(i,u1res[i]))
                    pymol.cmd.select(u1label,u1res[i])
                    #pymol.cmd.color('yellow',u1label)
                    pymol.cmd.set("cartoon_color",'yellow',u1label)
                #if u2res[i] != '(resi ) and %s'%(self.ID):
                if u2res[i] != 'model %s and ()'%(self.IDs[0][0]):
                    logInfo('u2res[%i]: %s'%(i,u2res[i]))
                    pymol.cmd.select(u2label,u2res[i])
                    #pymol.cmd.color('green',u2label)
                    pymol.cmd.set("cartoon_color", 'green', u2label)
        pymol.cmd.deselect()

    def setIDs(self,IDs):
        self.IDs = IDs

def get_flag_residues(flag,boundaries):
    '''gets the residue labeling for given flag (iflag,u1flag,u2flag)'''
    residues = []
    for bound in boundaries:
        tmpres = []
        start,stop = (int(boundaries[bound][0]),int(boundaries[bound][1]))
        logInfo('start: %i\nstop: %i'%(start,stop))
        for i in range(0,len(flag)):
            if flag[i] != '.' and i in range(start,stop+1):
                tmpres.append(str(i+1))
        if tmpres != []:
            residues.append('chain %s and (resi %s)'%(bound,','.join(tmpres)))
            logInfo(residues)
    return ' | '.join(residues)
    #return residues

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

    def setIDs(self,IDs):
        self.IDs = IDs

    def show(self,intermediates,pymol,boundaries):
        '''Displays the transition state on the pymol viewer'''
        #import pymol
        if self.u2 == 0:
            u2 = 0
        for intermediate in intermediates:
            if intermediate.number == self.f:
                f = intermediate
            if intermediate.number == self.u1:
                u1 = intermediate
            if intermediate.number == self.u2:
                u2 = intermediate
        #u1res = '(%s) and %s'%(get_flag_residues(u1.iflag),self.ID)
        u1res = 'model %s and (%s)'%(self.IDs[0][0],get_flag_residues(u1.iflag,boundaries))
        logInfo('u1res: %s'%(u1res))
        if u2 != 0:
            #u2res = '(%s) and %s'%(get_flag_residues(u2.iflag),self.ID)
            u2res = 'model %s and (%s)'%(self.IDs[0][0],get_flag_residues(u2.iflag,boundaries))
            logInfo('u2res: %s'%(u2res))
        #pymol.cmd.select('f',fres)
        f.show(pymol,boundaries)
        pymol.cmd.select('u1',u1res)
        if u2 != 0:
            pymol.cmd.select('u2',u2res)
            pymol.cmd.set("cartoon_color",'red','u1')
            pymol.cmd.set("cartoon_color",'blue','u2')
            pymol.cmd.deselect()
        #is seam
        else:
            pymol.cmd.hide('cartoon',self.IDs[0][0])
            pymol.cmd.hide('ribbon',self.IDs[0][0])
            u1.show(pymol,boundaries)


def parseImgMap(mapFile,dag='',IDs=[]):
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
                        tmpIntermediate = GFIntermediate(number,center,radius,dagfile)
                        tmpIntermediate.setIDs(IDs)
                        intermediates.append(tmpIntermediate)
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
                        tmpTransition = GFTransition(number,coords,dagfile)
                        tmpTransition.setIDs(IDs)
                        transitions.append(tmpTransition)
    return (intermediates,transitions)

def GetChainBoundaries(intermediates):
    '''Given a set of IDs find Intermediate 1 and extrapolate the boundaries
    from it'''
    boundaries = {}
    for intermediate in intermediates:
        if intermediate.number == 1:
            iflag = intermediate.iflag
            prev = ''
            first = 0
            second = 0
            for i in range(0,len(iflag)):
                if iflag[i] != prev:
                    if prev != '':
                        second = i-1
                        boundaries[prev] = (first,second)
                    prev = iflag[i]
                    first = i
                if i == len(iflag)-1:
                    second = i
                    boundaries[prev] = (first,second)
            logInfo(boundaries)
            return boundaries


class dagPanel(wx.lib.scrolledpanel.ScrolledPanel):
    def __init__(self,parent, dagImg,dagMap,dagFile,IDs=[]):
        wx.lib.scrolledpanel.ScrolledPanel.__init__(self,parent,-1)
        self.intermediates,self.transitions = parseImgMap(dagMap,dagFile,IDs)
        self.boundaries = GetChainBoundaries(self.intermediates)
        print self.boundaries
        vbox = wx.BoxSizer(wx.VERTICAL)
        img = wx.StaticBitmap(self, -1, wx.Bitmap(dagImg, wx.BITMAP_TYPE_ANY))
        vbox.Add(img)

        img.Bind(wx.EVT_LEFT_UP,self.onClick)
        if sys.platform == 'Darwin':
            img.Bind(wx.EVT_MOUSE_EVENTS,self.osxClick)
        self.SetSizer(vbox)
        self.SetupScrolling()


    def setPyMOL(self,pymol):
        self.pymol = pymol
        self.cmd = pymol.cmd
        self.stored = pymol.stored

    def osxClick(self,event):
        '''OSX doesn't recognize a click like linux does apparently'''
        if event.GetClickCount() == 1 and event.ButtonUp():
            self.onClick(event)

    def onClick(self,event):
        (x,y) = event.GetPosition()
        if platform.system() != 'Darwin' and platform.system() != 'Windows':
            (x,y) = self.CalcUnscrolledPosition(x,y)
        notFound = True
        for intermediate in self.intermediates:
            if intermediate.contains_point((x,y)):
                logInfo("intermediate %d"%(intermediate.number))
                intermediate.show(self.pymol,self.boundaries)
                notFound = False
                break
        if notFound:
            for transition in self.transitions:
                if transition.contains_point((x,y)):
                    logInfo("transition %d"%(transition.number))
                    transition.show(self.intermediates,self.pymol,self.boundaries)
                    notFound = False
                    break
        if notFound:
            logInfo("notFound")

class DagViewPanel(wx.lib.scrolledpanel.ScrolledPanel):
    def __init__(self,parent,W,H):

        #ScrolledPanel initialization
        winh = H-330
        wx.lib.scrolledpanel.ScrolledPanel.__init__(self,parent,id=-1,pos=(10,60),size=(340,winh),name="Pathway Visualization")
        self.SetBackgroundColour("#333333")
        self.parent = parent
        logInfo('385: ScrolledPanel initialized!')

        #Title labeling
        if platform.system() == 'Windows':
            self.lblDagView = wx.StaticText(self,-1,'Pathway Visualization',(25,15),(270,25),style=wx.ALIGN_CENTRE)
            self.lblDagView.SetFont(wx.Font(12,wx.DEFAULT,wx.ITALIC,wx.BOLD))
        #elif platform.system() == 'Darwin':
         #   self.lblDagView = wx.StaticBitmap(self,-1,wx.Image(self.parent.parent.scriptdir+"/images/osx/dagview/lblDagView.png",wx.BITMAP_TYPE_PNG).ConvertToBitmap(),pos=(25,15),size=(270,25))
        else:
            self.lblDagView = wx.StaticText(self,-1,'Pathway Visualization',pos=(90,15),style = wx.ALIGN_CENTRE)
            self.lblDagView.SetFont(wx.Font(12,wx.DEFAULT,wx.ITALIC,wx.BOLD))
        self.lblDagView.SetForegroundColour("#FFFFFF")
        logInfo('397: Title label set!')

        #Help Button
        if platform.system() == 'Darwin':
            self.HelpBtn = wx.BitmapButton(self,id=-1,bitmap=wx.Image(self.parent.parent.scriptdir+'/images/osx/HelpBtn.png',wx.BITMAP_TYPE_PNG).ConvertToBitmap(),pos=(295,10),size=(25,25))
        else:
            self.HelpBtn = wx.Button(self, id=-1,label='?',pos=(295,10),size=(25,25))
            self.HelpBtn.SetForegroundColour("#0000FF")
            self.HelpBtn.SetFont(wx.Font(10,wx.DEFAULT,wx.NORMAL,wx.BOLD))
        self.HelpBtn.Bind(wx.EVT_BUTTON,self.showHelp)
        self.HelpBtn.SetToolTipString("Display the help file for this window")
        logInfo('408: Help button set')

        #Subtile text
        if platform.system() == 'Windows':
            self.lblInst = wx.StaticText(self,-1,'View GeoFold pathways',(0,45),(320,25),wx.ALIGN_CENTRE)
            self.lblInst.SetFont(wx.Font(10,wx.DEFAULT,wx.ITALIC,wx.NORMAL))
        #elif platform.system() == 'Darwin':
        #    self.lblInst = wx.StaticBitmap(self,-1,wx.image(self.parent.parent.scriptdir+'/images/osx/dagview/lblInstDagView.png',wx.BITMAP_TYPE_PNG).ConvertToBitmap(),pos=(0,45),size=(320,25))
        else:
            self.lblInst = wx.StaticText(self,-1,'View GeoFold pathways',pos=(20,45),style=wx.ALIGN_CENTRE)
            self.lblInst.SetFont(wx.Font(10,wx.DEFAULT,wx.ITALIC,wx.NORMAL))
            resizeTextControlForUNIX(self.lblInst,0,self.GetSize()[0])
        self.lblInst.SetForegroundColour("#FFFFFF")
        logInfo('421: Subtitle set!')


        #PDB button
        self.lblPDB = wx.StaticText(self,-1,"None Uploaded", pos=(10,103),size=(180,25),style=wx.ALIGN_CENTRE)
        self.lblPDB.SetFont(wx.Font(10,wx.DEFAULT,wx.NORMAL,wx.BOLD))
        if platform.system() == 'Linux':
            resizeTextControlForUNIX(self.lblPDB,10,180)
        self.lblPDB.SetForegroundColour("#FFFFFF")
        #if platform.system() == 'Darwin':
        #    self.btnLoad = wx.BitmapButton(self,id=-1,bitmap=wx.Image(self.parent.parent.scriptdir+'/images/osx/dagview/btnLoadPDB.png',wx.BITMAP_TYPE_PNG).ConvertToBitmap(),pos=(200,100),size=(110,25))
        #else:
        self.btnLoad = wx.Button(self,id=-1,label='Load zip file',pos=(200,100),size=(110,25))
        self.btnLoad.SetForegroundColour("#000000")
        self.btnLoad.SetFont(wx.Font(10,wx.DEFAULT,wx.NORMAL,wx.BOLD))

        self.btnLoad.Bind(wx.EVT_BUTTON,self.loadZip)
        self.btnLoad.SetToolTipString('Load the zip file containing the GeoFold output')
        logInfo('437: Zip button set!')

        #combobox
        self.dagMenu = wx.ComboBox(self, pos=(10,138), size=(110, 25), choices=[], style=wx.CB_READONLY)
        #GO! Button ViewDag
        #ypos = self.btnDagPng.GetPosition()[1]+self.btnDagPng.GetSize()[1]+10
        ypos = self.dagMenu.GetPosition()[1]+self.dagMenu.GetSize()[1]+10
        #if platform.system() == 'Darwin':
        #    self.btnViewDag = wx.BitmapButton(self,id=-1,bitmap=wx.Image(self.parent.parent.scriptdir+'/images/osx/dagview/btnViewDag.png',wx.BITMAP_TYPE_PNG).ConvertToBitmap(),pos=(80,ypos),size=(150,25))
        #else:
        self.btnViewDag = wx.Button(self,id=-1,label="View Pathway!",pos=(80,ypos),size=(150,25))
        self.btnViewDag.SetForegroundColour("#000000")
        self.btnViewDag.SetFont(wx.Font(10,wx.DEFAULT,wx.ITALIC,wx.BOLD))

        self.btnViewDag.Bind(wx.EVT_BUTTON,self.ViewDagClick)
        logInfo('496: View Dag Button set')

        #Scrolling set up
        logInfo('499: Setting scrolling...')
        logInfo('501: GetSize = %d'%(self.btnViewDag.GetSize()[1]))
        self.scrollh = self.btnViewDag.GetPosition()[1] + self.btnViewDag.GetSize()[1] + 5
        logInfo('503: scrollh set to %d'%(self.scrollh))
        self.SetScrollbars(1,1,320,self.scrollh)
        logInfo('505: Scrollbars set!')
        self.winscrollpos = 0
        self.Bind(wx.EVT_SCROLLWIN, self.scrolled)
        logInfo('508: Scrolling set.')
        logInfo('509: Initialization complete!')

    def loadZip(self,event):
        '''opens a file dialog to open the zip file.  Loads the PDB and populates
        the dagMenu'''
        #create file dialog
        logInfo("Clicked Load Zip button")
        dlg = wx.FileDialog(self, message = 'Choose a File',defaultDir=self.seqWin.cwd,defaultFile='',
	    wildcard="Zip Files (*.zip)|*.zip",style=wx.OPEN | wx.CHANGE_DIR)
        if dlg.ShowModal() == wx.ID_OK:
            paths = dlg.GetPaths()
            #Change cwd to the last opened file
            if platform.system()=='Windows':
                lastDirIndx = paths[len(paths)-1].rfind('\\')
            else:
                lastDirIndx = paths[len(paths)-1].rfind('/')

            self.seqWin.cwd = str(paths[len(paths)-1][0:lastDirIndx])
            filename = str(paths[0])
            logInfo('filename = %s'%(filename))
            if platform.system() == 'Windows':
                localFile = filename.split('\\')
            else:
                localFile = filename.split('/')
            localFile = localFile[len(localFile)-1]
            logInfo('localFile: %s'%(localFile))
            self.lblPDB.SetLabel(localFile)
            self.lblPDB.SetForegroundColour('#FFFFFF')
            if platform.system() == 'Linux':
                resizeTextControlForUNIX(self.lblPDB,10,180)

        #run findFiles on the item
        status,dags = self.findFiles(filename)
        logInfo(dags)

        #Error handling
        logInfo(status)
        if status != 0:
            msgs = {-1:'The zip file selected is invalid.\nPlease try again',-2:'There was an error unzipping the file',-3:'No valid output was found in the zip file',-4:'PDB file could not be loaded from zip file'}
            logInfo(msgs[status])
            wx.MessageBox(msgs[status], "", wx.OK|wx.ICON_EXCLAMATION)
            return -1

        #Process dags to just show the base name
        newdags = []
        logInfo('newdags:')
        for dag in dags:
            dag = localFile.split('.zip')[0]+'_'+dag.split('_')[len(dag.split('_'))-1]
            newdags.append(dag)
            logInfo(dag)
        #Populate dagMenu
        self.dagMenu.Clear()
        self.dagMenu.AppendItems(newdags)
        return 0

    def setSeqWin(self, seqWin):
        self.seqWin = seqWin

    def showHelp(self, event):
        '''Open the help page'''
        if platform.system() == 'Darwin':
            try:
                browser = webbrowser.get('Safari')
            except Exception as e:
                print 'Could not load Safari!  The help files are located at %s/help'%(self.parent.parent.scriptdir)
                return
            browser.open(self.parent.parent.scriptdir+'/help/dagview.html')
        else:
            webbrowser.open(self.parent.parent.scriptdir+'/help/dagview.html')

    def scrolled(self,event):
        self.winscrollpos = self.GetScrollPos(wx.VERTICAL)
        event.Skip()

    def setPyMOL(self,pymol):
        '''Sets PyMOL to be used for this class'''
        self.pymol = pymol
        self.cmd = pymol.cmd
        self.stored = pymol.stored


    def activate(self):
        self.Scroll(0, self.winscrollpos)

    def loadPDB(self,event):
        '''Select PDB file to load'''
        logInfo("Clicked Load PDB button")
        dlg = wx.FileDialog(self, message = 'Choose a File',defaultDir=self.seqWin.cwd,defaultFile='',style=wx.OPEN | wx.CHANGE_DIR)
        if dlg.ShowModal() == wx.ID_OK:
            paths = dlg.GetPaths()
            #Change cwd to the last opened file
            if platform.system()=='Windows':
                lastDirIndx = paths[len(paths)-1].rfind('\\')
            else:
                lastDirIndx = paths[len(paths)-1].rfind('/')
            self.seqWin.cwd = str(paths[len(paths)-1][0:lastDirIndx])
            filename = str(paths[0])
            self.loadedPdb = filename
            localPdb = filename[lastDirIndx+1:]
            goToSandbox()
            try:
                shutil.copy(filename,'params.pdb')
            except:
                logInfo('File could not be copied')
            #Delete a file if we're loading a new one
            try:
                self.cmd.remove('params')
                self.cmd.delete('params')
            except:
                pass
            try:
                self.cmd.load('params.pdb','params')
            except:
                wx.MessageBox('The file %s could not be read!'%(filename),'File cannot be read', wx.OK|wx.ICON_EXCLAMATION)
                return
            logInfo('PDB file loaded',filename)
            self.cmd.select('paramssele','model params')
            self.cmd.hide('everything','paramssele')
            self.cmd.delete('paramssele')
            self.lblPDB.SetLabel(localPdb)
            self.lblPDB.SetForegroundColour('#FFFFFF')
            if platform.system() == 'Linux':
                resizeTextControlForUNIX(self.lblPDB,10,180)


    def loadDagOut(self,event):
        '''Load .dag.out file'''
        logInfo("Clicked Load DagOut button")
        dlg = wx.FileDialog(self, message = 'Choose a File',defaultDir=self.seqWin.cwd,defaultFile='',style=wx.OPEN | wx.CHANGE_DIR)
        if dlg.ShowModal() == wx.ID_OK:
            paths = dlg.GetPaths()
            #Change cwd to the last opened file
            if platform.system()=='Windows':
                lastDirIndx = paths[len(paths)-1].rfind('\\')
            else:
                lastDirIndx = paths[len(paths)-1].rfind('/')
            self.seqWin.cwd = str(paths[len(paths)-1][0:lastDirIndx])
            filename = str(paths[0])
            self.loadedDagOut = filename
            localDagOut = filename[lastDirIndx+1:]
            goToSandbox()
            try:
                shutil.copy(filename,'params.dag.out')
            except:
                pass
            logInfo('dag.out file loaded',filename)
            self.lblDagOut.SetLabel(localDagOut)
            self.lblDagOut.SetForegroundColour('#FFFFFF')
            if platform.system() == 'Linux':
                resizeTextControlForUNIX(self.lblDagOut,10,180)

    def loadDagHtml(self,event):
        '''Load .dag.html file'''
        logInfo("Clicked Load DagHtml button")
        dlg = wx.FileDialog(self, message = 'Choose a File',defaultDir=self.seqWin.cwd,defaultFile='',style=wx.OPEN | wx.CHANGE_DIR)
        if dlg.ShowModal() == wx.ID_OK:
            paths = dlg.GetPaths()
            #Change cwd to the last opened file
            if platform.system()=='Windows':
                lastDirIndx = paths[len(paths)-1].rfind('\\')
            else:
                lastDirIndx = paths[len(paths)-1].rfind('/')
            self.seqWin.cwd = str(paths[len(paths)-1][0:lastDirIndx])
            filename = str(paths[0])
            self.loadedDagHtml = filename
            localDagHtml = filename[lastDirIndx+1:]
            goToSandbox()
            try:
                shutil.copy(filename,'params.dag.html')
            except:
                pass
            logInfo('dag.html file loaded',filename)
            self.lblDagHtml.SetLabel(localDagHtml)
            self.lblDagHtml.SetForegroundColour('#FFFFFF')
            if platform.system() == 'Linux':
                resizeTextControlForUNIX(self.lblDagHtml,10,180)

    def loadDagPng(self,event):
        '''Load .dag.png file'''
        logInfo("Clicked Load DagPng button")
        dlg = wx.FileDialog(self, message = 'Choose a File',defaultDir=self.seqWin.cwd,defaultFile='',style=wx.OPEN | wx.CHANGE_DIR)
        if dlg.ShowModal() == wx.ID_OK:
            paths = dlg.GetPaths()
            #Change cwd to the last opened file
            if platform.system()=='Windows':
                lastDirIndx = paths[len(paths)-1].rfind('\\')
            else:
                lastDirIndx = paths[len(paths)-1].rfind('/')
            self.seqWin.cwd = str(paths[len(paths)-1][0:lastDirIndx])
            filename = str(paths[0])
            self.loadedDagPng = filename
            localDagPng = filename[lastDirIndx+1:]
            goToSandbox()
            try:
                shutil.copy(filename,'params.dag.png')
            except:
                pass
            logInfo('dag.png file loaded',filename)
            self.lblDagPng.SetLabel(localDagPng)
            self.lblDagPng.SetForegroundColour('#FFFFFF')
            if platform.system() == 'Linux':
                resizeTextControlForUNIX(self.lblDagPng,10,180)

    def ViewDagClick(self,event):
        logInfo('View Dag Button Clicked')
        self.cmd.show_as('cartoon','Native')
        #self.cmd.color('purple',self.ID)
        self.cmd.set("cartoon_color",'purple','Native')
        try:
            self.frame.Destroy()
        except:
            pass
        dagbase = self.dagMenu.GetValue()
        if platform.system() == 'Windows':
            self.loadedDagHtml = '%s\\%s.dag.html'%(self.cwd,dagbase)
            self.loadedDagOut = '%s\\%s.dag.out'%(self.cwd,dagbase)
            self.loadedDagPng = '%s\\%s.dag.png'%(self.cwd,dagbase)
        else:
            self.loadedDagHtml = '%s/%s.dag.html'%(self.cwd,dagbase)
            self.loadedDagOut = '%s/%s.dag.out'%(self.cwd,dagbase)
            self.loadedDagPng = '%s/%s.dag.png'%(self.cwd,dagbase)
        #self.intermediates,self.transitions = parseImgMap(self.loadedDagHtml,self.loadedDagOut,self.ID)
        self.frame = wx.Frame(None,-1)
        self.DagPanel = dagPanel(self.frame,self.loadedDagPng,self.loadedDagHtml,self.loadedDagOut,self.IDs)
        self.DagPanel.setPyMOL(self.pymol)
        self.frame.Show()

    def findFiles(self,zipDir):
        '''Takes a given zip file.  extracts it in the sandbox and picks out all
        the files able to be viewed.  It outputs a list to be put in a ComboBox
        to allow the user to select which one to view.  If an error occurs, outputs
        a negative number used to identify the error and handle it'''
        import zipfile
        logInfo('Calling findFiles')
        output = []
        #Check if selected file is valid
        if not zipfile.is_zipfile(zipDir):
            return -1, []
        #Unzip the file in the sandbox
        unzipped = zipfile.ZipFile(zipDir)
        info = unzipped.infolist()
        filename = info[0].filename[:len(info[0].filename)-1]
        goToSandbox()
        try:
            unzipped.extractall()
        except: #failed to unzip
            return -2, []
        #use glob to get all dag.out files
        if platform.system() == 'Windows':
            globDir = '%s\\%s\\*.dag.out'%(os.getcwd(),filename)
        else:
            globDir = '%s/%s/*.dag.out'%(os.getcwd(),filename)
            logInfo('globDir: %s'%(globDir))
        dagOuts = glob.glob(globDir)
        #find pdb file
        if platform.system() == 'Windows':
            self.cwd = '%s\\%s'%(os.getcwd(),filename)
            pdb = glob.glob('%s\\%s\\%s.pdb'%(os.getcwd(),filename,filename))
        else:
            self.cwd = '%s/%s'%(os.getcwd(),info[0].filename)
            pdb = glob.glob('%s/%s/%s.pdb'%(os.getcwd(),filename,filename))
            logInfo('pdb: %s'%(pdb))
        if len(pdb) == 0:
            return -4, []
        #for each file in dagOuts
        for dag in dagOuts:
            #get base filename, looks complicated in case '.dag.out'
            #is present elsewhere is file path
            base = '.dag.out'.join(dag.split('.dag.out')[:len(dag.split('.dag.out'))-1])
            logInfo('base: %s'%(base))
            #is dag.png there?
            dagPng = len(glob.glob(base+'.dag.png'))==1
            #is dag.html there?
            dagHtml = len(glob.glob(base+'.dag.html'))==1
            #if so, append to output
            if dagPng and dagHtml:
                output.append(base)
        #no valid output
        logInfo('len(output): %i'%(len(output)))
        logInfo(output)
        if len(output) == 0:
            return -3, []
        #everything worked!
        oldIDind = len(self.seqWin.IDs)
        self.seqWin.PyMOLPDBLoad(1, pdb[0], "Show")
        newIDs = self.seqWin.IDs[oldIDind:]
        self.IDs = [] #tuples of the form (model,chain)
        for ID in newIDs:
            logInfo('newID: %s'%(ID))
            ID = (ID[:len(ID)-2],ID[len(ID)-1])
            self.IDs.append(ID)
            logInfo('ID added: (%s,%s)'%ID)
        logInfo('self.IDs')
        logInfo(self.IDs)
        logInfo('end self.IDs')
        native = ''
        for ID in self.IDs:
            native += 'model %s and chain %s+'%ID
        native = native[:len(native)-1]
        logInfo('native: %s'%(native))
        self.cmd.select('Native',native)
        self.cmd.show_as('cartoon','Native')
        self.cmd.set('cartoon_color','purple','Native')
        self.cmd.deselect()
        return 0,output


def startPyMOL(pdb):
    '''starts PyMOL for us.  Only for testing.  PyMOL should already be opened
    by InteractiveROSETTA'''
    import __main__
    __main__.pymol_argv = ["pymol", "-qhxi"]
    import pymol
    pymol.finish_launching()
    pymol.cmd.load(pdb)
    pymol.cmd.show_as('cartoon')
    pymol.cmd.color('purple')
    return pymol


if __name__ == '__main__':
    intermediates,transitions = parseImgMap('2b3p_florynewtest.21846_1.dag.html','2b3p_florynewtest.21846_1.dag.out')
    pymol = startPyMOL('2b3p_florynewtest.21846.pdb')
    '''for intermediate in intermediates:
        intermediate.show(pymol)
        time.sleep(0.1)
    transitions[0].show(intermediates)
    for transition in transitions:
        transition.show(intermediates)
        time.sleep(0.1)
    for intermediate in intermediates:
        if intermediate.number == 302:
            intermediate.show(pymol)'''
    app = wx.App(0)
    frame = wx.Frame(None,-1)
    testPanel = dagPanel(frame,'2b3p_florynewtest.21846_1.dag.png','2b3p_florynewtest.21846_1.dag.html','2b3p_florynewtest.21846_1.dag.out')
    testPanel.setPyMOL(pymol)
    frame.Show()
    app.MainLoop()
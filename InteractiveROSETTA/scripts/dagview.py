# -*- coding: utf-8 -*-
import sys
import math
import time
import wx
import shutil
import wx.lib.scrolledpanel
import platform
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
            
    
    def show(self,pymol):
        ''' Will display the intermediate in the pymol window'''
        #import pymol 
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
        pymol.cmd.show_as('ribbon')
        pymol.cmd.color('white')
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
        
    def show(self,intermediates,pymol):
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
        u1res = get_flag_residues(u1.iflag)
        if u2 != 0:
            u2res = get_flag_residues(u2.iflag)
        #pymol.cmd.select('f',fres)
        f.show(pymol)
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
            u1.show(pymol)
        

def parseImgMap(mapFile,dag=''):
    """This function takes the html imagemap file generated by GeoFold
    and uses it to create a list of Intermediate and transition states"""
    
    transitions = []
    intermediates = []
    readMap = open(mapFile,'r')
    print("reading mapfile.... %s"%(mapFile))
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

'''def findFiles(folderPath):
    None    '''
    
    
class dagPanel(wx.lib.scrolledpanel.ScrolledPanel):
    def __init__(self,parent, dagImg,dagMap,dagFile):        
        wx.lib.scrolledpanel.ScrolledPanel.__init__(self,parent,-1)
        self.intermediates,self.transitions = parseImgMap(dagMap,dagFile)
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
        if platform.system() != 'Darwin':
            (x,y) = self.CalcUnscrolledPosition(x,y)
        notFound = True
        for intermediate in self.intermediates:
            if intermediate.contains_point((x,y)):
                print "intermediate %d"%(intermediate.number)
                intermediate.show(self.pymol)
                notFound = False
                break
        if notFound:
            for transition in self.transitions:
                if transition.contains_point((x,y)):
                    print "transition %d"%(transition.number)
                    transition.show(self.intermediates,self.pymol)
                    notFound = False
                    break
        if notFound:
            print "notFound"

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
            self.lblDagView = wx.StaticText(self,-1,'Pathway Visualization',(25,15),(270,25),style=wx.AALIGN_CENTRE)
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
        self.btnLoad = wx.Button(self,id=-1,label='Load PDB',pos=(200,100),size=(110,25))
        self.btnLoad.SetForegroundColour("#000000")
        self.btnLoad.SetFont(wx.Font(10,wx.DEFAULT,wx.NORMAL,wx.BOLD))
            
        self.btnLoad.Bind(wx.EVT_BUTTON,self.loadPDB)
        self.btnLoad.SetToolTipString('Load the PDB file used for GeoFold job')
        logInfo('437: PDB button set!')
        
        #Dag.out button
        self.lblDagOut = wx.StaticText(self,-1,"None Uploaded", pos=(10,138),size=(180,25),style=wx.ALIGN_CENTRE)
        self.lblDagOut.SetFont(wx.Font(10,wx.DEFAULT,wx.NORMAL,wx.BOLD))
        if platform.system() == 'Linux':
            resizeTextControlForUNIX(self.lblDagOut,10,180)
        self.lblDagOut.SetForegroundColour("#FFFFFF")
        #if platform.system() == 'Darwin':
        #    self.btnDagOut = wx.BitmapButton(self,id=-1,bitmap=wx.Image(self.parent.parent.scriptdir+'/images/osx/dagview/btnDagOut.png',wx.BITMAP_TYPE_PNG).ConvertToBitmap(),pos=(200,135),size=(110,25))
        #else:
        self.btnDagOut = wx.Button(self,id=-1,label='Load .dag.out',pos=(200,135),size=(110,25))
        self.btnDagOut.SetForegroundColour("#000000")
        self.btnDagOut.SetFont(wx.Font(10,wx.DEFAULT,wx.NORMAL,wx.BOLD))
        
        self.btnDagOut.Bind(wx.EVT_BUTTON,self.loadDagOut)
        self.btnDagOut.SetToolTipString('Load the .dag.out file generated by GeoFold')
        logInfo('453: DagOut Button Set')
        
        #Dag.html button
        self.lblDagHtml = wx.StaticText(self,-1,"None Uploaded", pos=(10,173),size=(180,25),style=wx.ALIGN_CENTRE)
        self.lblDagHtml.SetFont(wx.Font(10,wx.DEFAULT,wx.NORMAL,wx.BOLD))
        if platform.system() == 'Linux':
            resizeTextControlForUNIX(self.lblDagHtml,10,180)
        self.lblDagHtml.SetForegroundColour("#FFFFFF")
        #if platform.system() == 'Darwin':
        #    self.btnDagHtml = wx.BitmapButton(self,id=-1,bitmap=wx.Image(self.parent.parent.scriptdir+'/images/osx/dagview/btnDagHtml.png',wx.BITMAP_TYPE_PNG).ConvertToBitmap(),pos=(200,170),size=(110,25))
        #else:
        self.btnDagHtml = wx.Button(self,id=-1,label='Load .dag.html',pos=(200,170),size=(110,25))
        self.btnDagHtml.SetForegroundColour("#000000")
        self.btnDagHtml.SetFont(wx.Font(10,wx.DEFAULT,wx.NORMAL,wx.BOLD))
        
        self.btnDagHtml.Bind(wx.EVT_BUTTON,self.loadDagHtml)
        self.btnDagHtml.SetToolTipString('Load the .dag.html file generated by GeoFold')
        logInfo('469: DagHtml Button set!')
        
        #Dag.png button
        self.lblDagPng = wx.StaticText(self,-1,"None Uploaded", pos=(10,208),size=(180,25),style=wx.ALIGN_CENTRE)
        self.lblDagPng.SetFont(wx.Font(10,wx.DEFAULT,wx.NORMAL,wx.BOLD))
        if platform.system() == 'Linux':
            resizeTextControlForUNIX(self.lblDagPng,10,180)
        self.lblDagPng.SetForegroundColour("#FFFFFF")
        #if platform.system() == 'Darwin':
        #    self.btnDagPng = wx.BitmapButton(self,id=-1,bitmap=wx.Image(self.parent.parent.scriptdir+'/images/osx/dagview/btnDagPng.png',wx.BITMAP_TYPE_PNG).ConvertToBitmap(),pos=(200,205),size=(110,25))
        #else:
        self.btnDagPng = wx.Button(self,id=-1,label='Load .dag.png',pos=(200,205),size=(110,25))
        self.btnDagPng.SetForegroundColour("#000000")
        self.btnDagPng.SetFont(wx.Font(10,wx.DEFAULT,wx.NORMAL,wx.BOLD))
        
        self.btnDagPng.Bind(wx.EVT_BUTTON,self.loadDagPng)
        self.btnDagPng.SetToolTipString('Load the .dag.png file generated by GeoFold')
        logInfo('485: DagPng Button Set!')        
        
        #GO! Button ViewDag
        ypos = self.btnDagPng.GetPosition()[1]+self.btnDagPng.GetSize()[1]+10
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
        logInfo('500: GetPosition = %d'%(self.btnViewDag.GetPosition()[1]))
        logInfo('501: GetSize = %d'%(self.btnViewDag.GetSize()[1]))
        self.scrollh = self.btnViewDag.GetPosition()[1] + self.btnViewDag.GetSize()[1] + 5
        logInfo('503: scrollh set to %d'%(self.scrollh))
        self.SetScrollbars(1,1,320,self.scrollh)
        logInfo('505: Scrollbars set!')
        self.winscrollpos = 0
        self.Bind(wx.EVT_SCROLLWIN, self.scrolled)
        logInfo('508: Scrolling set.')
        logInfo('509: Initialization complete!')
        
    def setSeqWin(self, seqWin):
        self.seqWin = seqWin
    
    def showHelp(self, event):
        '''Open the help page'''
        if platform.system() == 'Darwin':
            try:
                browser = webbrowser.get('Safari')
            except:
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
        self.cmd.show_as('cartoon')
        self.cmd.color('purple')
        self.intermediates,self.transitions = parseImgMap(self.loadedDagHtml,self.loadedDagOut)
        self.frame = wx.Frame(None,-1)
        self.DagPanel = dagPanel(self.frame,self.loadedDagPng,self.loadedDagHtml,self.loadedDagOut)
        self.DagPanel.setPyMOL(self.pymol)
        self.frame.Show()
        
    def findFiles(self,zipDir):
        '''Takes a given zip file.  extracts it in the sandbox and picks out all
        the files able to be viewed.  It outputs a list to be put in a ComboBox
        to allow the user to select which one to view.  If an error occurs, outputs
        a negative number used to identify the error and handle it'''
        import zipfile
        output = []
        #Check if selected file is valid
        if not zipfile.is_zipfile(zipDir):
            return -1, []
        #Unzip the file in the sandbox
        unzipped = zipfile.ZipFile(zipDir)
        info = unzipped.infolist()
        goToSandbox()
        try:
            unzipped.extractall()
        except: #failed to unzip
            return -2, []
        #use glob to get all dag.out files
        globDir = '%s/%s*.dag.out'%(os.getcwd(),info[0].filename)
        dagOuts = glob.glob(globDir)
        #for each file in dagOuts
        for dag in dagOuts:
            #get base filename, looks complicated in case '.dag.out' 
            #is present elsewhere is file path
            base = '.dag.out'.join(dag.split('.dag.out')[:len(dag.split('.dag.out'))-1])
            #is dag.png there?
            dagPng = len(glob.glob(base+'.dag.png'))==1
            #is dag.html there?
            dagHtml = len(glob.glob(base+'.dag.html'))==1
            #if so, append to output
            if dagPng and dagHtml:
                output.append(base)
        #no valid output
        if len(output) == 0:
            return -3, []
        #everything worked!
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
    print 'parsed map!'
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
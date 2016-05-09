import wx
import os
import os.path
import platform
import subprocess
import multiprocessing
import psutil
from daemon import daemonLoop
from tools import resizeTextControlForUNIX
from tools import logInfo
from tools import icon
from selection import SelectPanel
from superimposition import SuperimpositionPanel
from minimization import MinimizationPanel
from ensemblebrowser import EnsembleBrowserPanel
from fixbb import FixbbPanel
try:
    # This uses openbabel to convert PDBs to MOL2 files, which makes it way easier
    # to add new parameters files since we can just read the PDB that the user is trying
    # to load directly
    from advresiduecreator import ResidueCreatorPanel
except:
    # Maybe the user couldn't compile openbabel, so we'll default to the original version
    # that requires the user to generate their own MOL2 files
    print "A Python OpenBabel installation was not detected on your system."
    print "Although OpenBabel is not required, it can greatly simplify parameterizing new residues"
    print "On Windows, you need to install the main OpenBabel as well: http://openbabel.org/wiki/Category:Installation"
    print "On Debian, Python OpenBabel can be installed using apt-get: sudo apt-get install python-openbabel"
    print "On Mac and RedHat, you need to compile from source: http://open-babel.readthedocs.org/en/latest/Installation/install.html#compiling-open-babel"
    from residuecreator import ResidueCreatorPanel
from pointmutations import PointMutationsPanel
from kic import KICPanel
from docking import DockingPanel
from msd import MSDPanel
from pmutscan import PointMutantScanPanel 
from surfaces import SurfacesPanel
from ensemblegen import EnsembleGenPanel

class ProtocolsPanel(wx.Panel):
    def __init__(self, parent, W, H):
        wx.Panel.__init__(self, parent, id=-1, pos=(0, 0), size=(350, H-265), name="ProtocolsPanel")
        self.W = W
        self.H = H
        self.parent = parent
        self.SetBackgroundColour("#333333")
        self.Show()
        
        if (platform.system() == "Windows"):
            self.label = wx.StaticText(self, -1, "Protocols", (5, 5), (340, 25), wx.ALIGN_CENTRE)
            self.label.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.label = wx.StaticBitmap(self, -1, wx.Image(self.parent.scriptdir + "/images/osx/lblProtocols.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(5, 5), size=(340, 25))
        else:
            self.label = wx.StaticText(self, -1, "Protocols", style=wx.ALIGN_CENTRE)
            self.label.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
            resizeTextControlForUNIX(self.label, 0, self.GetSize()[0])
        self.label.SetForegroundColour("#FFFFFF")
        
        self.protocols = ["Docking", 
                   "Energy Minimization", 
                   "Ensemble Browser", 
                   "Ensemble Generation",
                   "Loop Modeling (KIC)",
                   "Molecular Surfaces", 
                   "Point Mutant Scan", 
                   "Point Mutations", 
                   "Protein Design (Fixbb)", 
                   "Protein Design (MSD)", 
                   "Residue/Ligand Creator", 
                   "Superimposition"]
        if (platform.system() == "Darwin"):
            self.protMenu = wx.ComboBox(self, pos=(5, 30), size=(230, 25), choices=self.protocols, style=wx.CB_READONLY)
        else:
            self.protMenu = wx.ComboBox(self, pos=(5, 30), size=(230, 25), choices=self.protocols, style=wx.CB_READONLY | wx.CB_SORT)
        self.protMenu.SetSelection(self.protocols.index("Superimposition"))
        self.protMenu.SetToolTipString("List of currently available protocols")
        
        if (platform.system() == "Darwin"):
            self.GoBtn = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/GoBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(240, 30), size=(100, 25))
        else:
            self.GoBtn = wx.Button(self, id=-1, label="Go!", pos=(240, 30), size=(100, 25))
            self.GoBtn.SetForegroundColour("#000000")
            self.GoBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.GoBtn.Bind(wx.EVT_BUTTON, self.changeProtocol)
        self.GoBtn.SetToolTipString("Change to the selected protocol")
        self.currentProtocol = "Superimposition"
        
        self.protPanel = SuperimpositionPanel(self, W, H)
        
    def changeProtocol(self, event):
        logInfo("Go button clicked")
        selectedProtocol = self.protMenu.GetStringSelection()
        if (selectedProtocol != self.currentProtocol and selectedProtocol != ""):
            # Destroy the old panel so we can free memory up and create the new protocol panel
            if (self.currentProtocol == "Superimposition"):
                self.protPanel.Destroy()
                del self.protPanel
            elif (self.currentProtocol == "Energy Minimization"):
                # Check to see if the user has accepted the minimization
                # If not, ask if they really want to proceed
                if (self.protPanel.buttonState == "Finalize!"):
                    dlg = wx.MessageDialog(self, "You have not accepted your minimization and the results will be lost if you proceed.  Proceed anyway?", "Minimization Not Accepted", wx.YES_NO | wx.ICON_EXCLAMATION | wx.CENTRE)
                    if (dlg.ShowModal() == wx.ID_NO):
                        logInfo("Go cancelled due to unaccepted minimization")
                        dlg.Destroy()
                        return
                    dlg.Destroy()
                    logInfo("Minimization job rejected")
                # Try to delete the minimized model if the user doesn't want to accept it
                try:
                    self.cmd.remove("minimized_view")
                    self.cmd.delete("minimized_view")
                except:
                    pass
                self.cmd.label("all", "")
                self.protPanel.Destroy()
                del self.protPanel
            elif (self.currentProtocol == "Ensemble Browser"):
                self.protPanel.Destroy()
                del self.protPanel
            elif (self.currentProtocol == "Ensemble Generation"):
                # Check to see if the user has accepted the design
                # If not, ask if they really want to proceed
                if (self.protPanel.buttonState == "Cancel!"):
                    dlg = wx.MessageDialog(self, "Your ensemble is not finished and the results will be lost if you proceed.  Proceed anyway?", "Design Not Accepted", wx.YES_NO | wx.ICON_EXCLAMATION | wx.CENTRE)
                    if (dlg.ShowModal() == wx.ID_NO):
                        logInfo("Go cancelled due to unfinished ensemblegen job")
                        dlg.Destroy()
                        return
                    logInfo("Ensemblegen job rejected")
                    dlg.Destroy()
                self.protPanel.Destroy()
                del self.protPanel
            elif (self.currentProtocol == "Protein Design (Fixbb)"):
                # Check to see if the user has accepted the design
                # If not, ask if they really want to proceed
                if (self.protPanel.buttonState == "Finalize!"):
                    dlg = wx.MessageDialog(self, "You have not accepted your design and the results will be lost if you proceed.  Proceed anyway?", "Design Not Accepted", wx.YES_NO | wx.ICON_EXCLAMATION | wx.CENTRE)
                    if (dlg.ShowModal() == wx.ID_NO):
                        logInfo("Go cancelled due to unaccepted fixbb job")
                        dlg.Destroy()
                        return
                    logInfo("Fixbb job rejected")
                    dlg.Destroy()
                # Try to delete the design model if the user doesn't want to accept it
                try:
                    self.cmd.remove("designed_view")
                    self.cmd.delete("designed_view")
                except:
                    pass
                self.cmd.label("all", "")
                self.protPanel.Destroy()
                del self.protPanel
            elif (self.currentProtocol == "Residue/Ligand Creator"):
                try:
                    self.cmd.remove("params")
                    self.cmd.delete("params")
                except:
                    pass
                self.protPanel.Destroy()
                del self.protPanel
            elif (self.currentProtocol == "Point Mutations"):
                try:
                    self.cmd.remove("rotamer_view")
                    self.cmd.delete("rotamer_view")
                except:
                    pass
                self.protPanel.Destroy()
                del self.protPanel
            elif (self.currentProtocol == "Point Mutant Scan"):
                self.protPanel.Destroy()
                del self.protPanel
            elif (self.currentProtocol == "Loop Modeling (KIC)"):
                # Check to see if the user has accepted the loop model
                # If not, ask if they really want to proceed
                if (self.protPanel.buttonState == "Finalize!"):
                    dlg = wx.MessageDialog(self, "You have not accepted your loop model and the results will be lost if you proceed.  Proceed anyway?", "KIC Not Accepted", wx.YES_NO | wx.ICON_EXCLAMATION | wx.CENTRE)
                    if (dlg.ShowModal() == wx.ID_NO):
                        logInfo("Go cancelled due to unaccepted KIC job")
                        dlg.Destroy()
                        return
                    logInfo("KIC job rejected")
                    dlg.Destroy()
                # Try to delete the design model if the user doesn't want to accept it
                try:
                    self.cmd.remove("kic_view")
                    self.cmd.delete("kic_view")
                except:
                    pass
                self.cmd.label("all", "")
                self.protPanel.Destroy()
            elif (self.currentProtocol == "Docking"):
                # Check to see if the user has accepted the docking model
                # If not, ask if they really want to proceed
                if (self.protPanel.buttonState != "Dock!"):
                    dlg = wx.MessageDialog(self, "You have not accepted your docking model and the results will be lost if you proceed.  Proceed anyway?", "Docking Not Accepted", wx.YES_NO | wx.ICON_EXCLAMATION | wx.CENTRE)
                    if (dlg.ShowModal() == wx.ID_NO):
                        logInfo("Go cancelled due to unaccepted docking job")
                        dlg.Destroy()
                        return
                    logInfo("Docking job rejected")
                    dlg.Destroy()
                # Try to delete the design model if the user doesn't want to accept it
                try:
                    self.cmd.remove("dock_view")
                    self.cmd.delete("dock_view")
                except:
                    pass
                self.cmd.label("all", "")
                self.protPanel.Destroy()
            elif (self.currentProtocol == "Structure Prediction (Comparative Modeling)"):
                # Check to see if the user has accepted the docking model
                # If not, ask if they really want to proceed
                if (self.protPanel.buttonState == "Finalize!"):
                    dlg = wx.MessageDialog(self, "You have not accepted your comparative modeling structure and the results will be lost if you proceed.  Proceed anyway?", "Docking Not Accepted", wx.YES_NO | wx.ICON_EXCLAMATION | wx.CENTRE)
                    if (dlg.ShowModal() == wx.ID_NO):
                        logInfo("Go cancelled due to unaccepted comparative modeling job")
                        dlg.Destroy()
                        return
                    logInfo("Comparative modeling job rejected")
                    dlg.Destroy()
                # Try to delete the design model if the user doesn't want to accept it
                try:
                    self.cmd.remove("thread_view")
                    self.cmd.delete("thread_view")
                except:
                    pass
                self.cmd.label("all", "")
                self.protPanel.Destroy()
            elif (self.currentProtocol == "Protein Design (MSD)"):
                self.protPanel.Destroy()
                del self.protPanel
            elif (self.currentProtocol == "Antibody Modeling"):
                self.protPanel.Destroy()
                del self.protPanel
            elif (self.currentProtocol == "Molecular Surfaces"):
                try:
                    self.cmd.delete("curr_surf_recp")
                except:
                    pass
                try:
                    self.cmd.delete("curr_surf_lig")
                except:
                    pass
                self.cmd.hide("surface", "all")
                if (self.parent.Selection.showSurf):
                    self.currentProtocol = "n/a"
                    self.parent.Selection.displaySurfaces()
                self.protPanel.Destroy()
                del self.protPanel
            elif (self.currentProtocol == "Flexible Peptide Docking"):
                self.protPanel.Destroy()
                del self.protPanel
            self.currentProtocol = selectedProtocol
            self.seqWin.cannotDelete = False
            # Restart the Rosetta daemon to clear its memory up
            self.parent.restartDaemon()
            logInfo("Changed protocol to " + selectedProtocol)
            if (selectedProtocol == "Superimposition"):
                self.protPanel = SuperimpositionPanel(self, self.W, self.H)
            elif (selectedProtocol == "Energy Minimization"):
                self.protPanel = MinimizationPanel(self, self.W, self.H)
                self.protPanel.setSelectWin(self.selectWin)
            elif (selectedProtocol == "Ensemble Browser"):
                self.protPanel = EnsembleBrowserPanel(self, self.W, self.H)
            elif (selectedProtocol == "Ensemble Generation"):
                self.protPanel = EnsembleGenPanel(self, self.W, self.H)
                self.protPanel.setSelectWin(self.selectWin)
            elif (selectedProtocol == "Protein Design (Fixbb)"):
                self.protPanel = FixbbPanel(self, self.W, self.H)
                self.protPanel.setSelectWin(self.selectWin)
            elif (selectedProtocol == "Residue/Ligand Creator"):
                self.protPanel = ResidueCreatorPanel(self, self.W, self.H)
            elif (selectedProtocol == "Point Mutations"):
                self.protPanel = PointMutationsPanel(self, self.W, self.H)
                self.protPanel.setSelectWin(self.selectWin)
            elif (selectedProtocol == "Point Mutant Scan"):
                self.protPanel = PointMutantScanPanel(self, self.W, self.H)
                self.protPanel.setSelectWin(self.selectWin)
            elif (selectedProtocol == "Loop Modeling (KIC)"):
                self.protPanel = KICPanel(self, self.W, self.H)
                self.protPanel.setSelectWin(self.selectWin)
            elif (selectedProtocol == "Docking"):
                self.protPanel = DockingPanel(self, self.W, self.H)
                self.protPanel.setSelectWin(self.selectWin)
            elif (selectedProtocol == "Structure Prediction (Comparative Modeling)"):
                self.protPanel = CompModelPanel(self, self.W, self.H)
                self.protPanel.setSelectWin(self.selectWin)
            elif (selectedProtocol == "Protein Design (MSD)"):
                self.protPanel = MSDPanel(self, self.W, self.H)
                self.protPanel.setSelectWin(self.selectWin)
            elif (selectedProtocol == "Antibody Modeling"):
                self.protPanel = AntibodyPanel(self, self.W, self.H)
            elif (selectedProtocol == "Molecular Surfaces"):
                self.protPanel = SurfacesPanel(self, self.W, self.H)
                self.cmd.hide("surface", "all")
            elif (selectedProtocol == "Flexible Peptide Docking"):
                self.protPanel = FlexPepDockPanel(self, self.W, self.H)
                self.protPanel.setSelectWin(self.selectWin)
            self.protPanel.setSeqWin(self.seqWin)
            self.protPanel.setPyMOL(self.pymol)
            self.protPanel.activate()
            self.seqWin.setProtocolPanel(self.protPanel)
    
    def setSeqWin(self, seqWin):
        self.seqWin = seqWin
        self.protPanel.setSeqWin(seqWin)
        self.seqWin.setProtocolPanel(self.protPanel)
        
    def setPyMOL(self, pymol):
        self.pymol = pymol
        self.cmd = pymol.cmd
        self.stored = pymol.stored
        self.protPanel.setPyMOL(pymol)
        
    def setSelectWin(self, selectWin):
        self.selectWin = selectWin
        
    def activate(self):
        self.cmd.enable("seqsele")
        self.protPanel.activate()

class ProtocolsWin(wx.Frame):
    def __init__(self, W, H, scriptdir):
        if (platform.system() == "Darwin"):
            self.stdwinx = 0; self.stdwiny = 24
        else:
            self.stdwinx = 0; self.stdwiny = 0
        self.stdwinw = 370; self.stdwinh = H-40
        self.screenH = H; self.screenW = W
        winx = self.stdwinx; winy = self.stdwiny
        winw = self.stdwinw; winh = self.stdwinh
        self.scriptdir = scriptdir
        homedir = os.path.expanduser("~")
        # Try to get the save values from the cfg file
        try:
            if (platform.system() == "Windows"):
                f = open(homedir + "\\InteractiveROSETTA\\protwindow.cfg", "r")
            else:
                f = open(homedir + "/InteractiveROSETTA/protwindow.cfg", "r")
            for aline in f:
                if (aline.find("[OFFSET X]") >= 0):
                    winx = winx + int(aline.split()[len(aline.split())-1])
                elif (aline.find("[OFFSET Y]") >= 0):
                    winy = winy + int(aline.split()[len(aline.split())-1])
                elif (aline.find("[OFFSET WIDTH]") >= 0):
                    winw = winw + int(aline.split()[len(aline.split())-1])
                elif (aline.find("[OFFSET HEIGHT]") >= 0):
                    winh = winh + int(aline.split()[len(aline.split())-1])
            f.close()
        except:
            pass
        if (winx > self.screenW - 100):
            winx = self.stdwinx
        if (winy > self.screenH - 100):
            winy = self.stdwiny
        # Maybe the screen resolution has changed and the saved dimensions put the windows in
        # weird places, so default them to better positions and the user can change them later
        #if (winw < 350):
        #    winw = 370
        #elif (winw > W):
        #    winw = 370
        #if (winx < 0):
        #    winx = 0
        #elif (winx > W-winw):
        #    winx = W-winw
        #if (winh > H - 40):
        #    winh = H - 40
        #if (winy < 0):
        #    winy = 0
        #elif (winy > H-winh):
        #    winh = H-40
        wx.Frame.__init__(self, None, -1, "InteractiveROSETTA - Protocols", size=(winw, winh))
        self.SetPosition((winx, winy))
        self.SetBackgroundColour("#333333")
        self.SetSizeHints(330, 560, 370, H)
        self.SetIcon(icon.GetIcon())
        self.Show()
        
        self.Protocols = ProtocolsPanel(self, winw, winh)
        self.Selection = SelectPanel(self, winw, winh)
        self.Protocols.setSelectWin(self.Selection)
        self.Selection.setProtPanel(self.Protocols)
        
        self.saveTimer = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, self.saveWindowData, self.saveTimer)
        self.Bind(wx.EVT_SIZE, self.windowGeometryChange)
        self.Bind(wx.EVT_MOTION, self.windowGeometryChange)
        self.Bind(wx.EVT_ACTIVATE, self.focusEvent)
        
        # Start the Rosetta daemon that will run in the background looking for job input files generated
        # by the main GUI
        # It could be the case that the user was in the middle of a protocol, then quits suddenly so the
        # daemon doesn't terminate itself because it's in the middle of a protocol
        # The the user restarts the GUI before the daemon has a chance to finish the protocol
        # This checks to see if the daemon is already active and doesn't spawn a new one if one is already
        # running
        stillrunning = False
        count = 0
        for proc in psutil.process_iter():
            try:
                if (platform.system() == "Windows"):
                    if (len(proc.cmdline()) >= 3 and proc.cmdline()[0].find("python") >= 0 and proc.cmdline()[2].find("from multiprocessing.forking") >= 0):
                        stillrunning = True
                        break
                else:
                    # On Unix systems you just have to make sure two instances of python are running
                    # because there isn't any forking information in the daemon's instance of python
                    if (len(proc.cmdline()) >= 2 and proc.cmdline()[0].find("python") >= 0 and proc.cmdline()[1].find("InteractiveROSETTA.py") >= 0):
                        count = count + 1
            except:
                # In Windows it will crash if you try to read process information for the Administrator
                # Doesn't matter though since InteractiveROSETTA is run by a non-Administrator
                # But we need to catch these errors since we don't know which processes are admin ones
                pass
        if (platform.system() != "Windows" and count == 2):
            stillrunning = True
        if (not(stillrunning)):
            print "Starting Rosetta protocol daemon..."
            #if (platform.system() == "Darwin"):
                #self.daemon_process = subprocess.Popen(args=["/usr/bin/python", self.scriptdir + "/scripts/daemon.py"], shell=False)
            #else:
            self.daemon_process = multiprocessing.Process(target=daemonLoop)
            self.daemon_process.start()
        
    def restartDaemon(self):
        # This function kills the current daemon process and starts a new one
        # This function is called whenever the protocol panel changes to a different protocol
        # This function is necessary because on Windows PyRosetta can start to use well over
        # 4GB of memory if the daemon has loaded memory from multiple protocols
        # Killing the daemon and restarting it clears up that memory so the user's computer
        # doesn't slow down as it moves a lot of stuff into swap space
        savedir = os.getcwd()
        os.chdir(self.scriptdir)
        self.daemon_process.terminate()
        print "Restarting Rosetta protocol daemon..."
        self.daemon_process = multiprocessing.Process(target=daemonLoop)
        self.daemon_process.start()
        os.chdir(savedir)
        
    def focusEvent(self, event):
        if (event.GetActive()):
            # Protocols read selection information from the sequence window, so update the sequence window
            # If we're going from PyMOL->Protocols, since updates usually happen on the sequence window focus event
            self.seqWin.selectUpdate(False) # Update PyMOL changes in sequence
            self.Protocols.activate()
        event.Skip()
    
    def setSeqWin(self, seqWin):
        self.seqWin = seqWin
        self.Protocols.setSeqWin(seqWin)
        self.Selection.setSeqWin(seqWin)
        
    def windowGeometryChange(self, event):
        # This function starts a timer that will write out the size and position of this window to a cfg file
        # so the orientation is saved and can be loaded the next time InteractiveROSETTA is started
        if (not(self.saveTimer.IsRunning())):
            self.saveTimer.Start(5000)
        # We have to do some finagling if this window gets resized because the minimum height
        # of the specific protocol panel needs to be at least 300
        # If the window is shrunk, the protocol panel will shrink down to 300 but then the select
        # panel needs to shrink
        (w, h) = self.GetSize()
        return
        if (h > 560):
            try:
                self.Protocols.protPanel.SetSize((w, h-330))
                self.Protocols.protPanel.SetScrollbars(1, 1, 320, 800)
            except:
                pass
            (x, y) = self.Selection.GetPosition()
            self.Selection.SetPosition((x, h-270))
            self.Selection.SetSize((w-20, self.Selection.GetSize()[1]))
        event.Skip()
        
    def saveWindowData(self, event):
        homedir = os.path.expanduser("~")
        data = []
        try:
            if (platform.system() == "Windows"):
                f = open(homedir + "\\InteractiveROSETTA\\protwindow.cfg", "r")
            else:
                f = open(homedir + "/InteractiveROSETTA/protwindow.cfg", "r")
            for aline in f:
                data.append(aline)
            f.close()
        except:
            pass
        if (platform.system() == "Windows"):
            f = open(homedir + "\\InteractiveROSETTA\\protwindow.cfg", "w")
        else:
            f = open(homedir + "/InteractiveROSETTA/protwindow.cfg", "w")
        itemsFound = [False, False, False, False] # [offX, offY, offW, offH]
        (x, y) = self.GetPosition()
        (w, h) = self.GetSize()
        for aline in data:
            if (aline.find("[OFFSET X]") >= 0):
                itemsFound[0] = True
                f.write("[OFFSET X] " + str(x-self.stdwinx) + "\n")
            elif (aline.find("[OFFSET Y]") >= 0):
                itemsFound[1] = True
                f.write("[OFFSET Y] " + str(y-self.stdwiny) + "\n")
            elif (aline.find("[OFFSET WIDTH]") >= 0):
                itemsFound[2] = True
                f.write("[OFFSET WIDTH] " + str(w-self.stdwinw) + "\n")
            elif (aline.find("[OFFSET HEIGHT]") >= 0):
                itemsFound[3] = True
                f.write("[OFFSET HEIGHT] " + str(h-self.stdwinh) + "\n")
            else:
                f.write(aline)
        for i in range(0, len(itemsFound)):
            if (not(itemsFound[i])):
                if (i == 0):
                    f.write("[OFFSET X] " + str(x-self.stdwinx) + "\n")
                elif (i == 1):
                    f.write("[OFFSET Y] " + str(y-self.stdwiny) + "\n")
                elif (i == 2):
                    f.write("[OFFSET WIDTH] " + str(w-self.stdwinw) + "\n")
                elif (i == 3):
                    f.write("[OFFSET HEIGHT] " + str(h-self.stdwinh) + "\n")
        f.close()
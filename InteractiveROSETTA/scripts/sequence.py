import wx
import wx.grid
import os
import os.path
import glob
import platform
import urllib
import urllib2
import gzip
import math
import ctypes
import webbrowser
from threading import Thread
import Bio.PDB
import Bio.PDB.DSSP
from Bio.Align.Applications import MuscleCommandline
from tools import *

# There is apparently a libc bug on UNIX that caches DNS names that urllib2 tries to open
# If the Internet connection is not valid the first time it tries to establish a connection, it caches this bad
# DNS and never is able to connect in this session, even if the Internet connection comes back
# You have to get a handle on "res_init" in libc and call it before attempting to download more things
if (platform.system() == "Linux"):
    libc = ctypes.cdll.LoadLibrary("libc.so.6")
    res_init = libc.__res_init

class ProteinDialog(wx.Dialog):
    def __init__(self, parent, PDB, chains):
	if (platform.system() != "Linux"):
	    wx.Dialog.__init__(self, parent, -1, "Protein Loader", size=(305, 330))
	else:
	    wx.Dialog.__init__(self, parent, -1, "Protein Loader", size=(300, 300))
	
	self.lblPDB = wx.StaticText(self, -1, PDB, (0, 10), (300, 30))
	self.lblPDB.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
	resizeTextControlForUNIX(self.lblPDB, 0, 300)
	midpoint = int(math.ceil(len(chains) / 2.0))
	
	ygap = int(160 / midpoint)
	self.chkChains = []
	for i in range(0, len(chains)):
	    if (i < midpoint):
		self.chkChains.append(wx.CheckBox(self, -1, "Chain " + chains[i], (20, 40+(ygap*i))))
	    else:
		self.chkChains.append(wx.CheckBox(self, -1, "Chain " + chains[i], (170, 40+(ygap*(i-midpoint)))))
	    self.chkChains[i].SetValue(True)
	
	self.btnWater = wx.Button(self, id=-1, label="Load Waters", pos=(5, 230), size=(140, 30))
	self.btnWater.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.btnWater.Bind(wx.EVT_BUTTON, self.waterToggle)
	self.btnWater.SetToolTipString("Load the protein with waters present")
	self.btnHETATM = wx.Button(self, id=-1, label="Load HETATMs", pos=(155, 230), size=(140, 30))
	self.btnHETATM.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.btnHETATM.Bind(wx.EVT_BUTTON, self.HETATMToggle)
	self.btnHETATM.SetToolTipString("Load the protein with heteroatoms present")
	
	self.btnOK = wx.Button(self, id=-1, label="OK", pos=(5, 265), size=(140, 30))
	self.btnOK.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.btnOK.Bind(wx.EVT_BUTTON, self.okDialog)
	self.btnOK.SetToolTipString("Confirm load operation")
	self.btnCancel = wx.Button(self, id=-1, label="Cancel", pos=(155, 265), size=(140, 30))
	self.btnCancel.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.btnCancel.Bind(wx.EVT_BUTTON, self.cancelDialog)
	self.btnCancel.SetToolTipString("Cancel this load operation")
	
	self.SetPosition((wx.GetDisplaySize()[0]/2-150, wx.GetDisplaySize()[1]/2-150))
	
    def waterToggle(self, event):
	if (self.btnWater.GetLabel() == "Load Waters"):
	    self.btnWater.SetLabel("No Waters")
	    self.btnWater.SetToolTipString("Load the protein while ignoring waters")
	else:
	    self.btnWater.SetLabel("Load Waters")
	    self.btnWater.SetToolTipString("Load the protein with waters present")
    
    def HETATMToggle(self, event):
	if (self.btnHETATM.GetLabel() == "Load HETATMs"):
	    self.btnHETATM.SetLabel("No HETATMs")
	    self.btnHETATM.SetToolTipString("Load the protein while ignoring heteroatoms")
	else:
	    self.btnHETATM.SetLabel("Load HETATMs")
	    self.btnHETATM.SetToolTipString("Load the protein with heteroatoms present")
	    
    def okDialog(self, event):
	self.SetReturnCode(wx.OK)
	self.EndModal(wx.OK)
	
    def cancelDialog(self, event):
	self.SetReturnCode(wx.CANCEL)
	self.EndModal(wx.CANCEL)

class SequenceWin(wx.Frame):
    def __init__(self, W, H, cwd, frozen, poses, sequences, IDs, scriptdir):
	self.stdwinx = 370; self.stdwiny = H-270
	self.stdwinw = W-370; self.stdwinh = 270
	self.screenH = H; self.screenW = W
	#if (platform.system() == "Linux"):
	#    winy = winy + 22
	#    winh = winh - 22
	#    winx = winx + 5
	#    winw = winw - 5
	winx = self.stdwinx; winy = self.stdwiny
	winw = self.stdwinw; winh = self.stdwinh
	homedir = os.path.expanduser("~")
	self.scriptdir = scriptdir
	self.cwd = cwd
	# Try to get the save values from the cfg file
	try:
	    if (platform.system() == "Windows"):
		f = open(homedir + "\\InteractiveROSETTA\\seqwindow.cfg", "r")
	    else:
		f = open(homedir + "/InteractiveROSETTA/seqwindow.cfg", "r")
	    for aline in f:
		if (aline.find("[OFFSET X]") >= 0):
		    winx = winx + int(aline.split()[len(aline.split())-1])
		elif (aline.find("[OFFSET Y]") >= 0):
		    winy = winy + int(aline.split()[len(aline.split())-1])
		elif (aline.find("[OFFSET WIDTH]") >= 0):
		    winw = winw + int(aline.split()[len(aline.split())-1])
		elif (aline.find("[OFFSET HEIGHT]") >= 0):
		    winh = winh + int(aline.split()[len(aline.split())-1])
		elif (aline.find("[CWD]") >= 0):
		    self.cwd = aline.split("\t")[1].strip()
		elif (aline.find("[SERVER]") >= 0):
		    setServerName(aline.split("\t")[1].strip())
	    f.close()
	except:
	    pass
	if (winx > self.screenW - 100):
	    winx = self.stdwinx
	if (winy > self.screenH - 100):
	    winy = self.stdwiny
	# Maybe the screen resolution has changed and the saved dimensions put the windows in
	# weird places, so default them to better positions and the user can change them later
	#if (winw > W-370):
	#    winw = W-370
	#if (winx < 0):
	#    winx = 370
	#elif (winx > W-winw):
	#    winx = 370
	#if (winh < 230):
	#    winh = 250
	#elif (winh > H):
	#    winh = 270
	#if (winy < 0):
	#    winy = H-310
	#elif (winy > H-winh-40):
	#    winy = H-310
	wx.Frame.__init__(self, None, -1, "InteractiveROSETTA - Sequence Viewer", size=(winw, winh))
	self.frozen = frozen
	self.poses = poses
	self.sequences = sequences
	self.indxToSeqPos = []
	self.IDs = IDs
	self.SetPosition((winx, winy))
	self.SetBackgroundColour("#333333")
	#self.SetBackgroundColour("#333333")
	self.SetIcon(icon.GetIcon())
	#self.sizer = wx.GridBagSizer(1, 1)
	
	#self.SequencePanel = wx.Panel(self, id=-1, pos=(10, 10), size=(winw-10, winh-10))
	
	self.scroll = wx.ScrolledWindow(self, -1)
	self.scroll.SetBackgroundColour("#333333")
	self.scroll.SetSize((winw, winh))
	
	if (platform.system() == "Darwin"):
	    self.LoadPDBsBtn = wx.BitmapButton(self.scroll, id=-1, bitmap=wx.Image("images/osx/LoadPDBsBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(10, 10), size=(100, 25))
	else:
	    self.LoadPDBsBtn = wx.Button(self.scroll, id=-1, label="Load PDBs", pos=(10, 10), size=(100, 25))
	    #self.LoadPDBsBtn.SetBackgroundColour("#000000")
	    self.LoadPDBsBtn.SetForegroundColour("#000000")
	    self.LoadPDBsBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.LoadPDBsBtn.Bind(wx.EVT_BUTTON, self.loadPDBsClick)
	self.LoadPDBsBtn.SetToolTipString("Load new PDBs into PyMOL")
	#self.sizer.Add(self.LoadPDBsBtn, (0, 0), span=(1, 1), flag = wx.ALIGN_CENTER | wx.EXPAND, border=5)
	
	if (platform.system() == "Darwin"):
	    self.labelRCSB = wx.StaticBitmap(self.scroll, -1, wx.Image("images/osx/labelRSCB.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(115, 13), size=(43, 25))
	else:
	    self.labelRCSB = wx.StaticText(self.scroll, -1, "RCSB:", (115, 13), (43, 25))
	    self.labelRCSB.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
	    self.labelRCSB.SetForegroundColour("#FFFFFF")
	self.RCSBTxt = wx.TextCtrl(self.scroll, -1, pos=(158, 10), size=(50, 25))
	self.RCSBTxt.SetValue("")
	self.RCSBTxt.SetToolTipString("Four letter PDB code to search for in the RCSB database")
	
	if (platform.system() == "Darwin"):
	    self.FetchBtn = wx.BitmapButton(self.scroll, id=-1, bitmap=wx.Image("images/osx/FetchPDBBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(210, 10), size=(100, 25))
	else:
	    self.FetchBtn = wx.Button(self.scroll, id=-1, label="Fetch PDB", pos=(210, 10), size=(100, 25))
	    self.FetchBtn.SetForegroundColour("#000000")
	    self.FetchBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.FetchBtn.Bind(wx.EVT_BUTTON, self.fetchClick)
	self.FetchBtn.SetToolTipString("Fetch PDB code from RCSB")
	
	if (platform.system() == "Darwin"):
	    self.CloseBtn = wx.BitmapButton(self.scroll, id=-1, bitmap=wx.Image("images/osx/CloseBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(310, 10), size=(100, 25))
	else:
	    self.CloseBtn = wx.Button(self.scroll, id=-1, label="Close", pos=(310, 10), size=(100, 25))
	    #self.CloseBtn.SetBackgroundColour("#000000")
	    self.CloseBtn.SetForegroundColour("#000000")
	    self.CloseBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.CloseBtn.Bind(wx.EVT_BUTTON, self.closeClick)
	self.CloseBtn.SetToolTipString("Close selected/all PDBs")
	#self.sizer.Add(self.CloseBtn, (0, 1), span=(1, 1), flag = wx.ALIGN_CENTER | wx.EXPAND, border=5)
	
	if (platform.system() == "Darwin"):
	    self.SaveBtn = wx.BitmapButton(self.scroll, id=-1, bitmap=wx.Image("images/osx/SaveBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(410, 10), size=(100, 25))
	else:
	    self.SaveBtn = wx.Button(self.scroll, id=-1, label="Save PDB", pos=(410, 10), size=(100, 25))
	    self.SaveBtn.SetForegroundColour("#000000")
	    self.SaveBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.SaveBtn.Bind(wx.EVT_BUTTON, self.saveClick)
	self.SaveBtn.SetToolTipString("Save selected/all loaded models")
	#self.sizer.Add(self.SaveBtn, (0, 2), span=(1, 1), flag = wx.ALIGN_CENTER | wx.EXPAND, border=5)
	
	if (platform.system() == "Darwin"):
	    self.SaveImageBtn = wx.BitmapButton(self.scroll, id=-1, bitmap=wx.Image("images/osx/SaveImageBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(510, 10), size=(100, 25))
	else:
	    self.SaveImageBtn = wx.Button(self.scroll, id=-1, label="Save Image", pos=(510, 10), size=(100, 25))
	    self.SaveImageBtn.SetForegroundColour("#000000")
	    self.SaveImageBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.SaveImageBtn.Bind(wx.EVT_BUTTON, self.saveImage)
	self.SaveImageBtn.SetToolTipString("Save the current PyMOL view to a PNG image")
	#self.sizer.Add(self.SaveImageBtn, (0, 3), span=(1, 1), flag = wx.ALIGN_CENTER | wx.EXPAND, border=5)
	
	if (platform.system() == "Darwin"):
	    self.JoinBtn = wx.BitmapButton(self.scroll, id=-1, bitmap=wx.Image("images/osx/JoinChainsBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(610, 10), size=(100, 25))
	else:
	    self.JoinBtn = wx.Button(self.scroll, id=-1, label="Join Chains", pos=(610, 10), size=(100, 25))
	    self.JoinBtn.SetForegroundColour("#000000")
	    self.JoinBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.JoinBtn.Bind(wx.EVT_BUTTON, self.joinChains)
	self.JoinBtn.SetToolTipString("Join the selected chains into a single chain")
	
	if (platform.system() == "Darwin"):
	    self.RenumberBtn = wx.BitmapButton(self.scroll, id=-1, bitmap=wx.Image("images/osx/RenumberBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(710, 10), size=(100, 25))
	else:
	    self.RenumberBtn = wx.Button(self.scroll, id=-1, label="Renumber", pos=(710, 10), size=(100, 25))
	    self.RenumberBtn.SetForegroundColour("#000000")
	    self.RenumberBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.RenumberBtn.Bind(wx.EVT_BUTTON, self.renumber)
	self.RenumberBtn.SetToolTipString("Renumber the selected chain from 1, using the selected residue as the new start")
	
	if (platform.system() == "Darwin"):
	    self.ServerBtn = wx.BitmapButton(self.scroll, id=-1, bitmap=wx.Image("images/osx/ServerBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(max(810, W-445), 10), size=(25, 25))
	else:
	    self.ServerBtn = wx.Button(self.scroll, id=-1, label="S", pos=(max(810, W-445), 10), size=(25, 25))
	    self.ServerBtn.SetForegroundColour("#000000")
	    self.ServerBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.ServerBtn.Bind(wx.EVT_BUTTON, self.configureServer)
	self.ServerBtn.SetToolTipString("Connect to a remote server running PyRosetta and Rosetta (this is required for some protocols)")
	if (platform.system() == "Darwin"):
	    self.HelpBtn = wx.BitmapButton(self.scroll, id=-1, bitmap=wx.Image("images/osx/HelpBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(self.ServerBtn.GetPosition()[0]+25, 10), size=(25, 25))
	else:
	    self.HelpBtn = wx.Button(self.scroll, id=-1, label="?", pos=(self.ServerBtn.GetPosition()[0]+25, 10), size=(25, 25))
	    self.HelpBtn.SetForegroundColour("#0000FF")
	    self.HelpBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.HelpBtn.Bind(wx.EVT_BUTTON, self.showHelp)
	self.HelpBtn.SetToolTipString("Display the help file for this window")
	
	self.SeqViewer = wx.grid.Grid(self.scroll)
	self.SeqViewer.CreateGrid(0, 0)
	self.SeqViewer.SetSize((W-480, 150))
	self.SeqViewer.SetPosition((10, 50))
	self.SeqViewer.SetLabelFont(wx.Font(8, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.SeqViewer.DisableDragColSize()
	self.SeqViewer.DisableDragRowSize()
	self.rowpos = -1
	self.colpos = -1
	self.selectLookup = {}
	#self.sizer.Add(self.SeqViewer, (1, 0), span=(1, 5), flag = wx.ALIGN_CENTER | wx.EXPAND, border=5)
	self.cannotDelete = False
	
	xpos = self.SeqViewer.GetPosition()[0] + self.SeqViewer.GetSize()[0] + 5
	self.BGRecolorBtn = wx.BitmapButton(self.scroll, -1, wx.Image("images/colorwheel.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), (xpos, 50), (70, 35))
	self.BGRecolorBtn.Bind(wx.EVT_BUTTON, self.bgRecolor)
	self.BGRecolorBtn.SetToolTipString("Change PyMOL background color")
	if (platform.system() == "Darwin"):
	    self.StereoBtn = wx.BitmapButton(self.scroll, -1, wx.Image("images/osx/StereoBtn_Mono.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), (xpos, 90), (70, 35))
	else:
	    self.StereoBtn = wx.Button(self.scroll, id=-1, label="Mono", pos=(xpos, 90), size=(70, 35))
	    self.StereoBtn.SetForegroundColour("#000000")
	    self.StereoBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.StereoBtn.Bind(wx.EVT_BUTTON, self.stereoToggle)
	self.StereoBtn.SetToolTipString("Toggle stereo view on/off")
	if (platform.system() == "Darwin"):
	    self.ColoringBtn = wx.BitmapButton(self.scroll, -1, wx.Image("images/osx/ColoringBtn_NoColor.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), (xpos, 130), (70, 35))
	else:
	    self.ColoringBtn = wx.Button(self.scroll, id=-1, label="No Color", pos=(xpos, 130), size=(70, 35))
	    self.ColoringBtn.SetForegroundColour("#000000")
	    self.ColoringBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.ColoringBtn.Bind(wx.EVT_BUTTON, self.recolorClick)
	self.ColoringBtn.SetToolTipString("Color primary sequence by secondary structure/B-factor or turn off coloring")
	if (platform.system() == "Darwin"):
	    self.AlignBtn = wx.BitmapButton(self.scroll, -1, wx.Image("images/osx/AlignBtn_From1.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), (xpos, 170), (70, 35))
	else:
	    self.AlignBtn = wx.Button(self.scroll, id=-1, label="From 1", pos=(xpos, 170), size=(70, 35))
	    self.AlignBtn.SetForegroundColour("#000000")
	    self.AlignBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.AlignBtn.Bind(wx.EVT_BUTTON, self.realignToggle)
	self.AlignBtn.SetToolTipString("PDB sequences renumbered from 1")
	self.noMUSCLEWarned = False
	self.viewMode = "Normal"
	self.colorMode = "No Color"
	self.alignType = "From 1"
	
	self.saveTimer = wx.Timer(self)
	self.PyMOLUpdateTimer = wx.Timer(self)
	self.DownloadTimer = wx.Timer(self)
	self.activeJobs = []
	self.Bind(wx.EVT_TIMER, self.saveWindowData, self.saveTimer)
	self.Bind(wx.EVT_TIMER, self.updatePyMOLSelection, self.PyMOLUpdateTimer)
	self.Bind(wx.EVT_TIMER, self.downloader, self.DownloadTimer)
	self.DownloadTimer.Start(60000)
	self.Bind(wx.EVT_ACTIVATE, self.focusEvent)
	self.Bind(wx.grid.EVT_GRID_LABEL_LEFT_CLICK, self.labelClick)
	self.SeqViewer.GetGridWindow().Bind(wx.EVT_LEFT_UP, self.leftRelease)
	self.SeqViewer.GetGridWindow().Bind(wx.EVT_MOTION, self.onMouseOver)
	#if (platform.system() == "Windows"):
	self.Bind(wx.grid.EVT_GRID_CELL_LEFT_CLICK, self.leftClick)
	self.Bind(wx.grid.EVT_GRID_CELL_RIGHT_CLICK, self.rightClick)
	self.Bind(wx.grid.EVT_GRID_LABEL_LEFT_DCLICK, self.labelDClick)
	self.SeqViewer.GetGridWindow().Bind(wx.EVT_KEY_DOWN, self.keyPress)
	self.Bind(wx.EVT_SIZE, self.windowGeometryChange)
	self.Bind(wx.EVT_MOTION, self.windowGeometryChange)
	self.selectedResidues = []
	
	self.labelMsg = wx.StaticText(self.scroll, -1, "", (10, 205), (winw-50, 25))
	self.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
	self.labelMsg.SetForegroundColour("#FFFFFF")
	#self.sizer.Add(self.labelMsg, (2, 0), span=(1, 4), flag = wx.ALIGN_CENTER | wx.EXPAND, border=5)
	self.msgQueue = []
	
	self.scroll.SetScrollbars(1, 1, max(self.StereoBtn.GetPosition()[0] + self.StereoBtn.GetSize()[0] + 5, self.HelpBtn.GetPosition()[0] + self.HelpBtn.GetSize()[0] + 5), winh-40)
	
	#self.sizer.AddGrowableRow(1)
	#self.sizer.AddGrowableCol(5)
	#self.SequencePanel.SetSizerAndFit(self.sizer)
	self.Show()
	
	self.pdbreader = Bio.PDB.PDBParser()
	self.pdbwriter = Bio.PDB.PDBIO()
	
	self.sscolormap = {}
	self.sscolormap["H"] = ("red", "white")
	self.sscolormap["S"] = ("yellow", "black")
	self.sscolormap["L"] = ("white", "black")
	self.sscolormap[""] = ("white", "black")
	self.sscolormap["T"] = ("blue", "white")
	self.sscolormap["B"] = ("blue", "white")
	self.sscolormap["G"] = ("orange", "black")
	self.sscolormap["I"] = ("orange", "black")
	
	# Try to use DSSP if it is available, to give more accurate secondary structure coloring
	try:
	    if (platform.system() == "Windows"):
		self.dsspexe = glob.glob(self.scriptdir + "\\bin\\dssp_win*")[0]
	    else:
		self.dsspexe = glob.glob(self.scriptdir + "/bin/dssp_unix*")[0]
	except:
	    self.dsspexe = "N/A"
	    print "Note: DSSP could be used to improve secondary structure predictions."
	    if (platform.system() == "Windows"):
		print "      Install it to " + self.scriptdir + "\\bin\\dssp_win to make it available."
	    else:
		print "      Install it to " + self.scriptdir + "/bin/dssp_unix to make it available."
	
    def setProtWin(self, protWin):
	self.protWin = protWin
	
    def setProtocolPanel(self, protPanel):
	self.protPanel = protPanel
    
    def showHelp(self, event):
	# Open the help page
	webbrowser.open(self.scriptdir + "/help/sequence.html")
	
    def configureServer(self, event):
	# This button allows the user to give a name for the remote server
	# After submission, a test will be made and the user will be notified if the test was sucessful
	dlg = wx.TextEntryDialog(self, "Enter the URL of the remote server:", "Configure Remote Server", "", style=wx.OK | wx.CANCEL)
        dlg.SetValue(getServerName())
        if (dlg.ShowModal() == wx.ID_OK):
	    setServerName(dlg.GetValue().strip())
	    # Send test input
	    try:
		sendToServer("testinput")
		dlg2 = wx.MessageDialog(self, "Server test succeeded!", "Server Successful", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
		dlg2.ShowModal()
		dlg2.Destroy()
		self.saveWindowData(event)
	    except:
		dlg2 = wx.MessageDialog(self, "Server test failed!  Either the server is not set up to run Rosetta, the server is behind a firewall, or you do not have a network connection.", "Server Failed", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
		dlg2.ShowModal()
		dlg2.Destroy()
        dlg.Destroy()
    
    def windowGeometryChange(self, event):
	# This function starts a timer that will write out the size and position of this window to a cfg file
	# so the orientation is saved and can be loaded the next time InteractiveROSETTA is started
	if (not(self.saveTimer.IsRunning())):
	    self.saveTimer.Start(5000)
	event.Skip()

    def realign(self):
	if (self.alignType == "Align"):
	    # Use MUSCLE to perform the multiple sequence alignment using BioPython's tools
	    # First, does MUSCLE exist?  If not, tell the user this feature is disabled until it is download
	    try:
		if (platform.system() == "Windows"):
		    muscle = glob.glob(self.scriptdir + "\\bin\\muscle_win*")[0]
		else:
		    muscle = glob.glob(self.scriptdir + "/bin/muscle_unix*")[0]
	    except:
		if (not(self.noMUSCLEWarned)):
		    if (platform.system() == "Windows"):
			msg = "MUSCLE installation not found.  If you want to use the alignment feature, then download a copy of MUSCLE to " + self.scriptdir + "\\bin\\muscle_win."
		    else:
			msg = "MUSCLE installation not found.  If you want to use the alignment feature, then download a copy of MUSCLE to " + self.scriptdir + "/bin/muscle_unix."
		    dlg = wx.MessageDialog(self, msg, "MUSCLE Not Found", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
		    dlg.ShowModal()
		    dlg.Destroy()
		    self.noMUSCLEWarned = True
		self.realignToggle(None) # To get back to normal #
		return
	    # Get all the models
	    models = []
	    for r in range(0, len(self.IDs)):
		if (not(self.IDs[r][0:len(self.IDs[r])-2] in models)):
		    models.append(self.IDs[r][0:len(self.IDs[r])-2])
	    # Get all the sequences
	    nseqs = 0
	    if (len(self.IDs) == 0):
		self.realignToggle(None) # To get back to normal #
		return
	    else:
		lastmodel = self.IDs[0][0:len(self.IDs[0])-2]
	    seq = ""
	    seqs = []
	    # Generate a FASTA file with each chain as a different entry
	    # They need to be marked with an index, because MUSCLE can shuffle their order and the order
	    # needs to be reconstructed so the sequences match up in the window
	    # This is not a totally ideal MSA because the chains are all treated separately
	    # I may add a chain joining button at some point to make this easier so I don't have to mess with the MSA
	    f = open("alignin.fasta", "w")
	    indx = 0
	    for r in range(0, len(self.IDs)):
		model =  self.IDs[r][0:len(self.IDs[r])-2]
		if (model != lastmodel):
		    nseqs = nseqs + 1
		    lastmodel = model
		f.write("> " + str(indx) + "\n")
		f.write(self.sequences[r].replace("-", "").strip() + "\n")
		indx = indx + 1
	    f.close()
	    nseqs = nseqs + 1
	    if (nseqs <= 1):
		self.realignToggle(None) # To get back to normal #
		return
	    # Do alignment with MUSCLE
	    muscle_cline = MuscleCommandline(muscle, input="alignin.fasta", out="alignout.fasta")
	    muscle_cline()
	    alignment = range(0, indx)
	    # Read the alignment into a list
	    f = open("alignout.fasta", "r")
	    seq = ""
	    first = True
	    for aline in f:
		if (">" in aline and first):
		    first = False
		    indx = int(aline[aline.index(">")+1:])
		elif (">" in aline):
		    alignment[indx] = seq
		    indx = int(aline[aline.index(">")+1:])
		    seq = ""
		else:
		    seq = seq + aline.strip()
	    alignment[indx] = seq
	    f.close()
	    os.remove("alignin.fasta")
	    os.remove("alignout.fasta")
	    # Now update the data structures
	    thisseq = alignment.pop(0)
	    offset = 0
	    for r in range(0, len(self.sequences)):
		# Find the length of the aligned sequence that corresponds to this chain
		lenchain = 0
		for char in self.sequences[r]:
		    if (char != "-" and char != " "):
			lenchain = lenchain + 1
		thislen = 0
		for i in range(0, len(thisseq)):
		    if (thisseq[i] != "-" and thisseq[i] != " "):
			thislen = thislen + 1
		    if (thislen == lenchain):
			break
		self.sequences[r] = ""
		for j in range(0, offset):
		    self.sequences[r] = self.sequences[r] + "-"
		self.sequences[r] = self.sequences[r] + thisseq[0:i+1]
		charsleft = 0
		for char in thisseq[i+1:].strip():
		    if (char != "-" and len(char.strip()) != 0):
			charsleft = charsleft + 1
		if (charsleft == 0):
		    try:
			thisseq = alignment.pop(0)
			offset = 0
		    except:
			# Last one, the loop will exit after this iteration anyway
			pass
		else:
		    offset = offset + len(self.sequences[r])
		    thisseq = thisseq[i+1:]
	    # Now we have to fix indxToSeqPos
	    # It should already be in PDB order since it was just in PDB # mode, so let's take all the
	    # dashes out, then iterate down each of the sequences and insert dashes when appropriate
	    for r in range(0, len(self.indxToSeqPos)):
		for c in range(len(self.indxToSeqPos[r])-1, -1, -1):
		    if (self.indxToSeqPos[r][c] == "-"):
			self.indxToSeqPos[r].pop(c)
	    for r in range(0, len(self.sequences)):
		for c in range(0, len(self.sequences[r])):
		    if (self.sequences[r][c] == "-"):
			self.indxToSeqPos[r].insert(c, "-")
	    # MUSCLE likes to rename all the HETATMs as "X", so change them back to Os and Zs
	    for r in range(0, len(self.sequences)):
		poseindx = self.getPoseIndex(r)
		for c in range(0, len(self.sequences[r])):
		    if (self.sequences[r][c] == "-"):
			continue
		    resID = self.indxToSeqPos[r][c]
		    chain = self.IDs[r][len(self.IDs[r])-1]
		    if (chain == "_"):
			chain = " "
		    if (self.poses[poseindx][0][chain][resID].resname == "HOH"):
			self.sequences[r] = self.sequences[r][0:c] + "O" + self.sequences[r][c+1:]
		    elif (self.sequences[r][c] == "X"):
			self.sequences[r] = self.sequences[r][0:c] + "Z" + self.sequences[r][c+1:]
	    # Now update the SequenceWindow
	    for r in range(0, len(self.sequences)):
		if (r < len(self.sequences) - 1):
		    self.updateSeqViewer(r, relabel=False)
		else:
		    self.updateSeqViewer(r, relabel=True)
	elif (self.alignType == "PDB #"):
	    # This shouldn't be too hard at first: all you have to do is find the max PDB number for each chain
	    # then create a sequence of the same number of - and then iterate down the residues and replace each
	    # position with the appropriate letter
	    ncols = 0
	    for r in range(0, len(self.sequences)):
		poseindx = self.getPoseIndex(r)
		chain = self.IDs[r][len(self.IDs[r])-1]
		if (chain == "_"):
		    chain = " "
		nres = 0
		for residue in self.poses[poseindx][0][chain]:
		    nres = max(nres, residue.id[1])
		self.sequences[r] = ""
		self.indxToSeqPos[r] = []
		for i in range(0, nres):
		    self.sequences[r] = self.sequences[r] + "-"
		    self.indxToSeqPos[r].append("-")
		for residue in self.poses[poseindx][0][chain]:
		    ires = residue.id[1] - 1
		    self.sequences[r] = self.sequences[r][0:ires] + AA3to1(residue.resname) + self.sequences[r][ires+1:]
		    self.indxToSeqPos[r][ires] = residue.id
		ncols = max(ncols, len(self.sequences[r]))
	    # Here's where things get dicey...sometimes the PDB numbering is such that there are very large ranges
	    # of empty space (like HOHs being at position 5,000+, which generates a ton of columns)
	    # So now go through the sequences and remove columns that are only filled with -
	    # But then we're going to have to keep track of these deletions so the column labels can reflect what
	    # is actually at the position whenever there is a break like this
	    self.columnlabels = []
	    columns_to_delete = []
	    for c in range(0, ncols):
		columnbad = True
		for r in range(0, len(self.sequences)):
		    if (c < len(self.sequences[r])):
			if (self.sequences[r][c] != "-"):
			    columnbad = False
			    break
		    if (not(columnbad)):
			break
		if (columnbad):
		    columns_to_delete.append(c)
		else:
		    self.columnlabels.append(c+1)
	    # Now take out these columns
	    for i in range(len(columns_to_delete)-1, -1, -1):
		col = columns_to_delete[i]
		for r in range(0, len(self.sequences)):
		    if (col < len(self.sequences[r])):
			self.indxToSeqPos[r].pop(col)
			self.sequences[r] = self.sequences[r][0:col] + self.sequences[r][col+1:].strip()
	    # Now update the SequenceWindow
	    for r in range(0, len(self.sequences)):
		if (r < len(self.sequences) - 1):
		    self.updateSeqViewer(r, relabel=False)
		else:
		    self.updateSeqViewer(r, relabel=True)
	else:
	    # Simple, just make everything numbered from one again 
	    r = 0
	    for pose in self.poses:
		if (pose):
		    for ch in pose[0]:
			self.indxToSeqPos[r] = []
			self.sequences[r] = ""
			for residue in ch:
			    self.indxToSeqPos[r].append(residue.id)
			    self.sequences[r] = self.sequences[r] + AA3to1(residue.resname)
			r = r + 1
	    # Now update the SequenceWindow
	    for r in range(0, len(self.sequences)):
		if (r < len(self.sequences) - 1):
		    self.updateSeqViewer(r, relabel=False)
		else:
		    self.updateSeqViewer(r, relabel=True)

    def realignToggle(self, event):
	if (self.AlignBtn.GetLabel() == "From 1"):
	    self.alignType = "PDB #"
	    if (platform.system() == "Darwin"):
		self.AlignBtn.SetBitmapLabel(bitmap=wx.Image("images/osx/AlignBtn_PDB.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
	    else:
		self.AlignBtn.SetLabel(self.alignType)
	    self.regenerateLookupTable()
	#elif (self.AlignBtn.GetLabel() == "PDB #"):
	#    self.AlignBtn.SetLabel("Align")
	#    self.regenerateLookupTable()
	else:
	    if (platform.system() == "Darwin"):
		self.AlignBtn.SetBitmapLabel(bitmap=wx.Image("images/osx/AlignBtn_From1.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
	    else:
		self.AlignBtn.SetLabel(self.alignType)
	    self.AlignBtn.SetLabel(self.alignType)
	    self.regenerateLookupTable()

    def recoverFromError(self, protocol):
	# This function tells the user what the error was and tries to revert the protocol
	# back to the pre-daemon state so the main GUI can continue to be used
	f = open("errreport", "r")
	errmsg = "An error was encountered during a submitted " + protocol + " protocol:\n\n"
	for aline in f:
	    errmsg = errmsg + aline.strip()
	f.close()
	errmsg = str(errmsg)
	logInfo("Error encountered")
	logInfo(errmsg)
	if (platform.system() == "Windows"):
	    sessioninfo = os.path.expanduser("~") + "\\InteractiveRosetta\\sessionlog"
	else:
	    sessioninfo = os.path.expanduser("~") + "/InteractiveRosetta/sessionlog"
	errmsg = errmsg + "\n\nIf you don't know what caused this, send the file " + sessioninfo + " to a developer along with an explanation of what you did."
	# You have to use a MessageDialog because the MessageBox doesn't always work for some reason
	dlg = wx.MessageDialog(self, errmsg, "Error Encountered", wx.OK|wx.ICON_EXCLAMATION)
	dlg.ShowModal()
	dlg.Destroy()
	os.remove("errreport")

    def downloader(self, event):
	self.DownloadTimer.Stop()
	# You have to do this on Linux due to the libc bug
	if (platform.system() == "Linux"):
	    res_init()
	# Read the jobs that are active
	self.activeJobs = []
	removeJob = []
	changed = False
	goToSandbox()
	try:
	    f = open("downloadwatch", "r")
	except:
	    return
	for aline in f:
	    self.activeJobs.append(aline.strip())
	    removeJob.append(False)
	f.close()
	for i in range(0, len(self.activeJobs)):
	    job = self.activeJobs[i]
	    if (job[0:8] == "FRAGMENT"):
		ID = job.split("\t")[1].strip()
		URL = "http://www.robetta.org/downloads/fragments/" + ID
		downloadfiles = []
		try:
		    downloadpage = urllib2.urlopen(URL, timeout=1)
		    for aline in downloadpage:
			# Look for the links
			indx = aline.find("<a href")
			if (indx < 0):
			    continue
			linkbegin = aline.find(">", indx+1) + 1
			linkend = aline.find("</a>", indx+1)
			if (not("Parent" in aline[linkbegin:linkend]) and not("Name" in aline[linkbegin:linkend])):
			    if (aline[linkbegin:linkend].endswith(".200_v1_3")):
				downloadfiles.append(aline[linkbegin:linkend])
			    elif (aline[linkbegin:linkend].endswith(".fasta")):
				downloadfiles.append(aline[linkbegin:linkend])
			    elif (aline[linkbegin:linkend].endswith(".psipred")):
				downloadfiles.append(aline[linkbegin:linkend])
			    elif (aline[linkbegin:linkend].endswith(".psipred_ss2")):
				downloadfiles.append(aline[linkbegin:linkend])
		    downloadpage.close()
		    dlg = wx.MessageDialog(self, "Your fragments package for job ID " + ID + " is ready.", "Fragments Download Ready", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
		    dlg.ShowModal()
		    dlg.Destroy()
		    while (True):
			dlg = wx.FileDialog(
			    self, message="Save the Frag File",
			    defaultDir=self.cwd,
			    defaultFile=ID,
			    wildcard="Frag Files (*.frag)|*.frag",
			    style=wx.SAVE | wx.CHANGE_DIR)
			if (dlg.ShowModal() == wx.ID_OK):
			    paths = dlg.GetPaths()
			    # Change cwd to the last opened file
			    if (platform.system() == "Windows"):
				lastDirIndx = paths[len(paths)-1].rfind("\\")
			    else:
				lastDirIndx = paths[len(paths)-1].rfind("/")
			    self.cwd = str(paths[len(paths)-1][0:lastDirIndx])
			    self.saveWindowData(None)
			    # Load the PDBs into PyMOL
			    filename = str(paths[0]).split(".frag")[0] + ".frag"
			    # Does it exist already?  If so, ask if the user really wants to overwrite it
			    if (os.path.isfile(filename)):
				dlg2 = wx.MessageDialog(self, "The file " + filename + " already exists.  Overwrite it?", "Filename Already Exists", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
				if (dlg2.ShowModal() == wx.ID_NO):
				    dlg2.Destroy()
				    logInfo("Cancelled save operation due to filename already existing")
				    continue
				dlg2.Destroy() 
			    break
		    busyDlg = wx.BusyInfo("Downloading fragment file, please wait...")
		    gzipfile = gzip.open(str(filename), "wb")
		    for downloadfile in downloadfiles:
			serverdata = urllib2.urlopen(URL + "/" + downloadfile)
			gzipfile.write("BEGIN\t" + downloadfile + "\n")
			gzipfile.writelines(serverdata)
			gzipfile.write("END\t" + downloadfile + "\n")
			serverdata.close()
		    gzipfile.close()
		    busyDlg = None
		    changed = True
		    removeJob[i] = True
		except:
		    # Not there yet
		    pass
	    elif (job[0:3] == "MSD"):
		ID = job.split("\t")[1].strip()
		URL = "http://bach1.bio.rpi.edu/schenc3/InteractiveROSETTA/results/" + ID + "/results.gz"
		try:
		    downloadpage = urllib2.urlopen(URL, timeout=1) # To make sure its there before display the dialog
		    downloadpage.close()
		    dlg = wx.MessageDialog(self, "Your MSD package job ID " + ID + " is ready.", "MSD Download Ready", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
		    dlg.ShowModal()
		    dlg.Destroy()
		    busyDlg = wx.BusyInfo("Downloading MSD files, please wait...")
		    (oldfilename, info) = urllib.urlretrieve(URL, "results.gz")
		    busyDlg.Destroy()
		    del busyDlg
		    while (True):
			dlg = wx.FileDialog(
			    self, message="Save the MSD package",
			    defaultDir=self.cwd,
			    defaultFile=ID,
			    wildcard="Compressed Files (*.gz)|*.gz",
			    style=wx.SAVE | wx.CHANGE_DIR)
			if (dlg.ShowModal() == wx.ID_OK):
			    paths = dlg.GetPaths()
			    # Change cwd to the last opened file
			    if (platform.system() == "Windows"):
				lastDirIndx = paths[len(paths)-1].rfind("\\")
			    else:
				lastDirIndx = paths[len(paths)-1].rfind("/")
			    self.cwd = str(paths[len(paths)-1][0:lastDirIndx])
			    self.saveWindowData(None)
			    filename = str(paths[0]).split(".gz")[0] + ".gz"
			    # Does it exist already?  If so, ask if the user really wants to overwrite it
			    if (os.path.isfile(filename)):
				dlg2 = wx.MessageDialog(self, "The file " + filename + " already exists.  Overwrite it?", "Filename Already Exists", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
				if (dlg2.ShowModal() == wx.ID_NO):
				    dlg2.Destroy()
				    logInfo("Cancelled save operation due to filename already existing")
				    continue
				else:
				    os.remove(str(filename))
				dlg2.Destroy() 
			    break
		    goToSandbox()
		    os.rename("results.gz", str(filename))
		    # Unpackage the files
		    gzipfile = gzip.open(str(filename), "rb")
		    readingData = False
		    for aline in gzipfile:
			if (aline.startswith("BEGIN PDB")):
			    if (platform.system() == "Windows"):
				f = open(self.cwd + "\\" + aline.split()[len(aline.split())-1].strip(), "w")
			    else:
				f = open(self.cwd + "/" + aline.split()[len(aline.split())-1].strip(), "w")
			    readingData = True
			elif (aline.startswith("END PDB")):
			    f.close()
			    readingData = False
			elif (readingData):
			    f.write(aline.strip() + "\n")
		    gzipfile.close()
		    changed = True
		    removeJob[i] = True
		except:
		    # Not there yet
		    pass
		URL = "http://bach1.bio.rpi.edu/schenc3/InteractiveROSETTA/results/" + ID + "/errreport"
		try:
		    # Look for an error file
		    downloadpage = urllib2.urlopen(URL, timeout=1) # To make sure its there before display the dialog
		    f = open("errreport", "w")
		    for aline in downloadpage:
			f.write(aline.strip() + "\n")
		    downloadpage.close()
		    f.close()
		    self.recoverFromError("MSD")
		    changed = True
		    removeJob[i] = True
		except:
		    pass
	# If a change happened, update the "downloadwatch" file to take the completed jobs out
	if (changed):
	    goToSandbox()
	    try:
		os.remove("downloadwatch")
	    except:
		pass
	    f = open("downloadwatch", "w")
	    for i in range(0, len(self.activeJobs)):
		if (not(removeJob[i])):
		    f.write(self.activeJobs[i].strip() + "\n")
	    f.close()
	self.DownloadTimer.Start(60000)

    def onMouseOver(self, event):
        # Use CalcUnscrolledPosition() to get the mouse position
        # within the
        # entire grid including what's offscreen
        x, y = self.SeqViewer.CalcUnscrolledPosition(event.GetX(),event.GetY())
        coords = self.SeqViewer.XYToCell(x, y)
        # you only need these if you need the value in the cell
        row = coords[0]
        col = coords[1]
        if (row == self.rowpos and col == self.colpos):
	    event.Skip()
	    return
	else:
	    self.rowpos = row
	    self.colpos = col
	if (row < 0 or col < 0):
	    event.GetEventObject().SetToolTipString("")
	    event.Skip()
	    return
	try:
	    self.SeqViewer.GetCellValue(row, col)
	    resID = self.indxToSeqPos[row][col]
	    chain = self.IDs[row][len(self.IDs[row])-1]
	    if (chain == "_"):
		chain = " "
	    AA3 = self.poses[self.getPoseIndex(row)][0][chain][resID].resname.strip()
	    event.GetEventObject().SetToolTipString(AA3 + str(self.indxToSeqPos[row][col][1]))
	except:
	    event.GetEventObject().SetToolTipString("")
	event.Skip()

    def saveWindowData(self, event):
	self.saveTimer.Stop()
	homedir = os.path.expanduser("~")
	data = []
	try:
	    if (platform.system() == "Windows"):
		f = open(homedir + "\\InteractiveROSETTA\\seqwindow.cfg", "r")
	    else:
		f = open(homedir + "/InteractiveROSETTA/seqwindow.cfg", "r")
	    for aline in f:
		data.append(aline)
	    f.close()
	except:
	    pass
	if (platform.system() == "Windows"):
	    f = open(homedir + "\\InteractiveROSETTA\\seqwindow.cfg", "w")
	else:
	    f = open(homedir + "/InteractiveROSETTA/seqwindow.cfg", "w")
	itemsFound = [False, False, False, False, False, False, False, False] # [offX, offY, offW, offH, offpw, offph cwd, serverName]
	(x, y) = self.GetPosition()
	(w, h) = self.GetSize()
	(pw, ph) = self.cmd.get_viewport()
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
	    elif (aline.find("[OFFSET PWIDTH]") >= 0):
		itemsFound[4] = True
		f.write("[OFFSET PWIDTH] " + str(pw-self.stdwinw) + "\n")
	    elif (aline.find("[OFFSET PHEIGHT]") >= 0):
		itemsFound[5] = True
		f.write("[OFFSET PHEIGHT] " + str(ph-(self.screenH-340)) + "\n")
	    elif (aline.find("[CWD]") >= 0):
		itemsFound[6] = True
		f.write("[CWD]\t" + str(self.cwd) + "\n")
	    elif (aline.find("[SERVER]") >= 0):
		itemsFound[7] = True
		f.write("[SERVER]\t" + getServerName() + "\n")
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
		elif (i == 4):
		    f.write("[OFFSET PWIDTH] " + str(pw-self.stdwinw) + "\n")
		elif (i == 5):
		    f.write("[OFFSET PHEIGHT] " + str(ph-self.stdwinh) + "\n")
		elif (i == 6):
		    f.write("[CWD]\t" + str(self.cwd) + "\n")
		else:
		    f.write("[SERVER]\t" + getServerName() + "\n")
	f.close()
	
    def setPyMOL(self, pymol):
	self.pymol = pymol
	self.cmd = pymol.cmd
	self.stored = pymol.stored
    
    def focusEvent(self, event):
	self.inFocus = event.GetActive()
	if (self.inFocus):
	    try:
		self.cmd.enable("sele")
		self.selectUpdate(updatePyMOL=False)
	    except:
		# Immediately after the window is created this function gets called, yet self.cmd is still not defined
		pass
	else:
	    # When they are in PyMOL, everything has to belong to the "sele" selection
	    return
	    try:
		self.cmd.select("sele", "seqsele or sele")
		self.cmd.delete("seqsele")
		self.cmd.enable("sele")
	    except:
		try:
		    self.cmd.select("sele", "seqsele")
		    self.cmd.delete("seqsele")
		    self.cmd.enable("sele")
		except:
		    self.cmd.enable("sele")
	    return
	    # We're leaving the sequence window so update PyMOL with the current selections in the SeqViewer
	    self.inFocus = True
	    self.selectUpdate(self.inFocus)
	    self.inFocus = False
	    #self.selectTimer.Stop()

    def renumber(self, event):
	# Get the selected elements and make sure there is only one
	topLefts = self.SeqViewer.GetSelectionBlockTopLeft()
	bottomRights = self.SeqViewer.GetSelectionBlockBottomRight()
	cells = []
	for i in range(0, len(topLefts)):
	    for r in range(topLefts[i][0], bottomRights[i][0]+1):
		for c in range(topLefts[i][1], bottomRights[i][1]+1):
		    cells.append((r, c))
	if (len(cells) == 0):
	    return
	elif (len(cells) > 1):
	    dlg = wx.MessageDialog(self, "Please select only one residue to indicate the chain's new starting residue.", "Select Only One Residue", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
	    dlg.ShowModal()
	    dlg.Destroy()
	    return
	(r, c) = cells[0]
	poseindx = self.getPoseIndex(r)
	model = self.getModelForChain(r)
	ires = 1
	chain = self.IDs[r][len(self.IDs[r])-1]
	if (chain == "_"):
	    chain = " "
	leadingres = []
	# Take out the residues leading up to the new start
	for residue in self.poses[poseindx][0][chain]:
	    if (residue.id[1] == self.indxToSeqPos[r][c][1]):
		break
	    leadingres.append(residue)
	for residue in leadingres:
	    self.poses[poseindx][0][chain].detach_child(residue.id)
	# Add back the original leading residues at the end
	for residue in leadingres:
	    #thisID = residue.id
	    #residue.id = (thisID[0], ires, thisID[2])
	    #ires = ires + 1
	    self.poses[poseindx][0][chain].add(residue)
	# Renumber everything still in the chain
	for residue in self.poses[poseindx][0][chain]:
	    thisID = residue.id
	    residue.id = (thisID[0], ires, thisID[2])
	    ires = ires + 1
	# Change back to numbering from one so we can easily get back the proper sequences and indexes
	if (self.AlignBtn.GetLabel() != "From 1"):
	    self.realignToggle(None)
	else:
	    # Regenerating the lookup table calls the realigner which will regenerate the sequences properly
	    self.regenerateLookupTable()
	# Now we have to reload the structure in PyMOL
	self.pdbwriter.set_structure(self.poses[poseindx])
	self.pdbwriter.save(model + ".pdb")
	self.poses[poseindx] = self.pdbreader.get_structure(model, model + ".pdb")
	self.cmd.remove(model)
	self.cmd.delete(model)
	self.cmd.load(model + ".pdb", model)
	# Redo DSSP if available
	if (self.dsspexe != "N/A"):
	    try:
		dsspdata = Bio.PDB.DSSP(self.poses[poseindx][0], model + ".pdb", dssp=self.dsspexe)
	    except:
		self.dsspexe = "N/A"
		print "WARNING: DSSP could not be run.  Is its binary executable?"
	    for c in self.poses[poseindx][0]:
		for r in c:
		    if (self.dsspexe != "N/A"):
			# Alter the states in PyMOL with the DSSP secondary structure predictions
			try:
			    ss = dsspdata[(c.id, r.id[1])][1]
			    # Not all the DSSP codes match up with PyMOL, so make them match
			    if (ss == "E"):
				# Beta sheet
				ss = "S"
			    elif (ss == "S"):
				# Bend
				ss = "B"
			    elif (ss == "-"):
				# Loop
				ss = "L"
			    if (c.id != " "):
				self.cmd.alter("model " + model + " and chain " + c.id + " and resi " + str(r.id[1]), "ss=\"" + ss + "\"")
			    else:
				self.cmd.alter("model " + model + " and resi " + str(r.id[1]), "ss=\"" + ss + "\"")
			except:
			    ss = "S"
	defaultPyMOLView(self.cmd, model)
	# Clear the selection otherwise deleted rows are still considered selected
	self.SeqViewer.ClearSelection()

    def joinChains(self, event):
	# Do we have more than one chain selected, if not, do nothing
	selectedrows = self.getSelectedChains()
	if (len(selectedrows) <= 1):
	    return
	# Make sure the user really wants to do this
	dlg = wx.MessageDialog(self, "Are you sure you want to join these chains together?", "Confirm Chain Join", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
	if (dlg.ShowModal() == wx.ID_NO):
	    dlg.Destroy()
	    return
	dlg.Destroy()
	# Are multiple models selected?  If so, tell the user and ask if they really want to move a chain to a different model
	if (len(self.getSelectedModels()) > 1):
	    dlg = wx.MessageDialog(self, "You have multiple models selected.  The chain from the lower model will be moved to the higher model.  Is this okay?", "Confirm Chain Join", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
	    if (dlg.ShowModal() == wx.ID_NO):
		dlg.Destroy()
		return
	    dlg.Destroy()
	base_r = selectedrows.pop(0)
	basechain = self.IDs[base_r][len(self.IDs[base_r])-1]
	baseposeindx = self.getPoseIndex(base_r)
	basemodel = self.getModelForChain(base_r)
	# Renumber the base model from 1 to make sure there won't be any numbering conflicts
	ires = 1
	for residue in self.poses[baseposeindx][0][basechain]:
	    thisID = residue.id
	    residue.id = (thisID[0], ires, thisID[2])
	    ires = ires + 1
	nrowsdelete = len(selectedrows)
	for r in selectedrows:
	    chain = self.IDs[r][len(self.IDs[r])-1]
	    poseindx = self.getPoseIndex(r)
	    model = self.getModelForChain(r)
	    for residue in self.poses[poseindx][0][chain]:
		thisID = residue.id
		residue.id = (thisID[0], ires, thisID[2])
		ires = ires + 1
		residue.detach_parent()
		self.poses[baseposeindx][0][basechain].add(residue)
	    self.poses[poseindx][0].detach_child(chain)
	    self.sequences[r] = "RemoveMe"
	    if (r != len(self.sequences) - 1 and not(self.poses[r+1])):
		self.poses[r+1] = self.poses[r]
	    # Take it out of PyMOL
	    if (len(chain.strip()) > 0):
		self.cmd.remove("model " + model + " and chain " + chain)
	    else:
		self.cmd.remove("model " + model)
	# Delete the unoccupied rows
	for i in range(len(selectedrows)-1, -1, -1):
	    r = selectedrows[i]
	    self.SeqViewer.DeleteRows(pos=r)
	    self.sequences.pop(r)
	    self.indxToSeqPos.pop(r)
	    self.poses.pop(r)
	    self.IDs.pop(r)
	# Change back to numbering from one so we can easily get back the proper sequences and indexes
	if (self.AlignBtn.GetLabel() != "From 1"):
	    self.realignToggle(None)
	else:
	    # Regenerating the lookup table calls the realigner which will regenerate the sequences properly
	    self.regenerateLookupTable()
	# Now we have to reload the structure in PyMOL
	self.pdbwriter.set_structure(self.poses[baseposeindx])
	self.pdbwriter.save(basemodel + ".pdb")
	self.poses[baseposeindx] = self.pdbreader.get_structure(basemodel, basemodel + ".pdb")
	self.cmd.remove(basemodel)
	self.cmd.delete(basemodel)
	self.cmd.load(basemodel + ".pdb", basemodel)
	# Redo DSSP if available
	if (self.dsspexe != "N/A"):
	    try:
		dsspdata = Bio.PDB.DSSP(self.poses[baseposeindx][0], basemodel + ".pdb", dssp=self.dsspexe)
	    except:
		self.dsspexe = "N/A"
		print "WARNING: DSSP could not be run.  Is its binary executable?"
	    for c in self.poses[baseposeindx][0]:
		for r in c:
		    if (self.dsspexe != "N/A"):
			# Alter the states in PyMOL with the DSSP secondary structure predictions
			try:
			    ss = dsspdata[(c.id, r.id[1])][1]
			    # Not all the DSSP codes match up with PyMOL, so make them match
			    if (ss == "E"):
				# Beta sheet
				ss = "S"
			    elif (ss == "S"):
				# Bend
				ss = "B"
			    elif (ss == "-"):
				# Loop
				ss = "L"
			    if (c.id != " "):
				self.cmd.alter("model " + basemodel + " and chain " + c.id + " and resi " + str(r.id[1]), "ss=\"" + ss + "\"")
			    else:
				self.cmd.alter("model " + basemodel + " and resi " + str(r.id[1]), "ss=\"" + ss + "\"")
			except:
			    ss = "S"
	defaultPyMOLView(self.cmd, basemodel)
	# Clear the selection otherwise deleted rows are still considered selected
	self.SeqViewer.ClearSelection()

    def doesResidueExist(self, model, chain, seqpos):
	# Useful function for determining whether a residue exists on a model
	# You can easily figure this out by trying to access it in the lookup table
	# Return its coordinates so resfile data can be updated easily
	try:
	    return self.selectLookup[(model + "|" + chain, int(seqpos))]
	except:
	    return False

    def regenerateLookupTable(self):
	# Whenever residues/chains are deleted, the selectLookup table gets all messed up because the rows and columns
	# don't refer to the same things anymore, so we have to recalculate it after a deletion
	# First do a realignment if necessary to get everything in the right places in the grid
	self.realign()
	self.selectLookup = {}
	for r in range(0, self.SeqViewer.NumberRows):
	    for c in range(0, len(self.indxToSeqPos[r])):
		# Don't add dashes
		if (self.indxToSeqPos[r][c] != "-"):
		    self.selectLookup[(self.IDs[r], int(self.indxToSeqPos[r][c][1]))] = (r, c)
	# Recolor the residues
	self.recolorResidues()

    def keyPress(self, event):
	if (int(event.GetKeyCode()) == wx.WXK_DELETE):
	    if (self.cannotDelete):
		wx.MessageBox("You cannot perform deletions while a protocol is active!", "Cannot Delete During Protocol", wx.OK|wx.ICON_EXCLAMATION)
		return
	    logInfo("Pressed the DELETE key")
	    msg = "Are you sure you want to delete the selected residues?"
	    dlg = wx.MessageDialog(self, msg, "Delete Residues", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
	    if (dlg.ShowModal() == wx.ID_YES):
		# Get the selected residue ranges
		topLefts = self.SeqViewer.GetSelectionBlockTopLeft()
		bottomRights = self.SeqViewer.GetSelectionBlockBottomRight()
		goToSandbox()
		# Get the (r, c) cells to be deleted
		# This list needs to be sorted such that the columns are always decreasing, otherwise the deletions
		# will screw up what the c value refers to since it will alter the later c values
		# if you delete entries in front of them (i.e. if you delete position 6 before 7,
		# then position 7 becomes 6 and old position 8 will be deleted by accident)
		cells = []
		for i in range(0, len(topLefts)):
		    for r in range(topLefts[i][0], bottomRights[i][0]+1):
			for c in range(topLefts[i][1], bottomRights[i][1]+1):
			    cells.append((r, c))
		for i in range(0, len(cells)-1):
		    highest_c = i
		    for j in range(i+1, len(cells)):
			if (cells[j][1] > cells[highest_c][1]):
			    highest_c = j
			elif (cells[j][1] == cells[highest_c][1]):
			    if (cells[j][0] > cells[highest_c][0]):
				highest_c = j
		    temp = cells[i]
		    cells[i] = cells[highest_c]
		    cells[highest_c] = temp
		# First look to see if whole chains have been selected, because then we should use "deleteChain" on these to make it faster
		chain_deletes = []
		for row in range(0, self.SeqViewer.NumberRows):
		    count = 0
		    for (r, c) in cells:
			if (r == row):
			    count = count + 1
		    if (count == len(self.sequences[row])):
			chain_deletes.append(row)
			# Take these out of cells so we don't delete them one-by-one
			for i in range(len(cells)-1, -1, -1):
			    if (cells[i][0] == row):
				cells.pop(i)
		for (r, c) in cells:
		    if (c >= len(self.indxToSeqPos[r])):
			continue
		    logInfo("Deleting " + self.IDs[r] + " position " + str(c+1))
		    # If this is the last residue in a row, then delete the whole row
		    if (len(self.sequences[r]) == 1):
			self.deleteChain(r)
		    else:
			# Remove this residue from the sequence and pop the indxToSeqPos element
			if (self.AlignBtn.GetLabel() == "From 1"):
			    self.sequences[r] = self.sequences[r][0:c] + self.sequences[r][c+1:]
			    resID = self.indxToSeqPos[r].pop(c)
			else:
			    self.sequences[r] = self.sequences[r][0:c] + "-" + self.sequences[r][c+1:]
			    resID = self.indxToSeqPos[r][c]
			    self.indxToSeqPos[r][c] = "-"
			pyMOLPos = resID[1]
			self.selectLookup.pop((self.IDs[r], int(pyMOLPos)))
			# Find the pose that this residue belongs to
			poseloc = r
			while (not(self.poses[poseloc])):
			    poseloc = poseloc - 1
			# Find the real index of the residue in this pose by summing sequence lengths
			indx = 0
			for i in range(poseloc, r):
			    indx = indx + len(self.sequences[i])
			indx = indx + c + 1
			chain = self.IDs[r][len(self.IDs[r])-1]
			if (chain == "_"):
			    chain = " "
			self.poses[r][0][chain].detach_child(resID) # You have to give it the tuple otherwise it will screw up on HETATMs
			# IMPORTANT: You have to replace the model in the sandbox with the new model
			currID = self.getModelForChain(poseloc)
			self.pdbwriter.set_structure(self.poses[poseloc])
			self.pdbwriter.save(currID + ".pdb")
			# Update the SeqViewer to reflect this deletion
			for i in range(c, len(self.sequences[r])):
			    self.SeqViewer.SetCellValue(r, i, self.SeqViewer.GetCellValue(r, i+1))
			self.SeqViewer.SetCellValue(r, len(self.sequences[r]), "")
			# Remove this residue from PyMOL
			fields = self.IDs[r].split("|")
			currID = ""
			for i in range(0, len(fields)-1):
			    currID = currID + fields[i] + "|"
			currID = currID[0:len(currID)-1]
			chainID = fields[len(fields)-1]
			if (chainID != "_"):
			    self.cmd.remove("model " + currID + " and chain " + chainID + " and resi " + str(pyMOLPos))
			else:
			    self.cmd.remove("model " + currID + " and resi " + str(pyMOLPos))
		# Now delete whole chains that we may have identified
		for i in range(len(chain_deletes)-1, -1, -1):
		    self.deleteChain(chain_deletes[i])
		#self.recolorResidues()
		self.regenerateLookupTable()
	    else:
		logInfo("Canceled DELETE key press")
	    dlg.Destroy()
	    # Shave off any extra columns
	    maxseqlen = 0
	    for i in range(0, len(self.sequences)):
		if (len(self.sequences[i]) > maxseqlen):
		    maxseqlen = len(self.sequences[i])
	    for c in range(self.SeqViewer.NumberCols, maxseqlen, -1):
		self.SeqViewer.DeleteCols(c-1)
	    self.protPanel.activate()
	elif (int(event.GetKeyCode()) == wx.WXK_ESCAPE):
	    self.SeqViewer.ClearSelection()
	    self.PyMOLUpdateTimer.Stop()
	    self.PyMOLUpdateTimer.Start(100)
	event.Skip()

    def updatePyMOLSelection(self, event):
	self.PyMOLUpdateTimer.Stop()
	self.selectUpdate(updatePyMOL=True)

    def leftRelease(self, event):
	# User released the mouse on the grid, probably because they were selecting things
	# Update PyMOL
	# A timer is needed because the selection data is not available at the time this function executes
	# If this was a click and release on the same element, set that element to be the only thing selected
	# This is needed because otherwise you have to click and hold and drag a little to only select on box
	try:
	    if (self.clickc == self.colpos and self.clickr == self.rowpos):
		if (not(event.ControlDown())):
		    self.SeqViewer.ClearSelection()
		self.SeqViewer.SelectBlock(self.clickr, self.clickc, self.clickr, self.clickc, True)
	    self.PyMOLUpdateTimer.Stop()
	    self.PyMOLUpdateTimer.Start(100)
	except:
	    # This first click was on the grid, but not on a cell
	    pass
	event.Skip()

    def rightClick(self, event):
	# Use this to unselect everything because now clicking a single element selects that one element
	# so we need a way to unselect everything without having to do it in PyMOL
	self.SeqViewer.ClearSelection()
	self.PyMOLUpdateTimer.Stop()
	self.PyMOLUpdateTimer.Start(100)
	event.Skip()

    def leftClick(self, event):
	# This only is used on Windows
	# For some reason CTRL selections on Windows highlight the boxes but don't register as selections so 
	# I have to do a selection from the code
	# If CTRL is not held down then the event should be skipped so weird behavior doesn't ensue for normal
	# selections and SHIFT selections
	c = event.GetCol()
	r = event.GetRow()
	self.clickc = c
	self.clickr = r
	if (platform.system() != "Windows"):
	    event.Skip()
	    return
	if (event.ControlDown()):
	    # Find out if this cell is selected already or not
	    selected = False
	    topLefts = self.SeqViewer.GetSelectionBlockTopLeft()
	    bottomRights = self.SeqViewer.GetSelectionBlockBottomRight()
	    for i in range(0, len(topLefts)):
		if (r >= topLefts[i][0] and r <= bottomRights[i][0] and c >= topLefts[i][1] and c <= bottomRights[i][1]):
		    selected = True
		    break
	    if (not(selected)):
		self.SeqViewer.SelectBlock(r, c, r, c, True)
	    else:
		self.SeqViewer.DeselectCell(r, c)
	else:
	    event.Skip()
	
    def labelClick(self, event):
	# This is a dummy function so single clicks on the labels do nothing because I am not skipping the event
	pass
    
    def labelDClick(self, event):
	# For some odd reason, normal label clicks, although they look like they select the row/column, do not
	# actually select the row/column in the code so the Timer never picks up the selection
	# Here is an attempt to do it manually
	r = event.GetRow()
	c = event.GetCol()
	if (r >= 0):
	    for c in range(0, self.SeqViewer.NumberCols):
		self.SeqViewer.SelectBlock(r, c, r, c, True)
	elif (c >= 0):
	    for r in range(0, self.SeqViewer.NumberRows):
		self.SeqViewer.SelectBlock(r, c, r, c, True)
	self.selectUpdate(updatePyMOL=True)
	# Do not skip the event, because I don't want anything strange happening
	
    def getSelectedResidues(self):
	# Useful function for getting selected information
	# First we have to update the selection in case the user was fiddling in PyMOL and
	# then when straight to a protocol.  The sequence window doesn't know about the PyMOL
	# changes yet but the protocols read selections from it
	self.selectUpdate(None)
	topLefts = self.SeqViewer.GetSelectionBlockTopLeft()
	bottomRights = self.SeqViewer.GetSelectionBlockBottomRight()
	selectionkey = []
	for i in range(0, len(topLefts)):
	    for r in range(topLefts[i][0], bottomRights[i][0]+1):
		if (r >= len(self.poses)):
		    continue
		offset = 0
		for j in range(r, -1, -1):
		    if (self.poses[j]):
			poseindx = j
			break
		    offset = offset + len(self.sequences[j-1])
		for c in range(topLefts[i][1], bottomRights[i][1]+1):
		    if (c >= len(self.indxToSeqPos[r])):
			continue
		    if (self.indxToSeqPos[r][c] != "-"):
			selectionkey.append(str(r) + ":" + str(poseindx) + ":" + str(c+offset))
	# Sort
	selectionkey.sort()
	# Now change it back to meaningful integers
	lastr = -1
	selection = []
	for key in selectionkey:
	    r = int(key.split(":")[0])
	    poseindx = int(key.split(":")[1])
	    resi = int(key.split(":")[2])
	    if (lastr != r):
		selection.append([r, poseindx, []])
	    selection[len(selection)-1][2].append(resi)
	return selection
    
    def numChains(self):
	return self.SeqViewer.NumberRows
    
    def getChainPose(self, row):
	# Gets a pose object for only one chain given either by the row in the sequence window or the ID
	if (not("int" in str(type(row)))):
	    row = self.IDs.index(row)
	poseindx = self.getPoseIndex(row)
	chain = self.IDs[row][len(self.IDs[row])-1]
	if (chain == "_"):
	    chain = " "
	return self.poses[poseindx][0][chain]
    
    def getModelForChain(self, r):
	# Useful function to figure out what the name of the model is from an ID in the SeqViewer
	ID = self.IDs[r]
	model = ""
	fields = ID.split("|")
	for field in fields[0:len(fields)-1]:
	    model = model + field + "|"
	model = model[0:len(model)-1]
	return model
    
    def getPoseIndex(self, r):
	# Given a row in the SeqViewer, find out the index of the pose for this sequence
	poseindx = -1
	for i in range(r, -1, -1):
	    if (self.poses[i]):
		poseindx = i
		break
	return poseindx
    
    def getPoseIndexForModel(self, model):
	# Given a model, find the poseindx
	poseindx = -1
	for i in range(0, self.SeqViewer.NumberRows):
	    if (self.poses[i] and model.strip() == self.IDs[i][0:len(self.IDs[i].strip())-2]):
		poseindx = i
		break
	return poseindx
    
    def getRosettaIndex(self, model, chain, seqpos):
	# Useful function for finding what the Rosetta index of a model/chain/seqpos is
	poseindx = self.getPoseIndexForModel(model)
	ires = 1
	for ch in self.poses[poseindx][0]:
	    for residue in ch:
		thischain = ch.id
		thisseqpos = residue.id[1]
		if (chain == "_" and (thischain == "_" or len(thischain.strip()) == 0) and int(thisseqpos) == int(seqpos)):
		    return ires
		elif (chain == thischain and int(thisseqpos) == int(seqpos)):
		    return ires
		ires = ires + 1
	return -1
    
    def getResidueInfo(self, model, rosettaindx):
	# Useful function for getting the chain and seqpos given a model and rosetta index
	ires = 1
	poseindx = self.getPoseIndexForModel(model)
	for ch in self.poses[poseindx][0]:
	    for residue in ch:
		if (ires == rosettaindx):
		    return (ch.id, residue.id[1])
		ires = ires + 1
	return None
    
    def getIsCanonicalAA(self, r, c):
	# Find out if this is a canonical CAA or not
	poseindx = self.getPoseIndex(r)
	offset = 1
	for i in range(poseindx, r):
	    offset = offset + len(self.sequences[i])
	chain = self.IDs[r][len(self.IDs[r])-1]
	if (chain == "_"):
	    chain = " "
	restype = self.poses[poseindx][0][chain][self.indxToSeqPos[r][c]].resname
	if (not(restype in "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR ")):
	    return False
	return True
    
    def getSelectedModels(self):
	# Gets the names of the models that are currently selected for easy PyMOL access
	topLefts = self.SeqViewer.GetSelectionBlockTopLeft()
	bottomRights = self.SeqViewer.GetSelectionBlockBottomRight()
	models = []
	rowsSelected = []
	for i in range(0, len(topLefts)):
	    for r in range(topLefts[i][0], bottomRights[i][0]+1):
		if (not(r in rowsSelected)):
		    rowsSelected.append(r)
	rowsSelected.sort()
	for r in rowsSelected:
	    newmodel = self.getModelForChain(r)
	    if (not(newmodel in models)):
		models.append(newmodel)
	return models
    
    def getSelectedChains(self):
	# Gets all the selected chains
	topLefts = self.SeqViewer.GetSelectionBlockTopLeft()
	bottomRights = self.SeqViewer.GetSelectionBlockBottomRight()
	models = []
	rowsSelected = []
	for i in range(0, len(topLefts)):
	    for r in range(topLefts[i][0], bottomRights[i][0]+1):
		if (not(r in rowsSelected)):
		    rowsSelected.append(r)
	rowsSelected.sort()
	return rowsSelected
    
    def areResiduesSelected(self):
	# Used by the selection panel to see if residues are selected
	# If not, the button usually defaults to applying the operation to everything
	return (len(self.SeqViewer.GetSelectionBlockTopLeft()) > 0)
    
    def recolorClick(self, event):
	if (self.colorMode == "No Color"):
	    self.colorMode = "SS"
	    if (platform.system() == "Darwin"):
		self.ColoringBtn.SetBitmapLabel(bitmap=wx.Image("images/osx/ColoringBtn_SS.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
	    else:
		self.ColoringBtn.SetLabel(self.colorMode)
	    logInfo("Colored residues by SS")
	elif (self.colorMode == "SS"):
	    self.colorMode = "BFactor"
	    if (platform.system() == "Darwin"):
		self.ColoringBtn.SetBitmapLabel(bitmap=wx.Image("images/osx/ColoringBtn_BFactor.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
	    else:
		self.ColoringBtn.SetLabel(self.colorMode)
	    logInfo("Colored residues by B-factor")
	else:
	    self.colorMode = "No Color"
	    if (platform.system() == "Darwin"):
		self.ColoringBtn.SetBitmapLabel(bitmap=wx.Image("images/osx/ColoringBtn_NoColor.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
	    else:
		self.ColoringBtn.SetLabel(self.colorMode)
	    logInfo("Turned residue coloring off")
	self.recolorResidues()
    
    def recolorResidues(self):
	# This function will recolor the cells in the SeqViewer to either a default white mode,
	# a secondary structure coloring, or a B-factor coloring
	if (self.colorMode == "SS"):
	    # Default everything to white
	    for r in range(0, self.SeqViewer.NumberRows):
		for c in range(0, self.SeqViewer.NumberCols):
		    self.SeqViewer.SetCellBackgroundColour(r, c, "white")
		    self.SeqViewer.SetCellTextColour(r, c, "black")
	    # For each row in the SeqViewer, select the helices and B-sheets, get the selection
	    # in here, and color appropriately
	    for r in range(0, self.SeqViewer.NumberRows):
		model = self.getModelForChain(r)
		chain = self.IDs[r][len(self.IDs[r])-1]
		if (chain != "_" and chain != " " and chain != ""):
		    self.cmd.select("colorsele", "model " + model + " and chain " + chain)
		else:
		    self.cmd.select("colorsele", "model " + model)
		self.stored.selected = []
		self.cmd.iterate_state(1, "colorsele and name ca", "stored.selected.append(str(resi)+\"|\"+ss)")
		#self.stored.selected = list(set(self.stored.selected))
		for code in self.stored.selected:
		    resi = int(code.split("|")[0])
		    ss = code.split("|")[1].upper()
		    (color, tcolor) = self.sscolormap[ss]
		    for x in range(0, len(self.indxToSeqPos[r])):
			if (self.indxToSeqPos[r][x] != "-" and self.indxToSeqPos[r][x][1] == int(resi)):
			    c = x
			    break
		    self.SeqViewer.SetCellBackgroundColour(r, c, color)
		    self.SeqViewer.SetCellTextColour(r, c, tcolor)
	    self.cmd.delete("colorsele")
	elif (self.colorMode == "BFactor"):
	    # A B-Factor of 0 is blue, 100 is red, and everything else is a gradient in between
	    # Default everything to white
	    for r in range(0, self.SeqViewer.NumberRows):
		for c in range(0, self.SeqViewer.NumberCols):
		    self.SeqViewer.SetCellBackgroundColour(r, c, "white")
		    self.SeqViewer.SetCellTextColour(r, c, "black")
	    r = 0
	    while (r < self.SeqViewer.NumberRows):
		if (self.poses[r]):
		    # Skip over the blank entries
		    for ch in self.poses[r][0]:
			c = 0
			for residue in ch:
			    try:
				bfactor = residue["N"].get_bfactor()
			    except:
				# This NCAA does not have an "N" named atom
				continue
			    if (bfactor < 0):
				bfactor = 0.0
			    elif (bfactor > 100):
				bfactor = 100.0
			    red = int(2.55 * bfactor)
			    blue = 255 - red
			    color = "#%02x%02x%02x" % (red, 0, blue)
			    (row, col) = self.selectLookup[(self.IDs[r], residue.id[1])]
			    self.SeqViewer.SetCellBackgroundColour(row, col, color.upper())
			    self.SeqViewer.SetCellTextColour(row, col, "white")
			    c = c + 1
			r = r + 1
		else:
		    r = r + 1
	else:
	    # Default everything to white
	    for r in range(0, self.SeqViewer.NumberRows):
		for c in range(0, self.SeqViewer.NumberCols):
		    self.SeqViewer.SetCellBackgroundColour(r, c, "white")
		    self.SeqViewer.SetCellTextColour(r, c, "black")
	self.cmd.enable("sele")
	self.SeqViewer.Refresh() # Needed to repaint otherwise the changes don't occur
    
    def stereoToggle(self, event):
	if (self.viewMode == "Mono"):
	    self.viewMode = "Crosseye"
	    if (platform.system() == "Darwin"):
		self.StereoBtn.SetBitmapLabel(bitmap=wx.Image("images/osx/StereoBtn_Crosseye.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
	    else:
		self.StereoBtn.SetLabel(self.viewMode)
	    self.cmd.stereo("crosseye")
	    logInfo("Turned stereo crosseye on")
	elif (self.viewMode == "Crosseye"):
	    self.viewMode = "Walleye"
	    if (platform.system() == "Darwin"):
		self.StereoBtn.SetBitmapLabel(bitmap=wx.Image("images/osx/StereoBtn_Walleye.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
	    else:
		self.StereoBtn.SetLabel(self.viewMode)
	    self.cmd.stereo("walleye")
	    logInfo("Turned stereo walleye on")
	elif (self.viewMode == "Walleye"):
	    self.viewMode = "Quadbuffer"
	    if (platform.system() == "Darwin"):
		self.StereoBtn.SetBitmapLabel(bitmap=wx.Image("images/osx/StereoBtn_Quadbuffer.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
	    else:
		self.StereoBtn.SetLabel(self.viewMode)
	    self.cmd.stereo("quadbuffer")
	    logInfo("Turned stereo quadbuffer on")
	else:
	    self.viewMode = "Mono"
	    if (platform.system() == "Darwin"):
		self.StereoBtn.SetBitmapLabel(bitmap=wx.Image("images/osx/StereoBtn_Mono.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
	    else:
		self.StereoBtn.SetLabel(self.viewMode)
	    self.cmd.stereo("off")
	    logInfo("Turned stereo off")
    
    def bgRecolor(self, event):
	logInfo("Clicked on the background color button")
	dlg = wx.ColourDialog(self)
	dlg.GetColourData().SetChooseFull(True)
	if (dlg.ShowModal() == wx.ID_OK):
	    data = dlg.GetColourData()
	    mycolor = "0x%02x%02x%02x" % data.GetColour().Get()
	    self.cmd.bg_color(mycolor)
	    logInfo("Set background color to " + mycolor)
	else:
	    logInfo("Cancelled out of background color selection")
	dlg.Destroy()

    def selectInPyMOL(self, r, c, model, chain, res, resend, first):
	self.selectedResidues.append([r, c])
	# If this the first selection, we have to create "sele"
	if (first):
	    # We have to select by model, chain, and residue index using a bunch of "and" operators
	    if ("ABCDEFGHIJKLMNOPQRSTUVWXYZ".find(chain) < 0):
		self.cmd.select("temp", "model " + str(model) + " and resi " + str(res) + "-" + str(resend))
	    else:
		self.cmd.select("temp", "model " + str(model) + " and chain " + str(chain) + " and resi " + str(res) + "-" + str(resend))
	else:
	    # To get PyMOL to add to the existing selection "sele", we used the "or" operator to union "sele"
	    # with the new atoms
	    if ("ABCDEFGHIJKLMNOPQRSTUVWXYZ".find(chain) < 0):
		self.cmd.select("temp", "temp or (model " + str(model) + " and resi " + str(res) + "-" + str(resend) + ")")
	    else:
		self.cmd.select("temp", "temp or (model " + str(model) + " and chain " + str(chain) + " and resi " + str(res) + "-" + str(resend) + ")")
    
    def selectUpdate(self, updatePyMOL=False):
	# Find out which residues are selected in PyMOL and highlight the corresponding positions in the viewer
	try:
	    self.stored.selected = []
	except:
	    return
	try:
	    # There's this really annoying feature where if you click in the background in PyMOL, it looks like
	    # it deselected everything, but in reality it only disabled the selection, so when you come back here
	    # it sees there's still a selection and reselects everything again
	    # So first make sure the "seqsele" selection is enabled, if it's not then assume nothing is selected
	    # IMPORTANT NOTE: The "enabled" feature does not seem to be recognized properly if the selection is
	    # "sele", which is the default selection name.  You have to name it something besides "sele" to get
	    # this trick to work, hence I used "seqsele"
	    # It's just annoying because the default selection is "sele" so you have to make sure everything ends up
	    # in "seqsele" and not "sele"
	    # If there are issues, try getting rid of the "-i" option in InteractiveROSETTA.py for pymol initialization
	    # because then that internal side window showing the selections will be visible enabling you to see what
	    # is happening with the selections more easily
	    # 
	    # We have to try to get everything into the "seqsele" name, but sometimes if the user clicks on things in PyMOL
	    # they end up in sele, such that there is sele and seqsele
	    # If you don't combine the two, then everything in sele ends up getting dropped
	    try:
		# This will throw an error if there is nothing in seqsele, so maybe only sele is defined
		self.cmd.select("seqsele", "seqsele")
		try:
		    # This will throw an error if sele is not defined, so this determines whether it is [seqsele or seqsele/sele]
		    self.cmd.select("seqsele", "sele or seqsele")
		    defined = True
		except:
		    # This means only seqsele was defined
		    pass
	    except:
		try:
		    # Seqsele was not defined, so maybe only sele is
		    self.cmd.select("seqsele", "sele")
		except:
		    # Neither was defined, there's nothing selected at all
		    pass
	    defined = True
	    # At this point, everything should be in seqsele unless nothing is selected and seqsele should be enabled
	    # If it is not enabled, then "seqsele" is still defined but nothing is highlighted, in which case nothing is desired
	    # to be selected so we should remove the seqsele name so iterate_state doesn't see these selections
	    # Delete sele first otherwise it will keep coming back
	    try:
		self.cmd.delete("sele")
	    except:
		pass
	    if (not("seqsele" in self.cmd.get_names("selections", enabled_only=1))):
		try:
		    self.cmd.delete("seqsele")
		except:
		    pass
		defined = False
	    if (defined):
		try:
		    self.cmd.iterate_state(1, "seqsele", "stored.selected.append([model, chain, resi, resn])")
		    # If nothing is in stored but yet "seqsele" is still defined, delete it because some of the
		    # selecting options assume that if seqsele is empty it doesn't exist
		    if (len(self.stored.selected) == 0):
			self.cmd.delete("seqsele")
		except:
		    pass
	except:
	    pass
	#self.stored.selected = list(set(self.stored.selected)) # Removes all duplicate elements
	#self.stored.selected.sort()
	first = True
	amISelected = []
	for r in range(0, self.SeqViewer.NumberRows):
	    amISelected.append([])
	    for c in range(0, len(self.sequences[r])):
		amISelected[r].append(False)
	# If the SequenceWindow is in focus, then the user is fiddling with selecting residues in this window
	# so save the current selected residues BEFORE we add the ones from PyMOL (because otherwise PyMOL will
	# overwrite the selections in this window)
	if (updatePyMOL):
	    self.cmd.delete("seqsele")
	    self.cmd.delete("temp")
	    topLefts = self.SeqViewer.GetSelectionBlockTopLeft()
	    bottomRights = self.SeqViewer.GetSelectionBlockBottomRight()
	    for i in range(0, len(topLefts)):
		for x in range(topLefts[i][0], bottomRights[i][0]+1):
		    for y in range(topLefts[i][1], bottomRights[i][1]+1):
			# If this is a blank space, skip it
			if (y >= len(self.indxToSeqPos[x])):
			    continue
			model = ""
			# The chainID is appended to the ID via a "|"
			# This code recovers the model name if the PDB happened to have "|" as part of the filename
			fields = self.IDs[x].split("|")
			for field in fields[0:len(fields)-1]:
			    model = model + field + "|"
			model = model[0:len(model)-1] # Trim off the last "|"
			chain = fields[len(fields)-1]			    
			# res needs to be the actual sequence index from the PDB, not the renumbered-from-1 index
			#res = self.indxToSeqPos[x][y]
			# Save this as a selected residue
			amISelected[x][y] = True
			#self.selectCell(x, y, model, chain, res, first)
	    for r in range(0, len(amISelected)):
		inblock = False
		for c in range(0, len(amISelected[r])):
		    if (self.indxToSeqPos[r][c] == "-"):
			# Skip these, everything on both sides of the dashes should be continuous
			continue
		    elif (amISelected[r][c] and not(inblock)):
			inblock = True
			start = self.indxToSeqPos[r][c][1]
			lastc = c
		    elif (inblock and int(self.indxToSeqPos[r][c][1]) - int(self.indxToSeqPos[r][lastc][1]) < 0):
			# The residues aren't always arranged sequentially in the PDB ordering, so we have to
			# watch out for this
			model = self.IDs[r][0:len(self.IDs[r])-2]
			chain = self.IDs[r][len(self.IDs[r])-1]
			self.selectInPyMOL(r, c, model, chain, start, str(self.indxToSeqPos[r][c-1][1]), first)
			start = self.indxToSeqPos[r][c][1]
			first = False
		    elif (not(amISelected[r][c]) and inblock):
			inblock = False
			model = self.IDs[r][0:len(self.IDs[r])-2]
			chain = self.IDs[r][len(self.IDs[r])-1]
			self.selectInPyMOL(r, c, model, chain, start, str(self.indxToSeqPos[r][lastc][1]), first)
			first = False
		    elif (inblock and amISelected[r][c]):
			lastc = c
		if (inblock):
		    if (self.sequences[r][len(self.sequences[r])-1] != "-"):
			lastc = len(self.sequences[r])-1
		    model = self.IDs[r][0:len(self.IDs[r])-2]
		    chain = self.IDs[r][len(self.IDs[r])-1]
		    self.selectInPyMOL(r, c, model, chain, start, str(self.indxToSeqPos[r][lastc][1]), first)
		    first = False
	    if (not(first)):
		self.cmd.select("seqsele", "temp")
		self.cmd.delete("temp")
		self.cmd.delete("sele")
	else:
	    # Take the selection data from PyMOL and select those residues in the SequenceWindow
	    for atom in self.stored.selected:
		res = atom[2]
		chain = atom[1]
		if (chain == ""):
		    chain = "_"
		model = atom[0]
		if (model == "minimized_view"):
		    # We're selecting a PyMOL_Mover minimize preview so figure out what the real model is to select
		    # the real residues in the Sequence Window
		    model = self.protPanel.selectedModel
		elif (model == "designed_view"):
		    # We're selecting a PyMOL_Mover fixbb preview so figure out what the real model is to select
		    # the real residues in the Sequence Window
		    model = self.protPanel.selectedModel
		try:
		    (r, c) = self.selectLookup[(model + "|" + chain, int(res))]
		except:
		    # Some unrecognized view from a protocol
		    return
		#c = int(res)
		# The last Boolean value is set to "True" if you are adding to the current selection, "False" starts
		# a new selection, so we only want it to be "False" the first time through
		# There doesn't seem to be a SelectCell function that works, so we used SelectBlock instead
		amISelected[r][c] = True
		#self.SeqViewer.SelectBlock(r, c, r, c, not(first))
		#first = False
	    for r in range(0, len(amISelected)):
		inblock = False
		for c in range(0, len(amISelected[r])):
		    if (amISelected[r][c] and not(inblock)):
			inblock = True
			start = c
		    elif (not(amISelected[r][c]) and inblock):
			inblock = False
			self.SeqViewer.SelectBlock(r, start, r, c-1, not(first))
			first = False
		if (inblock):
		    self.SeqViewer.SelectBlock(r, start, r, len(amISelected[r])-1, not(first))
		    first = False
	    if (first):
		# Nothing selected, deselect everything
		for r in range(0, self.SeqViewer.NumberRows):
		    self.SeqViewer.DeselectRow(r)
	self.logSelection()
	self.cmd.enable("seqsele")
    
    def logSelection(self):
	# Periodically log which residues are selected
	topLefts = self.SeqViewer.GetSelectionBlockTopLeft()
	bottomRights = self.SeqViewer.GetSelectionBlockBottomRight()
	for i in range(0, len(topLefts)):
	    if (topLefts[i][0] >= len(self.IDs)):
		continue
	    topLeftPos = self.IDs[topLefts[i][0]] + " position " + str(topLefts[i][1]+1)
	    bottomRightPos = self.IDs[bottomRights[i][0]] + " position " + str(bottomRights[i][1]+1)
	    logInfo("New Selection Entry: ")
	    logInfo("Selected " + topLeftPos + " to " + bottomRightPos)
    
    def updateSeqViewer(self, rowToUpdate=-1, relabel=False):
	newnrows = len(self.sequences)
	newncols = 0
	for i in range(0, newnrows):
	    if (len(self.sequences[i]) > newncols):
		newncols = len(self.sequences[i])
	if (newncols > self.SeqViewer.NumberCols or relabel):
	    if (newncols > self.SeqViewer.NumberCols):
		self.SeqViewer.AppendCols(newncols - self.SeqViewer.NumberCols)
	    for i in range(0, self.SeqViewer.NumberCols):
		self.SeqViewer.SetColLabelValue(i, "A")
	    self.SeqViewer.AutoSizeColumns(True)
	    self.SeqViewer.SetColLabelValue(0, "")
	    for i in range(0, self.SeqViewer.NumberCols):
		if (self.AlignBtn.GetLabel() != "PDB #"):
		    if ((i+1) % 5 == 0):
			if (i+1 < 100):
			    self.SeqViewer.SetColLabelValue(i, str(i+1))
			else:
			    val = (i+1) % 100
			    if (val == 0):
				val = "00"
			    elif (val < 10):
				val = "0" + str(val)
			    else:
				val = str(val)
			    val2 = str((i+1) / 100)
			    self.SeqViewer.SetColLabelValue(i, val)
			    self.SeqViewer.SetColLabelValue(i-1, val2)
		    else:
			self.SeqViewer.SetColLabelValue(i, "")
		else:
		    # We have to read the labels from the list generated by realign
		    indx = self.columnlabels[i]
		    if (i > 0 and indx - self.columnlabels[i-1] != 1):
			foundbreak = True
		    elif (i == 0 and indx != 1):
			foundbreak = True
		    else:
			foundbreak = False
		    if (indx % 5 == 0 or foundbreak):
			# If the position right after a break would overwrite the entry before it, skip this one
			if (i > 0 and self.SeqViewer.GetColLabelValue(i-1).strip() != "" and indx >= 100):
			    self.SeqViewer.SetColLabelValue(i, "")
			    continue
			if (indx < 100):
			    self.SeqViewer.SetColLabelValue(i, str(indx))
			else:
			    val = indx % 100
			    if (val == 0):
				val = "00"
			    elif (val < 10):
				val = "0" + str(val)
			    else:
				val = str(val)
			    val2 = str(indx / 100)
			    self.SeqViewer.SetColLabelValue(i, val)
			    self.SeqViewer.SetColLabelValue(i-1, val2)
		    else:
			self.SeqViewer.SetColLabelValue(i, "")
	elif (newncols < self.SeqViewer.NumberCols):
	    self.SeqViewer.DeleteCols(newncols, self.SeqViewer.NumberCols - newncols)
	if (newnrows > self.SeqViewer.NumberRows and rowToUpdate < 0):
	    self.SeqViewer.AppendRows(newnrows - self.SeqViewer.NumberRows)
	    self.SeqViewer.SetRowLabelValue(newnrows-1, self.IDs[newnrows-1])
	    for i in range(0, len(self.sequences[newnrows-1])):
		self.SeqViewer.SetCellValue(newnrows-1, i, self.sequences[newnrows-1][i])
	elif (rowToUpdate >= 0):
	    for i in range(0, self.SeqViewer.NumberCols):
		self.SeqViewer.SetCellValue(rowToUpdate, i, "")
		if (i < len(self.sequences[rowToUpdate])):
		    self.SeqViewer.SetCellValue(rowToUpdate, i, self.sequences[rowToUpdate][i])
	    self.SeqViewer.SetRowLabelValue(rowToUpdate, self.IDs[rowToUpdate])
	readOnly = wx.grid.GridCellAttr()
	readOnly.SetReadOnly(True)
	for i in range(newnrows-1, newnrows):
	    self.SeqViewer.SetRowAttr(i, readOnly)
	# Resize the label width to fit long PDB filenames
	font = self.SeqViewer.GetFont()
	dc = wx.WindowDC(self.SeqViewer)
	dc.SetFont(font)
	maxwidth = 40
	for i in range(0, newnrows):
	    (w, h) = dc.GetTextExtent(self.SeqViewer.GetRowLabelValue(i))
	    if (w > maxwidth):
		maxwidth = w
	self.SeqViewer.DisableDragColSize()
	self.SeqViewer.DisableDragRowSize()
	self.SeqViewer.SetRowLabelSize(maxwidth+10)
	
    def rerunDSSPForModel(self, model, pdbmodel="None"):
	try:
	    if (self.dsspexe == "N/A"):
		return
	    if (pdbmodel == "None"):
		pdbmodel = model + ".pdb"
	    poseindx = self.getPoseIndexForModel(model)
	    dsspdata = Bio.PDB.DSSP(self.poses[poseindx][0], pdbmodel, dssp=self.dsspexe)
	    for c in self.poses[poseindx][0]:
		chain = c.id
		if (len(chain.strip()) == 0):
		    chain = "_"
		for r in c:
		    # Alter the states in PyMOL with the DSSP secondary structure predictions
		    try:
			ss = dsspdata[(c.id, r.id[1])][1]
			# Not all the DSSP codes match up with PyMOL, so make them match
			if (ss == "E"):
			    # Beta sheet
			    ss = "S"
			elif (ss == "S"):
			    # Bend
			    ss = "B"
			elif (ss == "-"):
			    # Loop
			    ss = "L"
			if (c.id != " "):
			    self.cmd.alter("model " + model + " and chain " + c.id + " and resi " + str(r.id[1]), "ss=\"" + ss + "\"")
			else:
			    self.cmd.alter("model " + model + " and resi " + str(r.id[1]), "ss=\"" + ss + "\"")
		    except:
			ss = "L"
	except:
	    # In case something goes wrong, just use PyMOL's SS coloring and don't crash
	    pass
	
    def reloadPose(self, poseindx, model, newpdb):
	# This function reloads a pose from a protocol and reloads the sequence, which may
	# have changed
	# Create a copy of the inputted pose
	self.poses[poseindx] = self.pdbreader.get_structure(model, newpdb)
	# Now reload the sequences
	i = 0
	for ch in self.poses[poseindx][0]:
	    r = poseindx + i
	    seq = ""
	    self.indxToSeqPos[r] = []
	    for residue in ch:
		seq = seq + AA3to1(residue.resname)
		self.indxToSeqPos[r].append(residue.id)
	    self.sequences[r] = seq
	    # Some of the protocols change the chain IDs, so we have to update the IDs and labels
	    chain = ch.id
	    if (len(chain.strip()) == 0):
		chain = "_"
	    self.IDs[r] = model + "|" + chain
	    self.updateSeqViewer(r)
	    i = i + 1
	# Rerun DSSP
	self.rerunDSSPForModel(model, newpdb)
	self.regenerateLookupTable()
	
    def reloadFromPyMOL(self, poseindx):
	# This function reloads the BioPython data dumped from PyMOL
	# You need this for docking in case the user was rotating the binding partners in PyMOL before the
	# docking simulation
	model = self.IDs[poseindx][0:len(self.IDs[poseindx])-2]
	self.cmd.save(model + ".pdb", model)
	fixPyMOLSave(model + ".pdb")
	self.poses[poseindx] = self.pdbreader.get_structure(model, model + ".pdb")
    
    def fetchClick(self, event):
	# You have to do this on Linux due to the libc bug
	if (platform.system() != "Windows"):
	    res_init()
	# Attempt to get a PDB directly from RCSB so the user doesn't have to download it manually
	pdbCode = str(self.RCSBTxt.GetValue().strip()[0:4])
	pdbUrl = 'http://www.rcsb.org/pdb/downloadFile.do?fileFormat=pdb&compression=NO&structureId='+pdbCode
	try:
	    logInfo("Attempting to fetch PDB " + pdbCode + " from RCSB")
	    try:
		pdbFile = urllib2.urlopen(pdbUrl, timeout=1)
	    except urllib2.URLError as e:
		# Probably no Internet connection
		raise Exception("Could not connect to RSCB.  Do you have an Internet connection?")
	    goToSandbox()
	    validfile = False
	    f = open(pdbCode + ".pdb", "w")
	    for aline in pdbFile:
		if (aline[0:4] == "ATOM"):
		    # Because a failed result returns HTML code, which will crash Rosetta if you try
		    # to load HTML into a code
		    validfile = True
		f.write(aline.strip() + "\n")
	    f.close()
	    if (validfile):
		self.labelMsg.SetLabel("Fetching PDB from RCSB, please be patient...")
		self.msgQueue.append("Fetching PDB from RCSB, please be patient...")
		self.Disable()
		self.protWin.Disable()
		# Ask if the user would like to save a copy before loading into Rosetta
		dlg = wx.MessageDialog(self, "Save a copy of this PDB before loading into Rosetta?", "Save RCSB PDB", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
		if (dlg.ShowModal() == wx.ID_YES):
		    while (True):
			dlg2 = wx.FileDialog(
			    self, message="Save a PDB File",
			    defaultDir=self.cwd,
			    defaultFile=pdbCode,
			    wildcard="PDB Files (*.pdb)|*.pdb",
			    style=wx.SAVE | wx.CHANGE_DIR)
			if (dlg2.ShowModal() == wx.ID_OK):
			    paths = dlg2.GetPaths()
			    # Change cwd to the last opened file
			    lastDirIndx = paths[len(paths)-1].rfind("/")
			    self.cwd = str(paths[len(paths)-1][0:lastDirIndx])
			    self.saveWindowData(None)
			    # Load the PDBs into PyMOL
			    filename = str(paths[0]).split(".pdb")[0] + ".pdb"
			    # Does it exist already?  If so, ask if the user really wants to overwrite it
			    if (os.path.isfile(filename)):
				dlg3 = wx.MessageDialog(self, "The file " + filename + " already exists.  Overwrite it?", "Filename Already Exists", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
				if (dlg3.ShowModal() == wx.ID_NO):
				    dlg3.Destroy()
				    dlg2.Destroy()
				    logInfo("Cancelled save operation due to filename already existing")
				    continue
				dlg3.Destroy()
			    goToSandbox()
			    f = open(pdbCode + ".pdb", "r")
			    f2 = open(str(filename), "w")
			    for aline in f:
				f2.write(aline.strip() + "\n")
			    f2.close()
			    f.close()
			    break
			else:
			    # User cancelled probably, get out of the loop
			    break
		    dlg2.Destroy()
		dlg.Destroy()
		# We have to check to see if there are multiple models in this PDB
		# If there are, then we're going to break them up into separate PDB files
		# and read them all as distinct models in PyMOL
		nummodels = 0
		filename = pdbCode + ".pdb"
		f = open(filename, "r")
		for aline in f:
		    if (aline[0:5] == "MODEL"):
			nummodels = nummodels + 1
		f.close()
		if (nummodels <= 1):
		    filenames = [filename]
		else:
		    # First tell the user this file has multiple models and ask if they really want to load
		    # all of the structures, or only one
		    dlg = wx.MessageDialog(self, "This PDB contains multiple models.  Do you want to load all models? (Yes loads all models, No loads only the first model)", "Structural Ensemble Detected", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
		    if (dlg.ShowModal() == wx.ID_YES):
			loadAllModels = True
		    else:
			loadAllModels = False
		    dlg.Destroy()
		    goToSandbox()
		    filenames = []
		    f = open(filename, "r")
		    if (platform.system() == "Windows"):
			indx = filename.rfind("\\")
		    else:
			indx = filename.rfind("/")
		    if (indx >= 0):
			modelbasename = filename[indx+1:].split(".pdb")[0]
		    else:
			modelbasename = filename.split(".pdb")[0]
		    writing = False
		    aline = "Begin"
		    i = 0
		    while (aline):
			aline = f.readline()
			if (aline[0:5] == "MODEL"):
			    writing = True
			    i = i + 1
			    f2 = open(modelbasename + "_S" + str(i) + ".pdb", "w")
			elif (aline[0:6] == "ENDMDL"):
			    writing = False
			    f2.close()
			    filenames.append(modelbasename + "_S" + str(i) + ".pdb")
			    # Break if only loading one model
			    if (not(loadAllModels)):
				break
			elif (writing):
			    f2.write(aline.strip() + "\n")
		    f.close()
		first = True
		for filename in filenames:
		    if (first):
			first = False
			self.PyMOLPDBLoad(1, filename, "Show")
		    else:
			self.PyMOLPDBLoad(1, filename, "Reuse")
		#self.PyMOLPDBLoad(1, pdbCode + ".pdb")
		logInfo("Fetched PDB file", pdbCode + ".pdb")
		# Pop this message out of the queue
		for i in range(0, len(self.msgQueue)):
		    if (self.msgQueue[i].find("Loading PDB") >= 0):
			self.msgQueue.pop(i)
			break
		for i in range(0, len(self.msgQueue)):
		    if (self.msgQueue[i].find("Fetching PDB") >= 0):
			self.msgQueue.pop(i)
			break
		if (len(self.msgQueue) > 0):
		    self.labelMsg.SetLabel(self.msgQueue[len(self.msgQueue)-1])
		else:
		    self.labelMsg.SetLabel("")
		self.Enable()
		self.protWin.Enable()
		self.recolorResidues()
		# Some of the protocols recognize changes in the sequence window, so let the protocol panel know something changed
		self.protPanel.activate()
	    else:
		raise Exception("The PDB code " + pdbCode + " could not be found in RCSB!")
	except Exception as e:
	    wx.MessageBox(e.message, "File Cannot Be Read", wx.OK|wx.ICON_EXCLAMATION)
	    try:
		os.remove(pdbCode + ".pdb")
	    except:
		pass
    
    def loadPDBsClick(self, event):
	logInfo("Clicked Load PDB button")
	self.labelMsg.SetLabel("Loading PDB into Rosetta, please be patient...")
	self.msgQueue.append("Loading PDB into Rosetta, please be patient...")
	dlg = wx.FileDialog(
	    self, message="Choose a File",
	    defaultDir=self.cwd,
	    defaultFile="",
	    wildcard="PDB Files (*.pdb)|*.pdb",
	    style=wx.OPEN | wx.CHANGE_DIR | wx.MULTIPLE)
	if (dlg.ShowModal() == wx.ID_OK):
	    paths = dlg.GetPaths()
	    # Change cwd to the last opened file
	    if (platform.system() == "Windows"):
		lastDirIndx = paths[len(paths)-1].rfind("\\")
	    else:
		lastDirIndx = paths[len(paths)-1].rfind("/")
	    self.cwd = str(paths[len(paths)-1][0:lastDirIndx])
	    self.saveWindowData(None)
	    # Load the PDBs into PyMOL
	    for i in range(0, len(paths)):
		filename = str(paths[i])
		filevalid = True
		try:
		    f = open(filename, "r")
		    f.close()
		except:
		    msg = "The file " + filename.strip() + " is invalid."
		    wx.MessageBox(msg, "File Cannot Be Read", wx.OK|wx.ICON_EXCLAMATION)
		    filevalid = False
		if (filevalid):
		    self.Disable()
		    self.protWin.Disable()
		    # We have to check to see if there are multiple models in this PDB
		    # If there are, then we're going to break them up into separate PDB files
		    # and read them all as distinct models in PyMOL
		    nummodels = 0
		    f = open(filename, "r")
		    for aline in f:
			if (aline[0:5] == "MODEL"):
			    nummodels = nummodels + 1
		    f.close()
		    if (nummodels <= 1):
			filenames = [filename]
		    else:
			# First tell the user this file has multiple models and ask if they really want to load
			# all of the structures, or only one
			dlg2 = wx.MessageDialog(self, "This PDB contains multiple models.  Do you want to load all models? (Yes loads all models, No loads only the first model)", "Structural Ensemble Detected", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
			if (dlg2.ShowModal() == wx.ID_YES):
			    loadAllModels = True
			else:
			    loadAllModels = False
			dlg2.Destroy()
			goToSandbox()
			filenames = []
			f = open(filename, "r")
			if (platform.system() == "Windows"):
			    indx = filename.rfind("\\")
			else:
			    indx = filename.rfind("/")
			if (indx >= 0):
			    modelbasename = filename[indx+1:].split(".pdb")[0]
			else:
			    modelbasename = filename.split(".pdb")[0]
			writing = False
			aline = "Begin"
			i = 0
			while (aline):
			    aline = f.readline()
			    if (aline[0:5] == "MODEL"):
				writing = True
				i = i + 1
				f2 = open(modelbasename + "_S" + str(i) + ".pdb", "w")
			    elif (aline[0:6] == "ENDMDL"):
				writing = False
				f2.close()
				filenames.append(modelbasename + "_S" + str(i) + ".pdb")
				# Break if only loading one model
				if (not(loadAllModels)):
				    break
			    elif (writing):
				f2.write(aline.strip() + "\n")
			f.close()
		    first = True
		    for filename in filenames:
			if (first):
			    first = False
			    self.PyMOLPDBLoad(1, filename, "Show")
			else:
			    self.PyMOLPDBLoad(1, filename, "Reuse")
		    logInfo("Loaded PDB file", filename)
		    # Pop this message out of the queue
		    for i in range(0, len(self.msgQueue)):
			if (self.msgQueue[i].find("Loading PDB") >= 0):
			    self.msgQueue.pop(i)
			    break
		    for i in range(0, len(self.msgQueue)):
			if (self.msgQueue[i].find("Fetching PDB") >= 0):
			    self.msgQueue.pop(i)
			    break
		    if (len(self.msgQueue) > 0):
			self.labelMsg.SetLabel(self.msgQueue[len(self.msgQueue)-1])
		    else:
			self.labelMsg.SetLabel("")
		    self.Enable()
		    self.protWin.Enable()
		    self.recolorResidues()
		    # Some of the protocols recognize changes in the sequence window, so let the protocol panel know something changed
		    self.protPanel.activate()
	else:
	    # Pop this message out of the queue
	    for i in range(0, len(self.msgQueue)):
		if (self.msgQueue[i].find("Loading PDB") >= 0):
		    self.msgQueue.pop(i)
		    break
	    if (len(self.msgQueue) > 0):
		self.labelMsg.SetLabel(self.msgQueue[len(self.msgQueue)-1])
	    else:
		self.labelMsg.SetLabel("")
	    logInfo("Cancelled out of Load PDB")
    
    def deleteChain(self, currRow):
	if (self.cannotDelete):
	    wx.MessageBox("You cannot perform deletions while a protocol is active!", "Cannot Delete During Protocol", wx.OK|wx.ICON_EXCLAMATION)
	    return
	logInfo("Deleted row " + self.IDs[currRow] + " in sequence viewer")
	goToSandbox()
	self.SeqViewer.DeleteRows(currRow)
	# Pop the sequences and IDs out that have been deleted
	seq = self.sequences.pop(currRow)
	ID = self.IDs.pop(currRow)
	indexes = self.indxToSeqPos.pop(currRow)
	# Unload the model from PyMOL
	fields = ID.split("|")
	currID = ""
	for i in range(0, len(fields)-1):
	    currID = currID + fields[i] + "|"
	currID = currID[0:len(currID)-1]
	chainID = fields[len(fields)-1]
	# Update the dictionary
	for indx in indexes:
	    if (indx != "-"):
		self.selectLookup.pop((ID, int(indx[1])))
	if (chainID != "_"):
	    self.cmd.remove("model " + currID + " and chain " + chainID)
	else:
	    self.cmd.remove("model " + currID)
	# Update the labels of everything after the deleted row
	for i in range(0, self.SeqViewer.NumberRows):
	    self.SeqViewer.SetRowLabelValue(i, self.IDs[i])
	# Remove the residues from the PyRosetta pose
	# This is a little tricky since the whole pose is only stored in the first chain and there are Nones
	# for the other chains
	# So we have to find this pose, delete all the residues from the deleted chain, and make sure the
	# pose ends up in the right spot in the list
	if (self.poses[currRow] and (currRow == len(self.poses) - 1 or self.poses[currRow+1])):
	    # Easy, only one chain, just delete this one pose and we're done
	    self.poses.pop(currRow) # This is popping a real pose
	    # If this is the last chain, we need to delete the model name from PyMOL, otherwise if you try
	    # to reload this PDB PyMOL will not load it because there will be a namespace clash
	    self.cmd.delete(currID)
	elif (self.poses[currRow]):
	    # We have the structure, but we're deleting the first chain, so we have to shift where the structure
	    # is stored in self.poses
	    s = self.poses.pop(currRow)
	    chain = ID[len(ID)-1]
	    if (chain == "_"):
		chain = " "
	    s[0].detach_child(chain)
	    self.poses[currRow] = s # Place the structure in the spot for the new first chain
	    # IMPORTANT: You have to replace the model in the sandbox with the new model
	    self.pdbwriter.set_structure(s)
	    self.pdbwriter.save(currID + ".pdb")
	else:
	    # We're on a later chain, which means we have to backtrack to find the structure
	    poseloc = currRow - 1
	    while (not(self.poses[poseloc])):
		poseloc = poseloc - 1
	    # Find the true delete window by adding the lengths of the sequences leading up to this chain
	    start = 0
	    for i in range(poseloc, currRow):
		start = start + len(self.sequences[i])
	    end = start + len(seq)
	    chain = ID[len(ID)-1]
	    if (chain == "_"):
		chain = " "
	    self.poses[poseloc][0].detach_child(chain)
	    self.poses.pop(currRow) # This is popping a blank entry, not a real structure
	    # IMPORTANT: You have to replace the model in the sandbox with the new model
	    self.pdbwriter.set_structure(self.poses[poseloc])
	    self.pdbwriter.save(currID + ".pdb")	
    
    def closeClick(self, event):
	# Find out which rows have a selection
	logInfo("Clicked on Close PDB button")
	selectedRows = []
	topLefts = self.SeqViewer.GetSelectionBlockTopLeft()
	bottomRights = self.SeqViewer.GetSelectionBlockBottomRight()
	for i in range(0, len(topLefts)):
	    for j in range(topLefts[i][0], bottomRights[i][0]+1):
		if (not(j in selectedRows)):
		    selectedRows.append(j)
	selectedRows.sort()
	# If no rows were selected, default to closing everything
	if (len(selectedRows) == 0):
	    for i in range(0, self.SeqViewer.NumberRows):
		selectedRows.append(i)
	if (len(selectedRows) == self.SeqViewer.NumberRows):
	    msg = "Are you sure you want to close all chains?"
	else:
	    msg = "Are you sure you want to close the selected chains?"
	dlg = wx.MessageDialog(self, msg, "Close PDBs", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
	if (dlg.ShowModal() == wx.ID_YES):
	    # Take off the rows backwards so the renumbering doesn't screw things up
	    while (len(selectedRows) > 0):
		currRow = selectedRows.pop()
		self.deleteChain(currRow)
	    # Shave off any extra columns
	    maxseqlen = 0
	    for i in range(0, len(self.sequences)):
		if (len(self.sequences[i]) > maxseqlen):
		    maxseqlen = len(self.sequences[i])
	    for c in range(self.SeqViewer.NumberCols, maxseqlen, -1):
		self.SeqViewer.DeleteCols(c-1)
	    self.regenerateLookupTable()
	else:
	    logInfo("Cancelled out of Close PDB")
	dlg.Destroy()
	#self.recolorResidues()
	self.SeqViewer.ClearSelection()
	self.protPanel.activate()
	
    def saveClick(self, event):
	logInfo("Clicked on Save PDB button")
	models = self.getSelectedModels()
	# Default to save all of them if none selected
	if (len(models) == 0):
	    for r in range(0, self.SeqViewer.NumberRows):
		currID = self.getModelForChain(r)
		if (not(currID in models)):
		    models.append(currID)
	# Find the pose for this model
	for model in models:
	    poseindx = self.getPoseIndexForModel(model)
	    dlg = wx.FileDialog(
		self, message="Save a PDB File",
		defaultDir=self.cwd,
		defaultFile=model,
		wildcard="PDB Files (*.pdb)|*.pdb",
		style=wx.SAVE | wx.CHANGE_DIR)
	    if (dlg.ShowModal() == wx.ID_OK):
		paths = dlg.GetPaths()
		# Change cwd to the last opened file
		if (platform.system() == "Windows"):
		    lastDirIndx = paths[len(paths)-1].rfind("\\")
		else:
		    lastDirIndx = paths[len(paths)-1].rfind("/")
		self.cwd = str(paths[len(paths)-1][0:lastDirIndx])
		self.saveWindowData(None)
		# Load the PDBs into PyMOL
		filename = str(paths[0]).split(".pdb")[0] + ".pdb"
		# Does it exist already?  If so, ask if the user really wants to overwrite it
		if (os.path.isfile(filename)):
		    dlg2 = wx.MessageDialog(self, "The file " + filename + " already exists.  Overwrite it?", "Filename Already Exists", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
		    if (dlg2.ShowModal() == wx.ID_NO):
			dlg2.Destroy()
			logInfo("Cancelled save operation due to filename already existing")
			continue
		    dlg2.Destroy()    
		goToSandbox()
		# These shenanigans are necessary because Rosetta cannot dump PDBs to filenames with
		# whitespace in them apparently
		logInfo("Saved a PDB to " + filename.strip())
		#self.pdbwriter.set_structure(self.poses[poseindx])
		#self.pdbwriter.save(filename)
		# Use PyMOL to save it because the coordinates could have changed using the manipulation panel
		# and these changes will not be reflected in BioPython
		self.cmd.save(filename.strip(), "model " + model)
		# PyMOL can do weird things with adding TER lines around NCAAs, which crashes Rosetta
		# Take those out
		fixPyMOLSave(filename.strip())
		# Update the IDs of this model in both the sequence window and in PyMOL
		newmodel = str(paths[len(paths)-1][lastDirIndx+1:]).split(".pdb")[0]
		self.cmd.set_name(model, newmodel)
		# Now we have to save it in the sandbox otherwise the protocols will be looking for a file with the new
		# name and will not find it
		self.cmd.save(newmodel + ".pdb", "model " + newmodel)
		for i in range(0, len(self.IDs)):
		    if (self.IDs[i][0:len(self.IDs[i])-2] == model):
			self.IDs[i] = newmodel + "|" + self.IDs[i][len(self.IDs[i])-1]
			self.SeqViewer.SetRowLabelValue(i, self.IDs[i])
		self.regenerateLookupTable() # Have to do this because the old dictionary is based off the old IDs
	    else:
		logInfo("Cancelled out of Save PDB")
	    dlg.Destroy()
	
    def saveImage(self, event):
	logInfo("Clicked on Save Image button")
	dlg = wx.FileDialog(
	    self, message="Save an Image",
	    defaultDir=self.cwd,
	    defaultFile="",
	    wildcard="PNG Images (*.png)|*.png",
	    style=wx.SAVE | wx.CHANGE_DIR)
	if (dlg.ShowModal() == wx.ID_OK):
	    filename = str(dlg.GetPaths()[0])
	    if (platform.system() == "Windows"):
		# Have to double the number of "\" otherwise the filename can get garbled when
		# PyMOL interprets things like \t as tabs...
		for i in range(len(filename)-1, -1, -1):
		    if (filename[i] == "\\"):
			filename = filename[0:i] + "\\" + filename[i:]
	    self.cmd.disable("seqsele")
	    self.cmd.disable("sele")
	    # Does it exist already?  If so, ask if the user really wants to overwrite it
	    if (os.path.isfile(filename)):
		dlg2 = wx.MessageDialog(self, "The file " + filename + " already exists.  Overwrite it?", "Filename Already Exists", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
		if (dlg2.ShowModal() == wx.ID_NO):
		    dlg2.Destroy()
		    logInfo("Canceled save operation due to filename already existing")
		dlg2.Destroy()
	    self.cmd.png(filename)
	    self.cmd.enable("seqsele")
	    self.cmd.enable("sele")
	    logInfo("Saved an image to " + filename)
	else:
	    logInfo("Cancelled out of Save Image")
	
    def PyMOLPDBLoad(self, dummy, pdbfile, showDialog="NoShow"):
	pdbfile = str(pdbfile)
	# Check the PDB for duplicate atoms
	# We want to rename dupliates otherwise BioPython will drop them
	cleanPDB(pdbfile)
	# Set the ID of this PDB as the filename without the .pdb extension
	if (platform.system() == "Windows"):
	    startindx = pdbfile.rfind("\\") + 1
	else:
	    startindx = pdbfile.rfind("/") + 1
	endindx = pdbfile.rfind(".pdb")
	data = []
	# This is also really annoying...PyRosetta apparently cannot handle file paths with whitespace, so we have
	# to create a directory for InteractiveROSETTA in the user's home directory and cd there
	# All of the uploaded files will be stored here, so we can use relative filenames and avoid whitespace,
	# unless there's whitespace in the PDB name itself, in which case the user can be told to rename it
	goToSandbox()
	newID = pdbfile[startindx:endindx]
	# Special characters in the file name will mess up PyMOL, so take them out
	for i in range(0, len(newID)):
	    if (newID[i] == "@"):
		newID = newID[0:i] + "A" + newID[i+1:]
	    elif (newID[i] == "$"):
		newID = newID[0:i] + "S" + newID[i+1:]
	    elif (newID[i] == "&"):
		newID = newID[0:i] + "A" + newID[i+1:]
	    elif (newID[i] == "*"):
		newID = newID[0:i] + "o" + newID[i+1:]
	    elif (newID[i] == "("):
		newID = newID[0:i] + "_" + newID[i+1:]
	    elif (newID[i] == ")"):
		newID = newID[0:i] + "_" + newID[i+1:]
	pdbID = newID
	# The PyMOL ID is normally just the PDB filename without the ".pdb"
	# We have to check to make sure that the same PDB has not been loaded twice because if it is an extra
	# copy, the ID needs to be changed
	taken = False
	for ID in self.IDs:
	    # This is fancy stuff for parsing the actual model name out of the ID, since the ID is the model
	    # + "|" + the chainID
	    fields = ID.split("|")
	    currID = ""
	    for i in range(0, len(fields)-1):
		currID = currID + fields[i] + "|"
	    currID = currID[0:len(currID)-1]
	    if (currID == newID):
		taken = True
		break
	if (taken):
	    for i in range(2, 9999):
		taken = False
		for ID in self.IDs:
		    fields = ID.split("|")
		    currID = ""
		    for j in range(0, len(fields)-1):
			currID = currID + fields[j] + "|"
		    currID = currID[0:len(currID)-1]
		    if (currID == newID + "_" + str(i)):
			taken = True
			break
		if (not(taken)):
		    newID = newID + "_" + str(i)
		    break
	chainIDs = []
	prev_res = "0000"
	f = open(pdbfile, "r")
	for aline in f:
	    data.append(aline)
	f.close
	if (platform.system() == "Windows"):
	    f = open(os.path.expanduser("~") + "\\InteractiveROSETTA\\" + newID + ".pdb", "w")
	else:
	    f = open(os.path.expanduser("~") + "/InteractiveROSETTA/" + newID + ".pdb", "w")
	offset = 0
	for aline in data:
	    #if (aline[0:3] != "TER"):
	    # Take out the unrecognized residues
	    if ((aline[0:4] == "ATOM" or aline[0:6] == "HETATM") and not(aline[17:20].strip() in getRecognizedTypes())):
		offset = offset + 1
	    elif (aline[0:4] == "ATOM" or aline[0:6] == "HETATM" or aline[0:3] == "TER"):
		try:
		    atomno = int(aline[7:11])
		    atomno = atomno - offset
		    aline = aline[0:7] + ("%4i" % atomno) + aline[11:]
		except:
		    pass
		f.write(aline)
	    else:
		f.write(aline)
	f.close()
	# Separate out the chain sequences, but keep the whole pose intact for later Rosetta operations
	try:
	    biopdb = self.pdbreader.get_structure(newID, newID + ".pdb")
	except:
	    # Tell the user the file failed to load and undo the additions to self.indxToSeqPos
	    msg = "There was an error loading " + pdbfile + "\n\n"
	    msg = msg + "Things to check for:\n"
	    msg = msg + "Is the file a valid PDB file?\n"
	    msg = msg + "Are some residues missing a large number of atoms?"
	    dlg = wx.MessageDialog(self, msg, "Close PDBs", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
	    if (dlg.ShowModal() == wx.ID_OK):
		self.labelMsg.SetLabel("")
		self.Enable()
		self.protWin.Enable()
		return
	doNotLoad = False
	# Before showing the chain selection dialog, let's first see if there's either multiple chains
	# or there are HETATMs present, otherwise there's no need to show the dialog
	if (len(biopdb[0]) > 1):
	    # Multiple chains, show it
	    pass
	else:
	    # Are there HETATMs?
	    foundHETATMs = False
	    for ch in biopdb[0]:
		for residue in ch:
		    if (not(residue.resname in "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR ")):
			foundHETATMs = True
			break
		if (foundHETATMs):
		    break
	    if (not(foundHETATMs)):
		showDialog = "NoShow"
	if (showDialog != "NoShow"):
	    # Get the chains, then display a dialog that will let the user uncheck chains they do not want to load
	    # and possibly opt to not load waters/hetatms
	    chains = []
	    for ch in biopdb[0]:
		if (ch.id == " "):
		    chains.append("_")
		else:
		    chains.append(ch.id)
	    if (showDialog == "Show"): # Otherwise this is an NMR ensemble and the options from the first model should be applied to all
		dlg = ProteinDialog(self, pdbID + ".pdb", chains)
		if (dlg.ShowModal() != wx.OK):
		    # Canceled or closed, do nothing
		    doNotLoad = True
		# Now check the status of those checkboxes
		self.chainsToDelete = []
		for i in range(0, len(chains)):
		    if (not(dlg.chkChains[i].GetValue())):
			self.chainsToDelete.append(chains[i])
		self.keepWaters = dlg.btnWater.GetLabel() == "Load Waters"
		self.keepHETATMs = dlg.btnHETATM.GetLabel() == "Load HETATMs"
		dlg.Destroy()
	    if (len(chains) == len(self.chainsToDelete)):
		# Load nothing
		doNotLoad = True
	    if (not(doNotLoad)):
		# Get rid of the indicated chains
		for chain in self.chainsToDelete:
		    if (chain == "_"):
			biopdb[0].detach_child(" ")
		    else:
			biopdb[0].detach_child(chain)
		# Now take out HETATMs/waters if requested
		toRemove = []
		for ch in biopdb[0]:
		    for residue in ch:
			if (not(self.keepWaters) and residue.resname == "HOH"):
			    toRemove.append((ch.id, residue.id))
			elif (not(self.keepHETATMs) and not(residue.resname in "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR ")):
			    toRemove.append((ch.id, residue.id))
		for (ch, res) in toRemove:
		    biopdb[0][ch].detach_child(res)
		# Now we have to look for empty chains incase there were HETATM only chains and remove them
		for ch in biopdb[0]:
		    length = 0
		    for residue in ch:
			length = length + 1
		    if (length == 0):
			biopdb[0].detach_child(ch.id)
		# Save it for PyMOL
		self.pdbwriter.set_structure(biopdb)
		self.pdbwriter.save(newID + ".pdb")
	if (not(doNotLoad)):
	    self.poses.append(biopdb)
	    self.cmd.load(newID + ".pdb", newID)
	    # Use DSSP if available to characterize secondary structure
	    if (self.dsspexe != "N/A"):
		try:
		    dsspdata = Bio.PDB.DSSP(biopdb[0], newID + ".pdb", dssp=self.dsspexe)
		except:
		    self.dsspexe = "N/A"
		    print "WARNING: DSSP could not be run.  Is its binary executable?"
	    first = True
	    i = 0
	    # Things get messed up if loading a PDB when in PDB# or Align mode, so turn these off
	    # Turn off Align for good and require the user to use superimposition again to get it back
	    # For PDB#, turn it back on after loading this whole structure to get the right alignment
	    if (self.AlignBtn.GetLabel() == "Align"):
		retoggle = True
	    elif (self.AlignBtn.GetLabel() == "PDB #"):
		retoggle = True
		self.AlignBtn.SetLabel("From 1")
	    else:
		retoggle = False
	    for c in biopdb[0]:
		if (not(first)):
		    self.poses.append(None) # Just to keep all the lists the same size
		first = False
		chain = c.id
		if (len(chain.strip()) == 0):
		    chain = "_"
		seq = ""
		for r in c:
		    seq = seq + AA3to1(r.resname)
		self.sequences.append(seq)
		self.IDs.append(newID + "|" + str(chain))
		self.indxToSeqPos.append([])
		ires = 0
		for r in c:
		    self.indxToSeqPos[len(self.indxToSeqPos)-1].append(r.id) # This is actually a 3-tuple, the seqpos is element 1
		    self.selectLookup[(newID + "|" + chain, r.id[1])] = (len(self.sequences)-1, ires)
		    if (self.dsspexe != "N/A"):
			# Alter the states in PyMOL with the DSSP secondary structure predictions
			try:
			    ss = dsspdata[(c.id, r.id[1])][1]
			    # Not all the DSSP codes match up with PyMOL, so make them match
			    if (ss == "E"):
				# Beta sheet
				ss = "S"
			    elif (ss == "S"):
				# Bend
				ss = "B"
			    elif (ss == "-"):
				# Loop
				ss = "L"
			    if (c.id != " "):
				self.cmd.alter("model " + newID + " and chain " + c.id + " and resi " + str(r.id[1]), "ss=\"" + ss + "\"")
			    else:
				self.cmd.alter("model " + newID + " and resi " + str(r.id[1]), "ss=\"" + ss + "\"")
			except:
			    ss = "S"
		    ires = ires + 1
		i = i + 1
		self.updateSeqViewer()
	    if (retoggle):
		self.realignToggle(None)
	    # Have PyMOL load from this duplicated PDB file so when protocols are run on the PDB we can easily dump
	    # the PDB from Rosetta and then reload it in PyMOL without touching the original PDB unless the user
	    # explicitly tells us to save the new PDB
	    if (self.dsspexe != "N/A"):
		self.cmd.rebuild()
	    defaultPyMOLView(self.cmd, newID)
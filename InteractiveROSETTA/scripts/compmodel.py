import wx
import wx.grid
import wx.lib.scrolledpanel
import os
import os.path
import time
import platform
import multiprocessing
import datetime
import requests
import webbrowser
from threading import Thread
from tools import *

class CompModelPanel(wx.lib.scrolledpanel.ScrolledPanel):
    def __init__(self, parent, W, H):
	#if (platform.system() == "Windows"):
	wx.lib.scrolledpanel.ScrolledPanel.__init__(self, parent, id=-1, pos=(10, 60), size=(340, H-330), name="ProtFixbb")
	winh = H-330
	#else:
	    #wx.lib.scrolledpanel.ScrolledPanel.__init__(self, parent, id=-1, pos=(10, 60), size=(340, H-330), name="ProtMinimization")
	    #winh = H-290
	self.SetBackgroundColour("#333333")
	self.parent = parent
	
	if (platform.system() == "Windows"):
	    self.lblProt = wx.StaticText(self, -1, "Comparative Modeling", (25, 15), (270, 25), wx.ALIGN_CENTRE)
	    self.lblProt.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
	elif (platform.system() == "Darwin"):
	    self.lblProt = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/lblCompModel.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(25, 15), size=(270, 25))
	else:
	    self.lblProt = wx.StaticText(self, -1, "Comparative Modeling", (70, 15), style=wx.ALIGN_CENTRE)
	    self.lblProt.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
	    resizeTextControlForUNIX(self.lblProt, 0, self.GetSize()[0]-20)
	self.lblProt.SetForegroundColour("#FFFFFF")
	
	if (platform.system() == "Darwin"):
	    self.HelpBtn = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/HelpBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(295, 10), size=(25, 25))
	else:
	    self.HelpBtn = wx.Button(self, id=-1, label="?", pos=(295, 10), size=(25, 25))
	    self.HelpBtn.SetForegroundColour("#0000FF")
	    self.HelpBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.HelpBtn.Bind(wx.EVT_BUTTON, self.showHelp)
	self.HelpBtn.SetToolTipString("Display the help file for this window")
	
	if (platform.system() == "Windows"):
	    self.lblInst = wx.StaticText(self, -1, "Model a structure ab initio using an existing\nstructure as a template.\nFirst generate fragment files using Robetta.\nRegistration with Robetta is required to do this.", (0, 45), (320, 25), wx.ALIGN_CENTRE)
	    self.lblInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
	elif (platform.system() == "Darwin"):
	    self.lblInst = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/lblInstCompModel.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 45), size=(320, 80))
	else:
	    self.lblInst = wx.StaticText(self, -1, "Model a structure ab initio using an existing\nstructure as a template.\nFirst generate fragment files using Robetta.\nRegistration with Robetta is required to do this.", (5, 45), style=wx.ALIGN_CENTRE)
	    self.lblInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
	    resizeTextControlForUNIX(self.lblInst, 0, self.GetSize()[0]-20)
	self.lblInst.SetForegroundColour("#FFFFFF")
	
	if (platform.system() == "Windows"):
	    self.lblVisit = wx.StaticText(self, -1, "Visit Robetta", (0, 130), (155, 20), wx.ALIGN_CENTRE)
	    self.lblVisit.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	elif (platform.system() == "Darwin"):
	    self.lblVisit = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/lblVisit.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 130), size=(155, 20))
	else:
	    self.lblVisit = wx.StaticText(self, -1, "Visit Robetta", (0, 130), style=wx.ALIGN_CENTRE)
	    self.lblVisit.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	    resizeTextControlForUNIX(self.lblVisit, 0, 155)
	self.lblVisit.SetForegroundColour("#FFFFFF")
	if (platform.system() == "Darwin"):
	    self.btnFragGen = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/btnFragments.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(175, 125), size=(120, 25))
	else:
	    self.btnFragGen = wx.Button(self, id=-1, label="Get Fragments", pos=(175, 125), size=(120, 25))
	    self.btnFragGen.SetForegroundColour("#000000")
	    self.btnFragGen.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
	self.btnFragGen.Bind(wx.EVT_BUTTON, self.fragGen)
	self.btnFragGen.SetToolTipString("Visit Robetta to generate fragment files")
	
	if (platform.system() == "Windows"):
	    self.lblJobID = wx.StaticText(self, -1, "Job ID:", (0, 163), (100, 20), wx.ALIGN_CENTRE)
	    self.lblJobID.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	elif (platform.system() == "Darwin"):
	    self.lblJobID = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/lblJobID.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 163), size=(100, 20))
	else:
	    self.lblJobID = wx.StaticText(self, -1, "Job ID:", (0, 163), style=wx.ALIGN_CENTRE)
	    self.lblJobID.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	    resizeTextControlForUNIX(self.lblJobID, 0, 100)
	self.lblJobID.SetForegroundColour("#FFFFFF")
	self.txtJobID = wx.TextCtrl(self, -1, pos=(110, 160), size=(100, 25))
	self.txtJobID.SetValue("")
	self.txtJobID.SetToolTipString("Job ID for submitted Robetta fragment job.")
	if (platform.system() == "Darwin"):
	    self.btnWatch = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/btnWatch.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(220, 160), size=(100, 25))
	else:
	    self.btnWatch = wx.Button(self, id=-1, label="Watch", pos=(220, 160), size=(100, 25))
	    self.btnWatch.SetForegroundColour("#000000")
	    self.btnWatch.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
	self.btnWatch.Bind(wx.EVT_BUTTON, self.jobWatch)
	self.btnWatch.SetToolTipString("Start automatically watching for the specified job ID to complete")
	
	if (platform.system() == "Windows"):
	    self.lblModel = wx.StaticText(self, -1, "Templates", (0, 190), (155, 20), wx.ALIGN_CENTRE)
	    self.lblModel.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	elif (platform.system() == "Darwin"):
	    self.lblModel = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/lblModelCompModel.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 190), size=(155, 20))
	else:
	    self.lblModel = wx.StaticText(self, -1, "Templates", (0, 190), style=wx.ALIGN_CENTRE)
	    self.lblModel.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	    resizeTextControlForUNIX(self.lblModel, 0, 155)
	self.lblModel.SetForegroundColour("#FFFFFF")
	self.modelMenu = wx.ComboBox(self, pos=(0, 210), size=(155, 25), choices=[], style=wx.CB_READONLY)
	self.modelMenu.Bind(wx.EVT_COMBOBOX, self.modelMenuSelect)
	self.modelMenu.SetToolTipString("Select the models to serve as templates for structure prediction")
	self.templates = []
	
	if (platform.system() == "Darwin"):
	    self.btnAddTemplate = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/btnAddFragments.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(165, 210), size=(77, 25))
	else:
	    self.btnAddTemplate = wx.Button(self, id=-1, label="Add", pos=(165, 210), size=(77, 25))
	    self.btnAddTemplate.SetForegroundColour("#000000")
	    self.btnAddTemplate.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.btnAddTemplate.Bind(wx.EVT_BUTTON, self.addModel)
	self.btnAddTemplate.SetToolTipString("Add model to the list of templates")
	if (platform.system() == "Darwin"):
	    self.btnRemoveTemplate = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/btnRemoveFragments.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(242, 210), size=(78, 25))
	else:
	    self.btnRemoveTemplate = wx.Button(self, id=-1, label="Remove", pos=(242, 210), size=(78, 25))
	    self.btnRemoveTemplate.SetForegroundColour("#000000")
	    self.btnRemoveTemplate.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.btnRemoveTemplate.Bind(wx.EVT_BUTTON, self.removeModel)
	self.btnRemoveTemplate.SetToolTipString("Remove model from the list of templates")
	
	self.grdTemplates = wx.grid.Grid(self)
	self.grdTemplates.CreateGrid(0, 1)
	self.grdTemplates.SetSize((320, 150))
	self.grdTemplates.SetPosition((0, 240))
	self.grdTemplates.SetLabelFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.grdTemplates.DisableDragColSize()
	self.grdTemplates.DisableDragRowSize()
	self.grdTemplates.SetColLabelValue(0, "Templates")
	self.grdTemplates.SetRowLabelSize(160)
	self.grdTemplates.SetColSize(0, 160)
	
	self.lblFrag = wx.StaticText(self, -1, "None Uploaded", pos=(10, 403), size=(180, 25), style=wx.ALIGN_CENTRE)
	self.lblFrag.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	if (platform.system() == "Linux"):
	    resizeTextControlForUNIX(self.lblFrag, 10, 180)
	self.lblFrag.SetForegroundColour("#FFFFFF")
	if (platform.system() == "Darwin"):
	    self.btnLoad = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/btnLoadFrag.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(200, 400), size=(110, 25))
	else:
	    self.btnLoad = wx.Button(self, id=-1, label="Load Frag", pos=(200, 400), size=(110, 25))
	    self.btnLoad.SetForegroundColour("#000000")
	    self.btnLoad.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.btnLoad.Bind(wx.EVT_BUTTON, self.loadFrag)
	self.btnLoad.SetToolTipString("Load a frag GZ file")
	self.loadedfile = ""
	
	if (platform.system() == "Darwin"):
	    self.scoretypeMenu = wx.ComboBox(self, pos=(7, 430), size=(305, 25), choices=[], style=wx.CB_READONLY)
	else:
	    self.scoretypeMenu = wx.ComboBox(self, pos=(7, 430), size=(305, 25), choices=[], style=wx.CB_READONLY | wx.CB_SORT)
	self.scoretypeMenu.Bind(wx.EVT_COMBOBOX, self.scoretypeMenuSelect)
	self.scoretypeMenu.Disable() # Is only enabled after a design and before accepting it
	self.scoretypeMenu.SetToolTipString("Scoretype by which PyMOL residues will be colored")
	
	if (platform.system() == "Windows"):
	    self.lblModelView = wx.StaticText(self, -1, "View Structure:", (20, 463), (120, 20), wx.ALIGN_CENTRE)
	    self.lblModelView.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	elif (platform.system() == "Darwin"):
	    self.lblModelView = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/lblModelView.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(20, 463), size=(120, 20))
	else:
	    self.lblModelView = wx.StaticText(self, -1, "View Structure:", (20, 463), style=wx.ALIGN_CENTRE)
	    self.lblModelView.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	    resizeTextControlForUNIX(self.lblModelView, 20, 120)
	self.lblModelView.SetForegroundColour("#FFFFFF")
	self.viewMenu = wx.ComboBox(self, pos=(175, 460), size=(120, 25), choices=[], style=wx.CB_READONLY)
	self.viewMenu.Bind(wx.EVT_COMBOBOX, self.viewMenuSelect)
	self.viewMenu.Disable()
	self.viewMenu.SetToolTipString("Select model positions to view in PyMOL")
	
	if (platform.system() == "Darwin"):
	    self.btnThread = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/btnThread.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(110, 490), size=(100, 25))
	else:
	    self.btnThread = wx.Button(self, id=-1, label="Thread!", pos=(110, 490), size=(100, 25))
	    self.btnThread.SetForegroundColour("#000000")
	    self.btnThread.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
	self.btnThread.Bind(wx.EVT_BUTTON, self.threadClick)
	self.btnThread.SetToolTipString("Begin threading the sequence onto the template structures")
	self.buttonState = "Thread!"
	
	self.scrollh = self.btnThread.GetPosition()[1] + self.btnThread.GetSize()[1] + 5
	self.SetScrollbars(1, 1, 320, self.scrollh)
    
    def showHelp(self, event):
	# Open the help page
	if (platform.system() == "Darwin"):
	    try:
		browser = webbrowser.get("Safari")
	    except:
		print "Could not load Safari!  The help files are located at " + self.scriptdir + "/help"
		return
	    browser.open(self.parent.parent.scriptdir + "/help/comparative.html")
	else:
	    webbrowser.open(self.parent.parent.scriptdir + "/help/comparative.html")
    
    def setSeqWin(self, seqWin):
	self.seqWin = seqWin
	# So the sequence window knows about what model "designed_view" really is
	self.seqWin.setProtocolPanel(self)
	
    def setPyMOL(self, pymol):
	self.pymol = pymol
	self.cmd = pymol.cmd
	self.stored = pymol.stored
	
    def setSelectWin(self, selectWin):
	self.selectWin = selectWin
	self.selectWin.setProtPanel(self)
	
    def activate(self):
	if (self.seqWin):		
	    # Find out which models are selected and the locations of their poses
	    selection = self.seqWin.getSelectedResidues()
	    rowsSelected = []
	    for line in selection:
		rowsSelected.append(line[0])
	    if (len(rowsSelected) == 0):
		# If nothing selected, select all models
		if (self.seqWin.numChains() > 0):
		    rowsSelected = range(0, self.seqWin.numChains())
		else:
		    return
	    self.models = []
	    for r in rowsSelected:
		ID = self.seqWin.IDs[r]
		model = self.seqWin.getModelForChain(r)
		if (not(model in self.models)):
		    self.models.append(str(model))
	    # Update the combo box to have these models as selections
	    self.modelMenu.Clear()
	    self.modelMenu.AppendItems(self.models)
	    if (len(self.models) > 0):
		if (platform.system() == "Windows"):
		    # Doing this on Linux freezes up the ComboBox, so don't do it on Linux
		    self.modelMenu.SetSelection(0)
    
    def uploadFASTA(self, event):
	if (len(self.txtFASTA.GetValue().strip()) > 0):
	    dlg = wx.MessageDialog(self, "Your current FASTA sequence will be replaced with the loaded data.  Do you want to continue?", "Replace FASTA Data", wx.YES_NO | wx.ICON_EXCLAMATION | wx.CENTRE)
	    if (dlg.ShowModal() == wx.ID_YES):
		# Don't clear it out until the user actually selects a resfile (instead of cancelling
		# at the file selection dialog
		pass
	    else:
		logInfo("Upload FASTA operation canceled due to data already being present in the text box.")
		return
	dlg = wx.FileDialog(
	    self, message="Choose a File",
	    defaultDir=self.seqWin.cwd,
	    defaultFile="",
	    wildcard="All Files (*)|*",
	    style=wx.OPEN | wx.CHANGE_DIR)
	if (dlg.ShowModal() == wx.ID_OK):
	    paths = dlg.GetPaths()
	    # Change cwd to the last opened file
	    lastDirIndx = paths[len(paths)-1].rfind("/")
	    self.seqWin.cwd = str(paths[len(paths)-1][0:lastDirIndx])
	    self.seqWin.saveWindowData(None)
	    filename = str(paths[0])
	    # Does it open?  If yes, then erase the resfile data and continue
	    try:
		f = open(filename, "r")
	    except:
		wx.MessageBox("The file " + filename.strip() + " cannot be opened!", "File Cannot Be Read", wx.OK|wx.ICON_EXCLAMATION)
		return
	    logInfo("Loaded data from a FASTA file", filename)
	    fastastr = ""
	    for aline in f:
		fastastr = fastastr + aline.strip() + "\n"
	    self.txtFASTA.SetValue(fastastr)
	    f.close()
    
    def fragGen(self, event):
	webbrowser.open("http://robetta.bakerlab.org/fragmentsubmit.jsp")
	dlg = wx.MessageDialog(self, "After submitting a job, please make note of the job ID integer and specify it below so InteractiveROSETTA can watch for it to be completed.", "Specify Job ID", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
	dlg.ShowModal()
	dlg.Destroy()
	
    def jobWatch(self, event):
	# First ask if the ID is valid (it has to be an integer, no letters)
	if (len(self.txtJobID.GetValue().strip()) == 0):
	    dlg = wx.MessageDialog(self, "Please specify a job ID!", "Job ID Required", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
	    dlg.ShowModal()
	    dlg.Destroy()
	    return
	try:
	    int(self.txtJobID.GetValue())
	except:
	    dlg = wx.MessageDialog(self, "Please specify a valid job ID!  Valid Job IDs are integers.", "Job ID Invalid", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
	    dlg.ShowModal()
	    dlg.Destroy()
	    return
	# Now write this to a text file
	# A timer on the SequenceWindow will use this information to look for the files to appear on the Robetta server
	goToSandbox()
	# First make sure this isn't a duplicate
	alreadythere = False
	try:
	    f = open("downloadwatch", "r")
	    for aline in f:
		if (len(aline.split("\t")) >= 2 and aline.split("\t")[0] == "FRAGMENT" and aline.split("\t")[1] == str(txtJobID.GetValue().strip())):
		    alreadythere = True
		    break
	    f.close()
	except:
	    pass
	if (not(alreadythere)):
	    f = open("downloadwatch", "a")
	    f.write("FRAGMENT\t" + str(self.txtJobID.GetValue()) + "\t" + str(datetime.datetime.now()) + "\n")
	    f.close()
	dlg = wx.MessageDialog(self, "InteractiveROSETTA is now watching the Robetta server for job ID " + str(self.txtJobID.GetValue()) + ".  You will be notified when the package is available for download.", "Listening for Download", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
	dlg.ShowModal()
	dlg.Destroy()
	
    def showModel(self, model):
	defaultPyMOLView(self.cmd)
	self.cmd.select("sele", "model " + model)
	self.cmd.zoom("sele")
	self.cmd.hide("everything", "not sele")
	self.seqWin.selectUpdate(False)
    
    def modelMenuSelect(self, event):
	logInfo("Selected model " + self.modelMenu.GetValue())
	self.showModel(self.modelMenu.GetStringSelection())
    
    def updateGrid(self):
	if (self.grdTemplates.NumberRows > 0):
	    self.grdTemplates.DeleteRows(0, self.grdTemplates.NumberRows)
	for i in range(0, len(self.templates)):
	    self.grdTemplates.AppendRows(1)
	    readOnly = wx.grid.GridCellAttr()
	    readOnly.SetReadOnly(True)
	    self.grdTemplates.SetRowAttr(i, readOnly)
	    self.grdTemplates.SetCellAlignment(i, 0, wx.ALIGN_CENTRE, wx.ALIGN_CENTRE)
	    self.grdTemplates.SetRowLabelValue(i, "")
	    self.grdTemplates.SetCellValue(i, 0, self.templates[i])
    
    def addModel(self, event):
	# Add this to the list of template structures
	model = str(self.modelMenu.GetValue())
	if (not(model in self.templates)):
	    self.templates.append(model)
	self.templates.sort()
	# Update the grid with this new information
	self.updateGrid()
	logInfo("Added " + model + " to the list of model structures")
	
    def removeModel(self, event):
	# Just pop out the indicated chain
	model = self.modelMenu.GetValue()
	if (model in self.templates):
	    indx = self.templates.index(model)
	    self.templates.pop(indx)
	# Update the grid with this new information
	self.updateGrid()
	logInfo("Removed " + model + " from the list of static chains")
	
    def loadFrag(self, event):
	# Get the structure from a MOL2 file and load it into PyMOL
	logInfo("Load Frag button clicked")
	dlg = wx.FileDialog(
	    self, message="Choose a File",
	    defaultDir=self.seqWin.cwd,
	    defaultFile="",
	    wildcard="Frag Files (*.frag)|*.frag",
	    style=wx.OPEN | wx.CHANGE_DIR)
	if (dlg.ShowModal() == wx.ID_OK):
	    paths = dlg.GetPaths()
	    # Change cwd to the last opened file
	    if (platform.system() == "Windows"):
		lastDirIndx = paths[len(paths)-1].rfind("\\")
	    else:
		lastDirIndx = paths[len(paths)-1].rfind("/")
	    self.seqWin.cwd = str(paths[len(paths)-1][0:lastDirIndx])
	    self.seqWin.saveWindowData(None)
	    filename = str(paths[0])
	    self.loadedfile = filename
	    localfilename = filename[lastDirIndx+1:]
	    self.lblFrag.SetLabel(localfilename)
	    self.lblFrag.SetForegroundColour("#FFFFFF")
	else:
	    logInfo("Load PDB operation cancelled")
	
    def viewMenuSelect(self, event):
	try:
	    self.focusView(self.viewMenu.GetStringSelection(), self.selectedModel, "thread_view")
	    logInfo("Viewing " + self.viewMenu.GetStringSelection())
	except:
	    # Probably the user left the field blank, do nothing
	    pass
    
    def focusView(self, posID, origmodel, newmodel=None):
	if (posID != "Whole View"):
	    chain = posID[0]
	    seqpos = posID[3:].strip()
	else:
	    self.cmd.select("viewsele", "model " + origmodel)
	    self.cmd.hide("all")
	    self.cmd.show("cartoon", "viewsele")
	    self.cmd.zoom("viewsele")
	    return
	# Find the neighborhood view
	if (newmodel):
	    firstmodel = newmodel
	else:
	    firstmodel = origmodel
	self.cmd.hide("all")
	if (chain == " " or chain == "_"):
	    self.cmd.select("viewsele", "resi " + seqpos + " and model " + firstmodel)
	else:
	    self.cmd.select("viewsele", "resi " + seqpos + " and model " + firstmodel + " and chain " + chain)
	self.cmd.select("exviewsele", "model " + firstmodel + " within 12 of viewsele")
	self.cmd.show("cartoon", "exviewsele")
	self.cmd.hide("ribbon", "exviewsele")
	self.cmd.show("sticks", "exviewsele")
	self.cmd.set_bond("stick_radius", 0.1, "exviewsele")
	# Display energy labels for new structures
	if (newmodel):
	    relabelEnergies(self.threadView, self.residue_E, newmodel, self.scoretypeMenu.GetStringSelection(), self.cmd, seqpos)
	    self.cmd.label("not exviewsele", "")
	self.cmd.zoom("exviewsele")
	self.cmd.show("sticks", "viewsele")
	self.cmd.set_bond("stick_radius", 0.25, "viewsele")
	# Highlight this residue in PyMOL
	self.cmd.select("sele", "viewsele")
	self.cmd.enable("sele")
	self.cmd.delete("viewsele")
	self.cmd.delete("exviewsele")
    
    def scoretypeMenuSelect(self, event):
	# Make sure there is even a PyMOL_Mover pose loaded
	if (self.selectedModel == ""):
	    return
	logInfo("Changed scoretype view to " + self.scoretypeMenu.GetStringSelection())
	recolorEnergies(self.threadView, self.residue_E, "thread_view", self.scoretypeMenu.GetStringSelection(), self.cmd)
	self.viewMenuSelect(event) # To update all the labels
    
    def threadClick(self, event):
	# This is also the "Finalize!" button
	if (self.buttonState == "Thread!"):
	    # Make sure there is at least one template structure specified and that a frag gz file is specified
	    if (len(self.templates) == 0):
		wx.MessageBox("Please select at least one loaded model as a template for threading!", "Template Required", wx.OK|wx.ICON_EXCLAMATION)
		return
	    if (len(self.loadedfile.strip()) == 0):
		wx.MessageBox("Please indicate a FRAG file generated from the FASTA sequence that will be modeled.  Use Robetta to generate this if you have not done so already.", "Movable Chain Required", wx.OK|wx.ICON_EXCLAMATION)
		return
	    self.seqWin.labelMsg.SetLabel("Performing protein threading, please be patient...")
	    self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
	    self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")
	    self.seqWin.msgQueue.append("Performing protein threading, please be patient...")
	    self.seqWin.cannotDelete = True
	    self.parent.GoBtn.Disable()
	    self.btnAddTemplate.Disable()
	    self.btnRemoveTemplate.Disable()
	    self.btnLoad.Disable()
	    self.btnThread.Disable()
	    self.stage = 1
	    logInfo("Clicked the Thread button")
	    self.tmrThread = wx.Timer(self)
	    self.Bind(wx.EVT_TIMER, self.threadThread, self.tmrThread)
	    self.tmrThread.Start(1000)
	else:
	    # Finalize button, ask whether the changes will be accepted or rejected
	    dlg = wx.MessageDialog(self, "Do you want to accept the threaded model?", "Accept/Reject Model", wx.YES_NO | wx.CANCEL | wx.ICON_QUESTION | wx.CENTRE)
	    result = dlg.ShowModal()
	    if (result == wx.ID_YES):
		logInfo("Accepted threaded model")
		accept = True
	    elif (result == wx.ID_NO):
		logInfo("Rejected threaded model")
		accept = False
	    else:
		logInfo("Canceled Finalize operation")
		dlg.Destroy()
		return
	    dlg.Destroy()
	    self.viewMenu.Disable()
	    self.scoretypeMenu.Disable()
	    self.parent.GoBtn.Enable()
	    self.btnAddTemplate.Enable()
	    self.btnRemoveTemplate.Enable()
	    self.btnLoad.Enable()
	    if (platform.system() == "Darwin"):
		self.btnThread.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/btnThread.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
	    else:
		self.btnThread.SetLabel("Finalize!")
	    self.buttonState = "Thread!"
	    self.btnThread.SetToolTipString("Perform docking simulation with selected parameters")
	    self.cmd.label("all", "")
	    self.seqWin.cannotDelete = False
	    if (not(accept)):
		self.cmd.remove("thread_view")
		self.cmd.delete("thread_view")
		return
	    # Load the structure in PyMOL and the sequence viewer
	    if (platform.system() == "Windows"):
		newname = os.path.expanduser("~") + "\\thread_T.pdb"
	    else:
		newname = os.path.expanduser("~") + "/thread_T.pdb"
	    os.rename(self.selectedModel, newname)
	    try:
		self.cmd.remove("thread_view")
		self.cmd.delete("thread_view")
		self.cmd.load(newname, "thread_T")
		defaultPyMOLView(self.cmd, "thread_T")
		self.seqWin.PyMOLPDBLoad(0, newname)
		del self.threadView
	    except:
		# Some weird error happened, do nothing instead of crashing
		print "Bug at accept button click"
		pass	
    
    def recoverFromError(self, msg=""):
	# This function tells the user what the error was and tries to revert the protocol
	# back to the pre-daemon state so the main GUI can continue to be used
	if (len(msg) == 0):
	    f = open("errreport", "r")
	    errmsg = "An error was encountered during the protocol:\n\n"
	    for aline in f:
		errmsg = errmsg + str(aline)
	    f.close()
	    os.remove("errreport")
	else:
	    errmsg = msg
	logInfo("Error Encountered")
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
	self.seqWin.cannotDelete = False
	self.parent.GoBtn.Enable()
	self.btnAddTemplate.Enable()
	self.btnRemoveTemplate.Enable()
	self.btnLoad.Enable()
	self.btnThread.SetLabel("Thread!")
	self.btnThread.SetToolTipString("Begin threading the sequence onto the template structures")
	# Get rid of the messages
	for i in range(0, len(self.seqWin.msgQueue)):
	    if (self.seqWin.msgQueue[i].find("Performing protein threading") >= 0):
		self.seqWin.msgQueue.pop(i)
		break
	if (len(self.seqWin.msgQueue) > 0):
	    self.seqWin.labelMsg.SetLabel(self.seqWin.msgQueue[len(self.seqWin.msgQueue)-1])
	else:
	    self.seqWin.labelMsg.SetLabel("")
	self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
	self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")
    
    def threadThread(self, event):
	# Dump a file with the loop modeling parameters for the daemon to pick up
	goToSandbox()
	if (self.stage == 1):
	    self.tmrThread.Stop()
	    self.timeoutCount = 0
	    f = open("threadinputtemp", "w")
	    for pdb in self.templates:
		f.write("PDBFILE\t" + pdb.strip() + ".pdb\n")
		f2 = open(pdb.strip() + ".pdb", "r")
		f.write("BEGIN PDB DATA\n")
		for aline in f2:
		    f.write(aline.strip() + "\n")
		f.write("END PDB DATA\n")
		f2.close()
	    f.write("FRAG\t" + self.loadedfile + "\n")
	    f.close()
	    appendScorefxnParamsInfoToFile("threadinputtemp", self.selectWin.weightsfile)
	    if (useServer and False):
		try: 
		    self.ID = sendToServer("threadinput")
		    self.usingServer = True
		    logInfo("Threading input sent to server daemon with ID " + self.ID)
		except:
		    # Something failed, default to the local daemon
		    os.rename("threadinputtemp", "threadinput")
		    self.usingServer = False
		    logInfo("Server daemon not available, threading input uploaded at threadinput")
	    else:
		os.rename("threadinputtemp", "threadinput")
		self.usingServer = False
		logInfo("Threading input uploaded locally at threadinput")
	    self.stage = 2
	    self.tmrThread.Start(1000)
	elif (self.stage == 2):
	    if (self.usingServer):
		# See if the file has been uploaded yet and bring it here if so
		queryServerForResults("threadoutput-" + self.ID)
		self.timeoutCount = self.timeoutCount + 1
	    if (self.timeoutCount >= serverTimeout):
		self.tmrThread.Stop()
		# If this is taking too long, maybe there's something wrong with the server
		# Ask the user if they want to continue waiting or use the local daemon instead
		dlg = wx.MessageDialog(self, "The server is taking a long time to respond.  Continue to wait?  Pressing No will run the calculations locally.", "Delayed Server Response", wx.YES_NO | wx.ICON_EXCLAMATION | wx.CENTRE)
		if (dlg.ShowModal() == wx.ID_YES):
		    # Reset the counter
		    self.timeoutCount = 0
		else:
		    self.usingServer = False
		    self.timeoutCount = 0
		    os.rename("threadinputtemp", "threadinput")
		    logInfo("Server took too long to respond so the local daemon was used")
		dlg.Destroy()
		self.tmrThread.Start(1000)
	    # Read the output dumped by the child process (finally!)
	    if (os.path.isfile("threadoutput")):
		self.tmrThread.Stop()
		self.residue_E = []
		f = open("threadoutput", "r")
		for aline in f:
		    if (aline[0:6] == "OUTPUT"):
			pdbfile = aline.split("\t")[1].strip()
			self.threadmodel = pdbfile
			#self.dockView = pose_from_pdb(pdbfile)
		    elif (aline[0:6] == "ENERGY"):
			if (aline.split()[1] == "total_score"):
			    # This is the scoretype line, row 0 in residue_E
			    self.residue_E.append(aline.split()[1:])
			else:
			    self.residue_E.append([])
			    indx = len(self.residue_E) - 1
			    for E in aline.split()[1:]:
				self.residue_E[indx].append(float(E))
		f.close()
		self.threadView = self.seqWin.pdbreader.get_structure("thread_view", self.threadmodel)
		self.selectedModel = self.threadmodel
		logInfo("Found threading output at threadoutput")
		# Pop this message out of the queue
		for i in range(0, len(self.seqWin.msgQueue)):
		    if (self.seqWin.msgQueue[i].find("Performing protein threading") >= 0):
			self.seqWin.msgQueue.pop(i)
			break
		if (len(self.seqWin.msgQueue) > 0):
		    self.seqWin.labelMsg.SetLabel(self.seqWin.msgQueue[len(self.seqWin.msgQueue)-1])
		else:
		    self.seqWin.labelMsg.SetLabel("")
		self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
		self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")
		# Add these loop residues to the view menu so the user can look at the new loop
		viewoptions = ["Whole View"]
		ires = 1
		for ch in self.threadView[0]:
		    for residue in ch:
			chain = ch.id
			if (len(chain.strip()) == 0):
			    chain = "_"
			seqpos = str(residue.id[1])
			resn = AA3to1(residue.resname)
			viewoptions.append(chain + ":" + resn + seqpos)
			ires = ires + 1
		self.viewMenu.Clear()
		self.viewMenu.AppendItems(viewoptions)
		self.viewMenu.SetSelection(0)
		self.viewMenu.Enable()
		# Add the nonzero scoretypes to the energy viewing list from the current score function
		self.scoretypeMenu.Clear()
		for scoretype in self.residue_E[0]:
		    try:
			toAdd = scoretypes[str(scoretype)]
		    except:
			toAdd = str(scoretype)
		    self.scoretypeMenu.Append(toAdd)
		self.scoretypeMenu.Enable()
		self.parent.GoBtn.Enable()
		self.btnThread.Enable()
		if (platform.system() == "Darwin"):
		    self.btnThread.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/btnThread_Finalize.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
		else:
		    self.btnThread.SetLabel("Finalize!")
		self.buttonState = "Finalize!"
		self.btnThread.SetToolTipString("Accept or reject protocol results")
		os.remove("threadoutput")
		# Load the docked pose as the "dock_view" model so the user can look at the results
		self.cmd.load(self.threadmodel, "thread_view")
		self.cmd.hide("everything", "model thread_view")
		recolorEnergies(self.threadView, self.residue_E, "thread_view", self.scoretypeMenu.GetValue(), self.cmd)
		self.focusView(self.viewMenu.GetValue(), "thread_view")
	    elif (os.path.isfile("errreport")):
		# Something went wrong, tell the user about it
		self.tmrThread.Stop()
		self.recoverFromError()
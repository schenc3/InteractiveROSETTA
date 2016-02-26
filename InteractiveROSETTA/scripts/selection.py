import wx
import wx.lib.scrolledpanel
import os
import os.path
import sys
import platform
import math
import numpy
import webbrowser
from tools import resizeTextControlForUNIX
from tools import logInfo
from tools import defaultPyMOLView
from tools import getChainColor

class SelectPanel(wx.lib.scrolledpanel.ScrolledPanel):
    def __init__(self, parent, W, H):
	#if (platform.system() == "Windows"):
	wx.lib.scrolledpanel.ScrolledPanel.__init__(self, parent, id=-1, pos=(10, H-270), size=(340, 230), name="SelectPanel")
	#else:
	#    wx.lib.scrolledpanel.ScrolledPanel.__init__(self, parent, id=-1, pos=(10, H-270), size=(340, 230), name="SelectPanel")
	#self.scorefxn = create_score_function("talaris2013")
	self.SetBackgroundColour("#333333")
	self.Show()
	self.parent = parent

	# Find the database location if it is not currently set as an environment variable
	try:
	    os.environ["PYROSETTA_DATABASE"]
	except:
	    if (platform.system() == "Windows"):
		cfgfile = os.path.expanduser("~") + "/InteractiveROSETTA/seqwindow.cfg"
	    else:
		cfgfile = os.path.expanduser("~") + "/.InteractiveROSETTA/seqwindow.cfg"
	    f = open(cfgfile.strip(), "r")
	    rosettadb = "Not Found"
	    for aline in f:
		if ("[ROSETTADB]" in aline):
		    rosettadb = aline.split("\t")[1].strip()
	    f.close()
	    os.environ["PYROSETTA_DATABASE"] = rosettadb

	if (platform.system() == "Windows"):
	    self.weightsfile = os.getenv("PYROSETTA_DATABASE") + "\\scoring\\weights\\talaris2013.wts"
	    self.labelAtoms = wx.StaticText(self, -1, "Atoms", (5, 5), (70, 20), wx.ALIGN_CENTRE)
	    self.labelAtoms.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
	elif (platform.system() == "Darwin"):
	    self.weightsfile = os.getenv("PYROSETTA_DATABASE") + "/scoring/weights/talaris2013.wts"
	    self.labelAtoms = wx.StaticBitmap(self, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/labelAtoms.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(5, 5), size=(70, 20))
	else:
	    self.weightsfile = os.getenv("PYROSETTA_DATABASE") + "/scoring/weights/talaris2013.wts"
	    self.labelAtoms = wx.StaticText(self, -1, "Atoms", pos=(13, 5), style=wx.ALIGN_CENTRE)
	    self.labelAtoms.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
	    resizeTextControlForUNIX(self.labelAtoms, 5, 70)
	self.labelAtoms.SetForegroundColour("#FFFFFF")

	if (platform.system() == "Darwin"):
	    self.AtomOffBtn = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/AtomOffBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(5, 30), size=(35, 35))
	else:
	    self.AtomOffBtn = wx.Button(self, id=-1, label="X", pos=(5, 30), size=(35, 35))
	    #self.AtomOffBtn.SetBackgroundColour("#000000")
	    self.AtomOffBtn.SetForegroundColour("#000000")
	    self.AtomOffBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.AtomOffBtn.Bind(wx.EVT_BUTTON, self.atomOff)
	self.AtomOffBtn.SetToolTipString("Hide selected atoms")

	self.AtomLinesBtn = wx.BitmapButton(self, -1, wx.Image(self.parent.scriptdir + "/images/lines.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(40, 30), size=(35, 35))
	#self.AtomLinesBtn.SetBackgroundColour("#000000")
	self.AtomLinesBtn.Bind(wx.EVT_BUTTON, self.atomLines)
	self.AtomLinesBtn.SetToolTipString("Apply thin bond-line view to selection")

	self.AtomBallsLinesBtn = wx.BitmapButton(self, -1, wx.Image(self.parent.scriptdir + "/images/balls_lines.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(5, 65), size=(35, 35))
	#self.AtomVDWBtn.SetBackgroundColour("#000000")
	self.AtomBallsLinesBtn.Bind(wx.EVT_BUTTON, self.atomBallsLines)
	self.AtomBallsLinesBtn.SetToolTipString("Apply thin bond-line-with-spheres view to selection")

	self.AtomSticksBtn = wx.BitmapButton(self, -1, wx.Image(self.parent.scriptdir + "/images/sticks.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(40, 65), size=(35, 35))
	#self.AtomSticksBtn.SetBackgroundColour("#000000")
	self.AtomSticksBtn.Bind(wx.EVT_BUTTON, self.atomSticks)
	self.AtomSticksBtn.SetToolTipString("Apply thick bond-line view to selection")

	self.AtomBallsSticksBtn = wx.BitmapButton(self, -1, wx.Image(self.parent.scriptdir + "/images/balls_sticks.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(5, 100), size=(35, 35))
	#self.AtomSticksBtn.SetBackgroundColour("#000000")
	self.AtomBallsSticksBtn.Bind(wx.EVT_BUTTON, self.atomBallsSticks)
	self.AtomBallsSticksBtn.SetToolTipString("Apply thick bond-line-with-spheres view to selection")

	self.AtomVDWBtn = wx.BitmapButton(self, -1, wx.Image(self.parent.scriptdir + "/images/VDW.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(40, 100), size=(35, 35))
	#self.AtomVDWBtn.SetBackgroundColour("#000000")
	self.AtomVDWBtn.Bind(wx.EVT_BUTTON, self.atomVDW)
	self.AtomVDWBtn.SetToolTipString("Apply VDW spheres view to selection")

	self.AtomRecolorBtn = wx.BitmapButton(self, -1, wx.Image(self.parent.scriptdir + "/images/colorwheel.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), (5, 135), (70, 35))
	#self.AtomWhiteBtn.SetBackgroundColour("#000000")
	self.AtomRecolorBtn.Bind(wx.EVT_BUTTON, self.atomRecolor)
	self.AtomRecolorBtn.SetToolTipString("Change the color of the selected atoms")

	if (platform.system() == "Darwin"):
	    self.AtomStandardBtn = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/AtomStandardBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(5, 169), size=(70, 20))
	else:
	    self.AtomStandardBtn = wx.Button(self, id=-1, label="Default", pos=(5, 169), size=(70, 20))
	    #self.AtomStandardBtn.SetBackgroundColour("#000000")
	    self.AtomStandardBtn.SetForegroundColour("#000000")
	    self.AtomStandardBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.AtomStandardBtn.Bind(wx.EVT_BUTTON, self.atomStandardColor)
	self.AtomStandardBtn.SetToolTipString("Change to default colors for selected atoms")

	if (platform.system() == "Darwin"):
	    self.AtomChainBtn = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/AtomChainBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(5, 189), size=(35, 20))
	else:
	    self.AtomChainBtn = wx.Button(self, id=-1, label="Ch", pos=(5, 189), size=(35, 20))
	    self.AtomChainBtn.SetForegroundColour("#000000")
	    self.AtomChainBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.AtomChainBtn.Bind(wx.EVT_BUTTON, self.atomChainColor)
	self.AtomChainBtn.SetToolTipString("Color selected atoms by chain")
	if (platform.system() == "Darwin"):
	    self.AtomTermBtn = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/AtomTerminiBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(40, 189), size=(35, 20))
	else:
	    self.AtomTermBtn = wx.Button(self, id=-1, label="Ter", pos=(40, 189), size=(35, 20))
	    self.AtomTermBtn.SetForegroundColour("#000000")
	    self.AtomTermBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.AtomTermBtn.Bind(wx.EVT_BUTTON, self.atomTermini)
	self.AtomTermBtn.SetToolTipString("Color selected atoms by termini")

	if (platform.system() == "Darwin"):
	    self.LabelOffBtn = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/LabelOffBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(5, 209), size=(35, 20))
	else:
	    self.LabelOffBtn = wx.Button(self, id=-1, label="X", pos=(5, 209), size=(35, 20))
	    self.LabelOffBtn.SetForegroundColour("#000000")
	    self.LabelOffBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.LabelOffBtn.Bind(wx.EVT_BUTTON, self.labelsOff)
	self.LabelOffBtn.SetToolTipString("Turn labels off on the selected residues")
	if (platform.system() == "Darwin"):
	    self.LabelOnBtn = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/LabelOnBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(40, 209), size=(35, 20))
	else:
	    self.LabelOnBtn = wx.Button(self, id=-1, label="Lbl", pos=(40, 209), size=(35, 20))
	    self.LabelOnBtn.SetForegroundColour("#000000")
	    self.LabelOnBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.LabelOnBtn.Bind(wx.EVT_BUTTON, self.labelResidues)
	self.LabelOnBtn.SetToolTipString("Label the selected residues")

	if (platform.system() == "Windows"):
	    self.labelRibbons = wx.StaticText(self, -1, "Ribbons", (85, 5), (70, 20), wx.ALIGN_CENTRE)
	    self.labelRibbons.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
	elif (platform.system() == "Darwin"):
	    self.labelRibbons = wx.StaticBitmap(self, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/labelRibbons.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(85, 5), size=(70, 20))
	else:
	    self.labelRibbons = wx.StaticText(self, -1, "Ribbons", pos=(85, 5), style=wx.ALIGN_CENTRE)
	    self.labelRibbons.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
	    resizeTextControlForUNIX(self.labelRibbons, 85, 70)
	self.labelRibbons.SetForegroundColour("#FFFFFF")

	if (platform.system() == "Darwin"):
	    self.RibbonOffBtn = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/RibbonOffBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(85, 30), size=(70, 35))
	else:
	    self.RibbonOffBtn = wx.Button(self, id=-1, label="X", pos=(85, 30), size=(70, 35))
	    #self.RibbonOffBtn.SetBackgroundColour("#000000")
	    self.RibbonOffBtn.SetForegroundColour("#000000")
	    self.RibbonOffBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.RibbonOffBtn.Bind(wx.EVT_BUTTON, self.ribbonsOff)
	self.RibbonOffBtn.SetToolTipString("Hide ribbon view for selected atoms")

	self.RibbonRibbonBtn = wx.BitmapButton(self, -1, wx.Image(self.parent.scriptdir + "/images/ribbons.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(85, 65), size=(70, 35))
	#self.RibbonOnBtn.SetBackgroundColour("#000000")
	self.RibbonRibbonBtn.Bind(wx.EVT_BUTTON, self.ribbonsRibbon)
	self.RibbonRibbonBtn.SetToolTipString("Apply thin string view to selected ribbon")

	self.RibbonCartoonBtn = wx.BitmapButton(self, -1, wx.Image(self.parent.scriptdir + "/images/cartoons.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(85, 100), size=(70, 35))
	#self.RibbonOnBtn.SetBackgroundColour("#000000")
	self.RibbonCartoonBtn.Bind(wx.EVT_BUTTON, self.ribbonsCartoon)
	self.RibbonCartoonBtn.SetToolTipString("Apply cartoon view to selected ribbon")

	self.RibbonRecolorBtn = wx.BitmapButton(self, -1, wx.Image(self.parent.scriptdir + "/images/colorwheel.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), (85, 135), (70, 35))
	#self.RibbonWhiteBtn.SetBackgroundColour("#000000")
	self.RibbonRecolorBtn.Bind(wx.EVT_BUTTON, self.ribbonRecolor)
	self.RibbonRecolorBtn.SetToolTipString("Change the color of the selected ribbons")

	if (platform.system() == "Darwin"):
	    self.RibbonStandardBtn = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/AtomStandardBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(85, 169), size=(70, 20))
	else:
	    self.RibbonStandardBtn = wx.Button(self, id=-1, label="Default", pos=(85, 169), size=(70, 20))
	    #self.RibbonStandardBtn.SetBackgroundColour("#000000")
	    self.RibbonStandardBtn.SetForegroundColour("#000000")
	    self.RibbonStandardBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.RibbonStandardBtn.Bind(wx.EVT_BUTTON, self.ribbonStandardColor)
	self.RibbonStandardBtn.SetToolTipString("Color selected ribbons by standard secondary structure coloring")

	if (platform.system() == "Darwin"):
	    self.RibbonChainBtn = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/AtomChainBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(85, 189), size=(35, 20))
	else:
	    self.RibbonChainBtn = wx.Button(self, id=-1, label="Ch", pos=(85, 189), size=(35, 20))
	    #self.RibbonChainBtn.SetBackgroundColour("#000000")
	    self.RibbonChainBtn.SetForegroundColour("#000000")
	    self.RibbonChainBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.RibbonChainBtn.Bind(wx.EVT_BUTTON, self.ribbonChainColor)
	self.RibbonChainBtn.SetToolTipString("Color selected ribbons by chain")
	if (platform.system() == "Darwin"):
	    self.RibbonTermBtn = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/AtomTerminiBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(120, 189), size=(35, 20))
	else:
	    self.RibbonTermBtn = wx.Button(self, id=-1, label="Ter", pos=(120, 189), size=(35, 20))
	    #self.RibbonTermBtn.SetBackgroundColour("#000000")
	    self.RibbonTermBtn.SetForegroundColour("#000000")
	    self.RibbonTermBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.RibbonTermBtn.Bind(wx.EVT_BUTTON, self.ribbonTermini)
	self.RibbonTermBtn.SetToolTipString("Color selected ribbons by termini")

	if (platform.system() == "Darwin"):
	    self.ToggleSurfBtn = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/ToggleSurfBtn_Off.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(85, 209), size=(70, 20))
	else:
	    self.ToggleSurfBtn = wx.Button(self, id=-1, label="Surf Off", pos=(85, 209), size=(70, 20))
	    self.ToggleSurfBtn.SetForegroundColour("#000000")
	    self.ToggleSurfBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.ToggleSurfBtn.Bind(wx.EVT_BUTTON, self.toggleSurf)
	self.ToggleSurfBtn.SetToolTipString("Display configured surfaces")
	self.showSurf = False

	if (platform.system() == "Windows"):
	    self.labelSelection = wx.StaticText(self, -1, "Selection", (175, 5), (140, 20), wx.ALIGN_CENTRE)
	    self.labelSelection.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
	elif (platform.system() == "Darwin"):
	    self.labelSelection = wx.StaticBitmap(self, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/labelSelection.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(175, 5), size=(140, 20))
	else:
	    self.labelSelection = wx.StaticText(self, -1, "Selection", pos=(205, 5), style=wx.ALIGN_CENTRE)
	    self.labelSelection.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
	    resizeTextControlForUNIX(self.labelSelection, 165, 160)
	self.labelSelection.SetForegroundColour("#FFFFFF")

	if (platform.system() == "Darwin"):
	    self.HelpBtn = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/HelpBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(315, 0), size=(25, 25))
	else:
	    self.HelpBtn = wx.Button(self, id=-1, label="?", pos=(315, 0), size=(25, 25))
	    self.HelpBtn.SetForegroundColour("#0000FF")
	    self.HelpBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.HelpBtn.Bind(wx.EVT_BUTTON, self.showHelp)
	self.HelpBtn.SetToolTipString("Display the help file for this window")

	self.SelectPanel = wx.Panel(self, id=-1, pos=(160, 25), size=(170, 185))
	self.SelectPanel.SetBackgroundColour("#333333")
	if (platform.system() == "Darwin"):
	    self.SelectAllBtn = wx.BitmapButton(self.SelectPanel, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/SelectAllBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(5, 5), size=(32, 30))
	else:
	    self.SelectAllBtn = wx.Button(self.SelectPanel, id=-1, label="All", pos=(5, 5), size=(32, 30))
	    self.SelectAllBtn.SetForegroundColour("#000000")
	    self.SelectAllBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.SelectAllBtn.Bind(wx.EVT_BUTTON, self.selectAll)
	self.SelectAllBtn.SetToolTipString("Select all atoms")
	if (platform.system() == "Darwin"):
	    self.SelectInvertBtn = wx.BitmapButton(self.SelectPanel, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/SelectInvertBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(37, 5), size=(32, 30))
	else:
	    self.SelectInvertBtn = wx.Button(self.SelectPanel, id=-1, label="Inv", pos=(37, 5), size=(32, 30))
	    self.SelectInvertBtn.SetForegroundColour("#000000")
	    self.SelectInvertBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.SelectInvertBtn.Bind(wx.EVT_BUTTON, self.selectInvert)
	self.SelectInvertBtn.SetToolTipString("Select all atoms not currently selected")
	if (platform.system() == "Darwin"):
	    self.SelectVisibleBtn = wx.BitmapButton(self.SelectPanel, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/SelectVisibleBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(69, 5), size=(32, 30))
	else:
	    self.SelectVisibleBtn = wx.Button(self.SelectPanel, id=-1, label="Vis", pos=(69, 5), size=(32, 30))
	    self.SelectVisibleBtn.SetForegroundColour("#000000")
	    self.SelectVisibleBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.SelectVisibleBtn.Bind(wx.EVT_BUTTON, self.selectVisible)
	self.SelectVisibleBtn.SetToolTipString("Select all visible atoms")
	if (platform.system() == "Darwin"):
	    self.SelectBackboneBtn = wx.BitmapButton(self.SelectPanel, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/SelectBBBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(101, 5), size=(32, 30))
	else:
	    self.SelectBackboneBtn = wx.Button(self.SelectPanel, id=-1, label="BB", pos=(101, 5), size=(32, 30))
	    self.SelectBackboneBtn.SetForegroundColour("#000000")
	    self.SelectBackboneBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.SelectBackboneBtn.Bind(wx.EVT_BUTTON, self.selectBB)
	self.SelectBackboneBtn.SetToolTipString("Select all backbone atoms")
	if (platform.system() == "Darwin"):
	    self.SelectSidechainsBtn = wx.BitmapButton(self.SelectPanel, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/SelectSCBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(133, 5), size=(32, 30))
	else:
	    self.SelectSidechainsBtn = wx.Button(self.SelectPanel, id=-1, label="SC", pos=(133, 5), size=(32, 30))
	    self.SelectSidechainsBtn.SetForegroundColour("#000000")
	    self.SelectSidechainsBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.SelectSidechainsBtn.Bind(wx.EVT_BUTTON, self.selectSC)
	self.SelectSidechainsBtn.SetToolTipString("Select all sidechain atoms")

	if (platform.system() == "Darwin"):
	    self.SelectHBtn = wx.BitmapButton(self.SelectPanel, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/SelectHBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(5, 35), size=(32, 30))
	else:
	    self.SelectHBtn = wx.Button(self.SelectPanel, id=-1, label="H", pos=(5, 35), size=(32, 30))
	    self.SelectHBtn.SetForegroundColour("#000000")
	    self.SelectHBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.SelectHBtn.Bind(wx.EVT_BUTTON, self.selectH)
	self.SelectHBtn.SetToolTipString("Select all hydrogen atoms")
	if (platform.system() == "Darwin"):
	    self.SelectCBtn = wx.BitmapButton(self.SelectPanel, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/SelectCBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(37, 35), size=(32, 30))
	else:
	    self.SelectCBtn = wx.Button(self.SelectPanel, id=-1, label="C", pos=(37, 35), size=(32, 30))
	    self.SelectCBtn.SetForegroundColour("#000000")
	    self.SelectCBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.SelectCBtn.Bind(wx.EVT_BUTTON, self.selectC)
	self.SelectCBtn.SetToolTipString("Select all carbon atoms")
	if (platform.system() == "Darwin"):
	    self.SelectNBtn = wx.BitmapButton(self.SelectPanel, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/SelectNBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(69, 35), size=(32, 30))
	else:
	    self.SelectNBtn = wx.Button(self.SelectPanel, id=-1, label="N", pos=(69, 35), size=(32, 30))
	    self.SelectNBtn.SetForegroundColour("#000000")
	    self.SelectNBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.SelectNBtn.Bind(wx.EVT_BUTTON, self.selectN)
	self.SelectNBtn.SetToolTipString("Select all nitrogen atoms")
	if (platform.system() == "Darwin"):
	    self.SelectOBtn = wx.BitmapButton(self.SelectPanel, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/SelectOBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(101, 35), size=(32, 30))
	else:
	    self.SelectOBtn = wx.Button(self.SelectPanel, id=-1, label="O", pos=(101, 35), size=(32, 30))
	    self.SelectOBtn.SetForegroundColour("#000000")
	    self.SelectOBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.SelectOBtn.Bind(wx.EVT_BUTTON, self.selectO)
	self.SelectOBtn.SetToolTipString("Select all oxygen atoms")
	if (platform.system() == "Darwin"):
	    self.SelectSolventBtn = wx.BitmapButton(self.SelectPanel, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/SelectSolventBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(133, 35), size=(32, 30))
	else:
	    self.SelectSolventBtn = wx.Button(self.SelectPanel, id=-1, label="Sol", pos=(133, 35), size=(32, 30))
	    self.SelectSolventBtn.SetForegroundColour("#000000")
	    self.SelectSolventBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.SelectSolventBtn.Bind(wx.EVT_BUTTON, self.selectSolvent)
	self.SelectSolventBtn.SetToolTipString("Select all solvent atoms")

	if (platform.system() == "Darwin"):
	    self.SelectExtendBtn = wx.BitmapButton(self.SelectPanel, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/SelectExtendBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(5, 65), size=(110, 30))
	else:
	    self.SelectExtendBtn = wx.Button(self.SelectPanel, id=-1, label="Extend", pos=(5, 65), size=(110, 30))
	    #self.SelectExtendBtn.SetBackgroundColour("#000000")
	    self.SelectExtendBtn.SetForegroundColour("#000000")
	    self.SelectExtendBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.SelectExtendBtn.Bind(wx.EVT_BUTTON, self.selectExtend)
	self.SelectExtendBtn.SetToolTipString("Extend the current selection by a given Angstrom radius")

	self.ExtendValueTxt = wx.TextCtrl(self.SelectPanel, -1, pos=(115, 67), size=(25, 25))
	self.ExtendValueTxt.SetValue("8")
	self.ExtendValueTxt.SetToolTipString("Radius in Angstroms of selection expansion")

	if (platform.system() == "Darwin"):
	    self.labelExtendA = wx.StaticBitmap(self.SelectPanel, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/labelAngstrom.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(142, 70), size=(25, 25))
	else:
	    self.labelExtendA = wx.StaticText(self.SelectPanel, -1, "A", (140, 68), (25, 25))
	    self.labelExtendA.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
	    self.labelExtendA.SetForegroundColour("#FFFFFF")

	if (platform.system() == "Darwin"):
	    self.SelectZoomBtn = wx.BitmapButton(self.SelectPanel, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/SelectZoomBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(5, 95), size=(80, 30))
	else:
	    self.SelectZoomBtn = wx.Button(self.SelectPanel, id=-1, label="Zoom", pos=(5, 95), size=(80, 30))
	    self.SelectZoomBtn.SetForegroundColour("#000000")
	    self.SelectZoomBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.SelectZoomBtn.Bind(wx.EVT_BUTTON, self.zoomSelection)
	self.SelectZoomBtn.SetToolTipString("Zoom in on the current selection")
	if (platform.system() == "Darwin"):
	    self.CenterBtn = wx.BitmapButton(self.SelectPanel, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/CenterBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(85, 95), size=(80, 30))
	else:
	    self.CenterBtn = wx.Button(self.SelectPanel, id=-1, label="Center", pos=(85, 95), size=(80, 30))
	    self.CenterBtn.SetForegroundColour("#000000")
	    self.CenterBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.CenterBtn.Bind(wx.EVT_BUTTON, self.centerSelection)
	self.CenterBtn.SetToolTipString("Center the origin of rotation in current selection")

	if (platform.system() == "Darwin"):
	    self.SelectNeighborhoodBtn = wx.BitmapButton(self.SelectPanel, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/SelectNeighborhoodBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(5, 125), size=(110, 30))
	else:
	    self.SelectNeighborhoodBtn = wx.Button(self.SelectPanel, id=-1, label="Neighborhood", pos=(5, 125), size=(110, 30))
	    #self.SelectNeighborhoodBtn.SetBackgroundColour("#000000")
	    self.SelectNeighborhoodBtn.SetForegroundColour("#000000")
	    self.SelectNeighborhoodBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.SelectNeighborhoodBtn.Bind(wx.EVT_BUTTON, self.neighborhoodView)
	self.SelectNeighborhoodBtn.SetToolTipString("View structure within a given radius in Angstroms of the current selection")

	self.NeighborhoodValueTxt = wx.TextCtrl(self.SelectPanel, -1, pos=(115, 127), size=(25, 25))
	self.NeighborhoodValueTxt.SetValue("8")
	self.NeighborhoodValueTxt.SetToolTipString("Radius of neighborhood view")

	if (platform.system() == "Darwin"):
	    self.labelNeighborhoodA = wx.StaticBitmap(self.SelectPanel, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/labelAngstrom.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(142, 130), size=(25, 25))
	else:
	    self.labelNeighborhoodA = wx.StaticText(self.SelectPanel, -1, "A", (140, 128), (25, 25))
	    self.labelNeighborhoodA.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
	    self.labelNeighborhoodA.SetForegroundColour("#FFFFFF")

	if (platform.system() == "Darwin"):
	    self.labelSF = wx.StaticBitmap(self.SelectPanel, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/labelSF.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(5, 155), size=(50, 25))
	else:
	    self.labelSF = wx.StaticText(self.SelectPanel, -1, "SFXN:", (5, 159), (50, 25))
	    self.labelSF.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
	    self.labelSF.SetForegroundColour("#FFFFFF")

	if (platform.system() == "Darwin"):
	    self.SelectScorefxnBtn = wx.Button(self.SelectPanel, id=-1, label="talaris2013", pos=(55, 155), size=(80, 25))
	else:
	    self.SelectScorefxnBtn = wx.Button(self.SelectPanel, id=-1, label="talaris2013", pos=(55, 155), size=(85, 25))
	self.SelectScorefxnBtn.SetForegroundColour("#000000")
	self.SelectScorefxnBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
	self.SelectScorefxnBtn.Bind(wx.EVT_BUTTON, self.changeSF)
	self.SelectScorefxnBtn.SetToolTipString("Select the Rosetta scoring function")

	if (platform.system() == "Darwin"):
	    self.SelectPanelToggleBtn = wx.BitmapButton(self.SelectPanel, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/SelectPanelFlipBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(140, 155), size=(25, 25))
	else:
	    self.SelectPanelToggleBtn = wx.Button(self.SelectPanel, id=-1, label="+", pos=(140, 155), size=(25, 25))
	    self.SelectPanelToggleBtn.SetForegroundColour("#000000")
	    self.SelectPanelToggleBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.SelectPanelToggleBtn.SetToolTipString("Switch to manipulation mode")
	self.SelectPanelToggleBtn.Bind(wx.EVT_LEFT_DOWN, self.panelToggle)
	self.panelMode = "Selection"

	self.ManipulatePanel = wx.Panel(self, id=-1, pos=(160, 25), size=(180, 205))
	self.XYPlane = wx.StaticBitmap(self.ManipulatePanel, -1, wx.Image(self.parent.scriptdir + "/images/plane.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(15, 15), size=(150, 150))
	self.XYPlane.Bind(wx.EVT_MOTION, self.xyMotion)
	self.XYPlane.Bind(wx.EVT_LEFT_UP, self.xyUp)
	self.ManipulatePanel.Bind(wx.EVT_LEFT_UP, self.xyUp)
	self.Bind(wx.EVT_LEFT_UP, self.xyUp)
	self.XYPlane.Bind(wx.EVT_LEFT_DOWN, self.xyDown)
	self.XYPlane.Bind(wx.EVT_ACTIVATE, self.xyFocus)
	self.leftDown = False

	if (platform.system() == "Darwin"):
	    self.UpBtn = wx.BitmapButton(self.ManipulatePanel, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/UpBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(60, 0), size=(60, 15))
	else:
	    self.UpBtn = wx.Button(self.ManipulatePanel, id=-1, label="", pos=(60, 0), size=(60, 15))
	    self.UpBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.NORMAL))
	    self.UpBtn.SetForegroundColour("#000000")
	self.UpBtn.Bind(wx.EVT_BUTTON, self.nudgeUp)
	self.UpBtn.SetToolTipString("Nudge up")
	if (platform.system() == "Darwin"):
	    self.DownBtn = wx.BitmapButton(self.ManipulatePanel, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/UpBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(60, 160), size=(60, 15))
	else:
	    self.DownBtn = wx.Button(self.ManipulatePanel, id=-1, label="", pos=(60, 165), size=(60, 15))
	    self.DownBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.NORMAL))
	    self.DownBtn.SetForegroundColour("#000000")
	self.DownBtn.Bind(wx.EVT_BUTTON, self.nudgeDown)
	self.DownBtn.SetToolTipString("Nudge down")
	if (platform.system() == "Darwin"):
	    self.LeftBtn = wx.BitmapButton(self.ManipulatePanel, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/LeftBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 60), size=(15, 60))
	else:
	    self.LeftBtn = wx.Button(self.ManipulatePanel, id=-1, label="", pos=(0, 60), size=(15, 60))
	    self.LeftBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.NORMAL))
	    self.LeftBtn.SetForegroundColour("#000000")
	self.LeftBtn.Bind(wx.EVT_BUTTON, self.nudgeLeft)
	self.LeftBtn.SetToolTipString("Nudge left")
	if (platform.system() == "Darwin"):
	    self.RightBtn = wx.BitmapButton(self.ManipulatePanel, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/LeftBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(165, 60), size=(15, 60))
	else:
	    self.RightBtn = wx.Button(self.ManipulatePanel, id=-1, label="", pos=(165, 60), size=(15, 60))
	    self.RightBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.NORMAL))
	    self.RightBtn.SetForegroundColour("#000000")
	self.RightBtn.Bind(wx.EVT_BUTTON, self.nudgeRight)
	self.RightBtn.SetToolTipString("Nudge right")
	if (platform.system() == "Darwin"):
	    self.LeftUpBtn = wx.BitmapButton(self.ManipulatePanel, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/CornerBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 0), size=(30, 30))
	else:
	    self.LeftUpBtn = wx.Button(self.ManipulatePanel, id=-1, label="", pos=(0, 0), size=(30, 30))
	    self.LeftUpBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.NORMAL))
	    self.LeftUpBtn.SetForegroundColour("#000000")
	self.LeftUpBtn.Bind(wx.EVT_BUTTON, self.nudgeLeftUp)
	self.LeftUpBtn.SetToolTipString("Nudge left up")
	if (platform.system() == "Darwin"):
	    self.RightUpBtn = wx.BitmapButton(self.ManipulatePanel, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/CornerBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(150, 0), size=(30, 30))
	else:
	    self.RightUpBtn = wx.Button(self.ManipulatePanel, id=-1, label="", pos=(150, 0), size=(30, 30))
	    self.RightUpBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.NORMAL))
	    self.RightUpBtn.SetForegroundColour("#000000")
	self.RightUpBtn.Bind(wx.EVT_BUTTON, self.nudgeRightUp)
	self.RightUpBtn.SetToolTipString("Nudge right up")
	if (platform.system() == "Darwin"):
	    self.LeftDownBtn = wx.BitmapButton(self.ManipulatePanel, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/CornerBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 150), size=(30, 30))
	else:
	    self.LeftDownBtn = wx.Button(self.ManipulatePanel, id=-1, label="", pos=(0, 150), size=(30, 30))
	    self.LeftDownBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.NORMAL))
	    self.LeftDownBtn.SetForegroundColour("#000000")
	self.LeftDownBtn.Bind(wx.EVT_BUTTON, self.nudgeLeftDown)
	self.LeftDownBtn.SetToolTipString("Nudge left down")
	if (platform.system() == "Darwin"):
	    self.RightDownBtn = wx.BitmapButton(self.ManipulatePanel, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/CornerBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(150, 150), size=(30, 30))
	else:
	    self.RightDownBtn = wx.Button(self.ManipulatePanel, id=-1, label="", pos=(150, 150), size=(30, 30))
	    self.RightDownBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.NORMAL))
	    self.RightDownBtn.SetForegroundColour("#000000")
	self.RightDownBtn.Bind(wx.EVT_BUTTON, self.nudgeRightDown)
	self.RightDownBtn.SetToolTipString("Nudge right down")

	if (platform.system() == "Darwin"):
	    self.RotateBtn = wx.BitmapButton(self.ManipulatePanel, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/RotateBtn_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 180), size=(70, 25))
	else:
	    self.RotateBtn = wx.Button(self.ManipulatePanel, id=-1, label="Rotate", pos=(0, 180), size=(70, 25))
	    self.RotateBtn.SetForegroundColour("#FF0000")
	    self.RotateBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.NORMAL))
	self.RotateBtn.Bind(wx.EVT_BUTTON, self.toggleRotate)
	self.RotateBtn.SetToolTipString("Enter PyMOL object rotation mode")
	if (platform.system() == "Darwin"):
	    self.TranslateBtn = wx.BitmapButton(self.ManipulatePanel, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/TranslateBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(75, 180), size=(70, 25))
	else:
	    self.TranslateBtn = wx.Button(self.ManipulatePanel, id=-1, label="Translate", pos=(75, 180), size=(70, 25))
	    self.TranslateBtn.SetForegroundColour("#000000")
	    self.TranslateBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.NORMAL))
	self.TranslateBtn.Bind(wx.EVT_BUTTON, self.toggleTranslate)
	self.TranslateBtn.SetToolTipString("Enter PyMOL object translation mode")
	if (platform.system() == "Darwin"):
	    self.ManipulatePanelToggleBtn = wx.BitmapButton(self.ManipulatePanel, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/SelectPanelFlipBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(150, 180), size=(30, 25))
	else:
	    self.ManipulatePanelToggleBtn = wx.Button(self.ManipulatePanel, id=-1, label="+", pos=(150, 180), size=(30, 25))
	    self.ManipulatePanelToggleBtn.SetForegroundColour("#000000")
	    self.ManipulatePanelToggleBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.ManipulatePanelToggleBtn.SetToolTipString("Switch to selection mode")
	self.ManipulatePanelToggleBtn.Bind(wx.EVT_LEFT_DOWN, self.panelToggle)
	self.manipType = "Rotate"
	self.cen_x = 75
	self.cen_y = 75
	self.rotval = 5.0

	self.ManipulatePanel.SetBackgroundColour("#333333")
	self.ManipulatePanel.Hide()

	self.Bind(wx.EVT_ACTIVATE, self.focusEvent)
	self.SetupScrolling()

    def showHelp(self, event):
	# Open the help page
	if (platform.system() == "Darwin"):
	    try:
		browser = webbrowser.get("Safari")
	    except:
		print "Could not load Safari!  The help files are located at " + self.scriptdir + "/help"
		return
	    browser.open(self.parent.scriptdir + "/help/selection.html")
	else:
	    webbrowser.open(self.parent.scriptdir + "/help/selection.html")

    def focusEvent(self, event):
	self.inFocus = event.GetActive()
	if (self.inFocus):
	    self.cmd.enable("sele")

    def setPyMOL(self, pymol):
	self.pymol = pymol
	self.cmd = pymol.cmd

    def toggleRotate(self, event):
	self.manipType = "Rotate"
	if (platform.system() == "Darwin"):
	    self.RotateBtn.SetBitmapLabel(bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/RotateBtn_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
	    self.TranslateBtn.SetBitmapLabel(bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/TranslateBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
	else:
	    self.RotateBtn.SetForegroundColour("#FF0000")
	    self.TranslateBtn.SetForegroundColour("#000000")

    def toggleTranslate(self, event):
	self.manipType = "Translate"
	if (platform.system() == "Darwin"):
	    self.RotateBtn.SetBitmapLabel(bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/RotateBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
	    self.TranslateBtn.SetBitmapLabel(bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/TranslateBtn_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
	else:
	    self.RotateBtn.SetForegroundColour("#000000")
	    self.TranslateBtn.SetForegroundColour("#FF0000")

    def xyFocus(self, event):
	if (not(event.GetActive())):
	    self.leftDown = False
	event.Skip()

    def xyDown(self, event):
	self.manip_x = event.GetX()
	self.manip_y = event.GetY()
	# Convert to XY with origin
	self.manip_x = self.manip_x - 75
	self.manip_y = 75 - self.manip_y
	if (self.manip_x**2 + self.manip_y**2 > 75**2):
	    self.manip_x = -900
	    self.manip_y = -900
	self.cen_x = event.GetX()
	self.cen_y = event.GetY()
	self.rotval = 5.0
	self.leftDown = True
	if (self.manipType == "Rotate"):
	    if (self.manip_x >= -75):
		self.manip_z = math.sqrt(75**2 - self.manip_x**2 - self.manip_y**2)
	    try:
		self.cmd.select("sele", "seqsele")
		self.cmd.center("sele")
		self.cmd.disable("sele")
	    except:
		self.cmd.select("sele", "all")
		self.cmd.disable("sele")
	else:
	    try:
		self.cmd.select("sele", "seqsele")
		self.cmd.disable("sele")
	    except:
		self.cmd.select("sele", "all")
		self.cmd.disable("sele")
	event.Skip()

    def xyUp(self, event):
	self.leftDown = False
	event.Skip()

    def selectForNudge(self):
	if (self.manipType == "Rotate"):
	    try:
		self.cmd.select("sele", "seqsele")
		self.cmd.center("sele")
		self.cmd.disable("sele")
	    except:
		self.cmd.select("sele", "all")
		self.cmd.disable("sele")
	else:
	    try:
		self.cmd.select("sele", "seqsele")
		self.cmd.disable("sele")
	    except:
		self.cmd.select("sele", "all")
		self.cmd.disable("sele")

    def nudgeUp(self, event):
	self.selectForNudge()
	if (self.manipType == "Rotate"):
	    self.cmd.rotate([1, 0, 0], -5.0, "sele")
	else:
	    self.cmd.translate([0, 1, 0], "sele")
    def nudgeDown(self, event):
	self.selectForNudge()
	if (self.manipType == "Rotate"):
	    self.cmd.rotate([1, 0, 0], 5.0, "sele")
	else:
	    self.cmd.translate([0, -1, 0], "sele")
    def nudgeLeft(self, event):
	self.selectForNudge()
	if (self.manipType == "Rotate"):
	    self.cmd.rotate([0, 1, 0], -5.0, "sele")
	else:
	    self.cmd.translate([-1, 0, 0], "sele")
    def nudgeRight(self, event):
	self.selectForNudge()
	if (self.manipType == "Rotate"):
	    self.cmd.rotate([0, 1, 0], 5.0, "sele")
	else:
	    self.cmd.translate([1, 0, 0], "sele")
    def nudgeLeftUp(self, event):
	self.selectForNudge()
	if (self.manipType == "Rotate"):
	    self.cmd.rotate([1, 1, 0], -5.0, "sele")
	else:
	    self.cmd.translate([-1, 1, 0], "sele")
    def nudgeLeftDown(self, event):
	self.selectForNudge()
	if (self.manipType == "Rotate"):
	    self.cmd.rotate([-1, 1, 0], -5.0, "sele")
	else:
	    self.cmd.translate([-1, -1, 0], "sele")
    def nudgeRightUp(self, event):
	self.selectForNudge()
	if (self.manipType == "Rotate"):
	    self.cmd.rotate([-1, 1, 0], 5.0, "sele")
	else:
	    self.cmd.translate([1, 1, 0], "sele")
    def nudgeRightDown(self, event):
	self.selectForNudge()
	if (self.manipType == "Rotate"):
	    self.cmd.rotate([1, 1, 0], 5.0, "sele")
	else:
	    self.cmd.translate([1, -1, 0], "sele")

    def xyMotion(self, event):
	x = event.GetX()
	y = event.GetY()
	# Convert to XY with origin
	x = x - 75
	y = 75 - y
	if (self.leftDown and math.fabs(x-self.manip_x) + math.fabs(y-self.manip_y) >= 4):
	    if (self.manipType == "Rotate"):
		if (x**2 + y**2 > 75**2):
		    # Outside the circle, do nothing
		    event.Skip()
		    return
		elif (x**2 + y**2 <= 75**2 and self.manip_x < -200):
		    # The original click was outside of the circle, but now we're inside of it, so set the original
		    # click as the edge of the circle
		    self.manip_x = x
		    self.manip_y = y
		    self.manip_z = math.sqrt(75**2 - self.manip_x**2 - self.manip_y**2)
		    event.Skip()
		    return
		z = math.sqrt(75**2 - x**2 - y**2)
		#if (x-self.manip_x == 0):
		    # Rotate around x-axis
		#    vec = [1, 0, 0]
		#elif (y-self.manip_y == 0):
		    # Rotate around y-axis
		#    vec = [0, 1, 0]
		#else:
		#    m = float(y-self.manip_y) / float(x-self.manip_x)
		    #m = -1.0 / m
		#    vec = [1, m, 0]
		vec = list(numpy.cross([x, y, z], [self.manip_x, self.manip_y, self.manip_z]))
		if (x-self.manip_x >= 0 and y-self.manip_y >= 0):
		    self.rotval = 5.0
		elif (x-self.manip_x < 0 and y-self.manip_y >= 0):
		    self.rotval = 5.0
		elif (x-self.manip_x < 0 and y-self.manip_y < 0):
		    self.rotval = -5.0
		else:
		    self.rotval = -5.0
		#if ((x-self.cen_x)**2 + (y-self.cen_y)**2 + 5 < (self.manip_x-self.cen_x)**2 + (self.manip_y-self.cen_y)**2):
		    #self.rotval = self.rotval * -1.0
		    #self.cen_x = x
		    #self.cen_y = y
		self.cmd.rotate(vec, self.rotval, "sele")
		self.manip_z = z
	    else:
		vec = [(x-self.manip_x)/2.0, (y-self.manip_y)/2.0, 0]
		self.cmd.translate(vec, "sele")
	    self.manip_x = x
	    self.manip_y = y
	event.Skip()

    def panelToggle(self, event):
	if (self.panelMode == "Selection"):
	    if (platform.system() == "Darwin"):
		self.labelSelection.SetBitmap(bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/labelManipulation.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
	    else:
		self.labelSelection.SetLabel("Manipulation")
	    self.panelMode = "Manipulation"
	    self.SelectPanel.Hide()
	    self.ManipulatePanel.Show()
	else:
	    if (platform.system() == "Darwin"):
		self.labelSelection.SetBitmap(bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/labelSelection.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
	    else:
		self.labelSelection.SetLabel("Selection")
	    self.panelMode = "Selection"
	    self.SelectPanel.Show()
	    self.ManipulatePanel.Hide()
	if (platform.system() == "Linux"):
	    resizeTextControlForUNIX(self.labelSelection, 165, 160)

    def atomOff(self, event):
	try:
	    self.cmd.hide("sticks", "seqsele")
	    self.cmd.hide("spheres", "seqsele")
	    self.cmd.hide("lines", "seqsele")
	except:
	    self.cmd.hide("sticks", "all")
	    self.cmd.hide("spheres", "all")
	    self.cmd.hide("lines", "all")
	logInfo("Turned atom display off")

    def atomLines(self, event):
	if (self.seqWin.protocol_view_active):
	    bonus = " and protocol_view"
	else:
	    bonus = ""
	try:
	    self.cmd.hide("spheres", "seqsele" + bonus)
	    self.cmd.show("sticks", "seqsele" + bonus)
	    self.cmd.set_bond("stick_radius", 0.1, "seqsele" + bonus)
	    self.cmd.set_bond("stick_transparency", 0, "seqsele" + bonus)
	except:
	    self.cmd.hide("spheres", "all" + bonus)
	    self.cmd.show("sticks", "all" + bonus)
	    self.cmd.set_bond("stick_radius", 0.1, "all" + bonus)
	    self.cmd.set_bond("stick_transparency", 0, "all" + bonus)
	logInfo("Turned atom line display on")

    def atomSticks(self, event):
	if (self.seqWin.protocol_view_active):
	    bonus = " and protocol_view"
	else:
	    bonus = ""
	try:
	    self.cmd.show("sticks", "seqsele" + bonus)
	    self.cmd.set_bond("stick_radius", 0.25, "seqsele" + bonus)
	    self.cmd.set_bond("stick_transparency", 0, "seqsele" + bonus)
	    self.cmd.hide("spheres", "seqsele" + bonus)
	except:
	    self.cmd.show("sticks", "all" + bonus)
	    self.cmd.set_bond("stick_radius", 0.25, "all" + bonus)
	    self.cmd.set_bond("stick_transparency", 0, "all" + bonus)
	    self.cmd.hide("spheres", "all" + bonus)
	logInfo("Turned atom stick display on")
	#self.cmd.hide("lines", "sele")

    def atomBallsLines(self, event):
	if (self.seqWin.protocol_view_active):
	    bonus = " and protocol_view"
	else:
	    bonus = ""
	try:
	    self.cmd.show("spheres", "seqsele" + bonus)
	    self.cmd.show("sticks", "seqsele" + bonus)
	    self.cmd.set_bond("stick_radius", 0.1, "seqsele" + bonus)
	    self.cmd.set_bond("stick_transparency", 0, "seqsele" + bonus)
	    self.cmd.set("sphere_scale", 0.15, "seqsele" + bonus)
	    self.cmd.set("sphere_transparency", 0, "seqsele" + bonus)
	except:
	    self.cmd.show("spheres", "all" + bonus)
	    self.cmd.show("sticks", "all" + bonus)
	    self.cmd.set_bond("stick_radius", 0.1, "all" + bonus)
	    self.cmd.set_bond("stick_transparency", 0, "all" + bonus)
	    self.cmd.set("sphere_scale", 0.15, "all" + bonus)
	    self.cmd.set("sphere_transparency", 0, "all" + bonus)
	logInfo("Turned atom balls and lines display on")

    def atomBallsSticks(self, event):
	if (self.seqWin.protocol_view_active):
	    bonus = " and protocol_view"
	else:
	    bonus = ""
	try:
	    self.cmd.show("spheres", "seqsele" + bonus)
	    self.cmd.show("sticks", "seqsele" + bonus)
	    self.cmd.set_bond("stick_radius", 0.25, "seqsele" + bonus)
	    self.cmd.set_bond("stick_transparency", 0, "seqsele" + bonus)
	    self.cmd.set("sphere_scale", 0.3, "seqsele" + bonus)
	    self.cmd.set("sphere_transparency", 0, "seqsele" + bonus)
	except:
	    self.cmd.show("spheres", "all" + bonus)
	    self.cmd.show("sticks", "all" + bonus)
	    self.cmd.set_bond("stick_radius", 0.25, "all" + bonus)
	    self.cmd.set_bond("stick_transparency", 0, "all" + bonus)
	    self.cmd.set("sphere_scale", 0.3, "all" + bonus)
	    self.cmd.set("sphere_transparency", 0, "all" + bonus)
	logInfo("Turned atom balls and sticks display on")

    def atomVDW(self, event):
	if (self.seqWin.protocol_view_active):
	    bonus = " and protocol_view"
	else:
	    bonus = ""
	try:
	    self.cmd.hide("sticks", "seqsele" + bonus)
	    self.cmd.show("spheres", "seqsele" + bonus)
	    # Sphere transparency uses "set", stick transparency uses "set_bond" :/
	    self.cmd.set("sphere_transparency", 0, "seqsele" + bonus)
	    self.cmd.set("sphere_scale", 1, "seqsele" + bonus)
	except:
	    self.cmd.hide("sticks", "all" + bonus)
	    self.cmd.show("spheres", "all" + bonus)
	    # Sphere transparency uses "set", stick transparency uses "set_bond" :/
	    self.cmd.set("sphere_transparency", 0, "all" + bonus)
	    self.cmd.set("sphere_scale", 1, "all" + bonus)
	logInfo("Turned atom VDW display on")

    def atomRecolor(self, event):
	logInfo("Click on the atom recolor button")
	if (self.seqWin.protocol_view_active):
	    bonus = " and protocol_view"
	else:
	    bonus = ""
	dlg = wx.ColourDialog(self)
	dlg.GetColourData().SetChooseFull(True)
	if (dlg.ShowModal() == wx.ID_OK):
	    data = dlg.GetColourData()
	    mycolor = "0x%02x%02x%02x" % data.GetColour().Get()
	    try:
		self.cmd.color(mycolor, "seqsele and symbol c" + bonus)
		try:
		    # Save the fact that these residues are not colored by chain in case the colors get moved around later
		    # (i.e. chains were deleted and other chains were moved up)
		    self.cmd.select("atomcolorsele", "atomcolorsele and not seqsele" + bonus)
		except:
		    pass
	    except:
		self.cmd.color(mycolor, "all and symbol c" + bonus)
		self.cmd.delete("atomcolorsele")
	    logInfo("Set atom colors to " + mycolor)
	dlg.Destroy()

    def atomStandardColor(self, event):
	if (self.seqWin.protocol_view_active):
	    bonus = " and protocol_view"
	else:
	    bonus = ""
	try:
	    self.cmd.select("temp", "symbol c in seqsele" + bonus)
	    self.cmd.color("gray", "temp")
	    self.cmd.select("temp", "symbol n in seqsele" + bonus)
	    self.cmd.color("blue", "temp")
	    self.cmd.select("temp", "symbol o in seqsele" + bonus)
	    self.cmd.color("red", "temp")
	    self.cmd.select("temp", "symbol s in seqsele" + bonus)
	    self.cmd.color("yellow", "temp")
	    self.cmd.select("temp", "symbol p in seqsele" + bonus)
	    self.cmd.color("violet", "temp")
	    self.cmd.select("temp", "symbol h in seqsele" + bonus)
	    self.cmd.color("white", "temp")
	    self.cmd.delete("temp")
	    try:
		# Save the fact that these residues are not colored by chain in case the colors get moved around later
		# (i.e. chains were deleted and other chains were moved up)
		self.cmd.select("atomcolorsele", "atomcolorsele and not seqsele" + bonus)
	    except:
		pass
	    self.cmd.enable("seqsele")
	except:
	    self.cmd.select("temp", "symbol c" + bonus)
	    self.cmd.color("gray", "temp")
	    self.cmd.select("temp", "symbol n" + bonus)
	    self.cmd.color("blue", "temp")
	    self.cmd.select("temp", "symbol o" + bonus)
	    self.cmd.color("red", "temp")
	    self.cmd.select("temp", "symbol s" + bonus)
	    self.cmd.color("yellow", "temp")
	    self.cmd.select("temp", "symbol p" + bonus)
	    self.cmd.color("violet", "temp")
	    self.cmd.select("temp", "symbol h" + bonus)
	    self.cmd.color("white", "temp")
	    self.cmd.delete("temp")
	    self.cmd.delete("atomcolorsele")
	logInfo("Turned atom coloring to standard")

    def atomChainColor(self, event):
	if (self.seqWin.protocol_view_active):
	    bonus = " and protocol_view"
	else:
	    bonus = ""
	# Color the selected atoms according to chain where the color of the chain
	# is the result of a hash of the model+chain string
	self.pymol.stored.selected = []
	try:
	    self.cmd.select("temp", "seqsele")
	    self.cmd.iterate_state(1, "seqsele" + bonus, "stored.selected.append(model+\"|\"+chain)")
	    # Get the unique set of model|chain pairs
	    self.pymol.stored.selected.sort()
	    IDs = list(set(self.pymol.stored.selected))
	    selection = "seqsele"
	except:
	    IDs = []
	    for ID in self.seqWin.IDs:
		modelchain = ID[0:len(ID)-2] + " " + ID[len(ID)-1]
		if (not(modelchain) in IDs):
		    IDs.append(modelchain)
	    selection = "all"
	for modelchain in IDs:
	    if (modelchain[len(modelchain)-1] == " " or modelchain[len(modelchain)-1] == "|"):
		modelchain = modelchain.strip() + "_"
	    model = modelchain.split("|")[0]
	    try:
		chain = modelchain.split("|")[1]
	    except:
		chain = "_"
	    #x = len(model) / 3
	    #r = hash(model[0:x] + chain) % 256
	    #g = hash(model[x:2*x] + chain) % 256
	    #b = hash(model[2*x:] + chain) % 256
	    #color = "0x%02x%02x%02x" % (r, g, b)
	    row = self.seqWin.IDs.index(modelchain)
	    color = getChainColor(row)
	    if (chain != "_" and chain != " " and chain != ""):
		self.cmd.select("temp", "model " + model + " and chain " + chain + " and " + selection + " and symbol c" + bonus)
	    else:
		self.cmd.select("temp", "model " + model + " and " + selection + " and symbol c" + bonus)
	    self.cmd.color(color, "temp")
	try:
	    # Save the fact that these residues are colored by chain in case the colors get moved around later
	    # (i.e. chains were deleted and other chains were moved up)
	    self.cmd.select("atomcolorsele", selection + " or atomcolorsele")
	except:
	    self.cmd.select("atomcolorsele", selection)
	try:
	    self.cmd.delete("temp")
	except:
	    pass
	if (selection == "seqsele"):
	    self.cmd.enable("seqsele")
	logInfo("Colored atoms by chain")

    def atomTermini(self, event):
	if (self.seqWin.protocol_view_active):
	    bonus = " and protocol_view"
	else:
	    bonus = ""
	# Color the atoms such that residues become more blue or red as they approach the
	# N and C termini respectively
	topLefts = self.seqWin.SeqViewer.GetSelectionBlockTopLeft()
	bottomRights = self.seqWin.SeqViewer.GetSelectionBlockBottomRight()
	if (len(topLefts) == 0):
	    for r in range(0, len(self.seqWin.IDs)):
		for c in range(0, len(self.seqWin.sequences[r])):
		    if (self.seqWin.sequences[r][c] == "-"):
			continue
		    modelchain = self.seqWin.IDs[r]
		    model = modelchain[0:len(modelchain)-2]
		    chain = modelchain[len(modelchain)-1]
		    indx = str(self.seqWin.indxToSeqPos[r][c][1])
		    if (chain == "_"):
			sel = "model " + model + " and resi " + indx + " and symbol c" + bonus
		    else:
			sel = "model " + model + " and chain " + chain + " and resi " + indx + " and symbol c" + bonus
		    pos = float(c) / float(len(self.seqWin.sequences[r]))
		    if (pos < 0.5):
			blue = 255
			red = 510 * pos
			green = red
		    else:
			red = 255
			blue = 510 * (1.0 - pos)
			green = blue
		    color = "0x%02x%02x%02x" % (red, green, blue)
		    self.cmd.color(color, sel)
	else:
	    for i in range(0, len(topLefts)):
		for r in range(topLefts[i][0], bottomRights[i][0]+1):
		    for c in range(topLefts[i][1], bottomRights[i][1]+1):
			if (self.seqWin.sequences[r][c] == "-"):
			    continue
			modelchain = self.seqWin.IDs[r]
			model = modelchain[0:len(modelchain)-2]
			chain = modelchain[len(modelchain)-1]
			indx = str(self.seqWin.indxToSeqPos[r][c][1])
			if (chain == "_"):
			    sel = "model " + model + " and resi " + indx + bonus
			else:
			    sel = "model " + model + " and chain " + chain + " and resi " + indx + bonus
			pos = float(c) / float(len(self.seqWin.sequences[r]))
			if (pos < 0.5):
			    blue = 255
			    red = 510 * pos
			    green = red
			else:
			    red = 255
			    blue = 510 * (1.0 - pos)
			    green = blue
			color = "0x%02x%02x%02x" % (red, green, blue)
			self.cmd.color(color, sel)

    def labelsOff(self, event):
	try:
	    self.cmd.label("seqsele", "\"\"")
	except:
	    self.cmd.label("all", "\"\"")
	logInfo("Turned atom VDW display on")

    def labelResidues(self, event):
	if (self.seqWin.protocol_view_active):
	    bonus = " and protocol_view"
	else:
	    bonus = ""
	self.pymol.stored.selected = []
	try:
	    self.cmd.iterate_state(1, "seqsele" + bonus, "stored.selected.append(model+\" \"+chain+\" \"+resi+\" \"+resn+\" \"+name)")
	    # Get the unique set of model|chain pairs
	    self.pymol.stored.selected.sort()
	    IDs = list(set(self.pymol.stored.selected))
	    IDs.sort()
	    uniqIDs = []
	    # This is needed in case you are trying to label a residue that does not have a CA atom
	    lastmodel = ""
	    lastchain = ""
	    lastresi = ""
	    # What does this code do?
	    # It is going through all the possible atoms that can be labeled on each selected residue
	    # We only want to label one atom, and the CA is the most likely candidate
	    # But we may have NCAAs that do not have an atom called "CA"
	    # So we'll take whatever the alphabetically first atom is as the atom to label, and if we later
	    # encounter a CA atom on the same residue we'll use that instead
	    # IDs needs to be sorted before doing this because the algorithm assumes that all atoms for a single
	    # residue are clustered together
	    for i in range(0, len(IDs)):
		ID = IDs[i]
		if (len(ID.split()) == 4):
		    model = ID.split()[0]
		    resi = ID.split()[1]
		    resn = ID.split()[2]
		    name = ID.split()[3]
		    if (len(uniqIDs) == 0):
			uniqIDs.append(ID)
		    elif (lastmodel == model and lastresi == resi and name.upper() == "CA"):
			uniqIDs[len(uniqIDs)-1] = ID
		    elif (lastmodel != model or lastresi != resi):
			uniqIDs.append(ID)
		    lastmodel = model
		    lastresi = resi
		else:
		    model = ID.split()[0]
		    chain = ID.split()[1]
		    resi = ID.split()[2]
		    resn = ID.split()[3]
		    name = ID.split()[4]
		    if (len(uniqIDs) == 0):
			uniqIDs.append(ID)
		    elif (lastmodel == model and lastresi == resi and name.upper() == "CA"):
			uniqIDs[len(uniqIDs)-1] = ID
		    elif (lastmodel != model or lastresi != resi):
			uniqIDs.append(ID)
		    lastmodel = model
		    lastchain = chain
		    lastresi = resi
	    for ID in uniqIDs:
		if (len(ID.split()) == 4):
		    model = ID.split()[0]
		    resi = ID.split()[1]
		    resn = ID.split()[2]
		    name = ID.split()[3]
		    self.cmd.select("labelsele", "model " + model + " and resi " + str(resi) + " and name " + name)
		    self.cmd.label("labelsele", "\"" + resn + resi + "\"")
		else:
		    model = ID.split()[0]
		    chain = ID.split()[1]
		    resi = ID.split()[2]
		    resn = ID.split()[3]
		    name = ID.split()[4]
		    self.cmd.select("labelsele", "model " + model + " and chain " + str(chain) + " and resi " + str(resi) + " and name " + name)
		    self.cmd.label("labelsele", "\"" + resn + resi + "\"")
	    self.cmd.delete("labelsele")
	    self.cmd.enable("seqsele")
	except:
	    pass
	logInfo("Displayed labels on the current selection")

    def ribbonsOff(self, event):
	try:
	    self.cmd.hide("ribbon", "seqsele")
	    self.cmd.hide("cartoon", "seqsele")
	except:
	    self.cmd.hide("ribbon", "all")
	    self.cmd.hide("cartoon", "all")
	logInfo("Turned ribbon displays off")

    def ribbonsCartoon(self, event):
	if (self.seqWin.protocol_view_active):
	    bonus = " and protocol_view"
	else:
	    bonus = ""
	try:
	    self.cmd.hide("ribbon", "seqsele" + bonus)
	    self.cmd.show("cartoon", "seqsele" + bonus)
	except:
	    self.cmd.hide("ribbon", "all" + bonus)
	    self.cmd.show("cartoon", "all" + bonus)
	logInfo("Turned cartoon display on")

    def ribbonsRibbon(self, event):
	if (self.seqWin.protocol_view_active):
	    bonus = " and protocol_view"
	else:
	    bonus = ""
	try:
	    self.cmd.show("ribbon", "seqsele" + bonus)
	    self.cmd.hide("cartoon", "seqsele" + bonus)
	except:
	    self.cmd.show("ribbon", "all")
	    self.cmd.hide("cartoon", "all")
	logInfo("Turned string display on")

    def ribbonRecolor(self, event):
	logInfo("Clicked on the ribbon recoloring button")
	if (self.seqWin.protocol_view_active):
	    bonus = " and protocol_view"
	else:
	    bonus = ""
	dlg = wx.ColourDialog(self)
	dlg.GetColourData().SetChooseFull(True)
	if (dlg.ShowModal() == wx.ID_OK):
	    data = dlg.GetColourData()
	    mycolor = "0x%02x%02x%02x" % data.GetColour().Get()
	    try:
		self.cmd.set("ribbon_color", mycolor, "seqsele" + bonus)
		self.cmd.set("cartoon_color", mycolor, "seqsele" + bonus)
		try:
		    # Save the fact that these residues are not colored by chain in case the colors get moved around later
		    # (i.e. chains were deleted and other chains were moved up)
		    self.cmd.select("chaincolorsele", "chaincolorsele and not seqsele")
		except:
		    pass
	    except:
		self.cmd.set("ribbon_color", mycolor, "resi 1-9999")
		self.cmd.set("cartoon_color", mycolor, "resi 1-9999")
		try:
		    self.cmd.delete("chaincolorsele")
		except:
		    pass
	    logInfo("Set ribbon color to " + mycolor)
	dlg.Destroy()

    def ribbonStandardColor(self, event):
	if (self.seqWin.protocol_view_active):
	    bonus = " and protocol_view"
	else:
	    bonus = ""
	try:
	    temp = "seqsele" + bonus
	    #self.cmd.select("temp", "seqsele")
	    self.cmd.set("ribbon_color", "white", temp)
	    self.cmd.set("cartoon_color", "white", temp)
	    #self.cmd.select("temp", "ss s in seqsele")
	    temp = "ss s in seqsele" + bonus
	    self.cmd.set("ribbon_color", "yellow", temp)
	    self.cmd.set("cartoon_color", "yellow", temp)
	    #self.cmd.select("temp", "ss h in seqsele")
	    temp = "ss h in seqsele" + bonus
	    self.cmd.set("ribbon_color", "red", temp)
	    self.cmd.set("cartoon_color", "red", temp)
	    #self.cmd.select("temp", "(ss b or ss t) and seqsele")
	    temp = "(ss b or ss t) and seqsele" + bonus
	    self.cmd.set("ribbon_color", "blue", temp)
	    self.cmd.set("cartoon_color", "blue", temp)
	    #self.cmd.select("temp", "(ss g or ss i) and seqsele")
	    temp = "(ss g or ss i) and seqsele" + bonus
	    self.cmd.set("ribbon_color", "orange", temp)
	    self.cmd.set("cartoon_color", "orange", temp)
         #DNA
         #Adenosine
	    logInfo('selection.py 1200')
	    logInfo('Adenosine coloring')
	    temp = '(resn ADE) and seqsele'+bonus
	    self.cmd.set('ribbon_color','green',temp)
	    self.cmd.set('cartoon_color','green',temp)
	    self.cmd.color('green',temp)
         #Thymine
	    logInfo('Thymine coloring')
	    temp = '(resn THY) and seqsele'+bonus
	    self.cmd.set('ribbon_color','red',temp)
	    self.cmd.set('cartoon_color','red',temp)
	    self.cmd.color('red',temp)
         #Cytosine
	    logInfo('Cytosine coloring')
	    temp = '(resn CYT) and seqsele'+bonus
	    self.cmd.set('ribbon_color','blue',temp)
	    self.cmd.set('cartoon_color','blue',temp)
	    self.cmd.color('blue',temp)
         #Guanine
	    logInfo('Guanine coloring')
	    temp = '(resn GUA) and seqsele'+bonus
	    self.cmd.set('ribbon_color','gray',temp)
	    self.cmd.set('cartoon_color','gray',temp)
	    self.cmd.color('gray',temp)
	    try:
		# Save the fact that these residues are not colored by chain in case the colors get moved around later
		# (i.e. chains were deleted and other chains were moved up)
		self.cmd.select("chaincolorsele", "chaincolorsele and not seqsele")
	    except:
		pass
	    self.cmd.delete("temp")
	    self.cmd.select("seqsele", "seqsele")
	    self.cmd.enable("seqsele")
	except:
	    logInfo('defaultPyMOLView')
	    defaultPyMOLView(self.cmd)
	    try:
		self.cmd.delete("colorchainsele")
	    except:
		pass
	logInfo("Reset ribbon coloring to standard")

    def recolorSavedChainColors(self):
	# This function is called whenever a chain is deleted to update chain-colored residues with their new
	# colors according to the colored strip in the sequence viewer
	try:
	    self.cmd.select("temp", "seqsele")
	except:
	    pass
	try:
	    self.cmd.select("seqsele", "atomcolorsele")
	    self.atomChainColor(None)
	    self.cmd.disable("atomcolorsele")
	except:
	    pass
	try:
	    self.cmd.select("seqsele", "chaincolorsele")
	    self.ribbonChainColor(None)
	    self.cmd.disable("chaincolorsele")
	except:
	    pass
	self.cmd.disable("seqsele")

    def ribbonChainColor(self, event):
	if (self.seqWin.protocol_view_active):
	    bonus = " and protocol_view"
	else:
	    bonus = ""
	# Color the selected ribbons according to chain where the color of the chain
	# is the result of a hash of the model+chain string
	self.pymol.stored.selected = []
	try:
	    self.cmd.select("temp", "seqsele")
	    self.cmd.iterate_state(1, "seqsele" + bonus, "stored.selected.append(model+\"|\"+chain)")
	    # Get the unique set of model|chain pairs
	    self.pymol.stored.selected.sort()
	    IDs = list(set(self.pymol.stored.selected))
	    selection = "seqsele"
	except:
	    IDs = []
	    for ID in self.seqWin.IDs:
		modelchain = ID[0:len(ID)-2] + "|" + ID[len(ID)-1]
		if (not(modelchain) in IDs):
		    IDs.append(modelchain)
	    selection = "all" + bonus
	for modelchain in IDs:
	    if (modelchain[len(modelchain)-1] == " " or modelchain[len(modelchain)-1] == "|"):
		modelchain = modelchain.strip() + "_"
	    model = modelchain.split("|")[0]
	    try:
		chain = modelchain.split("|")[1]
	    except:
		chain = "_"
	    #x = len(model) / 3
	    ##h = hash(model + chain)
	    ##r = h & 0xFF0000 >> 16
	    ##g = h & 0x00FF00 >> 8
	    ##b = h & 0x0000FF
	    #r = hash(model[0:x] + chain) % 256
	    #g = hash(model[x:2*x] + chain) % 256
	    #b = hash(model[2*x:] + chain) % 256
	    #color = "0x%02x%02x%02x" % (r, g, b)
	    row = self.seqWin.IDs.index(modelchain)
	    color = getChainColor(row)
	    if (chain != "_" and chain != " " and chain != ""):
		self.cmd.select("temp", "model " + model + " and chain " + chain + " and " + selection + bonus)
	    else:
		self.cmd.select("temp", "model " + model + " and " + selection + bonus)
	    self.cmd.set("ribbon_color", color, "temp")
	    self.cmd.set("cartoon_color", color, "temp")
	try:
	    try:
		# Save the fact that these residues are colored by chain in case the colors get moved around later
		# (i.e. chains were deleted and other chains were moved up)
		self.cmd.select("chaincolorsele", selection + " or chaincolorsele")
	    except:
		self.cmd.select("chaincolorsele", selection)
	    self.cmd.disable("chaincolorsele")
	    self.cmd.delete("temp")
	except:
	    pass
	if (selection == "seqsele"):
	    self.cmd.enable("seqsele")
	logInfo("Colored ribbons by chain")

    def ribbonTermini(self, event):
	if (self.seqWin.protocol_view_active):
	    bonus = " and protocol_view"
	else:
	    bonus = ""
	# Color the ribbons such that residues become more blue or red as they approach the
	# N and C termini respectively
	topLefts = self.seqWin.SeqViewer.GetSelectionBlockTopLeft()
	bottomRights = self.seqWin.SeqViewer.GetSelectionBlockBottomRight()
	if (len(topLefts) == 0):
	    for r in range(0, len(self.seqWin.IDs)):
		for c in range(0, len(self.seqWin.sequences[r])):
		    if (self.seqWin.sequences[r][c] == "-"):
			continue
		    modelchain = self.seqWin.IDs[r]
		    model = modelchain[0:len(modelchain)-2]
		    chain = modelchain[len(modelchain)-1]
		    indx = str(self.seqWin.indxToSeqPos[r][c][1])
		    if (chain == "_"):
			sel = "model " + model + " and resi " + indx + bonus
		    else:
			sel = "model " + model + " and chain " + chain + " and resi " + indx + bonus
		    pos = float(c) / float(len(self.seqWin.sequences[r]))
		    if (pos < 0.5):
			blue = 255
			red = 510 * pos
			green = red
		    else:
			red = 255
			blue = 510 * (1.0 - pos)
			green = blue
		    color = "0x%02x%02x%02x" % (red, green, blue)
		    self.cmd.set("cartoon_color", color, sel)
	else:
	    for i in range(0, len(topLefts)):
		for r in range(topLefts[i][0], bottomRights[i][0]+1):
		    for c in range(topLefts[i][1], bottomRights[i][1]+1):
			if (self.seqWin.sequences[r][c] == "-"):
			    continue
			modelchain = self.seqWin.IDs[r]
			model = modelchain[0:len(modelchain)-2]
			chain = modelchain[len(modelchain)-1]
			indx = str(self.seqWin.indxToSeqPos[r][c][1])
			if (chain == "_"):
			    sel = "model " + model + " and resi " + indx + bonus
			else:
			    sel = "model " + model + " and chain " + chain + " and resi " + indx + bonus
			pos = float(c) / float(len(self.seqWin.sequences[r]))
			if (pos < 0.5):
			    blue = 255
			    red = 510 * pos
			    green = red
			else:
			    red = 255
			    blue = 510 * (1.0 - pos)
			    green = blue
			color = "0x%02x%02x%02x" % (red, green, blue)
			self.cmd.set("cartoon_color", color, sel)

    def displaySurfaces(self):
	# Use the selection information to display pre-configured surfaces from the molecular surfaces protocol
	if (self.parent.Protocols.currentProtocol == "Molecular Surfaces"):
	    return
	self.cmd.flag("ignore", "all", "clear")
	for name in self.cmd.get_names("selections"):
	    if (name.startswith("surf_recp_")):
		try:
		    self.cmd.flag("ignore", "surf_lig_" + name[10:], "set")
		except:
		    pass
		self.cmd.show("surface", "surf_recp_" + name[10:])

    def toggleSurf(self, event):
	if (self.showSurf):
	    self.showSurf = False
	    self.ToggleSurfBtn.SetToolTipString("Display configured surfaces")
	    if (platform.system() == "Darwin"):
		self.ToggleSurfBtn.SetBitmapLabel(bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/ToggleSurfBtn_Off.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
	    else:
		self.ToggleSurfBtn.SetLabel("Surf Off")
	    self.cmd.hide("surface", "all")
	else:
	    self.showSurf = True
	    self.ToggleSurfBtn.SetToolTipString("Hide configured surfaces")
	    if (platform.system() == "Darwin"):
		self.ToggleSurfBtn.SetBitmapLabel(bitmap=wx.Image(self.parent.scriptdir + "/images/osx/selection/ToggleSurfBtn_On.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
	    else:
		self.ToggleSurfBtn.SetLabel("Surf On")
	    self.displaySurfaces()

    def selectAll(self, event):
	try:
	    self.cmd.delete("sele")
	except:
	    pass
	self.cmd.select("seqsele", "all")
	self.seqWin.selectUpdate(False)
	logInfo("Selected all atoms")

    def selectInvert(self, event):
	self.cmd.select("seqsele", "not seqsele")
	self.cmd.enable("seqsele")
	self.seqWin.selectUpdate(False)
	logInfo("Inverted current selection")

    def selectVisible(self, event):
	self.cmd.select("seqsele", "rep lines or rep sticks or rep spheres")
	self.cmd.enable("seqsele")
	self.seqWin.selectUpdate(False)
	logInfo("Clicked on the visible select button")

    def selectH(self, event):
	try:
	    self.cmd.select("seqsele", "byres seqsele") # Expand out to the whole residues first
	    self.cmd.select("seqsele", "seqsele and symbol h")
	except:
	    self.cmd.delete("seqsele")
	    self.cmd.select("seqsele", "symbol h")
	self.seqWin.selectUpdate(False)
	logInfo("Selected all hydrogen atoms")

    def selectC(self, event):
	try:
	    self.cmd.select("seqsele", "byres seqsele") # Expand out to the whole residues first
	    self.cmd.select("seqsele", "seqsele and symbol c")
	except:
	    self.cmd.delete("seqsele")
	    self.cmd.select("seqsele", "symbol c")
	self.seqWin.selectUpdate(False)
	logInfo("Selected all carbon atoms")

    def selectN(self, event):
	try:
	    self.cmd.select("seqsele", "byres seqsele") # Expand out to the whole residues first
	    self.cmd.select("seqsele", "seqsele and symbol n")
	except:
	    self.cmd.delete("seqsele")
	    self.cmd.select("seqsele", "symbol n")
	self.seqWin.selectUpdate(False)
	logInfo("Selected all nitrogen atoms")

    def selectO(self, event):
	try:
	    self.cmd.select("seqsele", "byres seqsele") # Expand out to the whole residues first
	    self.cmd.select("seqsele", "seqsele and symbol o")
	except:
	    self.cmd.delete("seqsele")
	    self.cmd.select("seqsele", "symbol o")
	self.seqWin.selectUpdate(False)
	logInfo("Selected all oxygen atoms")

    def selectSolvent(self, event):
	try:
	    self.cmd.select("seqsele", "byres seqsele") # Expand out to the whole residues first
	    self.cmd.select("seqsele", "seqsele and solvent")
	except:
	    self.cmd.delete("seqsele")
	    self.cmd.select("seqsele", "solvent")
	self.seqWin.selectUpdate(False)
	logInfo("Selected all solvent atoms")

    def selectBB(self, event):
	try:
	    self.cmd.select("seqsele", "byres seqsele") # Expand out to the whole residues first
	    self.cmd.select("seqsele", "seqsele and backbone")
	except:
	    self.cmd.delete("seqsele")
	    self.cmd.select("seqsele", "backbone")
	self.seqWin.selectUpdate(False)
	logInfo("Selected all BB atoms")

    def selectSC(self, event):
	try:
	    self.cmd.select("seqsele", "byres seqsele") # Expand out to the whole residues first
	    self.cmd.select("seqsele", "seqsele and sidechain")
	except:
	    self.cmd.delete("seqsele")
	    self.cmd.select("seqsele", "sidechain")
	self.seqWin.selectUpdate(False)
	logInfo("Selected all sidechain atoms")

    def selectExtend(self, event):
	try:
	    radius = float(self.ExtendValueTxt.GetValue())
	except:
	    return
	if (radius <= 0):
	    return
	self.cmd.select("seqsele", "all within " + str(radius) + " of seqsele")
	self.seqWin.selectUpdate(False)
	logInfo("Extended current selection by a radius of " + str(radius))

    def zoomSelection(self, event):
	try:
	    self.cmd.zoom("seqsele", 0, 0, 0, 2)
	except:
	    self.cmd.zoom("all", 0, 0, 0, 2)
	logInfo("Clicked the zoom button")

    def centerSelection(self, event):
	try:
	    self.cmd.center("seqsele")
	except:
	    self.cmd.center("all")
	logInfo("Clicked the center button")

    def neighborhoodView(self, event):
	try:
	    radius = float(self.NeighborhoodValueTxt.GetValue())
	except:
	    return
	if (radius <= 0):
	    return
	self.cmd.select("orig", "all in seqsele")
	self.cmd.select("seqsele", "all within " + str(radius) + " of seqsele")
	self.cmd.hide("all")
	self.cmd.hide("spheres", "seqsele")
	self.cmd.show("sticks", "seqsele")
	self.cmd.set_bond("stick_radius", 0.1, "seqsele")
	self.cmd.set_bond("stick_transparency", 0, "seqsele")
	self.cmd.show("cartoon", "seqsele")
	self.cmd.hide("spheres", "orig")
	self.cmd.show("sticks", "orig")
	self.cmd.set_bond("stick_radius", 0.25, "orig")
	self.cmd.delete("orig")
	self.cmd.enable("seqsele") # So you can see which ones are selected clearly
	self.cmd.zoom("seqsele", 0, 0, 0, 2)
	self.seqWin.selectUpdate(False)
	logInfo("Clicked the neighborhood button with a radius of " + str(radius))

    def setProtPanel(self, protPanel):
	self.protPanel = protPanel

    def setSeqWin(self, seqWin):
	self.seqWin = seqWin

    def changeSF(self, event):
	logInfo("Clicked the change scorefxn button")
	if (platform.system() == "Windows"):
	    wtsdir = os.getenv("PYROSETTA_DATABASE") + "\\scoring\\weights"
	else:
	    wtsdir = os.getenv("PYROSETTA_DATABASE") + "/scoring/weights"
	dlg = wx.FileDialog(
	    self, message="Choose a Weights File",
	    defaultDir=wtsdir,
	    defaultFile=self.SelectScorefxnBtn.GetLabel() + ".wts",
	    wildcard="Weights Files (*.wts)|*.wts",
	    style=wx.OPEN)
	if (dlg.ShowModal() == wx.ID_OK):
	    paths = dlg.GetPaths()
	    filename = str(paths[0])
	    if (platform.system() == "Windows"):
		lastDirIndx = filename.rfind("\\")
	    else:
		lastDirIndx = filename.rfind("/")
	    newlabel = filename[lastDirIndx+1:]
	    newlabel = newlabel.split(".wts")[0]
	    # Try to load the scoring function
	    try:
		self.SelectScorefxnBtn.SetLabel(newlabel)
		#self.scorefxn = ScoreFunction()
		#self.scorefxn.add_weights_from_file(filename)
		self.weightsfile = filename
		logInfo("Loaded a new weights file with the following data", filename)
		#if (len(self.scorefxn.get_nonzero_weighted_scoretypes()) == 0):
		#    raise Exception
		#self.protPanel.activate()
	    except:
		msg = "The file " + filename.strip() + " is invalid."
		wx.MessageBox(msg, "File Cannot Be Read", wx.OK|wx.ICON_EXCLAMATION)
		#filevalid = False
	    #if (filevalid):
		# Thread this load operation because it can take some time
		# For some strange reason, putting a string argument first expands the number of arguments to
		# the length of the string, so give it a dummy integer first to prevent this
		#self.Disable()
		#self.protWin.Disable()
		#thrLoad = Thread(target=self.PyMOLPDBLoad, args=(1, filename))
		#thrLoad.start()
		#self.PyMOLPDBLoad(1, filename)
	else:
	    logInfo("Cancelled change scorefxn operation")
	dlg.Destroy()
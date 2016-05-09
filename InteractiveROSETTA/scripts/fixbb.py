import wx
import wx.grid
import wx.lib.scrolledpanel
import os
import os.path
import time
import platform
import multiprocessing
import Bio.PDB
import webbrowser
from threading import Thread
from tools import *

class FixbbPanel(wx.lib.scrolledpanel.ScrolledPanel):
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
            self.lblProt = wx.StaticText(self, -1, "Fixed Backbone Design", (25, 15), (270, 25), wx.ALIGN_CENTRE)
            self.lblProt.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblProt = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/lblFixbb.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(25, 15), size=(270, 25))
        else:
            self.lblProt = wx.StaticText(self, -1, "Fixed Backbone Design", (70, 15), style=wx.ALIGN_CENTRE)
            self.lblProt.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
            resizeTextControlForUNIX(self.lblProt, 0, self.GetSize()[0])
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
            self.lblInst = wx.StaticText(self, -1, "Highlight residues to add/remove to design", (0, 45), (320, 25), wx.ALIGN_CENTRE)
            self.lblInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
        elif (platform.system() == "Darwin"):
            self.lblInst = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/lblInstFixbb.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 45), size=(320, 25))
        else:
            self.lblInst = wx.StaticText(self, -1, "Highlight residues to add/remove to design", (5, 45), style=wx.ALIGN_CENTRE)
            self.lblInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
            resizeTextControlForUNIX(self.lblInst, 0, self.GetSize()[0])
        self.lblInst.SetForegroundColour("#FFFFFF")
        
        if (platform.system() == "Darwin"):
            self.btnAminoA = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnA.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(7, 70), size=(42, 25))
        else:
            self.btnAminoA = wx.Button(self, id=-1, label="A", pos=(7, 70), size=(42, 25))
            self.btnAminoA.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoA.Bind(wx.EVT_BUTTON, self.aminoA)
        self.btnAminoA.SetToolTipString("Add ALA to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoC = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnC.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(51, 70), size=(42, 25))
        else:
            self.btnAminoC = wx.Button(self, id=-1, label="C", pos=(51, 70), size=(42, 25))
            self.btnAminoC.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoC.Bind(wx.EVT_BUTTON, self.aminoC)
        self.btnAminoC.SetToolTipString("Add CYS to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoD = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnD.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(95, 70), size=(42, 25))
        else:
            self.btnAminoD = wx.Button(self, id=-1, label="D", pos=(95, 70), size=(42, 25))
            self.btnAminoD.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoD.Bind(wx.EVT_BUTTON, self.aminoD)
        self.btnAminoD.SetToolTipString("Add ASP to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoE = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnE.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(139, 70), size=(42, 25))
        else:
            self.btnAminoE = wx.Button(self, id=-1, label="E", pos=(139, 70), size=(42, 25))
            self.btnAminoE.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoE.Bind(wx.EVT_BUTTON, self.aminoE)
        self.btnAminoE.SetToolTipString("Add GLU to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoF = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnF.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(183, 70), size=(42, 25))
        else:
            self.btnAminoF = wx.Button(self, id=-1, label="F", pos=(183, 70), size=(42, 25))
            self.btnAminoF.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoF.Bind(wx.EVT_BUTTON, self.aminoF)
        self.btnAminoF.SetToolTipString("Add PHE to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoG = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnG.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(227, 70), size=(42, 25))
        else:
            self.btnAminoG = wx.Button(self, id=-1, label="G", pos=(227, 70), size=(42, 25))
            self.btnAminoG.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoG.Bind(wx.EVT_BUTTON, self.aminoG)
        self.btnAminoG.SetToolTipString("Add GLY to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoH = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnH.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(271, 70), size=(42, 25))
        else:
            self.btnAminoH = wx.Button(self, id=-1, label="H", pos=(271, 70), size=(42, 25))
            self.btnAminoH.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoH.Bind(wx.EVT_BUTTON, self.aminoH)
        self.btnAminoH.SetToolTipString("Add HIS to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoI = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnI.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(7, 100), size=(42, 25))
        else:
            self.btnAminoI = wx.Button(self, id=-1, label="I", pos=(7, 100), size=(42, 25))
            self.btnAminoI.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoI.Bind(wx.EVT_BUTTON, self.aminoI)
        self.btnAminoI.SetToolTipString("Add ILE to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoK = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnK.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(51, 100), size=(42, 25))
        else:
            self.btnAminoK = wx.Button(self, id=-1, label="K", pos=(51, 100), size=(42, 25))
            self.btnAminoK.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoK.Bind(wx.EVT_BUTTON, self.aminoK)
        self.btnAminoK.SetToolTipString("Add LYS to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoL = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnL.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(95, 100), size=(42, 25))
        else:
            self.btnAminoL = wx.Button(self, id=-1, label="L", pos=(95, 100), size=(42, 25))
            self.btnAminoL.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoL.Bind(wx.EVT_BUTTON, self.aminoL)
        self.btnAminoL.SetToolTipString("Add LEU to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoM = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnM.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(139, 100), size=(42, 25))
        else:
            self.btnAminoM = wx.Button(self, id=-1, label="M", pos=(139, 100), size=(42, 25))
            self.btnAminoM.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoM.Bind(wx.EVT_BUTTON, self.aminoM)
        self.btnAminoM.SetToolTipString("Add MET to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoN = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnN.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(183, 100), size=(42, 25))
        else:
            self.btnAminoN = wx.Button(self, id=-1, label="N", pos=(183, 100), size=(42, 25))
            self.btnAminoN.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoN.Bind(wx.EVT_BUTTON, self.aminoN)
        self.btnAminoN.SetToolTipString("Add ASN to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoP = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnP.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(227, 100), size=(42, 25))
        else:
            self.btnAminoP = wx.Button(self, id=-1, label="P", pos=(227, 100), size=(42, 25))
            self.btnAminoP.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoP.Bind(wx.EVT_BUTTON, self.aminoP)
        self.btnAminoP.SetToolTipString("Add PRO to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoQ = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnQ.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(271, 100), size=(42, 25))
        else:
            self.btnAminoQ = wx.Button(self, id=-1, label="Q", pos=(271, 100), size=(42, 25))
            self.btnAminoQ.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoQ.Bind(wx.EVT_BUTTON, self.aminoQ)
        self.btnAminoQ.SetToolTipString("Add GLN to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoR = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnR.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(7, 130), size=(42, 25))
        else:
            self.btnAminoR = wx.Button(self, id=-1, label="R", pos=(7, 130), size=(42, 25))
            self.btnAminoR.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoR.Bind(wx.EVT_BUTTON, self.aminoR)
        self.btnAminoR.SetToolTipString("Add ARG to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoS = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnS.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(51, 130), size=(42, 25))
        else:
            self.btnAminoS = wx.Button(self, id=-1, label="S", pos=(51, 130), size=(42, 25))
            self.btnAminoS.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoS.Bind(wx.EVT_BUTTON, self.aminoS)
        self.btnAminoS.SetToolTipString("Add SER to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoT = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnT.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(95, 130), size=(42, 25))
        else:
            self.btnAminoT = wx.Button(self, id=-1, label="T", pos=(95, 130), size=(42, 25))
            self.btnAminoT.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoT.Bind(wx.EVT_BUTTON, self.aminoT)
        self.btnAminoT.SetToolTipString("Add THR to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoV = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnV.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(139, 130), size=(42, 25))
        else:
            self.btnAminoV = wx.Button(self, id=-1, label="V", pos=(139, 130), size=(42, 25))
            self.btnAminoV.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoV.Bind(wx.EVT_BUTTON, self.aminoV)
        self.btnAminoV.SetToolTipString("Add VAL to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoW = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnW.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(183, 130), size=(42, 25))
        else:
            self.btnAminoW = wx.Button(self, id=-1, label="W", pos=(183, 130), size=(42, 25))
            self.btnAminoW.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoW.Bind(wx.EVT_BUTTON, self.aminoW)
        self.btnAminoW.SetToolTipString("Add TRP to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoY = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnY.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(227, 130), size=(42, 25))
        else:
            self.btnAminoY = wx.Button(self, id=-1, label="Y", pos=(227, 130), size=(42, 25))
            self.btnAminoY.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoY.Bind(wx.EVT_BUTTON, self.aminoY)
        self.btnAminoY.SetToolTipString("Add TYR to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoX = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnX.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(271, 130), size=(42, 25))
        else:
            self.btnAminoX = wx.Button(self, id=-1, label="X", pos=(271, 130), size=(42, 25))
            self.btnAminoX.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoX.Bind(wx.EVT_BUTTON, self.aminoX)
        self.btnAminoX.SetToolTipString("Add all amino acids to the design palette")
        self.palette = ""
        self.addType = "WT"
        
        if (platform.system() == "Darwin"):
            self.btnAdd = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnAdd.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(7, 160), size=(57, 25))
        else:
            self.btnAdd = wx.Button(self, id=-1, label="Add", pos=(7, 160), size=(57, 25))
            self.btnAdd.SetForegroundColour("#000000")
            self.btnAdd.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAdd.Bind(wx.EVT_BUTTON, self.add)
        self.btnAdd.SetToolTipString("Add selected residues to the resfile")
        if (platform.system() == "Darwin"):
            self.btnRemove = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnRemove.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(69, 160), size=(57, 25))
        else:
            self.btnRemove = wx.Button(self, id=-1, label="Remove", pos=(69, 160), size=(57, 25))
            self.btnRemove.SetForegroundColour("#000000")
            self.btnRemove.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnRemove.Bind(wx.EVT_BUTTON, self.remove)
        self.btnRemove.SetToolTipString("Remove selected residues from the resfile")
        if (platform.system() == "Darwin"):
            self.btnRestrict = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnRestrict.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(131, 160), size=(57, 25))
        else:
            self.btnRestrict = wx.Button(self, id=-1, label="Restrict", pos=(131, 160), size=(57, 25))
            self.btnRestrict.SetForegroundColour("#000000")
            self.btnRestrict.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnRestrict.Bind(wx.EVT_BUTTON, self.restrict)
        self.btnRestrict.SetToolTipString("Restrict the resfile contents to the selected residues")
        if (platform.system() == "Darwin"):
            self.btnAll = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnAll.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(193, 160), size=(57, 25))
        else:
            self.btnAll = wx.Button(self, id=-1, label="All", pos=(193, 160), size=(57, 25))
            self.btnAll.SetForegroundColour("#000000")
            self.btnAll.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAll.Bind(wx.EVT_BUTTON, self.addAll)
        self.btnAll.SetToolTipString("Add all residues from the selected model to the resfile")
        if (platform.system() == "Darwin"):
            self.btnClear = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnClear.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(255, 160), size=(57, 25))
        else:
            self.btnClear = wx.Button(self, id=-1, label="Clear", pos=(255, 160), size=(57, 25))
            self.btnClear.SetForegroundColour("#000000")
            self.btnClear.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnClear.Bind(wx.EVT_BUTTON, self.clear)
        self.btnClear.SetToolTipString("Clear the contents of the resfile")
        self.selectedData = []
        
        self.desMenu = wx.ComboBox(self, pos=(7, 190), size=(119, 25), choices=[], style=wx.CB_READONLY)
        self.desMenu.Bind(wx.EVT_COMBOBOX, self.desMenuSelect)
        self.desMenu.SetToolTipString("Select resfile entries to edit")
        self.designView = ""
        self.selectedModel = ""
        if (platform.system() == "Darwin"):
            self.btnApply = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnApply.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(131, 190), size=(181, 25))
        else:
            self.btnApply = wx.Button(self, id=-1, label="Apply Selection", pos=(131, 190), size=(181, 25))
            self.btnApply.SetForegroundColour("#000000")
            self.btnApply.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnApply.Bind(wx.EVT_BUTTON, self.applySelection)
        self.btnApply.SetToolTipString("Apply the default add type selection to the current resfile selection")
        
        self.AAlist = ["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"]
        self.AAlist.sort()
        
        if (platform.system() == "Darwin"):
            self.scoretypeMenu = wx.ComboBox(self, pos=(7, 220), size=(305, 25), choices=[], style=wx.CB_READONLY)
        else:
            self.scoretypeMenu = wx.ComboBox(self, pos=(7, 220), size=(305, 25), choices=[], style=wx.CB_READONLY | wx.CB_SORT)
        self.scoretypeMenu.Bind(wx.EVT_COMBOBOX, self.scoretypeMenuSelect)
        self.scoretypeMenu.Disable() # Is only enabled after a design and before accepting it
        self.scoretypeMenu.SetToolTipString("Set the scoretype by which the PyMOL residues will be colored")
        
        self.grdResfile = wx.grid.Grid(self)
        self.grdResfile.CreateGrid(0, 2)
        if (winh-235 > 200):
            self.grdResfile.SetSize((320, winh-235))
        else:
            self.grdResfile.SetSize((320, 200))
        self.grdResfile.SetPosition((0, 250))
        self.grdResfile.SetLabelFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.grdResfile.DisableDragColSize()
        self.grdResfile.DisableDragRowSize()
        self.grdResfile.SetColLabelValue(0, "Residue Options")
        self.grdResfile.SetColLabelValue(1, "Model")
        self.grdResfile.SetRowLabelSize(80)
        self.grdResfile.SetColSize(0, 150)
        self.grdResfile.SetColSize(1, 90)
        self.grdResfile.Bind(wx.grid.EVT_GRID_CELL_LEFT_CLICK, self.gridClick)
        self.selectedr = -1
        self.resfile = []
        # Switch for telling the resfile to look for the restype from the sequence window pre-design
        # or from the outputted designed before accepting it
        self.useDesignedSeq = False 
        
        ypos = self.grdResfile.GetPosition()[1] + self.grdResfile.GetSize()[1] + 10
        if (platform.system() == "Darwin"):
            self.btnLoadResfile = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnLoadResfile.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(20, ypos), size=(120, 25))
        else:
            self.btnLoadResfile = wx.Button(self, id=-1, label="Load Resfile", pos=(20, ypos), size=(120, 25))
            self.btnLoadResfile.SetForegroundColour("#000000")
            self.btnLoadResfile.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.btnLoadResfile.Bind(wx.EVT_BUTTON, self.loadResfile)
        self.btnLoadResfile.SetToolTipString("Load the data in a premade resfile")
        if (platform.system() == "Darwin"):
            self.btnSaveResfile = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnSaveResfile.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(175, ypos), size=(120, 25))
        else:
            self.btnSaveResfile = wx.Button(self, id=-1, label="Save Resfile", pos=(175, ypos), size=(120, 25))
            self.btnSaveResfile.SetForegroundColour("#000000")
            self.btnSaveResfile.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.btnSaveResfile.Bind(wx.EVT_BUTTON, self.saveResfile)
        self.btnSaveResfile.SetToolTipString("Save the current resfile data to a real Rosetta resfile")
        
        if (platform.system() == "Darwin"):
            self.btnWTType = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnWTType_NATRO.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(20, ypos+30), size=(120, 25))
        else:
            self.btnWTType = wx.Button(self, id=-1, label="NATRO", pos=(20, ypos+30), size=(120, 25))
            self.btnWTType.SetForegroundColour("#000000")
            self.btnWTType.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.btnWTType.Bind(wx.EVT_BUTTON, self.changeWTType)
        self.btnWTType.SetToolTipString("Unspecified residues will select the wildtype rotamer only")
        if (platform.system() == "Darwin"):
            self.btnDesign = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnDesign.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(175, ypos+30), size=(120, 25))
        else:
            self.btnDesign = wx.Button(self, id=-1, label="Design!", pos=(175, ypos+30), size=(120, 25))
            self.btnDesign.SetForegroundColour("#000000")
            self.btnDesign.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.btnDesign.Bind(wx.EVT_BUTTON, self.designClick)
        self.btnDesign.SetToolTipString("Perform fixed backbone design")
        self.buttonState = "Design!"
        self.WTType = "NATRO"
        
        self.scrollh = self.btnDesign.GetPosition()[1] + self.btnDesign.GetSize()[1] + 5
        self.SetScrollbars(1, 1, 320, self.scrollh)
        self.winscrollpos = 0
        self.Bind(wx.EVT_SCROLLWIN, self.scrolled)
    
    def showHelp(self, event):
        # Open the help page
        if (platform.system() == "Darwin"):
            try:
                browser = webbrowser.get("Safari")
            except:
                print "Could not load Safari!  The help files are located at " + self.scriptdir + "/help"
                return
            browser.open(self.parent.parent.scriptdir + "/help/fixbb.html")
        else:
            webbrowser.open(self.parent.parent.scriptdir + "/help/fixbb.html")
    
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

    def gridClick(self, event):
        self.selectedr = event.GetRow()
        self.desMenu.SetSelection(event.GetRow())
        self.desMenuSelect(event)
        event.Skip()

    def updateResfile(self):
        # This function redraws the resfile grid to reflect changes to self.resfile
        scrollpos = self.grdResfile.GetScrollPos(wx.VERTICAL)
        self.desMenu.Clear()
        self.selectedModel = ""
        #try:
        #    self.cmd.remove("designed_view")
        #    self.cmd.delete("designed_view")
        #except:
        #    pass
        if (self.grdResfile.NumberRows > 0):
            self.grdResfile.DeleteRows(0, self.grdResfile.NumberRows)
        row = 0
        model = ""
        for [indx, r, seqpos, poseindx, chainID, chainoffset, reslist] in self.resfile:
            self.grdResfile.AppendRows(1)
            ID = self.seqWin.IDs[r]
            model = ID[0:len(ID)-2]
            resn = self.seqWin.SeqViewer.GetCellValue(r, indx)
            if (self.useDesignedSeq):
                if (chainID == "_" or chainID == ""):
                    chain = " "
                else:
                    chain = chainID
                mut = AA3to1(self.designedView[0][chain][self.seqWin.indxToSeqPos[r][indx]].resname)
                label = str(row+1) + ": " + str(resn) + str(seqpos) + str(mut)
            else:
                label = str(row+1) + ": " + str(resn) + str(seqpos)
            self.grdResfile.SetRowLabelValue(row, label)
            self.grdResfile.SetCellValue(row, 0, reslist)
            self.grdResfile.SetCellValue(row, 1, ID)
            # This needs to happen before the readOnly attr is set otherwise it doesn't apply the center on the
            # last cell for some bizarre reason
            self.grdResfile.SetCellAlignment(row, 0, wx.ALIGN_CENTRE, wx.ALIGN_CENTRE)
            self.grdResfile.SetCellAlignment(row, 1, wx.ALIGN_CENTRE, wx.ALIGN_CENTRE)
            # Very important note: you actually need to create a new GridCellAttr for each new row
            # You cannot just declare it outside of the loop and use the same one for each row otherwise you
            # get some pretty nasty crashes when you delete rows
            readOnly = wx.grid.GridCellAttr()
            readOnly.SetReadOnly(True)
            self.grdResfile.SetRowAttr(row, readOnly)
            # Now update the drop down menu so the user can tweak the settings of individual residues in the minmap
            self.desMenu.AppendItems([label])
            row = row + 1
        self.selectedModel = model
        # Resize the label width to fit long PDB filenames
        font = self.grdResfile.GetFont()
        dc = wx.WindowDC(self.grdResfile)
        dc.SetFont(font)
        maxwidth = 80
        for i in range(0, self.grdResfile.NumberRows):
            (w, h) = dc.GetTextExtent(self.grdResfile.GetCellValue(i, 1))
            if (w > maxwidth):
                maxwidth = w
        self.grdResfile.SetColSize(1, maxwidth+10)
        # Resize columns if necessary
        fitGridColumn(self.grdResfile, 0, 150)
        fitGridColumn(self.grdResfile, 1, 90)
        self.grdResfile.Scroll(0, scrollpos)
        # Update the coloring if after a design
        if (self.buttonState != "Design!"):
            self.recolorGrid(self.designedView, self.residue_E, self.grdResfile, self.scoretypeMenu.GetStringSelection())
        
    def scrolled(self, event):
        self.winscrollpos = self.GetScrollPos(wx.VERTICAL)
        event.Skip()
        
    def activate(self):
        # It's possible that the user could have deleted chains/residues that are currently in the minmap
        # Let's first make sure everything in the minmap still exists
        redrawNeeded = False
        for i in range(len(self.resfile)-1, -1, -1):
            [indx, r, seqpos, p, chainID, offset, resstring] = self.resfile[i]
            ID = self.grdResfile.GetCellValue(i, 1)
            #ID = ID.split(":")[0].strip()
            if (r >= len(self.seqWin.IDs) or ID != self.seqWin.IDs[r]):
                # Chain was deleted, pop this item
                self.resfile.pop(i)
                redrawNeeded = True
            elif (indx >= len(self.seqWin.indxToSeqPos[r]) or self.seqWin.indxToSeqPos[r][indx][1] != int(seqpos)):
                # Residue was deleted, pop this item
                self.resfile.pop(i)
                redrawNeeded = True
        if (redrawNeeded):
            self.updateResfile()
        # Grab the current selection of residues for processing with the buttons
        topLefts = self.seqWin.SeqViewer.GetSelectionBlockTopLeft()
        bottomRights = self.seqWin.SeqViewer.GetSelectionBlockBottomRight()
        self.selectedData = []
        for i in range(0, len(topLefts)):
            for r in range(topLefts[i][0], bottomRights[i][0]+1):
                for c in range(topLefts[i][1], bottomRights[i][1]+1):
                    if (self.seqWin.indxToSeqPos[r][c] == "-"):
                        continue
                    seqpos = str(self.seqWin.indxToSeqPos[r][c][1])
                    indx = c
                    poseindx = r
                    while (not(self.seqWin.poses[poseindx])):
                        poseindx = poseindx - 1
                    chainoffset = r - poseindx
                    chainID = self.seqWin.IDs[r][len(self.seqWin.IDs[r])-1]
                    chain = self.seqWin.IDs[r][len(self.seqWin.IDs[r])-1]
                    if (chain == "_"):
                        chain = " "
                    # Don't add any NCAAs or HETATMs for now
                    if ("ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".find(self.seqWin.poses[poseindx][0][chain][self.seqWin.indxToSeqPos[r][c]].resname) >= 0):
                        self.selectedData.append([indx, r, seqpos, poseindx, chainID, chainoffset])
        self.Scroll(0, self.winscrollpos)
    
    def checkIfNewModel(self):
        # This function asks whether the user is attempting to add things to the resfile from another model
        # There should only be one model designed at a time so a warning message will be displayed
        # saying that if the user continues the old resfile will be lost
        totalModels = []
        abort = False
        for data in self.selectedData:
            # data[3] is the poseindx and there should only be one in the selection and the resfile
            if (not(data[3] in totalModels)):
                totalModels.append(data[3])
        for data in self.resfile:
            # data[3] is the poseindx and there should only be one in the selection and the resfile
            if (not(data[3] in totalModels)):
                totalModels.append(data[3])
        if (len(totalModels) > 1):
            dlg = wx.MessageDialog(self, "You may only design one protein at a time and are attempting to add residues from many proteins to the resfile.  If you continue, the current entries in the resfile will be replaced with the selection.  Continue?", "Multiple Proteins in Resfile", wx.YES_NO | wx.ICON_EXCLAMATION | wx.CENTRE)
            if (dlg.ShowModal() == wx.ID_YES):
                # Clear it out and the new selections will be added by the calling function
                self.resfile = []
            else:
                abort = True
            dlg.Destroy()
        return abort
    
    def WT(self, event):
        self.addType = "WT"
        self.btnWT.SetForegroundColour("#FF0000")
        self.btnHydrophilic.SetForegroundColour("#000000")
        self.btnHydrophobic.SetForegroundColour("#000000")
        self.btnAromatic.SetForegroundColour("#000000")
        self.btnAllAA.SetForegroundColour("#000000")
        logInfo("Add type set to WT")

    def hydrophilic(self, event):
        self.addType = "Hydrophilic"
        self.btnWT.SetForegroundColour("#000000")
        self.btnHydrophilic.SetForegroundColour("#FF0000")
        self.btnHydrophobic.SetForegroundColour("#000000")
        self.btnAromatic.SetForegroundColour("#000000")
        self.btnAllAA.SetForegroundColour("#000000")
        logInfo("Add type set to HPH")
        
    def hydrophobic(self, event):
        self.addType = "Hydrophobic"
        self.btnWT.SetForegroundColour("#000000")
        self.btnHydrophilic.SetForegroundColour("#000000")
        self.btnHydrophobic.SetForegroundColour("#FF0000")
        self.btnAromatic.SetForegroundColour("#000000")
        self.btnAllAA.SetForegroundColour("#000000")
        logInfo("Add type set to HPO")
        
    def aromatic(self, event):
        self.addType = "Aromatic"
        self.btnWT.SetForegroundColour("#000000")
        self.btnHydrophilic.SetForegroundColour("#000000")
        self.btnHydrophobic.SetForegroundColour("#000000")
        self.btnAromatic.SetForegroundColour("#FF0000")
        self.btnAllAA.SetForegroundColour("#000000")
        logInfo("Add type set to ARO")
        
    def allAA(self, event):
        self.addType = "AllAA"
        self.btnWT.SetForegroundColour("#000000")
        self.btnHydrophilic.SetForegroundColour("#000000")
        self.btnHydrophobic.SetForegroundColour("#000000")
        self.btnAromatic.SetForegroundColour("#000000")
        self.btnAllAA.SetForegroundColour("#FF0000")
        logInfo("Add type set to All")
        
    def addAAToPalette(self, AA):
        # If this amino acid type is in the palette, remove it
        # Otherwise, add it
        if (AA in self.palette):
            indx = self.palette.find(AA)
            self.palette = self.palette[0:indx] + self.palette[indx+1:]
        else:
            added = False
            for i in range(0, len(self.palette)):
                if (ord(AA) < ord(self.palette[i])):
                    self.palette = self.palette[0:i] + AA + self.palette[i:]
                    added = True
                    break
            if (not(added)):
                self.palette = self.palette + AA
        
    def aminoA(self, event):
        if ("A" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoA.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnA.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoA.SetForegroundColour("#000000")
            logInfo("Removed ALA from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoA.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnA_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoA.SetForegroundColour("#FF0000")
            logInfo("Added ALA to palette")
        self.addAAToPalette("A")
    
    def aminoC(self, event):
        if ("C" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoC.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnC.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoC.SetForegroundColour("#000000")
            logInfo("Removed CYS from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoC.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnC_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoC.SetForegroundColour("#FF0000")
            logInfo("Added CYS to palette")
        self.addAAToPalette("C")
    
    def aminoD(self, event):
        if ("D" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoD.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnD.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoD.SetForegroundColour("#000000")
            logInfo("Removed ASP from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoD.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnD_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoD.SetForegroundColour("#FF0000")
            logInfo("Added ASP to palette")
        self.addAAToPalette("D")
    
    def aminoE(self, event):
        if ("E" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoE.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnE.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoE.SetForegroundColour("#000000")
            logInfo("Removed GLU from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoE.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnE_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoE.SetForegroundColour("#FF0000")
            logInfo("Added GLU to palette")
        self.addAAToPalette("E")
        
    def aminoF(self, event):
        if ("F" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoF.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnF.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoF.SetForegroundColour("#000000")
            logInfo("Removed PHE from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoF.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnF_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoF.SetForegroundColour("#FF0000")
            logInfo("Added PHE to palette")
        self.addAAToPalette("F")
        
    def aminoG(self, event):
        if ("G" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoG.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnG.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoG.SetForegroundColour("#000000")
            logInfo("Removed GLY from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoG.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnG_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoG.SetForegroundColour("#FF0000")
            logInfo("Added GLY to palette")
        self.addAAToPalette("G")
        
    def aminoH(self, event):
        if ("H" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoH.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnH.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoH.SetForegroundColour("#000000")
            logInfo("Removed HIS from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoH.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnH_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoH.SetForegroundColour("#FF0000")
            logInfo("Added HIS to palette")
        self.addAAToPalette("H")
        
    def aminoI(self, event):
        if ("I" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoI.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnI.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoI.SetForegroundColour("#000000")
            logInfo("Removed ILE from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoI.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnI_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoI.SetForegroundColour("#FF0000")
            logInfo("Added ILE to palette")
        self.addAAToPalette("I")
        
    def aminoK(self, event):
        if ("K" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoK.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnK.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoK.SetForegroundColour("#000000")
            logInfo("Removed LYS from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoK.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnK_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoK.SetForegroundColour("#FF0000")
            logInfo("Added LYS to palette")
        self.addAAToPalette("K")
        
    def aminoL(self, event):
        if ("L" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoL.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnL.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoL.SetForegroundColour("#000000")
            logInfo("Removed LEU from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoL.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnL_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoL.SetForegroundColour("#FF0000")
            logInfo("Added LEU to palette")
        self.addAAToPalette("L")
        
    def aminoM(self, event):
        if ("M" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoM.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnM.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoM.SetForegroundColour("#000000")
            logInfo("Removed MET from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoM.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnM_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoM.SetForegroundColour("#FF0000")
            logInfo("Added MET to palette")
        self.addAAToPalette("M")
    
    def aminoN(self, event):
        if ("N" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoN.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnN.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoN.SetForegroundColour("#000000")
            logInfo("Removed ASN from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoN.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnN_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoN.SetForegroundColour("#FF0000")
            logInfo("Added ASN to palette")
        self.addAAToPalette("N")
    
    def aminoP(self, event):
        if ("P" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoP.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnP.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoP.SetForegroundColour("#000000")
            logInfo("Removed PRO from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoP.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnP_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoP.SetForegroundColour("#FF0000")
            logInfo("Added PRO to palette")
        self.addAAToPalette("P")
    
    def aminoQ(self, event):
        if ("Q" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoQ.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnQ.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoQ.SetForegroundColour("#000000")
            logInfo("Removed GLN from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoQ.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnQ_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoQ.SetForegroundColour("#FF0000")
            logInfo("Added GLN to palette")
        self.addAAToPalette("Q")
        
    def aminoR(self, event):
        if ("R" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoR.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnR.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoR.SetForegroundColour("#000000")
            logInfo("Removed ARG from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoR.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnR_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoR.SetForegroundColour("#FF0000")
            logInfo("Added ARG to palette")
        self.addAAToPalette("R")
        
    def aminoS(self, event):
        if ("S" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoS.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnS.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoS.SetForegroundColour("#000000")
            logInfo("Removed SER from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoS.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnS_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoS.SetForegroundColour("#FF0000")
            logInfo("Added SER to palette")
        self.addAAToPalette("S")
        
    def aminoT(self, event):
        if ("T" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoT.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnT.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoT.SetForegroundColour("#000000")
            logInfo("Removed THR from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoT.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnT_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoT.SetForegroundColour("#FF0000")
            logInfo("Added THR to palette")
        self.addAAToPalette("T")
        
    def aminoV(self, event):
        if ("V" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoV.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnV.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoV.SetForegroundColour("#000000")
            logInfo("Removed VAL from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoV.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnV_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoV.SetForegroundColour("#FF0000")
            logInfo("Added VAL to palette")
        self.addAAToPalette("V")
        
    def aminoW(self, event):
        if ("W" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoW.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnW.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoW.SetForegroundColour("#000000")
            logInfo("Removed TRP from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoW.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnW_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoW.SetForegroundColour("#FF0000")
            logInfo("Added TRP to palette")
        self.addAAToPalette("W")
        
    def aminoY(self, event):
        if ("Y" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoY.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnY.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoY.SetForegroundColour("#000000")
            logInfo("Removed TYR from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoY.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnY_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoY.SetForegroundColour("#FF0000")
            logInfo("Added TYR to palette")
        self.addAAToPalette("Y")
        
    def aminoX(self, event):
        # For this "all AA" button, if all residues are not selected, select them all
        # If they are all selected, select none of them
        if (len(self.palette) == 20):
            self.palette = ""
            if (platform.system() == "Darwin"):
                self.btnAminoA.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnA.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoC.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnC.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoD.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnD.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoE.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnE.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoF.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnF.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoG.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnG.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoH.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnH.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoI.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnI.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoK.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnK.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoL.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnL.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoM.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnM.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoN.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnN.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoP.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnP.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoQ.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnQ.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoR.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnR.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoS.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnS.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoT.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnT.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoV.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnV.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoW.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnW.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoY.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnY.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoA.SetForegroundColour("#000000")
                self.btnAminoC.SetForegroundColour("#000000")
                self.btnAminoD.SetForegroundColour("#000000")
                self.btnAminoE.SetForegroundColour("#000000")
                self.btnAminoF.SetForegroundColour("#000000")
                self.btnAminoG.SetForegroundColour("#000000")
                self.btnAminoH.SetForegroundColour("#000000")
                self.btnAminoI.SetForegroundColour("#000000")
                self.btnAminoK.SetForegroundColour("#000000")
                self.btnAminoL.SetForegroundColour("#000000")
                self.btnAminoM.SetForegroundColour("#000000")
                self.btnAminoN.SetForegroundColour("#000000")
                self.btnAminoP.SetForegroundColour("#000000")
                self.btnAminoQ.SetForegroundColour("#000000")
                self.btnAminoR.SetForegroundColour("#000000")
                self.btnAminoS.SetForegroundColour("#000000")
                self.btnAminoT.SetForegroundColour("#000000")
                self.btnAminoV.SetForegroundColour("#000000")
                self.btnAminoW.SetForegroundColour("#000000")
                self.btnAminoY.SetForegroundColour("#000000")
        else:
            self.palette = "ACDEFGHIKLMNPQRSTVWY"
            if (platform.system() == "Darwin"):
                self.btnAminoA.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnA_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoC.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnC_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoD.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnD_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoE.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnE_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoF.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnF_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoG.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnG_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoH.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnH_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoI.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnI_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoK.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnK_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoL.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnL_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoM.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnM_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoN.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnN_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoP.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnP_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoQ.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnQ_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoR.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnR_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoS.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnS_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoT.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnT_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoV.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnV_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoW.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnW_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoY.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnY_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoA.SetForegroundColour("#FF0000")
                self.btnAminoC.SetForegroundColour("#FF0000")
                self.btnAminoD.SetForegroundColour("#FF0000")
                self.btnAminoE.SetForegroundColour("#FF0000")
                self.btnAminoF.SetForegroundColour("#FF0000")
                self.btnAminoG.SetForegroundColour("#FF0000")
                self.btnAminoH.SetForegroundColour("#FF0000")
                self.btnAminoI.SetForegroundColour("#FF0000")
                self.btnAminoK.SetForegroundColour("#FF0000")
                self.btnAminoL.SetForegroundColour("#FF0000")
                self.btnAminoM.SetForegroundColour("#FF0000")
                self.btnAminoN.SetForegroundColour("#FF0000")
                self.btnAminoP.SetForegroundColour("#FF0000")
                self.btnAminoQ.SetForegroundColour("#FF0000")
                self.btnAminoR.SetForegroundColour("#FF0000")
                self.btnAminoS.SetForegroundColour("#FF0000")
                self.btnAminoT.SetForegroundColour("#FF0000")
                self.btnAminoV.SetForegroundColour("#FF0000")
                self.btnAminoW.SetForegroundColour("#FF0000")
                self.btnAminoY.SetForegroundColour("#FF0000")
    
    def add(self, event, updateSelection=True):
        if (updateSelection):
            self.activate()
        #logInfo("Add button clicked")
        if (len(self.selectedData) == 0):
            return
        if (self.checkIfNewModel()):
            # This means the user does not want to start designing a new protein and will not add
            # new residues from another protein, so we should abort this operation
            return
        # It's possible that the user selected multiple proteins in a selection and then decided to
        # add more to the resfile.  The previous conditional doesn't make sure that the selection
        # itself doesn't have multiple entries, so just take whatever the first model in the selection
        # is and ignore everything else
        if (len(self.selectedData) > 0):
            posetoadd = self.selectedData[0][3]
        # Now figure out where this belongs so the list stays sorted
        if (len(self.palette) > 0):
            resAddType = self.palette
            self.addType = "Normal"
        else:
            # Default to WT if nothing selected
            self.addType = "WT"
        # For each of the selected entries, first verify that this entry is not already in the resfile and if it
        # isn't then add it in
        for i in range(0, len(self.selectedData)):
            [indx, r, seqpos, poseindx, chainID, chainoffset] = self.selectedData[i]
            # Make sure this is a CAA
            if (not(self.seqWin.getIsCanonicalAA(r, indx))):
                continue
            if (self.addType == "WT"):
                resAddType = self.seqWin.SeqViewer.GetCellValue(r, indx)
            if (poseindx != posetoadd):
                continue
            alreadyIn = False
            for j in range(0, len(self.resfile)):
                [rindx, rr, rseqpos, rposeindx, rchainID, rchainoffset, rreslist] = self.resfile[j]
                if (r == rr and indx == rindx):
                    alreadyIn = True
                    # Just make sure all the selected residues are in the list at this pos
                    for AA in resAddType:
                        if (not(AA in rreslist)):
                            addedIt = False
                            for k in range(0, len(rreslist)):
                                if (AA < rreslist[k]):
                                    rreslist = rreslist[0:k] + AA + rreslist[k:]
                                    addedIt = True
                                    break
                            if (not(addedIt)):
                                rreslist = rreslist + AA
                    self.resfile[j][6] = rreslist
                    break
            if (not(alreadyIn)):
                if (len(self.resfile) == 0):
                    # List empty, add new element
                    self.resfile.append([indx, r, seqpos, poseindx, chainID, chainoffset, resAddType])
                elif (r < self.resfile[0][1] or (r == self.resfile[0][1] and indx < self.resfile[0][0])):
                    # Belongs first
                    self.resfile.insert(0, [indx, r, seqpos, poseindx, chainID, chainoffset, resAddType])
                else:
                    notInYet = True
                    # Maybe it belongs somewhere in the middle?
                    for i in range(0, len(self.resfile)-1):
                        [indx1, r1, seqpos1, poseindx1, chainID1, chainoffset1, type1] = self.resfile[i]
                        [indx2, r2, seqpos2, poseindx2, chainID2, chainoffset2, type2] = self.resfile[i+1]
                        if (r == r1 and r == r2 and indx > indx1 and indx < indx2):
                            notInYet = False
                            self.resfile.insert(i+1, [indx, r, seqpos, poseindx, chainID, chainoffset, resAddType])
                        elif (r == r1 and r < r2 and indx > indx1):
                            notInYet = False
                            self.resfile.insert(i+1, [indx, r, seqpos, poseindx, chainID, chainoffset, resAddType])
                    if (notInYet):
                        # Belongs at the end
                        self.resfile.append([indx, r, seqpos, poseindx, chainID, chainoffset, resAddType])
        self.updateResfile()
    
    def remove(self, event):
        self.activate()
        logInfo("Remove button clicked")
        if (self.selectedr >= 0 and self.selectedr < len(self.resfile)):
            self.resfile.pop(self.selectedr)
            self.selectedr = -1
        # For each of the selected entries, find out if it is in the minmap and remove it if it is
        #for i in range(0, len(self.selectedData)):
        #    [indx, r, seqpos, poseindx, chainID, chainoffset] = self.selectedData[i]
        #    for j in range(0, len(self.resfile)):
        #        [rindx, rr, rseqpos, rposeindx, rchainID, rchainoffset, rreslist] = self.resfile[j]
        #        if (r == rr and indx == rindx):
        #            self.resfile.pop(j)
        #            break
        self.updateResfile()
    
    def restrict(self, event):
        self.activate()
        logInfo("Restrict button clicked")
        # Remove everything and add only the selected residues
        self.resfile = []
        self.add(event)
    
    def addAll(self, event):
        logInfo("Add all button clicked")
        # Add everything that is in the sequence viewer
        # Here "all" refers to everything in the protein already in the resfile or, if nothing is in
        # the resfile, everything in the first protein
        if (len(self.resfile) > 0):
            posetoadd = self.resfile[0][3]
        else:
            posetoadd = 0
        self.resfile = []
        allData = []
        for r in range(0, self.seqWin.SeqViewer.NumberRows):
            for c in range(0, len(self.seqWin.sequences[r])):
                if (self.seqWin.indxToSeqPos[r][c] == "-"):
                    continue
                seqpos = str(self.seqWin.indxToSeqPos[r][c][1])
                indx = c
                poseindx = r
                while (not(self.seqWin.poses[poseindx])):
                    poseindx = poseindx - 1
                chainoffset = r - poseindx
                chainID = self.seqWin.IDs[r][len(self.seqWin.IDs[r])-1]
                chain = self.seqWin.IDs[r][len(self.seqWin.IDs[r])-1]
                if (chain == "_"):
                    chain = " "
                if (poseindx != posetoadd):
                    continue
                # Don't add any NCAAs or HETATMs for now
                if ("ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".find(self.seqWin.poses[poseindx][0][chain][self.seqWin.indxToSeqPos[r][c]].resname) >= 0):
                    allData.append([indx, r, seqpos, poseindx, chainID, chainoffset])
        saveit = self.selectedData
        self.selectedData = allData
        self.add(event, updateSelection=False)
        self.selectedData = saveit
    
    def clear(self, event):
        logInfo("Clear button clicked")
        # Remove everything
        self.resfile = []
        self.updateResfile()
    
    def desMenuSelect(self, event):
        # Find this item in the resfile
        selectedValue = self.desMenu.GetStringSelection()
        logInfo("Design view menu was set to " + selectedValue)
        row = 0
        for r in range(0, self.grdResfile.NumberRows):
            if (self.grdResfile.GetRowLabelValue(r) == selectedValue):
                row = r
                break
        self.currSelection = row
        # Set the selected residue's row to blue so it is easy to see what the selection is
        for r in range(0, self.grdResfile.NumberRows):
            if (r == row):
                for c in range(0, self.grdResfile.NumberCols):
                    self.grdResfile.SetCellBackgroundColour(r, c, "light blue")
                    #self.grdResfile.SetCellTextColour(r, c, "white")
            else:
                for c in range(0, self.grdResfile.NumberCols):
                    self.grdResfile.SetCellBackgroundColour(r, c, "white")
                    #self.grdResfile.SetCellTextColour(r, c, "black")
        self.grdResfile.Refresh()
        # Do a neighborhood view on the selected position 
        model = ""
        fields = self.grdResfile.GetCellValue(self.currSelection, 1).split("|")
        for field in fields[0:len(fields)-1]:
            model = model + field + "|"
        model = model[0:len(model)-1]
        chain = fields[len(fields)-1][0]
        if (model != self.selectedModel):
            self.selectedModel = model
            for i in range(0, len(self.seqWin.IDs)):
                if (self.seqWin.IDs[i].find(model) >= 0):
                    if (self.buttonState == "Design!"):
                        pass
                    else:
                        recolorEnergies(self.designedView, self.residue_E, "designed_view", self.scoretypeMenu.GetStringSelection(), self.cmd)
                    break
        # Find the neighborhood view
        seqpos = selectedValue.split(":")[len(selectedValue.split(":"))-1]
        seqpos = seqpos.strip()
        seqpos = seqpos[1:] # The AA is index 0
        if (self.useDesignedSeq):
            # The label has the mutation at the end after a design, so we have to take that out
            seqpos = seqpos[0:len(seqpos)-1]
        if (self.buttonState == "Design!"):
            firstmodel = model
        else:
            firstmodel = "designed_view"
        self.cmd.hide("all")
        if (chain == " " or chain == "_"):
            self.cmd.select("dessele", "resi " + seqpos + " and model " + firstmodel)
        else:
            self.cmd.select("dessele", "resi " + seqpos + " and model " + firstmodel + " and chain " + chain)
        self.cmd.select("dessele", "model " + firstmodel + " within 12 of dessele")
        self.cmd.show("cartoon", "dessele")
        self.cmd.hide("ribbon", "dessele")
        self.cmd.show("sticks", "dessele")
        self.cmd.set_bond("stick_radius", 0.1, "dessele")
        # Display energy labels for designed structures
        if (self.buttonState == "Finalize!"):
            relabelEnergies(self.designedView, self.residue_E, "designed_view", self.scoretypeMenu.GetStringSelection(), self.cmd, seqpos)
            self.cmd.label("not dessele", "")
        self.cmd.zoom("dessele")
        if (chain == " " or chain == "_"):
            self.cmd.select("dessele", "resi " + seqpos + " and model " + firstmodel)
        else:
            self.cmd.select("dessele", "resi " + seqpos + " and model " + firstmodel + " and chain " + chain)
        self.cmd.show("sticks", "dessele")
        self.cmd.set_bond("stick_radius", 0.25, "dessele")
        # Highlight this residue in PyMOL
        self.cmd.select("seqsele", "dessele")
        if (self.buttonState == "Finalize!"):
            # If this is after a design, also show the original structure in green for comparison
            self.cmd.select("dessele", "model " + self.selectedModel + " and symbol c")
            self.cmd.color("green", "dessele")
            self.cmd.set("cartoon_color", "green", "dessele")
            if (chain == " " or chain == "_"):
                self.cmd.select("dessele", "resi " + seqpos + " and model " + self.selectedModel)
            else:
                self.cmd.select("dessele", "resi " + seqpos + " and model " + self.selectedModel + " and chain " + chain)
            self.cmd.select("dessele", "model " + self.selectedModel + " within 12 of dessele")
            self.cmd.show("cartoon", "dessele")
            self.cmd.hide("ribbon", "dessele")
            self.cmd.show("sticks", "dessele")
            self.cmd.set_bond("stick_radius", 0.1, "dessele")
            self.cmd.zoom("dessele")
            if (chain == " " or chain == "_"):
                self.cmd.select("dessele", "resi " + seqpos + " and model " + self.selectedModel)
            else:
                self.cmd.select("dessele", "resi " + seqpos + " and model " + self.selectedModel + " and chain " + chain)
            self.cmd.show("sticks", "dessele")
            self.cmd.set_bond("stick_radius", 0.25, "dessele")
        self.cmd.enable("seqsele")
        self.cmd.delete("dessele")
        self.seqWin.selectUpdate(False)
    
    def applySelection(self, event):
        logInfo("Apply selection button clicked")
        if (self.desMenu.GetStringSelection() != ""):
            indx = self.resfile[self.currSelection][0]
            if (len(self.palette) == 0):
                resAddType = self.seqWin.SeqViewer.GetCellValue(self.currSelection, indx)
            else:
                resAddType = self.palette
            self.resfile[self.currSelection][6] = resAddType
            self.grdResfile.SetCellValue(self.currSelection, 0, self.resfile[self.currSelection][6])
    
    def addRes(self, event):
        # Add the chosen AA to the list of residue options at the selected position
        if (self.desMenu.GetStringSelection() != "" and self.resMenu.GetStringSelection() != ""):
            AA3 = self.resMenu.GetStringSelection()
            logInfo("Added " + AA3 + " to the selected position")
            AAindx = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR ".find(AA3)
            AA1 = "ACDEFGHIKLMNPQRSTVWY"[AAindx/4]
            reslist = self.resfile[self.currSelection][6]
            if (AA1 in reslist):
                # Don't add duplicates
                return
            # Insert in alphabetical order
            added = False
            for i in range(0, len(reslist)):
                if (AA1 < reslist[i]):
                    reslist = reslist[0:i] + AA1 + reslist[i:]
                    added = True
                    break
            if (not(added)):
                reslist = reslist + AA1
            self.resfile[self.currSelection][6] = reslist
            self.grdResfile.SetCellValue(self.currSelection, 0, self.resfile[self.currSelection][6])
    
    def removeRes(self, event):
        # Remove the chosen AA to the list of residue options at the selected position
        if (self.desMenu.GetStringSelection() != "" and self.resMenu.GetStringSelection() != ""):
            AA3 = self.resMenu.GetStringSelection()
            logInfo("Removed " + AA3 + " from the selected position")
            AAindx = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR ".find(AA3)
            AA1 = "ACDEFGHIKLMNPQRSTVWY"[AAindx/4]
            reslist = self.resfile[self.currSelection][6]
            todeleteindx = reslist.find(AA1)
            if (todeleteindx >= 0 and len(reslist) > 1):
                reslist = reslist[0:todeleteindx] + reslist[todeleteindx+1:]
                self.resfile[self.currSelection][6] = reslist
                self.grdResfile.SetCellValue(self.currSelection, 0, self.resfile[self.currSelection][6])
            elif (todeleteindx >= 0 and len(reslist) == 1):
                # This was the only residue at this position so we have to remove this position from
                # the resfile
                self.resfile.pop(self.currSelection)
                self.updateResfile()
    
    def recolorGrid(self, pose, allresidue_E, grid, selectedScoretype):
        # Useful function for recoloring the text color of a grid by energy so the user can see all the relative
        # energies just by looking at the colors in the grid
        found = False
        for (scoretypestr, name) in scoretypes.items():
            if (name == selectedScoretype or scoretypestr == selectedScoretype):
                scoretype = scoretypestr
                found = True
                break
        if (found):
            sindx = allresidue_E[0].index(scoretype)
        else:
            sindx = 0 # Default to total_score
        residue_E = []
        for i in range(1, len(allresidue_E)):
            residue_E.append(allresidue_E[i][sindx])
        residue_E = scale_list(residue_E)
        for row in range(0, grid.NumberRows):
            chain = grid.GetCellValue(row, 1)[len(grid.GetCellValue(row, 1))-1]
            if (chain == "_"):
                chain = " "
            seqpos = grid.GetRowLabelValue(row).split()[1]
            seqpos = int(seqpos[1:len(seqpos)-1])
            # Find the rosetta index
            ires = 0
            for ch in pose[0]:
                for residue in ch:
                    if (ch.id == chain and residue.id[1] == seqpos):
                        break
                    ires = ires + 1
            r = residue_E[ires]
            b = 255 - r
            g = 0
            for c in range(0, grid.NumberCols):
                grid.SetCellTextColour(row, c, (r, g, b))
        grid.Refresh()
    
    def scoretypeMenuSelect(self, event):
        # Make sure there is even a PyMOL_Mover pose loaded
        if (self.selectedModel == ""):
            return
        logInfo("The selected scoretype was set to " + self.scoretypeMenu.GetStringSelection())
        if (self.buttonState != "Design!"):
            recolorEnergies(self.designedView, self.residue_E, "designed_view", self.scoretypeMenu.GetStringSelection(), self.cmd)
            self.recolorGrid(self.designedView, self.residue_E, self.grdResfile, self.scoretypeMenu.GetStringSelection())
        self.desMenuSelect(event) # To update all the labels
    
    def loadResfile(self, event):
        logInfo("Load Resfile button clicked")
        # Load data from an existing resfile into the resfile window
        # If there is already data in the graph, notify the user that this data will be erased if they
        # proceed further
        if (len(self.resfile) > 0):
            dlg = wx.MessageDialog(self, "The data in your current workflow will be lost and replaced with a loaded resfile.  Are you sure you want to proceed?", "Current Resfile Data Will Be Lost", wx.YES_NO | wx.ICON_EXCLAMATION | wx.CENTRE)
            if (dlg.ShowModal() == wx.ID_YES):
                # Don't clear it out until the user actually selects a resfile (instead of cancelling
                # at the file selection dialog
                pass
            else:
                logInfo("Load Resfile operation cancelled due to data already being in the resfile")
                return
        dlg = wx.FileDialog(
            self, message="Choose a File",
            defaultDir=self.seqWin.cwd,
            defaultFile="",
            wildcard="Resfiles (*.resfile)|*.resfile",
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
            # Does it open?  If yes, then erase the resfile data and continue
            try:
                f = open(filename, "r")
            except:
                wx.MessageBox("The file " + filename.strip() + " cannot be opened!", "File Cannot Be Read", wx.OK|wx.ICON_EXCLAMATION)
                return
            logInfo("Loaded data from a resfile", filename)
            # Make sure a PDB is loaded otherwise there's nothing to do
            if (self.seqWin.SeqViewer.NumberRows == 0):
                wx.MessageBox("There is no protein loaded!  Load a protein first.", "No Proteins Loaded", wx.OK|wx.ICON_EXCLAMATION)
                return
            self.resfile = []
            modelindx = 0
            model = self.seqWin.getModelForChain(modelindx)
            toAdd = []
            alreadystartedwriting = False
            for aline in f:
                if (aline.find("#MODEL") >= 0 and not(alreadystartedwriting)):
                    try:
                        resfilemodel = aline.split()[1].strip()
                    except:
                        # This entry is in someway corrupt, skip it
                        continue
                    # Now see if a model by this name is already loaded
                    modelexists = False
                    for i in range(0, self.seqWin.SeqViewer.NumberRows):
                        if (self.seqWin.IDs[i].find(resfilemodel) >= 0):
                            modelexists = True
                            modelindx = i
                            break
                    if (not(modelexists)):
                        # Tell the user that this resfile is for a different model, but the resfile
                        # will be loaded anyway and invalid positions will be discarded
                        # It defaults to the first model in the SequenceViewer
                        wx.MessageBox("This resfile was created for " + resfilemodel + " which does not appear to be loaded.  The contents of this resfile will be attempted to be applied to the first model you have loaded.", "Resfile Model Not Found", wx.OK|wx.ICON_WARNING)
                    model = self.seqWin.getModelForChain(modelindx)
                elif (aline.find("start") >= 0):
                    alreadystartedwriting = True
                elif (aline.find("PIKAA") >= 0):
                    fields = aline.split()
                    seqpos = fields[0]
                    chainID = fields[1]
                    reslist = fields[3]
                    # Now we have to find the row/column of this entry from the SequenceWindow
                    found = False
                    for i in range(0, self.seqWin.SeqViewer.NumberRows):
                        if (self.seqWin.IDs[i].find(model) >= 0 and self.seqWin.IDs[i][len(self.seqWin.IDs[i])-1] == chainID):
                            if (chainID == "_"):
                                chain = " "
                            else:
                                chain = chainID
                            indx = 0
                            for residue in self.seqWin.poses[self.seqWin.getPoseIndex(i)][0][chain]:
                                if (residue.id[1] == int(seqpos)):
                                    break
                                indx = indx + 1
                            if (indx < len(self.seqWin.sequences[i])):
                                r = i
                                found = True
                                break
                    if (not(found)):
                        # The models didn't match and this position is not in the default model
                        continue
                    # Make sure this is actually a canonical AA, because we shouldn't be designing
                    # at NCAA positions right now
                    chain = self.seqWin.IDs[r][len(self.seqWin.IDs[r])-1]
                    if (chain == "_"):
                        chain = " "
                    poseindx = self.seqWin.getPoseIndex(r)
                    if ("ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR ".find(self.seqWin.poses[poseindx][0][chain][self.seqWin.indxToSeqPos[r][indx]].resname) < 0):
                        continue
                    co = r - modelindx
                    toAdd.append([indx, r, seqpos, modelindx, chainID, co, reslist])
            # Sort the contents of the resfile so they appear in the graph in sequential order
            for i in range(0, len(toAdd)-1):
                lowest = i
                for j in range(i+1, len(toAdd)):
                    if (toAdd[j][1] < toAdd[lowest][1]):
                        lowest = j
                    elif (toAdd[j][1] == toAdd[lowest][1] and toAdd[j][0] < toAdd[lowest][0]):
                        lowest = j
                temp = toAdd[lowest]
                toAdd[lowest] = toAdd[i]
                toAdd[i] = temp
            # Add them to the resfile grid
            for entry in toAdd:
                self.resfile.append(entry)
            self.updateResfile()
            #self.add(event)
            #self.selectedData = saveit
            f.close()
        else:
            logInfo("Load resfile operation cancelled at the file select window")
    
    def dumpResfile(self, filename):
        poseindx = self.resfile[0][3]
        model = self.seqWin.getModelForChain(poseindx)
        f = open(filename, "w")
        f.write("#MODEL " + model + "\n")
        f.write("USE_INPUT_SC\n")
        f.write(self.WTType + "\n")
        f.write("\n")
        f.write("start\n")
        f.write("\n")
        for [indx, r, seqpos, poseindx, chainID, co, reslist] in self.resfile:
            f.write(str(seqpos) + "\t\t" + chainID + "\tPIKAA " + reslist + " EX 1 EX 2\n")
        f.close()
        logInfo("Wrote a resfile to " + filename)
    
    def saveResfile(self, event):
        if (len(self.resfile) == 0):
            wx.MessageBox("There's nothing to save to a resfile!", "No Resfile", wx.OK|wx.ICON_EXCLAMATION)
            return
        logInfo("Clicked the Save Resfile button")
        # Save the data in the resfile graph as an actual resfile for the user's convenience
        dlg = wx.FileDialog(
            self, message="Save a Resfile",
            defaultDir=self.seqWin.cwd,
            defaultFile="",
            wildcard="Resfiles (*.resfile)|*.resfile",
            style=wx.SAVE | wx.CHANGE_DIR)
        if (dlg.ShowModal() == wx.ID_OK):
            '''paths = dlg.GetPaths(); print 'paths:',paths,dlg.GetFilename(),dlg.GetDirectory()
            # Change cwd to the last opened file
            if (platform.system() == "Windows"):
                lastDirIndx = paths[len(paths)-1].rfind("\\")
            else:
                          try:
                                        lastDirIndx = paths[len(paths)-1].rfind("/")
                          except Exception as e:
                                        print e.message,len(paths)
                                        import traceback,sys; traceback.print_tb(sys.exc_info()[2])
            self.seqWin.cwd = str(paths[len(paths)-1][0:lastDirIndx])
            self.seqWin.saveWindowData(None)
            filename = str(paths[0]).split(".resfile")[0] + ".resfile"'''
            if platform.system()=='Windows':
                    filename = '\\'.join([dlg.GetDirectory(),dlg.GetFilename()])
            else:
                    filename = '/'.join([dlg.GetDirectory(),dlg.GetFilename()])
            print filename
            # Does it exist already?  If so, ask if the user really wants to overwrite it
            if (os.path.isfile(filename)):
                dlg2 = wx.MessageDialog(self, "The file " + filename + " already exists.  Overwrite it?", "Filename Already Exists", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
                if (dlg2.ShowModal() == wx.ID_NO):
                    dlg2.Destroy()
                    logInfo("Canceled save operation due to filename already existing")
            self.dumpResfile(filename)
        else:
            logInfo("Cancelled save resfile operation")
    
    def changeWTType(self, event):
        if (self.WTType == "NATRO"):
            self.WTType = "NATAA"
            if (platform.system() == "Darwin"):
                self.btnWTType.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnWTType_NATAA.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnWTType.SetLabel(self.WTType)
            self.btnWTType.SetToolTipString("Unspecified residues will select rotamers from the wildtype amino acid")
        else:
            self.WTType = "NATRO"
            if (platform.system() == "Darwin"):
                self.btnWTType.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnWTType_NATRO.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnWTType.SetLabel(self.WTType)
            self.btnWTType.SetToolTipString("Unspecified residues will select the wildtype rotamer only")
        logInfo("Changed default behavior to " + self.WTType)
    
    def enableControls(self, enable=True):
        if (enable):
            self.btnAminoA.Enable()
            self.btnAminoC.Enable()
            self.btnAminoD.Enable()
            self.btnAminoE.Enable()
            self.btnAminoF.Enable()
            self.btnAminoG.Enable()
            self.btnAminoH.Enable()
            self.btnAminoI.Enable()
            self.btnAminoK.Enable()
            self.btnAminoL.Enable()
            self.btnAminoM.Enable()
            self.btnAminoN.Enable()
            self.btnAminoP.Enable()
            self.btnAminoQ.Enable()
            self.btnAminoR.Enable()
            self.btnAminoS.Enable()
            self.btnAminoT.Enable()
            self.btnAminoV.Enable()
            self.btnAminoW.Enable()
            self.btnAminoY.Enable()
            self.btnAminoX.Enable()
            self.btnAdd.Enable()
            self.btnRemove.Enable()
            self.btnRestrict.Enable()
            self.btnAll.Enable()
            self.btnClear.Enable()
            self.btnApply.Enable()
            self.btnLoadResfile.Enable()
            self.btnWTType.Enable()
        else:
            self.btnAminoA.Disable()
            self.btnAminoC.Disable()
            self.btnAminoD.Disable()
            self.btnAminoE.Disable()
            self.btnAminoF.Disable()
            self.btnAminoG.Disable()
            self.btnAminoH.Disable()
            self.btnAminoI.Disable()
            self.btnAminoK.Disable()
            self.btnAminoL.Disable()
            self.btnAminoM.Disable()
            self.btnAminoN.Disable()
            self.btnAminoP.Disable()
            self.btnAminoQ.Disable()
            self.btnAminoR.Disable()
            self.btnAminoS.Disable()
            self.btnAminoT.Disable()
            self.btnAminoV.Disable()
            self.btnAminoW.Disable()
            self.btnAminoY.Disable()
            self.btnAminoX.Disable()
            self.btnAdd.Disable()
            self.btnRemove.Disable()
            self.btnRestrict.Disable()
            self.btnAll.Disable()
            self.btnClear.Disable()
            self.btnApply.Disable()
            self.btnLoadResfile.Disable()
            self.btnWTType.Disable()
    
    def cancelDesign(self):
        logInfo("Canceled fixbb operation")
        try:
            os.remove("designinput")
        except:
            pass
        try:
            os.remove("designinputtemp")
        except:
            pass
        self.tmrDesign.Stop()
        self.seqWin.cannotDelete = False
        self.enableControls()
        if (platform.system() == "Darwin"):
            self.btnDesign.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnDesign.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
        else:
            self.btnDesign.SetLabel("Design!")
        self.btnDesign.SetToolTipString("Perform fixed backbone design")
        deleteInputFiles()
        self.parent.parent.restartDaemon()
        self.parent.GoBtn.Enable()
        # Get rid of the messages
        for i in range(0, len(self.seqWin.msgQueue)):
            if (self.seqWin.msgQueue[i].find("Performing protein design") >= 0):
                self.seqWin.msgQueue.pop(i)
                break
        if (len(self.seqWin.msgQueue) > 0):
            self.seqWin.labelMsg.SetLabel(self.seqWin.msgQueue[len(self.seqWin.msgQueue)-1])
        else:
            self.seqWin.labelMsg.SetLabel("")
        self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")
        self.buttonState = "Design!"
    
    def designClick(self, event):
        # This is also the "Finalize!" button
        logInfo("Design button clicked")
        if (self.buttonState == "Design!"):
            if (len(self.resfile) > 0):
                self.seqWin.labelMsg.SetLabel("Performing protein design, please be patient...")
                self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
                self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")
                self.seqWin.msgQueue.append("Performing protein design, please be patient...")
                self.parent.GoBtn.Disable()
                self.enableControls(False)
                self.seqWin.cannotDelete = True
                #thrDesign = Thread(target=self.threadDesign, args=())
                #thrDesign.start()
                self.stage = 1
                if (platform.system() == "Darwin"):
                    self.btnDesign.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnDesign_Cancel.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                else:
                    self.btnDesign.SetLabel("Cancel!")
                self.buttonState = "Cancel!"
                self.btnDesign.SetToolTipString("Cancel the fixbb simulation")
                self.tmrDesign = wx.Timer(self)
                self.Bind(wx.EVT_TIMER, self.threadDesign, self.tmrDesign)
                self.tmrDesign.Start(1000)
            else:
                wx.MessageBox("There's nothing to design!", "Nothing to Design", wx.OK|wx.ICON_EXCLAMATION)
        elif (self.buttonState == "Cancel!"):
            dlg = wx.MessageDialog(self, "Are you sure you want to cancel the protein design simulation?  All progress will be lost.", "Cancel Protein Design Simulation", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
            result = dlg.ShowModal()
            if (result == wx.ID_YES):
                self.cancelDesign()
            dlg.Destroy()
        else:
            # Finalize button, ask whether the changes will be accepted or rejected
            dlg = wx.MessageDialog(self, "Do you want to accept the results of this protein design?", "Accept/Reject Design", wx.YES_NO | wx.CANCEL | wx.ICON_QUESTION | wx.CENTRE)
            result = dlg.ShowModal()
            if (result == wx.ID_YES):
                accept = True
                logInfo("Design accepted")
            elif (result == wx.ID_NO):
                accept = False
                logInfo("Design rejected")
            else:
                dlg.Destroy()
                logInfo("Finalize operation cancelled")
                return
            dlg.Destroy()
            self.useDesignedSeq = False
            self.scoretypeMenu.Disable()
            if (platform.system() == "Darwin"):
                self.btnDesign.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnDesign.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnDesign.SetLabel("Design!")
            self.buttonState = "Design!"
            self.btnDesign.SetToolTipString("Perform fixed backbone design")
            self.cmd.label("all", "")
            self.seqWin.cannotDelete = False
            if (not(accept)):
                self.cmd.remove("designed_view")
                self.cmd.delete("designed_view")
                self.updateResfile() # To get the old amino acid identities into the resfile
                return
            # Get rid of the original pose, save the designed pose, and reload the structure in PyMOL
            poseindx = self.seqWin.getPoseIndexForModel(self.desmodel)
            try:
                self.cmd.remove(self.desmodel)
                self.cmd.delete(self.desmodel)
                self.cmd.remove("designed_view")
                self.cmd.delete("designed_view")
                self.cmd.load(self.desmodel + "_D.pdb", self.desmodel)
                #self.designedView.pdb_info().name(str(self.desmodel + ".pdb"))
                self.seqWin.reloadPose(poseindx, self.desmodel, self.desmodel + "_D.pdb")
                defaultPyMOLView(self.cmd, self.desmodel)
                self.updateResfile() # To get the new amino acid identities into the resfile
                del self.designedView
                # IMPORTANT: You have to replace the model in the sandbox with the new designed model
                os.remove(self.desmodel + ".pdb")
                os.rename(self.desmodel + "_D.pdb", self.desmodel + ".pdb")
            except:
                # Some weird error happened, do nothing instead of crashing
                print "Bug at accept button click"
                pass        
            self.seqWin.recolorResidues()
    
    def recoverFromError(self):
        # This function tells the user what the error was and tries to revert the protocol
        # back to the pre-daemon state so the main GUI can continue to be used
        f = open("errreport", "r")
        errmsg = "An error was encountered during the protocol:\n\n"
        for aline in f:
            errmsg = errmsg + aline.strip()
        f.close()
        errmsg = str(errmsg)
        logInfo("Error encountered")
        logInfo(errmsg)
        if (platform.system() == "Windows"):
            sessioninfo = os.path.expanduser("~") + "\\InteractiveRosetta\\sessionlog"
        else:
            sessioninfo = os.path.expanduser("~") + "/.InteractiveRosetta/sessionlog"
        errmsg = errmsg + "\n\nIf you don't know what caused this, send the file " + sessioninfo + " to a developer along with an explanation of what you did."
        # You have to use a MessageDialog because the MessageBox doesn't always work for some reason
        dlg = wx.MessageDialog(self, errmsg, "Error Encountered", wx.OK|wx.ICON_EXCLAMATION)
        dlg.ShowModal()
        dlg.Destroy()
        os.remove("errreport")
        self.seqWin.cannotDelete = False
        self.parent.GoBtn.Enable()
        self.enableControls(True)
        if (platform.system() == "Darwin"):
            self.btnDesign.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnDesign.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
        else:
            self.btnDesign.SetLabel("Design!")
        self.buttonState = "Design!"
        # Get rid of the messages
        for i in range(0, len(self.seqWin.msgQueue)):
            if (self.seqWin.msgQueue[i].find("Performing protein design") >= 0):
                self.seqWin.msgQueue.pop(i)
                break
        if (len(self.seqWin.msgQueue) > 0):
            self.seqWin.labelMsg.SetLabel(self.seqWin.msgQueue[len(self.seqWin.msgQueue)-1])
        else:
            self.seqWin.labelMsg.SetLabel("")
        self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")
    
    def threadDesign(self, event):
        # Why am I using a Timer?  See the explanation in kic.py
        # Save the resfile to a temporary location that the daemon can find
        goToSandbox()
        if (self.stage == 1):
            self.tmrDesign.Stop()
            self.timeoutCount = 0
            self.dumpResfile("fixbb.resfile")
            f = open("designinputtemp", "w")
            r = self.resfile[0][1]
            pdbfile = self.seqWin.getModelForChain(r) + ".pdb"
            # Dump the PDB from PyMOL first in case the coordinates were altered by the user
            dumpmodel = pdbfile.split(".pdb")[0]
            self.cmd.save(pdbfile.strip(), "model " + dumpmodel)
            fixPyMOLSave(pdbfile.strip())
            f.write("PDBFILE\t" + pdbfile.strip() + "\n")
            f2 = open(pdbfile, "r")
            f.write("BEGIN PDB DATA\n")
            for aline in f2:
                f.write(aline.strip() + "\n")
            f.write("END PDB DATA\n")
            f2.close()
            f.write("RESFILE\tfixbb.resfile\n")
            f2 = open("fixbb.resfile", "r")
            f.write("BEGIN RESFILE DATA\n")
            for aline in f2:
                f.write(aline.strip() + "\n")
            f.write("END RESFILE DATA\n")
            f2.close()
            f.close()
            appendScorefxnParamsInfoToFile("designinputtemp", self.selectWin.weightsfile)
            if (useServer):
                try: 
                    self.ID = sendToServer("designinput")
                    self.usingServer = True
                    logInfo("Design input sent to server daemon with ID " + self.ID)
                except:
                    # Something failed, default to the local daemon
                    os.rename("designinputtemp", "designinput")
                    self.usingServer = False
                    logInfo("Server daemon not available, design input uploaded at designinput")
            else:
                os.rename("designinputtemp", "designinput")
                self.usingServer = False
                logInfo("Design input uploaded locally at designinput")
            self.stage = 2
            self.tmrDesign.Start(1000)
        else:
            if (self.usingServer):
                # See if the file has been uploaded yet and bring it here if so
                queryServerForResults("designoutput-" + self.ID)
                self.timeoutCount = self.timeoutCount + 1
            if (self.timeoutCount >= serverTimeout):
                self.tmrDesign.Stop()
                # If this is taking too long, maybe there's something wrong with the server
                # Ask the user if they want to continue waiting or use the local daemon instead
                dlg = wx.MessageDialog(self, "The server is taking a long time to respond.  Continue to wait?  Pressing No will run the calculations locally.", "Delayed Server Response", wx.YES_NO | wx.ICON_EXCLAMATION | wx.CENTRE)
                if (dlg.ShowModal() == wx.ID_YES):
                    # Reset the counter
                    self.timeoutCount = 0
                else:
                    self.usingServer = False
                    self.timeoutCount = 0
                    os.rename("designinputtemp", "designinput")
                    logInfo("Server took too long to respond so the local daemon was used")
                dlg.Destroy()
                self.tmrDesign.Start(1000)
            # Read the output dumped by the child process
            if (os.path.isfile("designoutput")):
                self.tmrDesign.Stop()
                self.residue_E = []
                pdbreader = Bio.PDB.PDBParser()
                f = open("designoutput", "r")
                for aline in f:
                    if (aline[0:6] == "OUTPUT"):
                        pdbfile = aline.split("\t")[1].strip()
                        self.designedView = pdbreader.get_structure(pdbfile.split("_D.pdb")[0], pdbfile)
                        self.desmodel = pdbfile.split("_D.pdb")[0]
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
                logInfo("Found fixbb output at designoutput")
                # Add the nonzero scoretypes to the energy viewing list from the current score function
                self.scoretypeMenu.Clear()
                for scoretype in self.residue_E[0]:
                    try:
                        toAdd = scoretypes[str(scoretype)]
                    except:
                        toAdd = str(scoretype)
                    self.scoretypeMenu.Append(toAdd)
                self.scoretypeMenu.Enable()
                # Pop this message out of the queue
                for i in range(0, len(self.seqWin.msgQueue)):
                    if (self.seqWin.msgQueue[i].find("Performing protein design") >= 0):
                        self.seqWin.msgQueue.pop(i)
                        break
                if (len(self.seqWin.msgQueue) > 0):
                    self.seqWin.labelMsg.SetLabel(self.seqWin.msgQueue[len(self.seqWin.msgQueue)-1])
                else:
                    self.seqWin.labelMsg.SetLabel("")
                self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
                self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")
                self.parent.GoBtn.Enable()
                self.enableControls()
                self.selectedModel = ""
                if (platform.system() == "Darwin"):
                    self.btnDesign.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/fixbb/btnDesign_Finalize.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                else:
                    self.btnDesign.SetLabel("Finalize!")
                self.buttonState = "Finalize!"
                self.btnDesign.SetToolTipString("Accept or reject the results of this protocol")
                os.remove("designoutput")
                # Load the designed pose as the "designed_view" model so the user can look at the results
                self.cmd.load(pdbfile, "designed_view")
                self.cmd.hide("everything", "model designed_view")
                self.useDesignedSeq = True
                self.updateResfile() # To get the new amino acid identities into the resfile
                self.recolorGrid(self.designedView, self.residue_E, self.grdResfile, self.scoretypeMenu.GetStringSelection())
                # To get the energy values in the B-factors
                recolorEnergies(self.designedView, self.residue_E, "designed_view", "Total Energy", self.cmd)
                self.seqWin.pdbwriter.set_structure(self.designedView)
                self.seqWin.pdbwriter.save(self.desmodel + "_D.pdb") 
                self.grdResfile.Refresh()
            elif (os.path.isfile("errreport")):
                self.tmrDesign.Stop()
                self.recoverFromError()
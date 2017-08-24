import wx
import wx.grid
import wx.lib.scrolledpanel
import os
import os.path
import sys
import time
import platform
import multiprocessing
import Bio.PDB
import webbrowser
import datetime
from threading import Thread
from tools import *

class PointMutantScanPanel(wx.lib.scrolledpanel.ScrolledPanel):
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
            self.lblProt = wx.StaticText(self, -1, "Point Mutant Scan", (25, 15), (270, 25), wx.ALIGN_CENTRE)
            self.lblProt.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblProt = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/lblPointMutantScan.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 15), size=(320, 25))
        else:
            self.lblProt = wx.StaticText(self, -1, "Point Mutant Scan", (70, 15), style=wx.ALIGN_CENTRE)
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
            self.lblInst = wx.StaticText(self, -1, "Highlight residues to add/remove to the scan", (0, 45), (320, 25), wx.ALIGN_CENTRE)
            self.lblInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
        elif (platform.system() == "Darwin"):
            self.lblInst = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/lblInstPointMutantScan.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 45), size=(320, 25))
        else:
            self.lblInst = wx.StaticText(self, -1, "Highlight residues to add/remove to the scan", (5, 45), style=wx.ALIGN_CENTRE)
            self.lblInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
            resizeTextControlForUNIX(self.lblInst, 0, self.GetSize()[0])
        self.lblInst.SetForegroundColour("#FFFFFF")
        
        if (platform.system() == "Darwin"):
            self.btnAminoA = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnA.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(7, 70), size=(42, 25))
        else:
            self.btnAminoA = wx.Button(self, id=-1, label="A", pos=(7, 70), size=(42, 25))
            self.btnAminoA.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoA.Bind(wx.EVT_BUTTON, self.aminoA)
        self.btnAminoA.SetToolTipString("Add ALA to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoC = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnC.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(51, 70), size=(42, 25))
        else:
            self.btnAminoC = wx.Button(self, id=-1, label="C", pos=(51, 70), size=(42, 25))
            self.btnAminoC.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoC.Bind(wx.EVT_BUTTON, self.aminoC)
        self.btnAminoC.SetToolTipString("Add CYS to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoD = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnD.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(95, 70), size=(42, 25))
        else:
            self.btnAminoD = wx.Button(self, id=-1, label="D", pos=(95, 70), size=(42, 25))
            self.btnAminoD.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoD.Bind(wx.EVT_BUTTON, self.aminoD)
        self.btnAminoD.SetToolTipString("Add ASP to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoE = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnE.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(139, 70), size=(42, 25))
        else:
            self.btnAminoE = wx.Button(self, id=-1, label="E", pos=(139, 70), size=(42, 25))
            self.btnAminoE.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoE.Bind(wx.EVT_BUTTON, self.aminoE)
        self.btnAminoE.SetToolTipString("Add GLU to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoF = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnF.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(183, 70), size=(42, 25))
        else:
            self.btnAminoF = wx.Button(self, id=-1, label="F", pos=(183, 70), size=(42, 25))
            self.btnAminoF.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoF.Bind(wx.EVT_BUTTON, self.aminoF)
        self.btnAminoF.SetToolTipString("Add PHE to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoG = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnG.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(227, 70), size=(42, 25))
        else:
            self.btnAminoG = wx.Button(self, id=-1, label="G", pos=(227, 70), size=(42, 25))
            self.btnAminoG.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoG.Bind(wx.EVT_BUTTON, self.aminoG)
        self.btnAminoG.SetToolTipString("Add GLY to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoH = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnH.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(271, 70), size=(42, 25))
        else:
            self.btnAminoH = wx.Button(self, id=-1, label="H", pos=(271, 70), size=(42, 25))
            self.btnAminoH.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoH.Bind(wx.EVT_BUTTON, self.aminoH)
        self.btnAminoH.SetToolTipString("Add HIS to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoI = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnI.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(7, 100), size=(42, 25))
        else:
            self.btnAminoI = wx.Button(self, id=-1, label="I", pos=(7, 100), size=(42, 25))
            self.btnAminoI.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoI.Bind(wx.EVT_BUTTON, self.aminoI)
        self.btnAminoI.SetToolTipString("Add ILE to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoK = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnK.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(51, 100), size=(42, 25))
        else:
            self.btnAminoK = wx.Button(self, id=-1, label="K", pos=(51, 100), size=(42, 25))
            self.btnAminoK.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoK.Bind(wx.EVT_BUTTON, self.aminoK)
        self.btnAminoK.SetToolTipString("Add LYS to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoL = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnL.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(95, 100), size=(42, 25))
        else:
            self.btnAminoL = wx.Button(self, id=-1, label="L", pos=(95, 100), size=(42, 25))
            self.btnAminoL.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoL.Bind(wx.EVT_BUTTON, self.aminoL)
        self.btnAminoL.SetToolTipString("Add LEU to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoM = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnM.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(139, 100), size=(42, 25))
        else:
            self.btnAminoM = wx.Button(self, id=-1, label="M", pos=(139, 100), size=(42, 25))
            self.btnAminoM.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoM.Bind(wx.EVT_BUTTON, self.aminoM)
        self.btnAminoM.SetToolTipString("Add MET to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoN = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnN.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(183, 100), size=(42, 25))
        else:
            self.btnAminoN = wx.Button(self, id=-1, label="N", pos=(183, 100), size=(42, 25))
            self.btnAminoN.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoN.Bind(wx.EVT_BUTTON, self.aminoN)
        self.btnAminoN.SetToolTipString("Add ASN to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoP = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnP.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(227, 100), size=(42, 25))
        else:
            self.btnAminoP = wx.Button(self, id=-1, label="P", pos=(227, 100), size=(42, 25))
            self.btnAminoP.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoP.Bind(wx.EVT_BUTTON, self.aminoP)
        self.btnAminoP.SetToolTipString("Add PRO to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoQ = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnQ.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(271, 100), size=(42, 25))
        else:
            self.btnAminoQ = wx.Button(self, id=-1, label="Q", pos=(271, 100), size=(42, 25))
            self.btnAminoQ.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoQ.Bind(wx.EVT_BUTTON, self.aminoQ)
        self.btnAminoQ.SetToolTipString("Add GLN to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoR = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnR.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(7, 130), size=(42, 25))
        else:
            self.btnAminoR = wx.Button(self, id=-1, label="R", pos=(7, 130), size=(42, 25))
            self.btnAminoR.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoR.Bind(wx.EVT_BUTTON, self.aminoR)
        self.btnAminoR.SetToolTipString("Add ARG to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoS = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnS.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(51, 130), size=(42, 25))
        else:
            self.btnAminoS = wx.Button(self, id=-1, label="S", pos=(51, 130), size=(42, 25))
            self.btnAminoS.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoS.Bind(wx.EVT_BUTTON, self.aminoS)
        self.btnAminoS.SetToolTipString("Add SER to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoT = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnT.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(95, 130), size=(42, 25))
        else:
            self.btnAminoT = wx.Button(self, id=-1, label="T", pos=(95, 130), size=(42, 25))
            self.btnAminoT.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoT.Bind(wx.EVT_BUTTON, self.aminoT)
        self.btnAminoT.SetToolTipString("Add THR to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoV = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnV.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(139, 130), size=(42, 25))
        else:
            self.btnAminoV = wx.Button(self, id=-1, label="V", pos=(139, 130), size=(42, 25))
            self.btnAminoV.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoV.Bind(wx.EVT_BUTTON, self.aminoV)
        self.btnAminoV.SetToolTipString("Add VAL to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoW = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnW.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(183, 130), size=(42, 25))
        else:
            self.btnAminoW = wx.Button(self, id=-1, label="W", pos=(183, 130), size=(42, 25))
            self.btnAminoW.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoW.Bind(wx.EVT_BUTTON, self.aminoW)
        self.btnAminoW.SetToolTipString("Add TRP to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoY = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnY.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(227, 130), size=(42, 25))
        else:
            self.btnAminoY = wx.Button(self, id=-1, label="Y", pos=(227, 130), size=(42, 25))
            self.btnAminoY.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoY.Bind(wx.EVT_BUTTON, self.aminoY)
        self.btnAminoY.SetToolTipString("Add TYR to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoX = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnX.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(271, 130), size=(42, 25))
        else:
            self.btnAminoX = wx.Button(self, id=-1, label="X", pos=(271, 130), size=(42, 25))
            self.btnAminoX.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoX.Bind(wx.EVT_BUTTON, self.aminoX)
        self.btnAminoX.SetToolTipString("Add all amino acids to the design palette")
        self.palette = ""
        self.addType = "WT"
        
        if (platform.system() == "Darwin"):
            self.btnAdd = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnAdd.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(7, 160), size=(57, 25))
        else:
            self.btnAdd = wx.Button(self, id=-1, label="Add", pos=(7, 160), size=(57, 25))
            self.btnAdd.SetForegroundColour("#000000")
            self.btnAdd.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAdd.Bind(wx.EVT_BUTTON, self.add)
        self.btnAdd.SetToolTipString("Add selected residues to the resfile")
        if (platform.system() == "Darwin"):
            self.btnRemove = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnRemove.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(69, 160), size=(57, 25))
        else:
            self.btnRemove = wx.Button(self, id=-1, label="Remove", pos=(69, 160), size=(57, 25))
            self.btnRemove.SetForegroundColour("#000000")
            self.btnRemove.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnRemove.Bind(wx.EVT_BUTTON, self.remove)
        self.btnRemove.SetToolTipString("Remove selected residues from the resfile")
        if (platform.system() == "Darwin"):
            self.btnRestrict = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnRestrict.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(131, 160), size=(57, 25))
        else:
            self.btnRestrict = wx.Button(self, id=-1, label="Restrict", pos=(131, 160), size=(57, 25))
            self.btnRestrict.SetForegroundColour("#000000")
            self.btnRestrict.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnRestrict.Bind(wx.EVT_BUTTON, self.restrict)
        self.btnRestrict.SetToolTipString("Restrict the resfile contents to the selected residues")
        if (platform.system() == "Darwin"):
            self.btnAll = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnAll.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(193, 160), size=(57, 25))
        else:
            self.btnAll = wx.Button(self, id=-1, label="All", pos=(193, 160), size=(57, 25))
            self.btnAll.SetForegroundColour("#000000")
            self.btnAll.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAll.Bind(wx.EVT_BUTTON, self.addAll)
        self.btnAll.SetToolTipString("Add all residues from the selected model to the resfile")
        if (platform.system() == "Darwin"):
            self.btnClear = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnClear.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(255, 160), size=(57, 25))
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
            self.btnApply = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnApply.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(131, 190), size=(181, 25))
        else:
            self.btnApply = wx.Button(self, id=-1, label="Apply Selection", pos=(131, 190), size=(181, 25))
            self.btnApply.SetForegroundColour("#000000")
            self.btnApply.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnApply.Bind(wx.EVT_BUTTON, self.applySelection)
        self.btnApply.SetToolTipString("Apply the default add type selection to the current resfile selection")
        
        self.AAlist = ["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"]
        self.AAlist.sort()
        
        self.grdResfile = wx.grid.Grid(self)
        self.grdResfile.CreateGrid(0, 2)
        if (winh-235 > 200):
            self.grdResfile.SetSize((320, winh-235))
        else:
            self.grdResfile.SetSize((320, 200))
        self.grdResfile.SetPosition((0, 220))
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
        if (platform.system() == "Windows"):
            self.lblProtReport = wx.StaticText(self, -1, "Mutant Report", (25, ypos), (270, 25), wx.ALIGN_CENTRE)
            self.lblProtReport.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblProtReport = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/lblReport.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(25, ypos), size=(270, 25))
        else:
            self.lblProtReport = wx.StaticText(self, -1, "Mutant Report", (70, ypos), style=wx.ALIGN_CENTRE)
            self.lblProtReport.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
            resizeTextControlForUNIX(self.lblProtReport, 0, self.GetSize()[0])
        self.lblProtReport.SetForegroundColour("#FFFFFF")
        
        self.grdReport = wx.grid.Grid(self)
        self.grdReport.CreateGrid(0, 1)
        if (winh-235 > 200):
            self.grdReport.SetSize((320, winh-235))
        else:
            self.grdReport.SetSize((320, 200))
        self.grdReport.SetPosition((0, ypos+30))
        self.grdReport.SetLabelFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.grdReport.DisableDragColSize()
        self.grdReport.DisableDragRowSize()
        self.grdReport.SetColLabelValue(0, u"Predicted \u25b3\u25b3G")
        self.grdReport.SetRowLabelSize(160)
        self.grdReport.SetColSize(0, 160)
        self.grdReport.Bind(wx.grid.EVT_GRID_CELL_LEFT_CLICK, self.reportClick)
        
        ypos = self.grdReport.GetPosition()[1] + self.grdReport.GetSize()[1] + 10
        if (platform.system() == "Darwin"):
            self.btnLoadReport = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnLoadReport.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(40, ypos), size=(100, 25))
        else:
            self.btnLoadReport = wx.Button(self, id=-1, label="Load Report", pos=(40, ypos), size=(100, 25))
            self.btnLoadReport.SetForegroundColour("#000000")
            self.btnLoadReport.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnLoadReport.Bind(wx.EVT_BUTTON, self.loadReport)
        self.btnLoadReport.SetToolTipString("Load previously-saved point mutation scanning data")
        if (platform.system() == "Darwin"):
            self.btnSaveReport = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnSaveReport.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(180, ypos), size=(100, 25))
        else:
            self.btnSaveReport = wx.Button(self, id=-1, label="Save Report", pos=(180, ypos), size=(100, 25))
            self.btnSaveReport.SetForegroundColour("#000000")
            self.btnSaveReport.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnSaveReport.Bind(wx.EVT_BUTTON, self.saveReport)
        self.btnSaveReport.Disable()
        self.btnSaveReport.SetToolTipString("Save this point mutation scan report for later viewing")
        
        if (platform.system() == "Windows"):
            self.lblMutantType = wx.StaticText(self, -1, "Mutant Type", (0, ypos+33), (150, 20), wx.ALIGN_CENTRE)
            self.lblMutantType.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblMutantType = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/lblMutantType.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, ypos+33), size=(150, 20))
        else:
            self.lblMutantType = wx.StaticText(self, -1, "Mutant Type", (0, ypos+33), style=wx.ALIGN_CENTRE)
            self.lblMutantType.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblMutantType, 0, 150)
        self.lblMutantType.SetForegroundColour("#FFFFFF")
        self.mutantType = "Single"
        if (platform.system() == "Darwin"):
            self.btnMutantType = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnMutantType_Single.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(160, ypos+30), size=(120, 25))
        else:
            self.btnMutantType = wx.Button(self, id=-1, label="Single Mutants", pos=(160, ypos+30), size=(150, 25))
            self.btnMutantType.SetForegroundColour("#000000")
            self.btnMutantType.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnMutantType.Bind(wx.EVT_BUTTON, self.toggleMutantType)
        self.btnMutantType.SetToolTipString("Search for stabilizing single point mutations")
        
        if (platform.system() == "Darwin"):
            self.btnServerToggle = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnServer_Off.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(40, ypos+60), size=(100, 25))
        else:
            self.btnServerToggle = wx.Button(self, id=-1, label="Server Off", pos=(40, ypos+60), size=(100, 25))
            self.btnServerToggle.SetForegroundColour("#000000")
            self.btnServerToggle.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.btnServerToggle.Bind(wx.EVT_BUTTON, self.serverToggle)
        self.btnServerToggle.SetToolTipString("Perform point mutation scans locally")
        self.serverOn = False
        if (platform.system() == "Darwin"):
            self.btnScan = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnScan.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(175, ypos+60), size=(100, 25))
        else:
            self.btnScan = wx.Button(self, id=-1, label="Scan!", pos=(180, ypos+60), size=(100, 25))
            self.btnScan.SetForegroundColour("#000000")
            self.btnScan.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.btnScan.Bind(wx.EVT_BUTTON, self.scanClick)
        self.btnScan.SetToolTipString("Perform a point mutant scan")
        self.buttonState = "Scan!"
        
        self.scrollh = self.btnScan.GetPosition()[1] + self.btnScan.GetSize()[1] + 5
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
        
    def reportClick(self, event):
        # Find the neighborhood view
        selectedValue = str(self.grdReport.GetRowLabelValue(event.GetRow()))
        seqpos = selectedValue.split("|")[len(selectedValue.split("|"))-1]
        seqpos = seqpos.strip()
        seqpos = seqpos[1:len(seqpos)-1] # The AA is index 0 and the mutation is the last char
        firstmodel = self.selectedModel
        chain = selectedValue.split("|")[0]
        self.cmd.hide("all")
        if (chain == " " or chain == "_"):
            self.cmd.select("repsele", "resi " + seqpos + " and model " + firstmodel)
        else:
            self.cmd.select("repsele", "resi " + seqpos + " and model " + firstmodel + " and chain " + chain)
        self.cmd.select("repsele", "model " + firstmodel + " within 12 of repsele")
        self.cmd.show("cartoon", "repsele")
        self.cmd.hide("ribbon", "repsele")
        self.cmd.show("sticks", "repsele")
        self.cmd.set_bond("stick_radius", 0.1, "repsele")
        self.cmd.zoom("repsele")
        if (chain == " " or chain == "_"):
            self.cmd.select("repsele", "resi " + seqpos + " and model " + firstmodel)
        else:
            self.cmd.select("repsele", "resi " + seqpos + " and model " + firstmodel + " and chain " + chain)
        self.cmd.show("sticks", "repsele")
        self.cmd.set_bond("stick_radius", 0.25, "repsele")
        # Highlight this residue in PyMOL
        self.cmd.select("seqsele", "repsele")
        self.cmd.enable("seqsele")
        self.cmd.delete("repsele")
        self.seqWin.selectUpdate(False)
        for r in range(0, self.grdReport.NumberRows):
            if (r == event.GetRow()):
                for c in range(0, self.grdReport.NumberCols):
                    self.grdReport.SetCellBackgroundColour(r, c, "light blue")
            else:
                for c in range(0, self.grdReport.NumberCols):
                    self.grdReport.SetCellBackgroundColour(r, c, "white")
        self.grdReport.Refresh()
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
                self.btnAminoA.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnA.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoA.SetForegroundColour("#000000")
            logInfo("Removed ALA from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoA.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnA_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoA.SetForegroundColour("#FF0000")
            logInfo("Added ALA to palette")
        self.addAAToPalette("A")
    
    def aminoC(self, event):
        if ("C" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoC.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnC.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoC.SetForegroundColour("#000000")
            logInfo("Removed CYS from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoC.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnC_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoC.SetForegroundColour("#FF0000")
            logInfo("Added CYS to palette")
        self.addAAToPalette("C")
    
    def aminoD(self, event):
        if ("D" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoD.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnD.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoD.SetForegroundColour("#000000")
            logInfo("Removed ASP from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoD.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnD_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoD.SetForegroundColour("#FF0000")
            logInfo("Added ASP to palette")
        self.addAAToPalette("D")
    
    def aminoE(self, event):
        if ("E" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoE.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnE.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoE.SetForegroundColour("#000000")
            logInfo("Removed GLU from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoE.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnE_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoE.SetForegroundColour("#FF0000")
            logInfo("Added GLU to palette")
        self.addAAToPalette("E")
        
    def aminoF(self, event):
        if ("F" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoF.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnF.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoF.SetForegroundColour("#000000")
            logInfo("Removed PHE from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoF.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnF_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoF.SetForegroundColour("#FF0000")
            logInfo("Added PHE to palette")
        self.addAAToPalette("F")
        
    def aminoG(self, event):
        if ("G" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoG.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnG.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoG.SetForegroundColour("#000000")
            logInfo("Removed GLY from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoG.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnG_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoG.SetForegroundColour("#FF0000")
            logInfo("Added GLY to palette")
        self.addAAToPalette("G")
        
    def aminoH(self, event):
        if ("H" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoH.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnH.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoH.SetForegroundColour("#000000")
            logInfo("Removed HIS from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoH.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnH_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoH.SetForegroundColour("#FF0000")
            logInfo("Added HIS to palette")
        self.addAAToPalette("H")
        
    def aminoI(self, event):
        if ("I" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoI.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnI.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoI.SetForegroundColour("#000000")
            logInfo("Removed ILE from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoI.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnI_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoI.SetForegroundColour("#FF0000")
            logInfo("Added ILE to palette")
        self.addAAToPalette("I")
        
    def aminoK(self, event):
        if ("K" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoK.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnK.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoK.SetForegroundColour("#000000")
            logInfo("Removed LYS from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoK.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnK_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoK.SetForegroundColour("#FF0000")
            logInfo("Added LYS to palette")
        self.addAAToPalette("K")
        
    def aminoL(self, event):
        if ("L" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoL.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnL.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoL.SetForegroundColour("#000000")
            logInfo("Removed LEU from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoL.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnL_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoL.SetForegroundColour("#FF0000")
            logInfo("Added LEU to palette")
        self.addAAToPalette("L")
        
    def aminoM(self, event):
        if ("M" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoM.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnM.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoM.SetForegroundColour("#000000")
            logInfo("Removed MET from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoM.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnM_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoM.SetForegroundColour("#FF0000")
            logInfo("Added MET to palette")
        self.addAAToPalette("M")
    
    def aminoN(self, event):
        if ("N" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoN.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnN.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoN.SetForegroundColour("#000000")
            logInfo("Removed ASN from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoN.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnN_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoN.SetForegroundColour("#FF0000")
            logInfo("Added ASN to palette")
        self.addAAToPalette("N")
    
    def aminoP(self, event):
        if ("P" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoP.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnP.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoP.SetForegroundColour("#000000")
            logInfo("Removed PRO from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoP.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnP_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoP.SetForegroundColour("#FF0000")
            logInfo("Added PRO to palette")
        self.addAAToPalette("P")
    
    def aminoQ(self, event):
        if ("Q" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoQ.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnQ.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoQ.SetForegroundColour("#000000")
            logInfo("Removed GLN from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoQ.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnQ_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoQ.SetForegroundColour("#FF0000")
            logInfo("Added GLN to palette")
        self.addAAToPalette("Q")
        
    def aminoR(self, event):
        if ("R" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoR.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnR.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoR.SetForegroundColour("#000000")
            logInfo("Removed ARG from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoR.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnR_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoR.SetForegroundColour("#FF0000")
            logInfo("Added ARG to palette")
        self.addAAToPalette("R")
        
    def aminoS(self, event):
        if ("S" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoS.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnS.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoS.SetForegroundColour("#000000")
            logInfo("Removed SER from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoS.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnS_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoS.SetForegroundColour("#FF0000")
            logInfo("Added SER to palette")
        self.addAAToPalette("S")
        
    def aminoT(self, event):
        if ("T" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoT.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnT.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoT.SetForegroundColour("#000000")
            logInfo("Removed THR from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoT.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnT_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoT.SetForegroundColour("#FF0000")
            logInfo("Added THR to palette")
        self.addAAToPalette("T")
        
    def aminoV(self, event):
        if ("V" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoV.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnV.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoV.SetForegroundColour("#000000")
            logInfo("Removed VAL from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoV.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnV_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoV.SetForegroundColour("#FF0000")
            logInfo("Added VAL to palette")
        self.addAAToPalette("V")
        
    def aminoW(self, event):
        if ("W" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoW.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnW.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoW.SetForegroundColour("#000000")
            logInfo("Removed TRP from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoW.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnW_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoW.SetForegroundColour("#FF0000")
            logInfo("Added TRP to palette")
        self.addAAToPalette("W")
        
    def aminoY(self, event):
        if ("Y" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoY.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnY.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoY.SetForegroundColour("#000000")
            logInfo("Removed TYR from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoY.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnY_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
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
                self.btnAminoA.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnA.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoC.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnC.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoD.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnD.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoE.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnE.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoF.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnF.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoG.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnG.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoH.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnH.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoI.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnI.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoK.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnK.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoL.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnL.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoM.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnM.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoN.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnN.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoP.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnP.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoQ.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnQ.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoR.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnR.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoS.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnS.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoT.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnT.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoV.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnV.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoW.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnW.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoY.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnY.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
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
                self.btnAminoA.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnA_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoC.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnC_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoD.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnD_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoE.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnE_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoF.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnF_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoG.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnG_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoH.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnH_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoI.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnI_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoK.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnK_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoL.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnL_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoM.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnM_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoN.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnN_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoP.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnP_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoQ.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnQ_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoR.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnR_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoS.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnS_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoT.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnT_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoV.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnV_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoW.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnW_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoY.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnY_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
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
    
    def add(self, event):
        #self.activate()
        logInfo("Add button clicked")
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
        self.add(event)
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
        # Find the neighborhood view
        seqpos = selectedValue.split(":")[len(selectedValue.split(":"))-1]
        seqpos = seqpos.strip()
        seqpos = seqpos[1:] # The AA is index 0
        if (self.useDesignedSeq):
            # The label has the mutation at the end after a design, so we have to take that out
            seqpos = seqpos[0:len(seqpos)-1]
        firstmodel = model
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
        self.cmd.zoom("dessele")
        if (chain == " " or chain == "_"):
            self.cmd.select("dessele", "resi " + seqpos + " and model " + firstmodel)
        else:
            self.cmd.select("dessele", "resi " + seqpos + " and model " + firstmodel + " and chain " + chain)
        self.cmd.show("sticks", "dessele")
        self.cmd.set_bond("stick_radius", 0.25, "dessele")
        # Highlight this residue in PyMOL
        self.cmd.select("seqsele", "dessele")
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
    
    def toggleMutantType(self, event):
        if (self.mutantType == "Single"):
            self.mutantType = "Double"
            if (platform.system() == "Darwin"):
                self.btnMutantType.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnMutant_Double.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnMutantType.SetLabel("Double Mutants")
            self.btnMutantType.SetToolTipString("Search for stabilizing double point mutations")
        else:
            self.mutantType = "Single"
            if (platform.system() == "Darwin"):
                self.btnMutantType.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnMutant_Single.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnMutantType.SetLabel("Single Mutants")
            self.btnMutantType.SetToolTipString("Search for stabilizing single point mutations")
        logInfo("Changed default behavior to " + self.mutantType)
    
    def dumpMutantList(self, filename):
        poseindx = self.resfile[0][3]
        model = self.seqWin.getModelForChain(poseindx)
        f = open(filename, "w")
        f.write("#MODEL " + model + "\n")
        if (self.mutantType == "Single"):
            for i in range(0, len(self.resfile)):
                [indx, r, seqpos, poseindx, chainID, chainoffset, reslist] = self.resfile[i]
                if (len(chainID.strip()) == 0):
                    chainID = "_"
                currtype = self.seqWin.sequences[r][indx]
                for restype in reslist:
                    if (restype != currtype):
                        f.write(chainID + " " + currtype + " " + str(seqpos) + " " + restype + "\n")
        else:
            for i in range(0, len(self.resfile)-1):
                [indx, r, seqpos, poseindx, chainID, chainoffset, reslist] = self.resfile[i]
                if (len(chainID.strip()) == 0):
                    chainID = "_"
                currtype = self.seqWin.sequences[r][indx]
                for j in range(i+1, len(self.resfile)):
                    [indx2, r2, seqpos2, poseindx2, chainID2, chainoffset2, reslist2] = self.resfile[j]
                    if (len(chainID2.strip()) == 0):
                        chainID2 = "_"
                    currtype2 = self.seqWin.sequences[r2][indx2]
                    for restype in reslist:
                        if (restype != currtype):
                            for restype2 in reslist2:
                                if (restype2 != currtype2):
                                    f.write(chainID + " " + currtype + " " + str(seqpos) + " " + restype + " " + chainID2 + " " + currtype2 + " " + str(seqpos2) + " " + restype2 + "\n")
        f.close()
        logInfo("Wrote a mutant list file to " + filename)
    
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
            self.btnMutantType.Enable()
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
            self.btnMutantType.Disable()
    
    def saveReport(self, event):
        logInfo("Clicked the Save Report button")
        # Save the data in the resfile graph as an actual resfile for the user's convenience
        dlg = wx.FileDialog(
            self, message="Save a Scan File",
            defaultDir=self.seqWin.cwd,
            defaultFile="",
            wildcard="Point Mutant Scan (*.scan)|*.scan",
            style=wx.SAVE | wx.CHANGE_DIR)
        if (dlg.ShowModal() == wx.ID_OK):
            paths = dlg.GetPaths()
            # Change cwd to the last opened file
            if (platform.system() == "Windows"):
                lastDirIndx = paths[len(paths)-1].rfind("\\")
            else:
                lastDirIndx = paths[len(paths)-1].rfind("/")
            self.seqWin.cwd = str(paths[len(paths)-1][0:lastDirIndx])
            self.seqWin.saveWindowData(None)
            filename = str(paths[0]).split(".scan")[0] + ".scan"
            # Does it exist already?  If so, ask if the user really wants to overwrite it
            if (os.path.isfile(filename)):
                dlg2 = wx.MessageDialog(self, "The file " + filename + " already exists.  Overwrite it?", "Filename Already Exists", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
                if (dlg2.ShowModal() == wx.ID_NO):
                    dlg2.Destroy()
                    logInfo("Canceled save operation due to filename already existing")
            goToSandbox()
            fin = open("scanoutput", "r")
            fout = open(filename, "w")
            for aline in fin:
                fout.write(aline)
            fout.close()
            fin.close()
        else:
            logInfo("Cancelled save constraints operation")
    
    def loadReport(self, event):
        logInfo("Load Report button clicked")
        dlg = wx.FileDialog(
            self, message="Choose a File",
            defaultDir=self.seqWin.cwd,
            defaultFile="",
            wildcard="Point Mutant Scans (*.scan)|*.scan",
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
            logInfo("Loaded data from a point mutant scan file", filename)
            f = open(filename, "r")
            self.ddGData = []
            offset = 0
            model = "Unknown"
            for aline in f:
                if ("divide_up_pdbs" in aline):
                    offset = 1
                if ("go()" in aline or "mutation" in aline or "main()" in aline):
                    # Header lines, skip it
                    continue
                if (aline.startswith("#MODEL")):
                    model = aline.split("\t")[1].strip()
                try:
                    mutation = str(aline.split()[1+offset]).replace("-", "|")
                    ddG = float(aline.split()[3+offset])
                except:
                    continue
                self.ddGData.append((mutation, ddG))
            f.close()
            # Sort the data so the most optimal mutations are at the front
            for i in range(0, len(self.ddGData)-1):
                best = i
                for j in range(i+1, len(self.ddGData)):
                    if (self.ddGData[j][1] < self.ddGData[best][1]):
                        best = j
                temp = self.ddGData[best]
                self.ddGData[best] = self.ddGData[i]
                self.ddGData[i] = temp
            # Load the data into the report grid
            if (self.grdReport.NumberRows > 0):
                self.grdReport.DeleteRows(0, self.grdReport.NumberRows)
            self.grdReport.AppendRows(len(self.ddGData)+1)
            self.grdReport.SetRowLabelValue(0, "Model:")
            self.grdReport.SetCellValue(0, 0, model)
            self.grdReport.SetCellAlignment(0, 0, wx.ALIGN_CENTRE, wx.ALIGN_CENTRE)
            readOnly = wx.grid.GridCellAttr()
            readOnly.SetReadOnly(True)
            self.grdReport.SetRowAttr(0, readOnly)
            for i in range(0, len(self.ddGData)):
                self.grdReport.SetRowLabelValue(i+1, self.ddGData[i][0])
                self.grdReport.SetCellValue(i+1, 0, str(self.ddGData[i][1]))
                self.grdReport.SetCellAlignment(i+1, 0, wx.ALIGN_CENTRE, wx.ALIGN_CENTRE)
                readOnly = wx.grid.GridCellAttr()
                readOnly.SetReadOnly(True)
                self.grdReport.SetRowAttr(i+1, readOnly)
    
    def serverToggle(self, event):
        if (self.serverOn):
            self.serverOn = False
            if (platform.system() == "Darwin"):
                self.btnServerToggle.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnServer_Off.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnServerToggle.SetLabel("Server Off")
            self.btnServerToggle.SetToolTipString("Perform point mutation scanning locally")
            logInfo("Turned off point mutation scanning server usage")
        else:
            self.serverOn = True
            if (platform.system() == "Darwin"):
                self.btnServerToggle.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnServer_On.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnServerToggle.SetLabel("Server On")
            self.btnServerToggle.SetToolTipString("Perform point mutation scanning on a remote server")
            logInfo("Turned on point mutation scanning server usage")
    
    def cancelScan(self):
        logInfo("Canceled point mutant scan operation")
        try:
            os.remove("scaninput")
        except:
            pass
        try:
            os.remove("scaninputtemp")
        except:
            pass
        self.tmrScan.Stop()
        self.seqWin.cannotDelete = False
        self.enableControls()
        if (platform.system() == "Darwin"):
            self.btnScan.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnScan.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
        else:
            self.btnScan.SetLabel("Scan!")
        self.btnScan.SetToolTipString("Perform point mutant scan")
        deleteInputFiles()
        self.parent.parent.restartDaemon()
        self.parent.GoBtn.Enable()
        # Get rid of the messages
        for i in range(0, len(self.seqWin.msgQueue)):
            if (self.seqWin.msgQueue[i].find("Performing point mutant scan") >= 0):
                self.seqWin.msgQueue.pop(i)
                break
        if (len(self.seqWin.msgQueue) > 0):
            self.seqWin.labelMsg.SetLabel(self.seqWin.msgQueue[len(self.seqWin.msgQueue)-1])
        else:
            self.seqWin.labelMsg.SetLabel("")
        self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")
        self.buttonState = "Scan!"
    
    def scanClick(self, event):
        # This is also the "Finalize!" button
        logInfo("Scan button clicked")
        if (self.buttonState == "Scan!"):
            if (len(self.resfile) > 0):
                self.seqWin.labelMsg.SetLabel("Performing point mutant scan, please be patient...")
                self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
                self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")
                self.seqWin.msgQueue.append("Performing point mutant scan, please be patient...")
                self.parent.GoBtn.Disable()
                self.enableControls(False)
                self.seqWin.cannotDelete = True
                self.stage = 1
                if (platform.system() == "Darwin"):
                    self.btnScan.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnScan_Cancel.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                else:
                    self.btnScan.SetLabel("Cancel!")
                self.buttonState = "Cancel!"
                self.btnScan.SetToolTipString("Cancel the point mutant scan")
                self.tmrScan = wx.Timer(self)
                self.Bind(wx.EVT_TIMER, self.threadScan, self.tmrScan)
                self.tmrScan.Start(1000)
            else:
                wx.MessageBox("There's nothing to design!", "Nothing to Design", wx.OK|wx.ICON_EXCLAMATION)
        elif (self.buttonState == "Cancel!"):
            dlg = wx.MessageDialog(self, "Are you sure you want to cancel the point mutant scan?  All progress will be lost.", "Cancel Point Mutant Scan", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
            result = dlg.ShowModal()
            if (result == wx.ID_YES):
                self.cancelScan()
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
                self.btnDesign.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnDesign.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
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
            self.btnScan.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnScan.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
        else:
            self.btnScan.SetLabel("Scan!")
        self.buttonState = "Scan!"
        # Get rid of the messages
        for i in range(0, len(self.seqWin.msgQueue)):
            if (self.seqWin.msgQueue[i].find("Performing point mutant scan") >= 0):
                self.seqWin.msgQueue.pop(i)
                break
        if (len(self.seqWin.msgQueue) > 0):
            self.seqWin.labelMsg.SetLabel(self.seqWin.msgQueue[len(self.seqWin.msgQueue)-1])
        else:
            self.seqWin.labelMsg.SetLabel("")
        self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")
    
    def threadScan(self, event):
        # Why am I using a Timer?  See the explanation in kic.py
        # Save the resfile to a temporary location that the daemon can find
        goToSandbox()
        if (self.stage == 1):
            self.tmrScan.Stop()
            self.timeoutCount = 0
            self.dumpMutantList("mutant.list")
            f = open("scaninputtemp", "w")
            r = self.resfile[0][1]
            self.selectedModel = self.seqWin.getModelForChain(r)
            pdbfile = self.selectedModel + ".pdb"
            # Dump the PDB from PyMOL first in case the coordinates were altered by the user
            dumpmodel = pdbfile.split(".pdb")[0]
            self.cmd.save(pdbfile.strip(), "model " + dumpmodel)
            fixPyMOLSave(pdbfile.strip())
            # PMutScan does not accept empty chainIDs, so we have to assign one
            # Find taken chainIDs
            f2 = open(pdbfile.strip(), "r")
            pdbdata = []
            takenIDs = []
            for aline in f2:
                pdbdata.append(aline.strip())
                if (aline.startswith("ATOM") or aline.startswith("HETATM")):
                    if (not(aline[21].strip() in takenIDs)):
                        takenIDs.append(aline[21])
            f2.close()
            # Find a new one
            for chainID in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
                if (not(chainID in takenIDs)):
                    self.replaceID = chainID
                    break
            # Fix the PDB chainIDs
            f2 = open(pdbfile.strip(), "w")
            for aline in pdbdata:
                if ((aline.startswith("ATOM") or aline.startswith("HETATM")) and aline[21] == " "):
                    aline = aline[0:21] + self.replaceID + aline[22:]
                f2.write(aline.strip() + "\n")
            f2.close()
            # Fix the mutant.list file
            f2 = open("mutant.list", "r")
            mutantdata = []
            for aline in f2:
                mutantdata.append(aline.strip())
            f2.close()
            f2 = open("mutant.list", "w")
            for aline in mutantdata:
                if (aline[0] == "_"):
                    aline = self.replaceID + aline[1:]
                try:
                    indx = aline.index("_")
                    aline = aline[0:indx] + self.replaceID + aline[indx+1:]
                except:
                    pass
                f2.write(aline.strip() + "\n")
            f2.close()
            f.write("PDBFILE\t" + pdbfile.strip() + "\n")
            f2 = open(pdbfile, "r")
            f.write("BEGIN PDB DATA\n")
            for aline in f2:
                f.write(aline.strip() + "\n")
            f.write("END PDB DATA\n")
            f2.close()
            f.write("LIST\tmutant.list\n")
            f2 = open("mutant.list", "r")
            f.write("BEGIN LIST DATA\n")
            for aline in f2:
                f.write(aline.strip() + "\n")
            f.write("END LIST DATA\n")
            f2.close()
            f.write("MUTANTTYPE\t" + self.mutantType.upper() + "\n")
            f.close()
            appendScorefxnParamsInfoToFile("scaninputtemp", self.selectWin.weightsfile)
            if (self.serverOn):
                #try: 
                self.ID = sendToServer("scaninput")
                dlg = wx.TextEntryDialog(None, "Enter a description for this submission:", "Job Description", "")
                if (dlg.ShowModal() == wx.ID_OK):
                    desc = dlg.GetValue()
                    desc = desc.replace("\t", " ").replace("\n", " ").strip()
                else:
                    desc = self.ID
                # First make sure this isn't a duplicate
                alreadythere = False
                try:
                    f = open("downloadwatch", "r")
                    for aline in f:
                        if (len(aline.split("\t")) >= 2 and aline.split("\t")[0] == "PMUTSCAN" and aline.split("\t")[1] == self.ID.strip()):
                            alreadythere = True
                            break
                    f.close()
                except:
                    pass
                if (not(alreadythere)):
                    f = open("downloadwatch", "a")
                    f.write("PMUTSCAN\t" + self.ID.strip() + "\t" + str(datetime.datetime.now().strftime("%A, %B %d - %I:%M:%S %p")) + "\t" + getServerName() + "\t" + desc + "\n")
                    f.close()
                dlg = wx.MessageDialog(self, "InteractiveROSETTA is now watching the server for job ID " + desc.strip() + ".  You will be notified when the package is available for download.", "Listening for Download", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
                dlg.ShowModal()
                dlg.Destroy()
                # Re-enable everything since we're not waiting for the local daemon to do anything
                # Pop this message out of the queue
                for i in range(0, len(self.seqWin.msgQueue)):
                    if (self.seqWin.msgQueue[i].find("Performing point mutant scan") >= 0):
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
                if (platform.system() == "Darwin"):
                    self.btnScan.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnScan.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                else:
                    self.btnScan.SetLabel("Scan!")
                self.buttonState = "Scan!"
                self.btnScan.SetToolTipString("Perform a point mutant scan")
                logInfo("Point mutant scanning input sent to server daemon with ID " + self.ID)
                return
                #except:
                    #dlg = wx.MessageDialog(self, "The server could not be reached!  Ensure that you have specified a valid server and that you have an network connection.", "Server Could Not Be Reached", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
                    #dlg.ShowModal()
                    #dlg.Destroy()
                    #return
            else:
                os.rename("scaninputtemp", "scaninput")
                self.usingServer = False
                logInfo("Scanning input uploaded locally at scaninput")
            # Now keep a list of every mutant that will be attempted so we can easily update the progress bar
            self.mutants = []
            if (self.mutantType.upper() == "DOUBLE"):
                for i in range(0, self.grdResfile.NumberRows-1):
                    for j in range(i+1, self.grdResfile.NumberRows):
                        for resi in self.grdResfile.GetCellValue(i, 0):
                            for resj in self.grdResfile.GetCellValue(j, 0):
                                if (self.grdResfile.GetCellValue(i, 1).split("|")[1] == "_"):
                                    key1 = self.replaceID + "-" + str(self.grdResfile.GetRowLabelValue(i).split()[1]) + str(resi)
                                else:
                                    key1 = str(self.grdResfile.GetCellValue(i, 1).split("|")[1]) + "-" + str(self.grdResfile.GetRowLabelValue(i).split()[1]) + str(resi)
                                if (self.grdResfile.GetCellValue(j, 1).split("|")[1] == "_"):
                                    key2 = self.replaceID + "-" + str(self.grdResfile.GetRowLabelValue(j).split()[1]) + str(resj)
                                else:
                                    key2 = str(self.grdResfile.GetCellValue(j, 1).split("|")[1]) + "-" + str(self.grdResfile.GetRowLabelValue(j).split()[1]) + str(resj)
                                self.mutants.append(key1 + "," + key2)
            else:
                for i in range(0, self.grdResfile.NumberRows):
                    for resi in self.grdResfile.GetCellValue(i, 0):
                        if (self.grdResfile.GetCellValue(i, 1).split("|")[1] == "_"):
                            self.mutants.append(self.replaceID + "-" + str(self.grdResfile.GetRowLabelValue(i).split()[1]) + str(resi))
                        else:
                            self.mutants.append(str(self.grdResfile.GetCellValue(i, 1).split("|")[1]) + "-" + str(self.grdResfile.GetRowLabelValue(i).split()[1]) + str(resi))
            self.progress = wx.ProgressDialog("Point Mutation Scan", "Performing a point mutation scan...", len(self.mutants)+1, style=wx.PD_CAN_ABORT | wx.PD_APP_MODAL | wx.PD_ELAPSED_TIME | wx.PD_REMAINING_TIME)
            self.stage = 2
            self.tmrScan.Start(1000)
        else:
            if (self.usingServer):
                # See if the file has been uploaded yet and bring it here if so
                queryServerForResults("scanoutput-" + self.ID)
                self.timeoutCount = self.timeoutCount + 1
            if (self.timeoutCount >= serverTimeout):
                self.tmrScan.Stop()
                # If this is taking too long, maybe there's something wrong with the server
                # Ask the user if they want to continue waiting or use the local daemon instead
                dlg = wx.MessageDialog(self, "The server is taking a long time to respond.  Continue to wait?  Pressing No will run the calculations locally.", "Delayed Server Response", wx.YES_NO | wx.ICON_EXCLAMATION | wx.CENTRE)
                if (dlg.ShowModal() == wx.ID_YES):
                    # Reset the counter
                    self.timeoutCount = 0
                else:
                    self.usingServer = False
                    self.timeoutCount = 0
                    os.rename("scaninputtemp", "scaninput")
                    logInfo("Server took too long to respond so the local daemon was used")
                dlg.Destroy()
                self.tmrScan.Start(1000)
            # Read the output dumped by the child process
            if (os.path.isfile("scanoutput")):
                try:
                    self.progress.Destroy()
                except:
                    pass
                self.tmrScan.Stop()
                self.residue_E = []
                pdbreader = Bio.PDB.PDBParser()
                f = open("scanoutput", "r")
                self.ddGData = []
                offset = 0
                model = "Unknown"
                for aline in f:
                    if ("go()" in aline or "mutation" in aline or "main()" in aline):
                        # Header lines, skip it
                        continue
                    if ("divide_up_pdbs" in aline):
                        offset = 1
                    if (aline.startswith("#MODEL")):
                        model = aline.split("\t")[1].strip()
                    try:
                        mutation = str(aline.split()[1+offset]).replace("-", "|")
                        if (mutation[0] == self.replaceID):
                            mutation = "_" + mutation[1:] # Get blank chainIDs back
                        ddG = float(aline.split()[3+offset])
                    except:
                        continue
                    self.ddGData.append((mutation, ddG))
                f.close()
                # Sort the data so the most optimal mutations are at the front
                for i in range(0, len(self.ddGData)-1):
                    best = i
                    for j in range(i+1, len(self.ddGData)):
                        if (self.ddGData[j][1] < self.ddGData[best][1]):
                            best = j
                    temp = self.ddGData[best]
                    self.ddGData[best] = self.ddGData[i]
                    self.ddGData[i] = temp
                # Load the data into the report grid
                if (self.grdReport.NumberRows > 0):
                    self.grdReport.DeleteRows(0, self.grdReport.NumberRows)
                self.grdReport.AppendRows(len(self.ddGData)+1)
                self.grdReport.SetRowLabelValue(0, "Model:")
                self.grdReport.SetCellValue(0, 0, model)
                self.grdReport.SetCellAlignment(0, 0, wx.ALIGN_CENTRE, wx.ALIGN_CENTRE)
                readOnly = wx.grid.GridCellAttr()
                readOnly.SetReadOnly(True)
                self.grdReport.SetRowAttr(0, readOnly)
                for i in range(0, len(self.ddGData)):
                    self.grdReport.SetRowLabelValue(i+1, self.ddGData[i][0])
                    self.grdReport.SetCellValue(i+1, 0, str(self.ddGData[i][1]))
                    self.grdReport.SetCellAlignment(i+1, 0, wx.ALIGN_CENTRE, wx.ALIGN_CENTRE)
                    readOnly = wx.grid.GridCellAttr()
                    readOnly.SetReadOnly(True)
                    self.grdReport.SetRowAttr(i+1, readOnly)
                logInfo("Found point mutation scanning output at scanoutput")
                # Pop this message out of the queue
                for i in range(0, len(self.seqWin.msgQueue)):
                    if (self.seqWin.msgQueue[i].find("Performing point mutant scan") >= 0):
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
                self.btnSaveReport.Enable()
                if (platform.system() == "Darwin"):
                    self.btnScan.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pmutscan/btnScan.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                else:
                    self.btnScan.SetLabel("Scan!")
                self.buttonState = "Scan!"
                self.btnScan.SetToolTipString("Perform a point mutant scan")
            elif (os.path.isfile("errreport")):
                self.tmrScan.Stop()
                self.recoverFromError()
            elif (os.path.isfile("scanprogress")):
                f = open("scanprogress", "r")
                data = f.readlines()
                f.close()
                if (len(data) == 0):
                    return
                try:
                    lastline = data[len(data)-1].strip()
                    if (self.mutantType.upper() == "DOUBLE"):
                        curpos1 = lastline.split()[1].split(",")[0]
                        curpos2 = lastline.split()[1].split(",")[1]
                        curpos = curpos1 + "," + curpos2
                    else:
                        curpos = lastline.split()[1].strip()
                    indx = self.mutants.index(curpos)
                except:
                    return
                if (indx == len(self.mutants)+1):
                    try:
                        self.progress.Destroy()
                    except:
                        pass
                else:
                    (keepGoing, skip) = self.progress.Update(indx)
                    if (not(keepGoing)):
                        # User clicked "Cancel" on the progress bar
                        self.cancelScan()
                        self.progress.Destroy()
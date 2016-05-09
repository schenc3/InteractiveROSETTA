import wx
import wx.grid
import wx.lib.scrolledpanel
import wx.stc
import os
import os.path
import time
import platform
import multiprocessing
import datetime
import Bio.PDB
import webbrowser
from threading import Thread
from tools import *

class MSDPanel(wx.lib.scrolledpanel.ScrolledPanel):
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
            self.lblProt = wx.StaticText(self, -1, "Multi-State Design", (25, 15), (270, 25), wx.ALIGN_CENTRE)
            self.lblProt.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblProt = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/lblMSD.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(25, 15), size=(270, 25))
        else:
            self.lblProt = wx.StaticText(self, -1, "Multi-State Design", (70, 15), style=wx.ALIGN_CENTRE)
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
            self.lblInst = wx.StaticText(self, -1, "Generate MSD inputs to run on a server", (0, 45), (320, 25), wx.ALIGN_CENTRE)
            self.lblInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
        elif (platform.system() == "Darwin"):
            self.lblInst = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/lblInstMSD.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 45), size=(320, 25))
        else:
            self.lblInst = wx.StaticText(self, -1, "Generate MSD inputs to run on a server", (5, 45), style=wx.ALIGN_CENTRE)
            self.lblInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
            resizeTextControlForUNIX(self.lblInst, 0, self.GetSize()[0])
        self.lblInst.SetForegroundColour("#FFFFFF")
        
        if (platform.system() == "Windows"):
            self.lblEntity = wx.StaticText(self, -1, "Entity MSD Resfile", (0, 70), (320, 25), wx.ALIGN_CENTRE)
            self.lblEntity.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblEntity = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/lblEntity.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(25, 70), size=(270, 25))
        else:
            self.lblEntity = wx.StaticText(self, -1, "Entity MSD Resfile", (70, 70), style=wx.ALIGN_CENTRE)
            self.lblEntity.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
            resizeTextControlForUNIX(self.lblEntity, 0, self.GetSize()[0])
        self.lblEntity.SetForegroundColour("#FFFFFF")
        if (platform.system() == "Windows"):
            self.lblEntityInst = wx.StaticText(self, -1, "Specify which state positions are equivalent", (0, 100), (320, 25), wx.ALIGN_CENTRE)
            self.lblEntityInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
        elif (platform.system() == "Darwin"):
            self.lblEntityInst = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/lblInstEntity.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 100), size=(320, 25))
        else:
            self.lblEntityInst = wx.StaticText(self, -1, "Specify which state positions are equivalent", (5, 100), style=wx.ALIGN_CENTRE)
            self.lblEntityInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
            resizeTextControlForUNIX(self.lblEntityInst, 0, self.GetSize()[0])
        self.lblEntityInst.SetForegroundColour("#FFFFFF")
        
        if (platform.system() == "Darwin"):
            self.btnAminoA = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnA.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(7, 130), size=(42, 25))
        else:
            self.btnAminoA = wx.Button(self, id=-1, label="A", pos=(7, 130), size=(42, 25))
            self.btnAminoA.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoA.Bind(wx.EVT_BUTTON, self.aminoA)
        self.btnAminoA.SetToolTipString("Add ALA to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoC = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnC.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(51, 130), size=(42, 25))
        else:
            self.btnAminoC = wx.Button(self, id=-1, label="C", pos=(51, 130), size=(42, 25))
            self.btnAminoC.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoC.Bind(wx.EVT_BUTTON, self.aminoC)
        self.btnAminoC.SetToolTipString("Add CYS to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoD = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnD.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(95, 130), size=(42, 25))
        else:
            self.btnAminoD = wx.Button(self, id=-1, label="D", pos=(95, 130), size=(42, 25))
            self.btnAminoD.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoD.Bind(wx.EVT_BUTTON, self.aminoD)
        self.btnAminoD.SetToolTipString("Add ASP to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoE = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnE.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(139, 130), size=(42, 25))
        else:
            self.btnAminoE = wx.Button(self, id=-1, label="E", pos=(139, 130), size=(42, 25))
            self.btnAminoE.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoE.Bind(wx.EVT_BUTTON, self.aminoE)
        self.btnAminoE.SetToolTipString("Add GLU to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoF = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnF.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(183, 130), size=(42, 25))
        else:
            self.btnAminoF = wx.Button(self, id=-1, label="F", pos=(183, 130), size=(42, 25))
            self.btnAminoF.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoF.Bind(wx.EVT_BUTTON, self.aminoF)
        self.btnAminoF.SetToolTipString("Add PHE to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoG = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnG.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(227, 130), size=(42, 25))
        else:
            self.btnAminoG = wx.Button(self, id=-1, label="G", pos=(227, 130), size=(42, 25))
            self.btnAminoG.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoG.Bind(wx.EVT_BUTTON, self.aminoG)
        self.btnAminoG.SetToolTipString("Add GLY to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoH = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnH.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(271, 130), size=(42, 25))
        else:
            self.btnAminoH = wx.Button(self, id=-1, label="H", pos=(271, 130), size=(42, 25))
            self.btnAminoH.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoH.Bind(wx.EVT_BUTTON, self.aminoH)
        self.btnAminoH.SetToolTipString("Add HIS to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoI = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnI.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(7, 160), size=(42, 25))
        else:
            self.btnAminoI = wx.Button(self, id=-1, label="I", pos=(7, 160), size=(42, 25))
            self.btnAminoI.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoI.Bind(wx.EVT_BUTTON, self.aminoI)
        self.btnAminoI.SetToolTipString("Add ILE to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoK = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnK.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(51, 160), size=(42, 25))
        else:
            self.btnAminoK = wx.Button(self, id=-1, label="K", pos=(51, 160), size=(42, 25))
            self.btnAminoK.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoK.Bind(wx.EVT_BUTTON, self.aminoK)
        self.btnAminoK.SetToolTipString("Add LYS to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoL = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnL.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(95, 160), size=(42, 25))
        else:
            self.btnAminoL = wx.Button(self, id=-1, label="L", pos=(95, 160), size=(42, 25))
            self.btnAminoL.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoL.Bind(wx.EVT_BUTTON, self.aminoL)
        self.btnAminoL.SetToolTipString("Add LEU to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoM = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnM.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(139, 160), size=(42, 25))
        else:
            self.btnAminoM = wx.Button(self, id=-1, label="M", pos=(139, 160), size=(42, 25))
            self.btnAminoM.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoM.Bind(wx.EVT_BUTTON, self.aminoM)
        self.btnAminoM.SetToolTipString("Add MET to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoN = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnN.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(183, 160), size=(42, 25))
        else:
            self.btnAminoN = wx.Button(self, id=-1, label="N", pos=(183, 160), size=(42, 25))
            self.btnAminoN.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoN.Bind(wx.EVT_BUTTON, self.aminoN)
        self.btnAminoN.SetToolTipString("Add ASN to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoP = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnP.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(227, 160), size=(42, 25))
        else:
            self.btnAminoP = wx.Button(self, id=-1, label="P", pos=(227, 160), size=(42, 25))
            self.btnAminoP.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoP.Bind(wx.EVT_BUTTON, self.aminoP)
        self.btnAminoP.SetToolTipString("Add PRO to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoQ = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnQ.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(271, 160), size=(42, 25))
        else:
            self.btnAminoQ = wx.Button(self, id=-1, label="Q", pos=(271, 160), size=(42, 25))
            self.btnAminoQ.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoQ.Bind(wx.EVT_BUTTON, self.aminoQ)
        self.btnAminoQ.SetToolTipString("Add GLN to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoR = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnR.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(7, 190), size=(42, 25))
        else:
            self.btnAminoR = wx.Button(self, id=-1, label="R", pos=(7, 190), size=(42, 25))
            self.btnAminoR.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoR.Bind(wx.EVT_BUTTON, self.aminoR)
        self.btnAminoR.SetToolTipString("Add ARG to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoS = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnS.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(51, 190), size=(42, 25))
        else:
            self.btnAminoS = wx.Button(self, id=-1, label="S", pos=(51, 190), size=(42, 25))
            self.btnAminoS.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoS.Bind(wx.EVT_BUTTON, self.aminoS)
        self.btnAminoS.SetToolTipString("Add SER to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoT = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnT.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(95, 190), size=(42, 25))
        else:
            self.btnAminoT = wx.Button(self, id=-1, label="T", pos=(95, 190), size=(42, 25))
            self.btnAminoT.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoT.Bind(wx.EVT_BUTTON, self.aminoT)
        self.btnAminoT.SetToolTipString("Add THR to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoV = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnV.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(139, 190), size=(42, 25))
        else:
            self.btnAminoV = wx.Button(self, id=-1, label="V", pos=(139, 190), size=(42, 25))
            self.btnAminoV.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoV.Bind(wx.EVT_BUTTON, self.aminoV)
        self.btnAminoV.SetToolTipString("Add VAL to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoW = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnW.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(183, 190), size=(42, 25))
        else:
            self.btnAminoW = wx.Button(self, id=-1, label="W", pos=(183, 190), size=(42, 25))
            self.btnAminoW.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoW.Bind(wx.EVT_BUTTON, self.aminoW)
        self.btnAminoW.SetToolTipString("Add TRP to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoY = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnY.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(227, 190), size=(42, 25))
        else:
            self.btnAminoY = wx.Button(self, id=-1, label="Y", pos=(227, 190), size=(42, 25))
            self.btnAminoY.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoY.Bind(wx.EVT_BUTTON, self.aminoY)
        self.btnAminoY.SetToolTipString("Add TYR to the design palette")
        if (platform.system() == "Darwin"):
            self.btnAminoX = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnX.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(271, 190), size=(42, 25))
        else:
            self.btnAminoX = wx.Button(self, id=-1, label="X", pos=(271, 190), size=(42, 25))
            self.btnAminoX.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAminoX.Bind(wx.EVT_BUTTON, self.aminoX)
        self.btnAminoX.SetToolTipString("Add all amino acids to the design palette")
        self.palette = ""
        
        self.states = []
        if (platform.system() == "Darwin"):
            self.btnAddEntity = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnAddEntity.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 220), size=(100, 25))
        else:
            self.btnAddEntity = wx.Button(self, id=-1, label="Add Entry", pos=(0, 220), size=(100, 25))
            self.btnAddEntity.SetForegroundColour("#000000")
            self.btnAddEntity.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAddEntity.Bind(wx.EVT_BUTTON, self.addEntityPosition)
        self.btnAddEntity.SetToolTipString("Add selected residues to the resfile")
        self.entityMenu = wx.ComboBox(self, pos=(110, 220), size=(100, 25), choices=[], style=wx.CB_READONLY)
        self.entityMenu.SetToolTipString("Select entity resfile entries to remove")
        if (platform.system() == "Darwin"):
            self.btnRemoveEntity = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnRemoveEntity.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(220, 220), size=(100, 25))
        else:
            self.btnRemoveEntity = wx.Button(self, id=-1, label="Remove", pos=(220, 220), size=(100, 25))
            self.btnRemoveEntity.SetForegroundColour("#000000")
            self.btnRemoveEntity.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnRemoveEntity.Bind(wx.EVT_BUTTON, self.removeEntityPosition)
        self.btnRemoveEntity.SetToolTipString("Remove selected residues from the resfile")
        
        if (platform.system() != "Linux"):
            self.lblModel = wx.StaticText(self, -1, "None Selected", (0, 253), (160, 25), wx.ALIGN_CENTRE)
            self.lblModel.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
        else:
            self.lblModel = wx.StaticText(self, -1, "None Selected", (0, 253), style=wx.ALIGN_CENTRE)
            self.lblModel.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
            resizeTextControlForUNIX(self.lblModel, 0, 160)
        self.lblModel.SetForegroundColour("#FFFFFF")
        self.posMenu = wx.ComboBox(self, pos=(160, 250), size=(160, 25), choices=[], style=wx.CB_READONLY)
        self.posMenu.Bind(wx.EVT_COMBOBOX, self.entityPosMenuSelect)
        self.posMenu.SetToolTipString("Specify the residue for this entity state on the indicated model")
        self.selectedModel = ""
        
        self.grdEntity = wx.grid.Grid(self)
        self.grdEntity.CreateGrid(0, 1)
        self.grdEntity.SetSize((320, 250))
        self.grdEntity.SetPosition((0, 280))
        self.grdEntity.SetLabelFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.grdEntity.DisableDragColSize()
        self.grdEntity.DisableDragRowSize()
        self.grdEntity.SetColLabelValue(0, "Residue Options")
        self.grdEntity.SetRowLabelSize(40)
        self.grdEntity.SetColSize(0, 190)
        self.grdEntity.Bind(wx.grid.EVT_GRID_CELL_LEFT_CLICK, self.entityClick)
        self.grdEntity.Bind(wx.grid.EVT_GRID_CELL_RIGHT_CLICK, self.entityRClick)
        self.grdRow = -1
        self.grdCol = -1
        self.entities = []
        self.stateoptions = []
        
        self.AAlist = ["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"]
        self.AAlist.sort()
        if (platform.system() == "Darwin"):
            self.resEntityMenu = wx.ComboBox(self, pos=(7, 540), size=(119, 25), choices=self.AAlist, style=wx.CB_READONLY)
        else:
            self.resEntityMenu = wx.ComboBox(self, pos=(7, 540), size=(119, 25), choices=self.AAlist, style=wx.CB_READONLY | wx.CB_SORT)
        self.resEntityMenu.SetToolTipString("Amino acid type to add/remove from the selected entity resfile entry")
        if (platform.system() == "Darwin"):
            self.btnAddRes = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnAddRes.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(131, 540), size=(88, 25))
        else:
            self.btnAddRes = wx.Button(self, id=-1, label="Add", pos=(131, 540), size=(88, 25))
            self.btnAddRes.SetForegroundColour("#000000")
            self.btnAddRes.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAddRes.Bind(wx.EVT_BUTTON, self.addRes)
        self.btnAddRes.SetToolTipString("Add the selected amino acid to the selected entity resfile entry")
        if (platform.system() == "Darwin"):
            self.btnRemoveRes = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnRemoveRes.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(224, 540), size=(88, 25))
        else:
            self.btnRemoveRes = wx.Button(self, id=-1, label="Remove", pos=(224, 540), size=(88, 25))
            self.btnRemoveRes.SetForegroundColour("#000000")
            self.btnRemoveRes.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnRemoveRes.Bind(wx.EVT_BUTTON, self.removeRes)
        self.btnRemoveRes.SetToolTipString("Remove the selected amino acid from the selected entity resfile entry")
        
        if (platform.system() == "Windows"):
            self.lblSecondary = wx.StaticText(self, -1, "Secondary Resfiles (Optional)", (0, 570), (320, 25), wx.ALIGN_CENTRE)
            self.lblSecondary.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblSecondary = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/lblSecondary.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(25, 570), size=(270, 25))
        else:
            self.lblSecondary = wx.StaticText(self, -1, "Secondary Resfiles (Optional)", (70, 570), style=wx.ALIGN_CENTRE)
            self.lblSecondary.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
            resizeTextControlForUNIX(self.lblSecondary, 0, self.GetSize()[0])
        self.lblSecondary.SetForegroundColour("#FFFFFF")
        if (platform.system() == "Windows"):
            self.lblSecondaryInst = wx.StaticText(self, -1, "Specify if other residues should be repacked\non each PDB. Defaults to NATRO.", (0, 600), (320, 25), wx.ALIGN_CENTRE)
            self.lblSecondaryInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
        elif (platform.system() == "Darwin"):
            self.lblSecondaryInst = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/lblInstSecondary.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 600), size=(320, 35))
        else:
            self.lblSecondaryInst = wx.StaticText(self, -1, "Specify if other residues should be repacked\non each PDB. Defaults to NATRO.", (5, 600), style=wx.ALIGN_CENTRE)
            self.lblSecondaryInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
            resizeTextControlForUNIX(self.lblSecondaryInst, 0, self.GetSize()[0])
        self.lblSecondaryInst.SetForegroundColour("#FFFFFF")
        
        if (platform.system() != "Linux"):
            self.lblSecondaryModel = wx.StaticText(self, -1, "None Selected", (0, 638), (160, 25), wx.ALIGN_CENTRE)
            self.lblSecondaryModel.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
        else:
            self.lblSecondaryModel = wx.StaticText(self, -1, "None Selected", (0, 638), style=wx.ALIGN_CENTRE)
            self.lblSecondaryModel.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
            resizeTextControlForUNIX(self.lblSecondaryModel, 0, 160)
        self.lblSecondaryModel.SetForegroundColour("#FFFFFF")
        self.modelMenu = wx.ComboBox(self, pos=(160, 635), size=(160, 25), choices=[], style=wx.CB_READONLY)
        self.modelMenu.Bind(wx.EVT_COMBOBOX, self.modelMenuSelect)
        self.modelMenu.SetToolTipString("Choose a state whose secondary behavior will be edited")
        
        if (platform.system() == "Darwin"):
            self.btnAdd = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnAdd.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(7, 665), size=(57, 25))
        else:
            self.btnAdd = wx.Button(self, id=-1, label="Add", pos=(7, 665), size=(57, 25))
            self.btnAdd.SetForegroundColour("#000000")
            self.btnAdd.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAdd.Bind(wx.EVT_BUTTON, self.add)
        self.btnAdd.SetToolTipString("Add selected residues to the resfile")
        if (platform.system() == "Darwin"):
            self.btnRemove = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnRemove.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(69, 665), size=(57, 25))
        else:
            self.btnRemove = wx.Button(self, id=-1, label="Remove", pos=(69, 665), size=(57, 25))
            self.btnRemove.SetForegroundColour("#000000")
            self.btnRemove.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnRemove.Bind(wx.EVT_BUTTON, self.remove)
        self.btnRemove.SetToolTipString("Remove selected residues from the resfile")
        if (platform.system() == "Darwin"):
            self.btnRestrict = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnRestrict.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(131, 665), size=(57, 25))
        else:
            self.btnRestrict = wx.Button(self, id=-1, label="Restrict", pos=(131, 665), size=(57, 25))
            self.btnRestrict.SetForegroundColour("#000000")
            self.btnRestrict.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnRestrict.Bind(wx.EVT_BUTTON, self.restrict)
        self.btnRestrict.SetToolTipString("Restrict the resfile contents to the selected residues")
        if (platform.system() == "Darwin"):
            self.btnAll = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnAll.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(193, 665), size=(57, 25))
        else:
            self.btnAll = wx.Button(self, id=-1, label="All", pos=(193, 665), size=(57, 25))
            self.btnAll.SetForegroundColour("#000000")
            self.btnAll.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAll.Bind(wx.EVT_BUTTON, self.addAll)
        self.btnAll.SetToolTipString("Add all residues from the selected model to the resfile")
        if (platform.system() == "Darwin"):
            self.btnClear = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnClear.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(255, 665), size=(57, 25))
        else:
            self.btnClear = wx.Button(self, id=-1, label="Clear", pos=(255, 665), size=(57, 25))
            self.btnClear.SetForegroundColour("#000000")
            self.btnClear.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnClear.Bind(wx.EVT_BUTTON, self.clear)
        self.btnClear.SetToolTipString("Clear the contents of the resfile")
        self.selectedData = []
        
        self.grdResfile = wx.grid.Grid(self)
        self.grdResfile.CreateGrid(0, 1)
        self.grdResfile.SetSize((320, 200))
        self.grdResfile.SetPosition((0, 695))
        self.grdResfile.SetLabelFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.grdResfile.DisableDragColSize()
        self.grdResfile.DisableDragRowSize()
        self.grdResfile.SetColLabelValue(0, "Behavior")
        self.grdResfile.SetRowLabelSize(160)
        self.grdResfile.SetColSize(0, 160)
        self.grdResfile.Bind(wx.grid.EVT_GRID_CELL_LEFT_CLICK, self.secondaryClick)
        self.resfile = []
        # Switch for telling the resfile to look for the restype from the sequence window pre-design
        # or from the outputted designed before accepting it
        self.useDesignedSeq = False 
        
        ypos = self.grdResfile.GetPosition()[1] + self.grdResfile.GetSize()[1] + 10
        if (platform.system() == "Windows"):
            self.lblStates = wx.StaticText(self, -1, "Define States", (0, ypos), (320, 25), wx.ALIGN_CENTRE)
            self.lblStates.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblStates = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/lblStates.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(25, ypos), size=(270, 25))
        else:
            self.lblStates = wx.StaticText(self, -1, "Define States", (0, ypos), style=wx.ALIGN_CENTRE)
            self.lblStates.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
            resizeTextControlForUNIX(self.lblStates, 0, self.GetSize()[0])
        self.lblStates.SetForegroundColour("#FFFFFF")
        if (platform.system() == "Windows"):
            self.lblStatesInst = wx.StaticText(self, -1, "Group PDBs into distinct states.  Click on the\ngraph to change the state to which each PDB\nbelongs.  Multiple conformations of the same\nstructure should be in a single state.", (0, ypos+30), (320, 25), wx.ALIGN_CENTRE)
            self.lblStatesInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
        elif (platform.system() == "Darwin"):
            self.lblStatesInst = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/lblInstStates.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, ypos+30), size=(320, 70))
        else:
            self.lblStatesInst = wx.StaticText(self, -1, "Group PDBs into distinct states.  Click on the\ngraph to change the state to which each PDB\nbelongs.  Multiple conformations of the same\nstructure should be in a single state.", (5, ypos+30), style=wx.ALIGN_CENTRE)
            self.lblStatesInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
            resizeTextControlForUNIX(self.lblStatesInst, 0, self.GetSize()[0])
        self.lblStatesInst.SetForegroundColour("#FFFFFF")
        
        self.grdStates = wx.grid.Grid(self)
        self.grdStates.CreateGrid(0, 1)
        self.grdStates.SetSize((320, 150))
        self.grdStates.SetPosition((0, ypos+100))
        self.grdStates.SetLabelFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.grdStates.DisableDragColSize()
        self.grdStates.DisableDragRowSize()
        self.grdStates.SetColLabelValue(0, "State")
        self.grdStates.SetRowLabelSize(160)
        self.grdStates.SetColSize(0, 160)
        self.grdStates.Bind(wx.grid.EVT_GRID_CELL_LEFT_CLICK, self.stateClick)
        self.stateAssignments = []
        
        ypos = self.grdStates.GetPosition()[1] + self.grdStates.GetSize()[1] + 10
        if (platform.system() == "Windows"):
            self.lblVectorInst = wx.StaticText(self, -1, "If there are multiple conformations in a state,\nwill the score of the state be the minimum or\nmaximum score of the conformations?", (0, ypos+10), (320, 25), wx.ALIGN_CENTRE)
            self.lblVectorInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
        if (platform.system() == "Darwin"):
            self.lblVectorInst = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/lblInstVector.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, ypos+10), size=(320, 50))
        else:
            self.lblVectorInst = wx.StaticText(self, -1, "If there are multiple conformations in a state,\nwill the score of the state be the minimum or\nmaximum score of the conformations?", (5, ypos+10), style=wx.ALIGN_CENTRE)
            self.lblVectorInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
            resizeTextControlForUNIX(self.lblVectorInst, 0, self.GetSize()[0])
        self.lblVectorInst.SetForegroundColour("#FFFFFF")
        
        self.vectorMenu = wx.ComboBox(self, pos=(0, ypos+60), size=(160, 25), choices=[], style=wx.CB_READONLY)
        self.vectorMenu.Bind(wx.EVT_COMBOBOX, self.vectorMenuSelect)
        self.vectorMenu.SetToolTipString("Specify the state for which to alter the vector function")
        if (platform.system() == "Darwin"):
            self.btnMin = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnMin_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(165, ypos+60), size=(70, 25))
        else:
            self.btnMin = wx.Button(self, id=-1, label="MIN", pos=(165, ypos+60), size=(70, 25))
            self.btnMin.SetForegroundColour("#FF0000")
            self.btnMin.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnMin.Bind(wx.EVT_BUTTON, self.vectorMin)
        self.btnMin.SetToolTipString("Apply a minimization operation to this state")
        if (platform.system() == "Darwin"):
            self.btnMax = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnMax.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(245, ypos+60), size=(70, 25))
        else:
            self.btnMax = wx.Button(self, id=-1, label="MAX", pos=(245, ypos+60), size=(70, 25))
            self.btnMax.SetForegroundColour("#000000")
            self.btnMax.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnMax.Bind(wx.EVT_BUTTON, self.vectorMax)
        self.btnMax.SetToolTipString("Apply a maximization operation to this state")
        self.vectorFunctions = []
        
        if (platform.system() == "Windows"):
            self.lblFitness = wx.StaticText(self, -1, "Fitness Function", (0, ypos+90), (320, 25), wx.ALIGN_CENTRE)
            self.lblFitness.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblFitness = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/lblFitness.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(25, ypos+90), size=(270, 25))
        else:
            self.lblFitness = wx.StaticText(self, -1, "Fitness Function", (0, ypos+90), style=wx.ALIGN_CENTRE)
            self.lblFitness.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
            resizeTextControlForUNIX(self.lblFitness, 0, self.GetSize()[0])
        self.lblFitness.SetForegroundColour("#FFFFFF")
        if (platform.system() == "Windows"):
            self.lblFitnessInst = wx.StaticText(self, -1, "Define the mathematical relationship between\nthe scores of each of the states when\ndetermining fitness.", (0, ypos+120), (320, 25), wx.ALIGN_CENTRE)
            self.lblFitnessInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
        elif (platform.system() == "Darwin"):
            self.lblFitnessInst = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/lblInstFitness.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, ypos+120), size=(320, 50))
        else:
            self.lblFitnessInst = wx.StaticText(self, -1, "Define the mathematical relationship between\nthe scores of each of the states when\ndetermining fitness.", (5, ypos+120), style=wx.ALIGN_CENTRE)
            self.lblFitnessInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
            resizeTextControlForUNIX(self.lblFitnessInst, 0, self.GetSize()[0])
        self.lblFitnessInst.SetForegroundColour("#FFFFFF")
        
        self.txtFitness = wx.stc.StyledTextCtrl(self, -1, pos=(0, ypos+170), size=(320, 50))
        self.txtFitness.SetText("?")
        self.txtFitness.SetReadOnly(True)
        self.txtFitness.SetToolTipString("Fitness function for multistate design")
        self.txtFitness.Bind(wx.EVT_LEFT_UP, self.fitnessRelease)
        if (platform.system() == "Windows"):
            self.lblStateLiteral = wx.StaticText(self, -1, "State/Literal Value", (0, ypos+233), (160, 25), wx.ALIGN_CENTRE)
            self.lblStateLiteral.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
        elif (platform.system() == "Darwin"):
            self.lblStateLiteral = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/lblStateLiteral.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, ypos+233), size=(160, 25))
        else:
            self.lblStateLiteral = wx.StaticText(self, -1, "State/Literal Value", (5, ypos+233), style=wx.ALIGN_CENTRE)
            self.lblStateLiteral.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
            resizeTextControlForUNIX(self.lblStateLiteral, 0, 160)
        self.lblStateLiteral.SetForegroundColour("#FFFFFF")
        self.txtStateLiteral = wx.TextCtrl(self, -1, pos=(165, ypos+230), size=(150, 25))
        self.txtStateLiteral.SetValue("")
        self.txtStateLiteral.Disable()
        self.txtStateLiteral.SetToolTipString("State variable or literal at the selected position in the Fitness Function text box")
        self.txtStateLiteral.Bind(wx.EVT_TEXT, self.stateliteralUpdate)
        self.operationMenu = wx.ComboBox(self, pos=(0, ypos+260), size=(160, 25), choices=[], style=wx.CB_READONLY)
        self.operationMenu.Bind(wx.EVT_COMBOBOX, self.operationMenuSelect)
        self.operationMenu.SetToolTipString("Operations and functions available for the fitness function")
        self.operationMenu.AppendItems(["Addition (? + ?)", "Subtraction (? - ?)", "Multiplication (? * ?)",
                                 "Division (? / ?)", "Parentheses (())", "Absolute Value (abs(?))", "Exponential Base e (exp(?))",
                                 "Natural Logarithm (ln(?))", "Power Function (pow(?, ?))", "Square Root (sqrt(?))",
                                 "Boolean Equality (eq(?, ?))", "Boolean Greater Than (gt(?, ?))", "Boolean Greater Than Equal To (gte(?, ?))",
                                 "Boolean Less Than (lt(?, ?))", "Boolean Less Than Equal To (lte(?, ?))", "Boolean And (and(?, ?))",
                                 "Boolean Or (or(?, ?))", "Boolean Not (not(?))", "If-Then-Else (ite(?, ?, ?))"])
        self.numCommas = {"abs": 0, "exp": 0, "ln": 0, "pow": 1, "sqrt": 1, "eq": 1, "gt": 1, "gte": 1, "lt": 1,
                   "lte": 1, "and": 1, "or": 1, "not": 0, "ite": 2}
        if (platform.system() == "Darwin"):
            self.btnAddOp = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnAddOp.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(165, ypos+260), size=(70, 25))
        else:
            self.btnAddOp = wx.Button(self, id=-1, label="Add", pos=(165, ypos+260), size=(70, 25))
            self.btnAddOp.SetForegroundColour("#000000")
            self.btnAddOp.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAddOp.Bind(wx.EVT_BUTTON, self.addOperation)
        self.btnAddOp.SetToolTipString("Apply a minimization operation to this state")
        if (platform.system() == "Darwin"):
            self.btnRemoveOp = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnRemoveOp.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(245, ypos+260), size=(70, 25))
        else:
            self.btnRemoveOp = wx.Button(self, id=-1, label="Remove", pos=(245, ypos+260), size=(70, 25))
            self.btnRemoveOp.SetForegroundColour("#000000")
            self.btnRemoveOp.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnRemoveOp.Bind(wx.EVT_BUTTON, self.removeOperation)
        self.btnRemoveOp.SetToolTipString("Apply a maximization operation to this state")
        
        self.fitness = self.txtFitness.GetText()
        
        if (platform.system() == "Darwin"):
            self.btnLoadMSD = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnLoadMSD.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(20, ypos+290), size=(120, 25))
        else:
            self.btnLoadMSD = wx.Button(self, id=-1, label="Load MSD", pos=(20, ypos+290), size=(120, 25))
            self.btnLoadMSD.SetForegroundColour("#000000")
            self.btnLoadMSD.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.btnLoadMSD.Bind(wx.EVT_BUTTON, self.loadMSD)
        self.btnLoadMSD.SetToolTipString("Load data stored in an MSD file")
        if (platform.system() == "Darwin"):
            self.btnSaveMSD = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnSaveMSD.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(175, ypos+290), size=(120, 25))
        else:
            self.btnSaveMSD = wx.Button(self, id=-1, label="Save MSD", pos=(175, ypos+290), size=(120, 25))
            self.btnSaveMSD.SetForegroundColour("#000000")
            self.btnSaveMSD.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.btnSaveMSD.Bind(wx.EVT_BUTTON, self.saveMSD)
        self.btnSaveMSD.SetToolTipString("Save the MSD setup to a file for later use")
        
        if (platform.system() == "Darwin"):
            self.btnMSD = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnMSD.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(100, ypos+320), size=(120, 25))
        else:
            self.btnMSD = wx.Button(self, id=-1, label="Run MSD!", pos=(100, ypos+320), size=(120, 25))
            self.btnMSD.SetForegroundColour("#000000")
            self.btnMSD.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.btnMSD.Bind(wx.EVT_BUTTON, self.MSDClick)
        self.btnMSD.SetToolTipString("Send the multi-state design job to the server")
        self.WTType = "NATRO"
        
        self.scrollh = self.btnMSD.GetPosition()[1] + self.btnMSD.GetSize()[1] + 5
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
            browser.open(self.parent.parent.scriptdir + "/help/msd.html")
        else:
            webbrowser.open(self.parent.parent.scriptdir + "/help/msd.html")
    
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

    def updateResfile(self):
        scrollpos = self.grdResfile.GetScrollPos(wx.VERTICAL)
        if (self.grdResfile.NumberRows > 0):
            self.grdResfile.DeleteRows(0, self.grdResfile.NumberRows)
        row = 0
        for [indx, r, seqpos, poseindx, chainID, chainoffset, reslist] in self.resfile[self.resfilemodelindx]:
            self.grdResfile.AppendRows(1)
            ID = self.seqWin.IDs[r]
            resn = self.seqWin.SeqViewer.GetCellValue(r, indx)
            if (self.useDesignedSeq):
                if (chainID == "_" or chainID == ""):
                    chain = " "
                else:
                    chain = chainID
                mut = AA3to1(self.designedView[0][chain][self.seqWin.indxToSeqPos[r][indx]].resname)
                label = chainID + "|" + str(resn) + str(seqpos) + str(mut)
            else:
                label = chainID + "|" + str(resn) + str(seqpos)
            self.grdResfile.SetRowLabelValue(row, label)
            self.grdResfile.SetCellValue(row, 0, reslist)
            # This needs to happen before the readOnly attr is set otherwise it doesn't apply the center on the
            # last cell for some bizarre reason
            self.grdResfile.SetCellAlignment(row, 0, wx.ALIGN_CENTRE, wx.ALIGN_CENTRE)
            # Very important note: you actually need to create a new GridCellAttr for each new row
            # You cannot just declare it outside of the loop and use the same one for each row otherwise you
            # get some pretty nasty crashes when you delete rows
            readOnly = wx.grid.GridCellAttr()
            readOnly.SetReadOnly(True)
            self.grdResfile.SetRowAttr(row, readOnly)
            # Now update the drop down menu so the user can tweak the settings of individual residues in the minmap
            row = row + 1
        self.grdResfile.Scroll(0, scrollpos)
        
    def redrawEntity(self):
        # Function for redrawing the entity grid
        # Figure out the number of columns = 1 + #models (the extra 1 is for the residue option)
        scrollpos = self.grdEntity.GetScrollPos(wx.VERTICAL)
        ncols = 1 + len(self.states)
        if (self.grdEntity.NumberCols < ncols):
            self.grdEntity.AppendCols(ncols-self.grdEntity.NumberCols)
        elif (self.grdEntity.NumberCols > ncols):
            self.grdEntity.DeleteCols(numCols=self.grdEntity.NumberCols-ncols)
        # Relabel the columns
        self.grdEntity.SetColLabelValue(0, "Residue Options")
        for i in range(1, len(self.states)+1):
            self.grdEntity.SetColLabelValue(i, self.states[i-1][0])
        # See if the number of rows needs updating
        if (len(self.states) > 0):
            nrows = len(self.states[0]) - 1 # The first element is the model name
        else:
            nrows = 0
        if (self.grdEntity.NumberRows < nrows):
            self.grdEntity.AppendRows(nrows-self.grdEntity.NumberRows)
        elif (self.grdEntity.NumberRows > nrows):
            self.grdEntity.DeleteRows(numRows=self.grdEntity.NumberRows-nrows)
        # Update the labels and fill in residue data
        for i in range(0, nrows):
            self.grdEntity.SetRowLabelValue(i, str(i+1))
            self.grdEntity.SetCellValue(i, 0, self.stateoptions[i])
            # Fill in the correspondance data
            for c in range(1, ncols):
                self.grdEntity.SetCellValue(i, c, self.states[c-1][i+1])
            # Set as read only
            readOnly = wx.grid.GridCellAttr()
            readOnly.SetReadOnly(True)
            readOnly.SetAlignment(wx.ALIGN_CENTRE, wx.ALIGN_CENTRE)
            self.grdEntity.SetRowAttr(i, readOnly)
        # Fix column widths
        # Resize columns if necessary
        fitGridColumn(self.grdEntity, 0, 190)
        self.grdEntity.SetColSize(0, 190)
        for i in range(1, ncols):
            fitGridColumn(self.grdEntity, i, 90)
        self.grdEntity.Scroll(0, scrollpos)
        self.grdEntity.Refresh()
        
    def scrolled(self, event):
        self.winscrollpos = self.GetScrollPos(wx.VERTICAL)
        event.Skip()
        
    def activate(self):
        # This function assumes that everything currently loaded in the Sequence Window is eligible for MSD
        # If the user doesn't intend to use all of the loaded structures in MSD, then they just never specify
        # an entry for it in the entity resfile and it will not be included
        # First let's get an updated list of all of the models loaded in the Sequence Window
        models = []
        for ID in self.seqWin.IDs:
            #model = ID[0:len(ID)-2]
            if (not(ID in models)):
                models.append(ID)
        currmodels = []
        for i in range(1, self.grdEntity.NumberCols):
            currmodels.append(self.grdEntity.GetColLabelValue(i))
        # Now this is a little tricky, because if something new is added we want to retain all the information
        # for things that were already loaded, but the might have shifted positions relative to other chains
        # so we have to move this information around
        if (models != currmodels):
            # If the currently-viewed model was removed, take it out of the combobox
            if (not(self.selectedModel in models)):
                if (self.grdRow >= 0 and self.grdCol >= 0):
                    self.grdEntity.SetCellBackgroundColour(self.grdRow, self.grdCol, "white")
                self.grdRow = -1
                self.grdCol = -1
                self.posMenu.SetSelection(0) # To get the combobox to display "Nothing" on Linux
                self.posMenu.Clear()
                self.selectedModel = ""
                self.lblModel.SetLabel("None Selected")
                self.lblModel.SetForegroundColour("#FFFFFF")
                if (platform.system() == "Linux"):
                    resizeTextControlForUNIX(self.lblModel, 0, 160)
            newstates = []
            newresfile = []
            newassignments = []
            for model in models:
                newstates.append([model])
                newresfile.append([])
                newassignments.append("")
                for r in range(0, self.grdEntity.NumberRows):
                    newstates[len(newstates)-1].append("") # Blank entries for the correspondances
            for stateindx in range(0, len(self.states)):
                try:
                    modelindx = models.index(self.states[stateindx][0]) # Element 0 is the model name
                    newstates[modelindx] = self.states[stateindx][:]
                    newresfile[modelindx] = self.resfile[stateindx][:]
                    newassignments[modelindx] = self.stateAssignments[stateindx]
                except:
                    pass
            self.states = newstates
            self.resfile = newresfile
            self.stateAssignments = newassignments
            self.distinctmodels = []
            for ID in models:
                model = ID[0:len(ID)-2]
                if (not(model in self.distinctmodels)):
                    self.distinctmodels.append(model)
            # Now we have to give valid state assignments for the new models
            for model in self.distinctmodels:
                for i in range(0, len(self.stateAssignments)):
                    if (self.states[i][0][0:len(self.states[i][0])-2] == model):
                        if (len(self.stateAssignments[i]) == 0):
                            for ID in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
                                if (not(ID in self.stateAssignments)):
                                    break
                            self.stateAssignments[i] = ID
                        elif ("ABCDEFGHIJKLMNOPQRSTUVWXYZ".index(self.stateAssignments[i]) >= len(self.distinctmodels)):
                            self.stateAssignments[i] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[len(self.distinctmodels)-1]
                        break
            # Update the secondary resfile model menu
            self.modelMenu.Clear()
            self.modelMenu.AppendItems(self.distinctmodels)
            if (len(models) > 0):
                self.modelMenu.SetSelection(0) # To get the combobox to display the first element on Linux
                self.lblSecondaryModel.SetLabel(self.modelMenu.GetValue())
                self.lblSecondaryModel.SetForegroundColour("#FFFFFF")
                self.resfileModel = self.modelMenu.GetValue()
                self.resfilemodelindx = 0
            else:
                self.lblSecondaryModel.SetLabel("None Selected")
                self.resfileModel = ""
            self.grdEntity.ClearGrid()
            self.redrawEntity() # Redraw the grid
            # Redraw the states grid
            if (self.grdStates.NumberRows > len(self.distinctmodels)):
                self.grdStates.DeleteRows(self.grdStates.NumberRows - len(self.distinctmodels))
            elif (self.grdStates.NumberRows < len(self.distinctmodels)):
                self.grdStates.AppendRows(len(self.distinctmodels) - self.grdStates.NumberRows)
            j = 0
            for i in range(0, len(self.stateAssignments)):
                if (len(self.stateAssignments[i]) == 0):
                    continue
                self.grdStates.SetRowLabelValue(j, self.distinctmodels[j])
                self.grdStates.SetCellValue(j, 0, self.stateAssignments[i])
                readOnly = wx.grid.GridCellAttr()
                readOnly.SetReadOnly(True)
                readOnly.SetAlignment(wx.ALIGN_CENTRE, wx.ALIGN_CENTRE)
                self.grdStates.SetRowAttr(j, readOnly)
                j = j + 1
            # Regenerate the vector function menu
            if (len(self.distinctmodels) < len(self.vectorFunctions)):
                # Pop things off
                for i in range(0, len(self.vectorFunctions) - len(self.distinctmodels)):
                    self.vectorFunctions.pop()
            elif (len(self.distinctmodels) > len(self.vectorFunctions)):
                # Default the new states to MIN
                for i in range(0, len(self.distinctmodels) - len(self.vectorFunctions)):
                    self.vectorFunctions.append("MIN")
            self.vectorMenu.Clear()
            for i in range(0, len(self.distinctmodels)):
                self.vectorMenu.Append("ABCDEFGHIJKLMNOPQRSTUVWXYZ"[i])
            if (len(models) > 0):
                self.vectorMenu.SetSelection(0)
                if (self.vectorFunctions[0] == "MIN"):
                    self.btnMin.SetForegroundColour("#FF0000")
                    self.btnMax.SetForegroundColour("#000000")
                else:
                    self.btnMin.SetForegroundColour("#000000")
                    self.btnMax.SetForegroundColour("#FF0000")
            else:
                self.btnMin.SetForegroundColour("#FF0000")
                self.btnMax.SetForegroundColour("#000000")
            # Just erase the fitness function every time a new model gets added/removed because it's too
            # complicated to figure out how to remove the state while keeping the function valid
            # The user shouldn't be adding/removing things at that step anyway
            #self.txtFitness.SetValue("?")
            #self.fitness = "?"
        # Grab the current selection of residues for processing with the buttons
        topLefts = self.seqWin.SeqViewer.GetSelectionBlockTopLeft()
        bottomRights = self.seqWin.SeqViewer.GetSelectionBlockBottomRight()
        self.selectedData = []
        for i in range(0, len(topLefts)):
            for r in range(topLefts[i][0], bottomRights[i][0]+1):
                for c in range(topLefts[i][1], bottomRights[i][1]+1):
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
    
    def addEntityPosition(self, event):
        # Do nothing if nothing is in the palette (since we can't default to WT since there is no WT in MSD)
        if (len(self.palette) == 0):
            return
        # Just add another row of blank fields for each of the states and then update the grid
        for state in self.states:
            state.append("")
        self.stateoptions.append(self.palette) # Default to ALLAAxc
        self.entityMenu.Clear()
        for i in range(1, len(self.stateoptions)+1):
            self.entityMenu.AppendItems(str(i))
        self.redrawEntity()
    
    def removeEntityPosition(self, event):
        # Pop out the values for the indicated row
        try:
            popindx = int(self.entityMenu.GetValue())
            if (popindx > self.grdEntity.NumberRows):
                raise Exception()
        except:
            return
        if (self.grdRow >= 0 and self.grdCol >= 0):
            self.grdEntity.SetCellBackgroundColour(self.grdRow, self.grdCol, "white")
        self.grdRow = -1
        self.grdCol = -1
        self.posMenu.SetSelection(0)
        self.posMenu.Clear()
        self.lblModel.SetLabel("None Selected")
        if (platform.system() != "Windows"):
            resizeTextControlForUNIX(self.lblModel, 0, 160)
        # There might have been a selection on the row that was removed so unselect it
        for state in self.states:
            state.pop(popindx)
        self.stateoptions.pop(popindx-1) # There's no model name at the beginning of this one
        self.entityMenu.Clear()
        for i in range(1, len(self.stateoptions)+1):
            self.entityMenu.AppendItems(str(i))
        self.redrawEntity()
    
    def entityClick(self, event):
        # Set the background of the last selection back to white, and change this selection to light blue
        if (self.grdRow >= 0 and self.grdCol >= 0):
            self.grdEntity.SetCellBackgroundColour(self.grdRow, self.grdCol, "white")
            if (event.GetCol() <= 0):
                if (event.GetCol() < 0):
                    self.grdRow = -1
                    self.grdCol = -1
                self.posMenu.SetSelection(0)
                self.posMenu.Clear()
                self.lblModel.SetLabel("None Selected")
                if (platform.system() != "Windows"):
                    resizeTextControlForUNIX(self.lblModel, 0, 160)
        # Fill up the menu with positions from the model in this column
        self.grdRow = event.GetRow()
        self.grdCol = event.GetCol()
        self.entityMenu.SetSelection(self.grdRow)
        if (self.grdRow < 0 or self.grdCol < 1):
            if (self.grdCol == 0):
                self.grdEntity.SetCellBackgroundColour(self.grdRow, self.grdCol, "light blue")
            self.grdEntity.Refresh()
            event.Skip()
            return
        self.grdEntity.SetCellBackgroundColour(self.grdRow, self.grdCol, "light blue")
        # Get the structure
        self.selectedModel = self.grdEntity.GetColLabelValue(self.grdCol)
        model = self.selectedModel[0:len(self.selectedModel)-2]
        chain = self.selectedModel[len(self.selectedModel)-1]
        if (chain == "_"):
            chain = " "
        pose = self.seqWin.poses[self.seqWin.getPoseIndexForModel(model)]
        # Fill up the menu
        positions = ["Nothing"]
        ires = 1
        for residue in pose[0][chain]:
            # Skip NCAAs
            if ("ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR ".find(residue.resname) >= 0):
                chain2 = chain
                if (chain2 == " " or chain2 == ""):
                    chain2 = "_"
                label = chain2 + ":" + AA3to1(residue.resname) + str(residue.id[1])
                positions.append(label)
            ires = ires + 1
        self.posMenu.Clear()
        self.posMenu.AppendItems(positions)
        # Automatically select whatever is already in the grid
        try:
            indx = positions.index(self.grdEntity.GetCellValue(self.grdRow, self.grdCol))
        except:
            # Default to nothing
            indx = 0
        self.posMenu.SetSelection(indx)
        # Update the label
        self.lblModel.SetLabel(self.selectedModel)
        if (platform.system() != "Windows"):
            resizeTextControlForUNIX(self.lblModel, 0, 160)
        # View it in PyMOL
        self.focusView(self.posMenu.GetValue(), model)
        self.grdEntity.Refresh()
        event.Skip()
        
    def entityRClick(self, event):
        # This is a function for quickly filling up a column with sequential residues in a selection
        indx = 0
        r = event.GetRow()
        c = event.GetCol()
        while (self.selectedData):
            if (r >= self.grdEntity.NumberRows):
                break
            while (True):
                if (indx >= len(self.selectedData) or self.seqWin.IDs[self.selectedData[indx][1]] == self.grdEntity.GetColLabelValue(c)):
                    break
                indx = indx + 1
            if (indx >= len(self.selectedData)):
                break
            chain = self.selectedData[indx][4]
            if (not(self.seqWin.sequences[self.selectedData[indx][1]][self.selectedData[indx][0]] in "ACDEFGHIKLMNPQRSTVWY")):
                indx = indx + 1
                continue
            self.grdEntity.SetCellValue(r, c, chain + ":" + self.seqWin.sequences[self.selectedData[indx][1]][self.selectedData[indx][0]] + self.selectedData[indx][2])
            self.states[c-1][r+1] = self.grdEntity.GetCellValue(r, c)
            r = r + 1
            indx = indx + 1
        event.Skip()
        
    def focusView(self, focusstr, model):
        # Find the neighborhood view on this residue
        if (focusstr == "Nothing"):
            # Show the whole structure
            self.cmd.select("msdsele", "model " + model)
            self.cmd.hide("everything", "all")
            self.cmd.show("cartoon", "msdsele")
            self.cmd.zoom("msdsele")
            self.cmd.delete("msdsele")
            return
        chain = focusstr[0]
        seqpos = focusstr.strip()[3:] # indx 2 is the AA
        logInfo("Position " + focusstr + " was selected")
        self.cmd.hide("all")
        if (chain == " " or chain == "_"):
            self.cmd.select("msdsele", "resi " + seqpos + " and model " + model)
        else:
            self.cmd.select("msdsele", "resi " + seqpos + " and model " + model + " and chain " + chain)
        self.cmd.select("msdsele", "model " + model + " within 12 of msdsele")
        self.cmd.show("cartoon", "msdsele")
        self.cmd.hide("ribbon", "msdsele")
        self.cmd.show("sticks", "msdsele")
        self.cmd.set_bond("stick_radius", 0.1, "msdsele")
        self.cmd.zoom("msdsele")
        if (chain == " " or chain == "_"):
            self.cmd.select("msdsele", "resi " + seqpos + " and model " + model)
        else:
            self.cmd.select("msdsele", "resi " + seqpos + " and model " + model + " and chain " + chain)
        self.cmd.show("sticks", "msdsele")
        self.cmd.set_bond("stick_radius", 0.25, "msdsele")
        # Highlight this residue in PyMOL
        self.cmd.select("sele", "msdsele")
        self.cmd.enable("sele")
        self.cmd.delete("msdsele")
        self.seqWin.selectUpdate(False)
    
    def entityPosMenuSelect(self, event):
        # Change the current residue correspondance entry to the current selection
        if (self.posMenu.GetValue().strip() == "Nothing"):
            self.grdEntity.SetCellValue(self.grdRow, self.grdCol, "")
            self.states[self.grdCol-1][self.grdRow+1] = "" # Remember that in the grid we skip column 1 (residue options) and in self.states we skip row 1 (model name)
        else:
            self.grdEntity.SetCellValue(self.grdRow, self.grdCol, self.posMenu.GetValue().strip())
            self.states[self.grdCol-1][self.grdRow+1] = self.posMenu.GetValue().strip()
        self.focusView(self.posMenu.GetValue(), self.selectedModel[0:len(self.selectedModel)-2])
        event.Skip()
        
    def modelMenuSelect(self, event):
        # Update the model variable and change the label
        self.resfileModel = self.modelMenu.GetValue()
        for i in range(0, len(self.states)):
            if (self.states[i][0][0:len(self.states[i][0])-2] == self.resfileModel):
                self.resfilemodelindx = i
                break
        self.lblSecondaryModel.SetLabel(self.resfileModel)
        self.lblSecondaryModel.SetForegroundColour("#FFFFFF")
        if (platform.system() == "Linux"):
            resizeTextControlForUNIX(self.lblSecondaryModel, 0, 160)
        # Update the grid with the data for this particular resfile
        self.updateResfile()
        
    def secondaryClick(self, event):
        # Highlight the selected row
        for r in range(0, self.grdResfile.NumberRows):
            if (r == event.GetRow()):
                self.grdResfile.SetCellBackgroundColour(r, 0, "light blue")
            else:
                self.grdResfile.SetCellBackgroundColour(r, 0, "white")
        # Focus on it in PyMOL
        self.focusView(self.grdResfile.GetRowLabelValue(event.GetRow()), self.resfileModel)
        self.grdResfile.Refresh()
        event.Skip()
        
    def stateClick(self, event):
        indx = "ABCDEFGHIJKLMNOPQRSTUVWXYZ".index(self.grdStates.GetCellValue(event.GetRow(), 0))
        if (indx >= len(self.distinctmodels) - 1):
            indx = -1
        newID = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[indx+1]
        self.grdStates.SetCellValue(event.GetRow(), 0, newID)
        self.stateAssignments[event.GetRow()] = newID
        event.Skip()
        
    def vectorMenuSelect(self, event):
        if (self.vectorFunctions[self.vectorMenu.GetSelection()] == "MIN"):
            self.btnMin.SetForegroundColour("#FF0000")
            self.btnMax.SetForegroundColour("#000000")
        else:
            self.btnMin.SetForegroundColour("#000000")
            self.btnMax.SetForegroundColour("#FF0000")
        
    def vectorMin(self, event):
        if (platform.system() == "Darwin"):
            self.btnMin.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnMin_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            self.btnMax.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnMax.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
        else:
            self.btnMin.SetForegroundColour("#FF0000")
            self.btnMax.SetForegroundColour("#000000")
        self.vectorFunctions[self.vectorMenu.GetSelection()] = "MIN"
    
    def vectorMax(self, event):
        if (platform.system() == "Darwin"):
            self.btnMin.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnMin.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            self.btnMax.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnMax_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
        else:
            self.btnMin.SetForegroundColour("#000000")
            self.btnMax.SetForegroundColour("#FF0000")
        self.vectorFunctions[self.vectorMenu.GetSelection()] = "MAX"
    
    def operationMenuSelect(self, event):
        pass
    
    def addOperation(self, event):
        operation = self.operationMenu.GetValue()
        try:
            self.fitnessStart < self.fitnessEnd
        except:
            # Not defined yet, do nothing
            return
        # If there is an operator at the beginning or end of a selection, then simply change that
        # operator to the new one specified
        # If not, then add the operator plus a new blank entry
        if (operation.startswith("Addition")):
            if (self.fitness[self.fitnessStart] in "+-*/"):
                self.fitness = self.fitness[0:self.fitnessStart] + "+" + self.fitness[self.fitnessStart+1:]
            elif (self.fitness[self.fitnessEnd-1] in "+-*/"):
                self.fitness = self.fitness[0:self.fitnessEnd-1] + "+" + self.fitness[self.fitnessEnd:]
            else:
                self.fitness = self.fitness[0:self.fitnessEnd] + "+?" + self.fitness[self.fitnessEnd:]
        elif (operation.startswith("Subtraction")):
            if (self.fitness[self.fitnessStart] in "+-*/"):
                self.fitness = self.fitness[0:self.fitnessStart] + "-" + self.fitness[self.fitnessStart+1:]
            elif (self.fitness[self.fitnessEnd-1] in "+-*/"):
                self.fitness = self.fitness[0:self.fitnessEnd-1] + "-" + self.fitness[self.fitnessEnd:]
            else:
                self.fitness = self.fitness[0:self.fitnessEnd] + "-?" + self.fitness[self.fitnessEnd:]
        elif (operation.startswith("Multiplication")):
            # Here there is a special consideration, since * and / take precedence over + and -
            # So first we have to see if the selection is enclosed in ()
            # If yes, then nothing special needs to happen
            # If no, then if + or - is in the selection the selection needs to be enclosed in () first
            if (self.fitness[self.fitnessStart] in "+-*/"):
                self.fitness = self.fitness[0:self.fitnessStart] + "*" + self.fitness[self.fitnessStart+1:]
            elif (self.fitness[self.fitnessEnd-1] in "+-*/"):
                self.fitness = self.fitness[0:self.fitnessEnd-1] + "*" + self.fitness[self.fitnessEnd:]
            else:
                if (self.fitness[self.fitnessStart] != "(" or self.fitness[self.fitnessEnd-1] != ")"):
                    needParen = False
                    for i in range(self.fitnessStart, self.fitnessEnd):
                        if (self.fitness[i] in "+-"):
                            needParen = True
                            break
                    # Add the () in
                    if (needParen):
                        self.fitness = self.fitness[0:self.fitnessStart] + "(" + self.fitness[self.fitnessStart:self.fitnessEnd] + ")" + self.fitness[self.fitnessEnd:]
                        self.fitnessEnd = self.fitnessEnd + 2
                self.fitness = self.fitness[0:self.fitnessEnd] + "*?" + self.fitness[self.fitnessEnd:]
        elif (operation.startswith("Division")):
            # Here there is a special consideration, since * and / take precedence over + and -
            # So first we have to see if the selection is enclosed in ()
            # If yes, then nothing special needs to happen
            # If no, then if + or - is in the selection the selection needs to be enclosed in () first
            if (self.fitness[self.fitnessStart] in "+-*/"):
                self.fitness = self.fitness[0:self.fitnessStart] + "/" + self.fitness[self.fitnessStart+1:]
            elif (self.fitness[self.fitnessEnd-1] in "+-*/"):
                self.fitness = self.fitness[0:self.fitnessEnd-1] + "/" + self.fitness[self.fitnessEnd:]
            else:
                if (self.fitness[self.fitnessStart] != "(" or self.fitness[self.fitnessEnd-1] != ")"):
                    needParen = False
                    for i in range(self.fitnessStart, self.fitnessEnd):
                        if (self.fitness[i] in "+-"):
                            needParen = True
                            break
                    # Add the () in
                    if (needParen):
                        self.fitness = self.fitness[0:self.fitnessStart] + "(" + self.fitness[self.fitnessStart:self.fitnessEnd] + ")" + self.fitness[self.fitnessEnd:]
                        self.fitnessEnd = self.fitnessEnd + 2
                self.fitness = self.fitness[0:self.fitnessEnd] + "/?" + self.fitness[self.fitnessEnd:]
        elif (operation.startswith("Parentheses")):
            # The following functions are not compatible with trailing operators
            if (not(self.fitness[self.fitnessStart] in "+-*/") and not(self.fitness[self.fitnessEnd-1] in "+-*/")):
                self.fitness = self.fitness[0:self.fitnessStart] + "(" + self.fitness[self.fitnessStart:self.fitnessEnd] + ")" + self.fitness[self.fitnessEnd:]
                self.fitnessEnd = self.fitnessEnd + 2
        elif (operation.startswith("Absolute Value")):
            # The following functions are not compatible with trailing operators
            if (not(self.fitness[self.fitnessStart] in "+-*/") and not(self.fitness[self.fitnessEnd-1] in "+-*/")):
                self.fitness = self.fitness[0:self.fitnessStart] + "abs(" + self.fitness[self.fitnessStart:self.fitnessEnd] + ")" + self.fitness[self.fitnessEnd:]
                self.fitnessEnd = self.fitnessEnd + 5
        elif (operation.startswith("Exponential")):
            # The following functions are not compatible with trailing operators
            if (not(self.fitness[self.fitnessStart] in "+-*/") and not(self.fitness[self.fitnessEnd-1] in "+-*/")):
                self.fitness = self.fitness[0:self.fitnessStart] + "exp(" + self.fitness[self.fitnessStart:self.fitnessEnd] + ")" + self.fitness[self.fitnessEnd:]
                self.fitnessEnd = self.fitnessEnd + 5
        elif (operation.startswith("Natural Logarithm")):
            # The following functions are not compatible with trailing operators
            if (not(self.fitness[self.fitnessStart] in "+-*/") and not(self.fitness[self.fitnessEnd-1] in "+-*/")):
                self.fitness = self.fitness[0:self.fitnessStart] + "ln(" + self.fitness[self.fitnessStart:self.fitnessEnd] + ")" + self.fitness[self.fitnessEnd:]
                self.fitnessEnd = self.fitnessEnd + 4
        elif (operation.startswith("Power Function")):
            # The following functions are not compatible with trailing operators
            if (not(self.fitness[self.fitnessStart] in "+-*/") and not(self.fitness[self.fitnessEnd-1] in "+-*/")):
                self.fitness = self.fitness[0:self.fitnessStart] + "pow(" + self.fitness[self.fitnessStart:self.fitnessEnd] + ",?)" + self.fitness[self.fitnessEnd:]
                self.fitnessEnd = self.fitnessEnd + 7
        elif (operation.startswith("Square Root")):
            # The following functions are not compatible with trailing operators
            if (not(self.fitness[self.fitnessStart] in "+-*/") and not(self.fitness[self.fitnessEnd-1] in "+-*/")):
                self.fitness = self.fitness[0:self.fitnessStart] + "sqrt(" + self.fitness[self.fitnessStart:self.fitnessEnd] + ")" + self.fitness[self.fitnessEnd:]
                self.fitnessEnd = self.fitnessEnd + 6
        elif (operation.startswith("Boolean Equality")):
            # The following functions are not compatible with trailing operators
            if (not(self.fitness[self.fitnessStart] in "+-*/") and not(self.fitness[self.fitnessEnd-1] in "+-*/")):
                self.fitness = self.fitness[0:self.fitnessStart] + "eq(" + self.fitness[self.fitnessStart:self.fitnessEnd] + ",?)" + self.fitness[self.fitnessEnd:]
                self.fitnessEnd = self.fitnessEnd + 6
        elif (operation.startswith("Boolean Greater Than Equal")):
            # The following functions are not compatible with trailing operators
            if (not(self.fitness[self.fitnessStart] in "+-*/") and not(self.fitness[self.fitnessEnd-1] in "+-*/")):
                self.fitness = self.fitness[0:self.fitnessStart] + "gte(" + self.fitness[self.fitnessStart:self.fitnessEnd] + ",?)" + self.fitness[self.fitnessEnd:]
                self.fitnessEnd = self.fitnessEnd + 7
        elif (operation.startswith("Boolean Less Than Equal")):
            # The following functions are not compatible with trailing operators
            if (not(self.fitness[self.fitnessStart] in "+-*/") and not(self.fitness[self.fitnessEnd-1] in "+-*/")):
                self.fitness = self.fitness[0:self.fitnessStart] + "lte(" + self.fitness[self.fitnessStart:self.fitnessEnd] + ",?)" + self.fitness[self.fitnessEnd:]
                self.fitnessEnd = self.fitnessEnd + 7
        elif (operation.startswith("Boolean Greater Than")):
            # The following functions are not compatible with trailing operators
            if (not(self.fitness[self.fitnessStart] in "+-*/") and not(self.fitness[self.fitnessEnd-1] in "+-*/")):
                self.fitness = self.fitness[0:self.fitnessStart] + "gt(" + self.fitness[self.fitnessStart:self.fitnessEnd] + ",?)" + self.fitness[self.fitnessEnd:]
                self.fitnessEnd = self.fitnessEnd + 6
        elif (operation.startswith("Boolean Less Than")):
            # The following functions are not compatible with trailing operators
            if (not(self.fitness[self.fitnessStart] in "+-*/") and not(self.fitness[self.fitnessEnd-1] in "+-*/")):
                self.fitness = self.fitness[0:self.fitnessStart] + "lt(" + self.fitness[self.fitnessStart:self.fitnessEnd] + ",?)" + self.fitness[self.fitnessEnd:]
                self.fitnessEnd = self.fitnessEnd + 6
        elif (operation.startswith("Boolean Or")):
            # The following functions are not compatible with trailing operators
            if (not(self.fitness[self.fitnessStart] in "+-*/") and not(self.fitness[self.fitnessEnd-1] in "+-*/")):
                self.fitness = self.fitness[0:self.fitnessStart] + "or(" + self.fitness[self.fitnessStart:self.fitnessEnd] + ",?)" + self.fitness[self.fitnessEnd:]
                self.fitnessEnd = self.fitnessEnd + 6
        elif (operation.startswith("Boolean And")):
            # The following functions are not compatible with trailing operators
            if (not(self.fitness[self.fitnessStart] in "+-*/") and not(self.fitness[self.fitnessEnd-1] in "+-*/")):
                self.fitness = self.fitness[0:self.fitnessStart] + "and(" + self.fitness[self.fitnessStart:self.fitnessEnd] + ",?)" + self.fitness[self.fitnessEnd:]
                self.fitnessEnd = self.fitnessEnd + 7
        elif (operation.startswith("Boolean Not")):
            # The following functions are not compatible with trailing operators
            if (not(self.fitness[self.fitnessStart] in "+-*/") and not(self.fitness[self.fitnessEnd-1] in "+-*/")):
                self.fitness = self.fitness[0:self.fitnessStart] + "not(" + self.fitness[self.fitnessStart:self.fitnessEnd] + ")" + self.fitness[self.fitnessEnd:]
                self.fitnessEnd = self.fitnessEnd + 5
        elif (operation.startswith("If-Then-Else")):
            # The following functions are not compatible with trailing operators
            if (not(self.fitness[self.fitnessStart] in "+-*/") and not(self.fitness[self.fitnessEnd-1] in "+-*/")):
                self.fitness = self.fitness[0:self.fitnessStart] + "ite(" + self.fitness[self.fitnessStart:self.fitnessEnd] + ",?,?)" + self.fitness[self.fitnessEnd:]
                self.fitnessEnd = self.fitnessEnd + 9
        self.txtFitness.SetReadOnly(False) # It won't update if readonly is True
        self.txtFitness.SetText(self.fitness)
        self.txtFitness.SetReadOnly(True) # But we need readonly to be false so the user doesn't mess with this control
        self.txtFitness.SetSelection(self.fitnessStart, self.fitnessEnd)
        self.expandFitnessSelection()
    
    def removeOperation(self, event):
        # If the selection starts with an operator, remove the operator and the following selection
        if (self.fitness[self.fitnessStart] in "+-*/"):
            # If the last character is also an operator, ignore it but remove everything before
            if (self.fitness[self.fitnessEnd-1] in "+-*/"):
                self.fitness = self.fitness[:self.fitnessStart] + self.fitness[self.fitnessEnd-1:]
            else:
                self.fitness = self.fitness[:self.fitnessStart] + self.fitness[self.fitnessEnd:]
            self.fitnessEnd = self.fitnessStart
        elif (self.fitness[self.fitnessEnd-1] in "+-*/"):
            self.fitness = self.fitness[:self.fitnessStart] + self.fitness[self.fitnessEnd:]
            self.fitnessEnd = self.fitnessStart
        elif (self.fitness[self.fitnessStart] == "(" and self.fitness[self.fitnessEnd-1] == ")"):
            # Just take the parentheses away
            self.fitness = self.fitness[:self.fitnessStart] + self.fitness[self.fitnessStart+1:self.fitnessEnd-1] + self.fitness[self.fitnessEnd:]
            self.fitnessEnd = self.fitnessStart
        elif (self.fitness[self.fitnessStart] in "abcdefghijklmnopqrstuvwxyz"):
            # For removal of functions, turn the function into a ? so the user cannot turn the original ?
            # into a function, remove the function, and then have a blank textbox which will not allow you to 
            # do anything at all
            self.fitness = self.fitness[:self.fitnessStart] + "?" + self.fitness[self.fitnessEnd:]
            self.fitnessEnd = self.fitnessStart
        self.txtFitness.SetReadOnly(False) # It won't update if readonly is True
        self.txtFitness.SetText(self.fitness)
        self.txtFitness.SetReadOnly(True) # But we need readonly to be false so the user doesn't mess with this control
        self.txtFitness.SetSelection(self.fitnessStart, self.fitnessEnd)
        self.expandFitnessSelection()
    
    def expandFitnessSelection(self):
        # This function expands what the user selected in the fitness function window so it actually
        # is selecting something that is validly editable
        # If there is nothing selected, just select the first thing after the cursor, then expand
        # If it is at the end of the function, grab the element before
        self.txtFitness.StyleSetSpec(2, "back:#FFFFFF")
        self.txtFitness.StyleSetSpec(3, "back:#FFFF00")
        self.txtFitness.StartStyling(0, 0xffffff)
        self.txtFitness.SetStyling(len(self.txtFitness.GetText()), 2)
        (start, end) = self.txtFitness.GetSelection()
        if (start == end and start < len(self.txtFitness.GetText())):
            # Grab the next element
            end = end + 1
        elif (start == end and start == len(self.txtFitness.GetText())):
            # Grab the last element
            start = start - 1
        # First off, is the end of the selection in the middle of a function name?  If yes, extend it out to the (
        if (self.fitness[end-1] in "abcdefghijklmnopqrstuvwxyz"):
            for i in range(end, len(self.fitness)):
                if (self.fitness[i] == "("):
                    end = end + 1 # end is not inclusive
                    break
                end = end + 1
        # Now we have to match commas, which is tricky
        # First count all the commas, then start at the last comma and rewind, keeping track of
        # all the functions passed.  When you've covered all the commas, stop
        numcommas = 0
        for i in range(start, end):
            if (self.fitness[i] == ","):
                numcommas = numcommas + 1
        commaindx = self.fitness[start:end].rfind(",") + start
        if (commaindx >= start):
            # Find the first function name preceeding this comma
            functionindx = commaindx - 1
            while (functionindx >= 0):
                if (self.fitness[functionindx] in "abcdefghijklmnopqrstuvwxyz"):
                    # Get the whole function name
                    endfunction = functionindx + 1
                    while (functionindx >= 0 and self.fitness[functionindx] in "abcdefghijklmnopqrstuvwxyz"):
                        functionindx = functionindx - 1
                    # Get back to the start of the function
                    functionindx = functionindx + 1
                    # Subtract commas off of the total
                    numcommas = numcommas - self.numCommas[self.fitness[functionindx:endfunction]]
                    # If numcommas <= 0, get out of here
                    if (numcommas <= 0):
                        break
                elif (self.fitness[functionindx] == "," and functionindx < start):
                    numcommas = numcommas + 1
                functionindx = functionindx - 1
            if (functionindx < 0):
                start = 0
            elif (functionindx < start):
                start = functionindx
        # How many ) do we have at the beginning, because we need to match at least an equal amount of (
        parentheses = 0
        for i in range(start, end):
            if (self.fitness[i] == ")"):
                parentheses = parentheses + 1
            elif (self.fitness[i] == "("):
                parentheses = parentheses - 1
        # If we have a positive value for parentheses, then we need to at least match an equal number of (
        if (parentheses > 0):
            for i in range(start-1, -1, -1):
                if (self.fitness[i] == "("):
                    parentheses = parentheses - 1
                    if (parentheses == 0):
                        # If this is a function call, back up one so the next step knows to complete the function name
                        if (i > 0 and self.fitness[i-1] in "abcdefghijklmnopqrstuvwxyz"):
                            start = start - 1
                        break
                elif (self.fitness[i] == ")"):
                    parentheses = parentheses + 1
                start = start - 1
        # If this is a function call, back up one so the next step knows to complete the function name
        if (start > 0 and self.fitness[start-1] in "abcdefghijklmnopqrstuvwxyz"):
            start = start - 1
        # Now that we've matched all the parentheses, we're either at a (, in which case do nothing, or we are
        # in the middle of a literal or function name, so just go to the end of this name/literal
        if (self.fitness[start] in "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789?abcdefghijklmnopqrstuvwxyz"):
            for i in range(start-1, -1, -1):
                if (not(self.fitness[i] in "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789?abcdefghijklmnopqrstuvwxyz")):
                    break
                start = start - 1
        # The parentheses might not be matched, because we may have had more ( originally and need to match some )
        parentheses = 0
        for i in range(end-1, start-1, -1):
            if (self.fitness[i] == "("):
                parentheses = parentheses + 1
            elif (self.fitness[i] == ")"):
                parentheses = parentheses - 1
        # Now match all of the )
        if (parentheses > 0):
            for i in range(end, len(self.fitness)+1):
                if (i < len(self.fitness)):
                    if (self.fitness[i] == ")"):
                        parentheses = parentheses - 1
                        if (parentheses == 0):
                            end = end + 1 # Need to do this because end is not inclusive
                            break
                    elif (self.fitness[i] == "("):
                        parentheses = parentheses + 1
                end = end + 1
        # If we're at a literal, finish up the literal
        if (self.fitness[end-1] in "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789?"):
            for i in range(end, len(self.fitness)+1):
                if (i < len(self.fitness)):
                    if (not(self.fitness[i] in "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789?")):
                        break
                end = end + 1
        if (end > len(self.fitness)):
            end = len(self.fitness)
        # Highlight the selection
        self.txtFitness.StartStyling(start, 0xffffff)
        self.txtFitness.SetStyling(end-start, 3)
        #self.txtFitness.StartStyling(end-1, 0xffffff)
        #self.txtFitness.SetStyling(1, 2)
        self.fitnessStart = start
        self.fitnessEnd = end
        # Now figure out if the literal text box should be enabled
        # It will not be enabled if the first or last selected character is an operator
        if (self.fitness[start] in "+-*/"):
            self.txtStateLiteral.Disable()
        else:
            self.txtStateLiteral.Enable()
            # Give this textbox the focus so they can start editing without having to click it
            self.txtStateLiteral.SetFocus()
            self.txtStateLiteral.SetSelection(len(self.txtStateLiteral.GetValue()), len(self.txtStateLiteral.GetValue()))
    
    def fitnessRelease(self, event):
        self.expandFitnessSelection()
        event.Skip()
        
    def stateliteralUpdate(self, event):
        # Do not add the text if it is invalid
        text = self.txtStateLiteral.GetValue().upper()
        validstates = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[0:self.grdEntity.NumberCols-1]
        badtext = False
        decimals = 0
        for i in range(0, len(text)):
            if (text[i] == "."):
                decimals = decimals + 1
                if (decimals >= 2):
                    badtext = True
                    break
            elif (i == 0 and not(text[i] in validstates) and not(text[i] in "0123456789")):
                badtext = True
                break
            elif (i > 0 and not(text[i] in "0123456789")):
                badtext = True
                break
        if (len(text) == 0):
            # Empty, change this position back to a question mark
            self.fitness = self.fitness[0:self.fitnessStart] + "?" + self.fitness[self.fitnessStart+1:]
        elif (badtext):
            self.fitness = self.fitness[0:self.fitnessStart] + "?" + self.fitness[self.fitnessEnd:]
            self.fitnessEnd = self.fitnessStart + 1
        else:
            self.fitness = self.fitness[0:self.fitnessStart] + text + self.fitness[self.fitnessEnd:]
            self.fitnessEnd = self.fitnessStart + len(text)
        self.txtFitness.SetReadOnly(False) # It won't update if readonly is True
        self.txtFitness.SetText(self.fitness)
        self.txtFitness.SetReadOnly(True) # But we need readonly to be false so the user doesn't mess with this control
        self.txtFitness.SetSelection(self.fitnessStart, self.fitnessEnd)
        self.expandFitnessSelection()
        event.Skip()
    
    def loadMSD(self, event):
        logInfo("Load MSD button clicked")
        #try:
        # Load data from a previous MSD session
        if (len(self.stateoptions) > 0):
            dlg = wx.MessageDialog(self, "The data in your current workflow will be lost and replaced with a loaded session.  Are you sure you want to proceed?", "Current MSD Data Will Be Lost", wx.YES_NO | wx.ICON_EXCLAMATION | wx.CENTRE)
            if (dlg.ShowModal() == wx.ID_YES):
                # Don't clear it out until the user actually selects a resfile (instead of cancelling
                # at the file selection dialog
                pass
            else:
                logInfo("Load MSD operation cancelled due to data already being in the session")
                return
        dlg = wx.FileDialog(
            self, message="Choose a File",
            defaultDir=self.seqWin.cwd,
            defaultFile="",
            wildcard="MSD Sessions (*.msd)|*.msd",
            style=wx.OPEN | wx.CHANGE_DIR)
        if (dlg.ShowModal() == wx.ID_OK):
            dlg.Destroy()
            paths = dlg.GetPaths()
            # Change cwd to the last opened file
            lastDirIndx = paths[len(paths)-1].rfind("/")
            self.seqWin.cwd = str(paths[len(paths)-1][0:lastDirIndx])
            self.seqWin.saveWindowData(None)
            filename = str(paths[0])
            f = open(filename, "r")
            #try:
            models = []
            states = []
            resfile = []
            assignments = []
            vector = []
            for aline in f:
                if (aline[0:6] == "MODELS"):
                    for model in aline.split("\t")[1:]:
                        states.append([model.strip()])
                        resfile.append([])
                        models.append(model.strip())
                        assignments.append("")
                        vector.append("MIN")
                    # Now let's see if the user has all these models loaded
                    modelsnotfound = models[:]
                    for i in range(1, self.grdEntity.NumberCols):
                        try:
                            modelsnotfound.pop(modelsnotfound.index(self.grdEntity.GetColLabelValue(i)))
                        except:
                            # Not here, leave it alone
                            pass
                    if (len(modelsnotfound) > 0):
                        msg = "It seems that model"
                        if (len(modelsnotfound) == 1):
                            msg = msg + " " + modelsnotfound[0] + " is not loaded."
                        else:
                            msg = msg + "s "
                            for i in range(0, len(modelsnotfound)):
                                if (i < len(modelsnotfound)-1):
                                    msg = msg + modelsnotfound[i] + ", "
                                else:
                                    msg = msg + "and " + modelsnotfound[i] + " are not loaded."
                        msg = msg + "  The data can still be loaded, but data for the missing models will be left out.  Proceed?"
                        dlg = wx.MessageDialog(self, msg, "Models Missing", wx.YES_NO | wx.ICON_EXCLAMATION | wx.CENTRE)
                        if (dlg.ShowModal() == wx.ID_NO):
                            f.close()
                            dlg.Destroy()
                            return
                        dlg.Destroy()
                elif (aline[0:10] == "ENTVARIPOS"):
                    nEntries = int(aline.split("\t")[1])
                elif (aline[0:12] == "STATEOPTIONS"):
                    stateoptions = []
                    for option in aline.split("\t")[1:]:
                        stateoptions.append(option.strip())
                elif (aline[0:4] == "CORR"):
                    ID = aline.split("\t")[1].strip()
                    model = ID[0:len(ID)-2]
                    indx = models.index(ID)
                    for correspondence in aline.split("\t")[2:]:
                        if (len(correspondence) > 0):
                            chain = correspondence[0]
                            seqpos = correspondence[3:].strip()
                            # Does this residue exist on the model, because if not then we should just have a blank
                            # space (the user could have changed the model since this session was made)
                            if (self.seqWin.doesResidueExist(model, chain, seqpos)):
                                states[indx].append(correspondence.strip())
                            else:
                                states[indx].append("")
                        else:
                            # Blank entry
                            states[indx].append("")
                elif (aline[0:9] == "SECONDARY"):
                    ID = aline.split("\t")[1].strip()
                    model = ID[0:len(ID)-2]
                    for indx in range(0, len(models)):
                        if (models[indx][0:len(models[indx])-2] == model):
                            break
                    #indx = models.index(ID)
                    lineindx = 2
                    data = aline.split("\t")
                    while (lineindx < len(data)):
                        # Again, make sure this residue exists in the model or don't add it
                        chain = data[lineindx+4].strip()
                        seqpos = data[lineindx+2].strip()
                        if (self.seqWin.doesResidueExist(model, chain, seqpos)):
                            resfile[indx].append([int(data[lineindx]), int(data[lineindx+1]), data[lineindx+2].strip(), int(data[lineindx+3]), data[lineindx+4].strip(), int(data[lineindx+5]), data[lineindx+6].strip()])
                        lineindx = lineindx + 7
                elif (aline[0:5] == "STATE"):
                    model = aline.split("\t")[1].strip()
                    indx = models.index(model)
                    if (len(aline.split("\t")) > 2):
                        assignments[indx] = aline.split("\t")[2].strip()
                    else:
                        assignments[indx] = ""
                elif (aline[0:6] == "VECTOR"):
                    state = aline.split("\t")[1].strip()
                    indx = "ABCDEFGHIJKLMNOPQRSTUVWXYZ".index(state)
                    vector[indx] = aline.split("\t")[2].strip()
                elif (aline[0:7] == "FITNESS"):
                    #fitness = aline.split("\t")[1].strip()
                    # I'm not going to read the fitness function right now because if the order of things got
                    # changed then the letter assignments may not be the same which will ruin the equation
                    # Plus this forces the user to redo the fitness function which is good, because then they'll
                    # make sure all the rest of the information is correct
                    fitness = "?" 
            #except:
                # Some kind of error loading the data, maybe the file got corrupted and entries were either
                # missing or out-of-order
                #dlg = wx.MessageDialog(self, "An error was encountered loading this MSD session.  The file may be corrupt.", "Error Loading MSD Session", wx.OK | wx.ICON_ERROR | wx.CENTRE)
                #dlg.ShowModal()
                #dlg.Destroy()
            f.close()
            # Now we have to play with the states
            # If a model was missing, then the number of available letters is lower than in the original session
            # For instance, if originally we had 3 models, A, B, and C are available, but only two models are
            # actually loaded, so only A and B are available, but one of the loaded models was C originally
            # So we have to rename them but keep them in the same groups (the letters don't really matter, only
            # that the same groups are formed)
            # Easiest way to do this is group them based on the old names and then just assign valid letters to
            # all the models in the groups
            grouping = []
            for i in range(0, len(models)):
                grouping.append([])
            for i in range(0, len(models)):
                model = models[i]
                state = assignments[i]
                indx = "ABCDEFGHIJKLMNOPQRSTUVWXYZ".index(state)
                grouping[indx].append(model)
            # Now take out all the blank groups
            for i in range(len(grouping)-1, -1, -1):
                if (len(grouping[i]) == 0):
                    grouping.pop(i)
            # Now reassign proper letters based on what bins they are grouped in
            for i in range(0, len(models)):
                for j in range(0, len(grouping)):
                    try:
                        indx = grouping[j].index(models[i])
                        assignments[i] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[j]
                    except:
                        # Not here
                        pass
            # Since I fixed the issue with not allowing individual chains in the entities file, the assignments can get messed
            # up by the preceding code, so now I have to set all the second+ chains to have no state
            for i in range(1, len(assignments)):
                prevmodel = states[i-1][0]
                prevmodel = prevmodel[0:len(prevmodel)-2]
                thismodel = states[i][0]
                thismodel = thismodel[0:len(thismodel)-2]
                if (prevmodel == thismodel):
                    assignments[i] = ""
            # Now put this loaded data into the proper data structures
            for i in range(0, len(self.states)):
                model = self.states[i][0]
                try:
                    indx = models.index(model)
                    self.states[i] = states[indx][:]
                    self.resfile[i] = resfile[indx][:]
                    self.stateAssignments[i] = assignments[indx]
                    self.vectorFunctions[i] = vector[indx]
                except:
                    pass
            # Now we have to go through all the resfile data and make sure it is up to date
            # The reason for this is that the old resfile data had row and column information from the Sequence
            # Window, which will not be accurate if the orders of things is different
            for i in range(0, len(self.resfile)):
                ID = self.states[i][0]
                model = ID[0:len(ID)-2]
                for j in range(len(self.resfile[i])-1, -1, -1):
                    chain = self.resfile[i][j][4]
                    seqpos = self.resfile[i][j][2]
                    # Does it still exist, get it's coordinate
                    ret = self.seqWin.doesResidueExist(model, chain, seqpos)
                    if (ret):
                        self.resfile[i][j][0] = ret[1]
                        self.resfile[i][j][1] = ret[0]
                    else:
                        # Get rid of this entry
                        self.resfile[i].pop(j)
            # Now we have to make sure all of the entries in self.states are of the same lengths
            # They wouldn't be necessarily if there are extra models already loaded that are not in the MSD file
            maxlen = 0
            for i in range(0, len(self.states)):
                if (len(self.states[i]) > maxlen):
                    maxlen = len(self.states[i])
            for i in range(0, len(self.states)):
                if (len(self.states[i]) < maxlen):
                    for j in range(0, maxlen-len(self.states[i])):
                        self.states[i].append("")
            self.updateResfile()
            j = 0
            for i in range(0, len(self.stateAssignments)):
                if (len(self.stateAssignments[i].strip()) > 0):
                    self.grdStates.SetCellValue(j, 0, self.stateAssignments[i])
                    j = j + 1
            self.stateoptions = stateoptions[:]
            self.fitness = fitness
            self.redrawEntity()
            self.txtFitness.SetReadOnly(False) # It won't update if readonly is True
            self.txtFitness.SetText(self.fitness)
            self.txtFitness.SetReadOnly(True) # But we need readonly to be false so the user doesn't mess with this control
            self.expandFitnessSelection()
        else:
            dlg.Destroy()
            logInfo("Load MSD operation cancelled at the file select window")
        #except:
            #dlg = wx.MessageDialog(self, "An error was encountered attempting to import your data.  The save file may be corrupt a model in the save file does not match the structure of the model loaded.", "Error Encountered", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
            #dlg.ShowModal()
            #dlg.Destroy()
    
    def saveMSD(self, event):
        if (len(self.states) == 0):
            wx.MessageBox("There's nothing to save!", "No MSD Data", wx.OK|wx.ICON_EXCLAMATION)
            return
        logInfo("Clicked the Save MSD button")
        # Save the data in the MSD controls for easy access at a later time
        dlg = wx.FileDialog(
            self, message="Save an MSD session",
            defaultDir=self.seqWin.cwd,
            defaultFile="",
            wildcard="MSD Sessions (*.msd)|*.msd",
            style=wx.SAVE | wx.CHANGE_DIR)
        if (dlg.ShowModal() == wx.ID_OK):
            paths = dlg.GetPaths()
            # Change cwd to the last opened file
            lastDirIndx = paths[len(paths)-1].rfind("/")
            self.seqWin.cwd = str(paths[len(paths)-1][0:lastDirIndx])
            self.seqWin.saveWindowData(None)
            filename = str(paths[0]).split(".msd")[0] + ".msd"
            # Does it exist already?  If so, ask if the user really wants to overwrite it
            if (os.path.isfile(filename)):
                dlg2 = wx.MessageDialog(self, "The file " + filename + " already exists.  Overwrite it?", "Filename Already Exists", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
                if (dlg2.ShowModal() == wx.ID_NO):
                    dlg2.Destroy()
                    logInfo("Canceled save operation due to filename already existing")
            # Start saving the MSD session data
            # First let's save a list of all the models so when this session is loaded we can easily see
            # if all the same models are available in the SequenceWindow
            f = open(filename, "w")
            f.write("MODELS")
            for i in range(0, len(self.states)):
                f.write("\t" + self.states[i][0]) # Remember the first element of each sublist is the model|chain
            f.write("\n")
            # Now write how many entries are in the entity resfile
            f.write("ENTVARIPOS\t" + str(self.grdEntity.NumberRows) + "\n")
            # Now write out all the resfile strings for these positions
            f.write("STATEOPTIONS")
            for i in range(0, len(self.stateoptions)):
                f.write("\t" + self.stateoptions[i])
            f.write("\n")
            # Now write all the correspondences in the entity resfile
            # It is important that you delimit by tabs, since the loader will have to split fields by tabs
            # You cannot split by whitespace because some of the correspondences might be empty, so you need
            # to split those blank fields properly
            for i in range(0, len(self.states)):
                f.write("CORR\t" + self.states[i][0])
                for j in range(1, len(self.states[i])):
                    f.write("\t" + str(self.states[i][j]))
                f.write("\n")
            # Now write out all the secondary resfiles
            # There are potentially a lot of fields here, first field is SECONDARY, then the model name, then
            # groups of 7 for the resfile data
            for i in range(0, len(self.resfile)):
                if (len(self.resfile[i]) == 0):
                    continue
                f.write("SECONDARY\t" + self.states[i][0])
                for j in range(0, len(self.resfile[i])):
                    for data in self.resfile[i][j]:
                        f.write("\t" + str(data))
                f.write("\n")
            # Now write out the state assignments
            for i in range(0, len(self.stateAssignments)):
                f.write("STATE\t" + str(self.states[i][0]) + "\t" + self.stateAssignments[i] + "\n")
            # Now write out all the vector operations that will be applied to these states
            for i in range(0, len(self.vectorFunctions)):
                f.write("VECTOR\t" + "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[i] + "\t" + str(self.vectorFunctions[i]) + "\n")
            # Finally, write out the fitness function
            f.write("FITNESS\t" + self.fitness + "\n")
            f.close()
        else:
            logInfo("Cancelled save resfile operation")
    
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
                self.btnAminoA.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnA.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoA.SetForegroundColour("#000000")
            logInfo("Removed ALA from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoA.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnA_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoA.SetForegroundColour("#FF0000")
            logInfo("Added ALA to palette")
        self.addAAToPalette("A")
    
    def aminoC(self, event):
        if ("C" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoC.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnC.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoC.SetForegroundColour("#000000")
            logInfo("Removed CYS from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoC.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnC_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoC.SetForegroundColour("#FF0000")
            logInfo("Added CYS to palette")
        self.addAAToPalette("C")
    
    def aminoD(self, event):
        if ("D" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoD.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnD.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoD.SetForegroundColour("#000000")
            logInfo("Removed ASP from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoD.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnD_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoD.SetForegroundColour("#FF0000")
            logInfo("Added ASP to palette")
        self.addAAToPalette("D")
    
    def aminoE(self, event):
        if ("E" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoE.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnE.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoE.SetForegroundColour("#000000")
            logInfo("Removed GLU from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoE.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnE_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoE.SetForegroundColour("#FF0000")
            logInfo("Added GLU to palette")
        self.addAAToPalette("E")
        
    def aminoF(self, event):
        if ("F" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoF.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnF.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoF.SetForegroundColour("#000000")
            logInfo("Removed PHE from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoF.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnF_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoF.SetForegroundColour("#FF0000")
            logInfo("Added PHE to palette")
        self.addAAToPalette("F")
        
    def aminoG(self, event):
        if ("G" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoG.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnG.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoG.SetForegroundColour("#000000")
            logInfo("Removed GLY from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoG.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnG_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoG.SetForegroundColour("#FF0000")
            logInfo("Added GLY to palette")
        self.addAAToPalette("G")
        
    def aminoH(self, event):
        if ("H" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoH.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnH.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoH.SetForegroundColour("#000000")
            logInfo("Removed HIS from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoH.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnH_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoH.SetForegroundColour("#FF0000")
            logInfo("Added HIS to palette")
        self.addAAToPalette("H")
        
    def aminoI(self, event):
        if ("I" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoI.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnI.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoI.SetForegroundColour("#000000")
            logInfo("Removed ILE from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoI.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnI_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoI.SetForegroundColour("#FF0000")
            logInfo("Added ILE to palette")
        self.addAAToPalette("I")
        
    def aminoK(self, event):
        if ("K" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoK.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnK.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoK.SetForegroundColour("#000000")
            logInfo("Removed LYS from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoK.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnK_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoK.SetForegroundColour("#FF0000")
            logInfo("Added LYS to palette")
        self.addAAToPalette("K")
        
    def aminoL(self, event):
        if ("L" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoL.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnL.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoL.SetForegroundColour("#000000")
            logInfo("Removed LEU from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoL.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnL_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoL.SetForegroundColour("#FF0000")
            logInfo("Added LEU to palette")
        self.addAAToPalette("L")
        
    def aminoM(self, event):
        if ("M" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoM.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnM.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoM.SetForegroundColour("#000000")
            logInfo("Removed MET from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoM.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnM_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoM.SetForegroundColour("#FF0000")
            logInfo("Added MET to palette")
        self.addAAToPalette("M")
    
    def aminoN(self, event):
        if ("N" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoN.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnN.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoN.SetForegroundColour("#000000")
            logInfo("Removed ASN from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoN.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnN_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoN.SetForegroundColour("#FF0000")
            logInfo("Added ASN to palette")
        self.addAAToPalette("N")
    
    def aminoP(self, event):
        if ("P" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoP.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnP.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoP.SetForegroundColour("#000000")
            logInfo("Removed PRO from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoP.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnP_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoP.SetForegroundColour("#FF0000")
            logInfo("Added PRO to palette")
        self.addAAToPalette("P")
    
    def aminoQ(self, event):
        if ("Q" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoQ.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnQ.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoQ.SetForegroundColour("#000000")
            logInfo("Removed GLN from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoQ.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnQ_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoQ.SetForegroundColour("#FF0000")
            logInfo("Added GLN to palette")
        self.addAAToPalette("Q")
        
    def aminoR(self, event):
        if ("R" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoR.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnR.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoR.SetForegroundColour("#000000")
            logInfo("Removed ARG from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoR.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnR_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoR.SetForegroundColour("#FF0000")
            logInfo("Added ARG to palette")
        self.addAAToPalette("R")
        
    def aminoS(self, event):
        if ("S" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoS.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnS.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoS.SetForegroundColour("#000000")
            logInfo("Removed SER from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoS.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnS_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoS.SetForegroundColour("#FF0000")
            logInfo("Added SER to palette")
        self.addAAToPalette("S")
        
    def aminoT(self, event):
        if ("T" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoT.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnT.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoT.SetForegroundColour("#000000")
            logInfo("Removed THR from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoT.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnT_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoT.SetForegroundColour("#FF0000")
            logInfo("Added THR to palette")
        self.addAAToPalette("T")
        
    def aminoV(self, event):
        if ("V" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoV.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnV.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoV.SetForegroundColour("#000000")
            logInfo("Removed VAL from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoV.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnV_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoV.SetForegroundColour("#FF0000")
            logInfo("Added VAL to palette")
        self.addAAToPalette("V")
        
    def aminoW(self, event):
        if ("W" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoW.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnW.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoW.SetForegroundColour("#000000")
            logInfo("Removed TRP from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoW.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnW_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoW.SetForegroundColour("#FF0000")
            logInfo("Added TRP to palette")
        self.addAAToPalette("W")
        
    def aminoY(self, event):
        if ("Y" in self.palette):
            if (platform.system() == "Darwin"):
                self.btnAminoY.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnY.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnAminoY.SetForegroundColour("#000000")
            logInfo("Removed TYR from palette")
        else:
            if (platform.system() == "Darwin"):
                self.btnAminoY.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnY_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
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
                self.btnAminoA.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnA.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoC.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnC.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoD.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnD.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoE.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnE.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoF.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnF.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoG.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnG.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoH.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnH.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoI.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnI.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoK.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnK.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoL.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnL.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoM.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnM.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoN.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnN.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoP.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnP.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoQ.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnQ.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoR.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnR.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoS.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnS.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoT.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnT.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoV.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnV.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoW.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnW.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoY.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnY.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
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
                self.btnAminoA.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnA_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoC.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnC_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoD.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnD_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoE.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnE_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoF.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnF_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoG.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnG_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoH.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnH_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoI.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnI_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoK.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnK_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoL.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnL_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoM.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnM_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoN.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnN_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoP.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnP_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoQ.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnQ_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoR.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnR_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoS.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnS_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoT.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnT_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoV.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnV_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoW.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnW_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnAminoY.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/msd/btnY_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
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
        logInfo("Add button clicked")
        # For each of the selected entries, first verify that this entry is not already in the resfile and if it
        # isn't then add it in
        for i in range(0, len(self.selectedData)):
            [indx, r, seqpos, poseindx, chainID, chainoffset] = self.selectedData[i]
            if (self.seqWin.getModelForChain(r) != self.resfileModel):
                # Ignore selections from models not currently selected
                continue
            # Make sure this is a CAA
            if (not(self.seqWin.getIsCanonicalAA(r, indx))):
                continue
            alreadyIn = False
            for j in range(0, len(self.resfile[self.resfilemodelindx])):
                [rindx, rr, rseqpos, rposeindx, rchainID, rchainoffset, rreslist] = self.resfile[self.resfilemodelindx][j]
                if (r == rr and indx == rindx):
                    alreadyIn = True
                    break
            if (not(alreadyIn)):
                if (len(self.resfile[self.resfilemodelindx]) == 0):
                    # List empty, add new element
                    self.resfile[self.resfilemodelindx].append([indx, r, seqpos, poseindx, chainID, chainoffset, "NATAA"])
                elif (r < self.resfile[self.resfilemodelindx][0][1] or (r == self.resfile[self.resfilemodelindx][0][1] and indx < self.resfile[self.resfilemodelindx][0][0])):
                    # Belongs first
                    self.resfile[self.resfilemodelindx].insert(0, [indx, r, seqpos, poseindx, chainID, chainoffset, "NATAA"])
                else:
                    notInYet = True
                    # Maybe it belongs somewhere in the middle?
                    for i in range(0, len(self.resfile[self.resfilemodelindx])-1):
                        [indx1, r1, seqpos1, poseindx1, chainID1, chainoffset1, type1] = self.resfile[self.resfilemodelindx][i]
                        [indx2, r2, seqpos2, poseindx2, chainID2, chainoffset2, type2] = self.resfile[self.resfilemodelindx][i+1]
                        if (r == r1 and r == r2 and indx > indx1 and indx < indx2):
                            notInYet = False
                            self.resfile[self.resfilemodelindx].insert(i+1, [indx, r, seqpos, poseindx, chainID, chainoffset, "NATAA"])
                        elif (r == r1 and r < r2 and indx > indx1):
                            notInYet = False
                            self.resfile[self.resfilemodelindx].insert(i+1, [indx, r, seqpos, poseindx, chainID, chainoffset, "NATAA"])
                    if (notInYet):
                        # Belongs at the end
                        self.resfile[self.resfilemodelindx].append([indx, r, seqpos, poseindx, chainID, chainoffset, "NATAA"])
        self.updateResfile()
    
    def remove(self, event):
        self.activate()
        logInfo("Remove button clicked")
        # For each of the selected entries, find out if it is in the minmap and remove it if it is
        for i in range(0, len(self.selectedData)):
            [indx, r, seqpos, poseindx, chainID, chainoffset] = self.selectedData[i]
            for j in range(0, len(self.resfile[self.resfilemodelindx])):
                [rindx, rr, rseqpos, rposeindx, rchainID, rchainoffset, rreslist] = self.resfile[self.resfilemodelindx][j]
                if (r == rr and indx == rindx):
                    self.resfile[self.resfilemodelindx].pop(j)
                    break
        self.updateResfile()
    
    def restrict(self, event):
        self.activate()
        logInfo("Restrict button clicked")
        # Remove everything and add only the selected residues
        self.resfile[self.resfilemodelindx] = []
        self.add(event)
    
    def addAll(self, event):
        logInfo("Add all button clicked")
        # Add everything that is in the sequence viewer
        # Here "all" refers to everything in the protein already in the resfile or, if nothing is in
        # the resfile, everything in the first protein
        posetoadd = self.seqWin.getPoseIndexForModel(self.resfileModel)
        self.resfile[self.resfilemodelindx] = []
        allData = []
        for r in range(0, self.seqWin.SeqViewer.NumberRows):
            for c in range(0, len(self.seqWin.sequences[r])):
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
        self.resfile[self.resfilemodelindx] = []
        self.updateResfile()
    
    def applySelection(self, event):
        logInfo("Apply selection button clicked")
        if (self.desMenu.GetStringSelection() != ""):
            indx = self.resfile[self.currSelection][0]
            if (self.addType == "WT"):
                resAddType = self.seqWin.SeqViewer.GetCellValue(self.currSelection, indx)
            elif (self.addType == "Hydrophobic"):
                resAddType = "AFILMVWY"
            elif (self.addType == "Hydrophilic"):
                resAddType = "DEHKNQRST"
            elif (self.addType == "Aromatic"):
                resAddType = "FWY"
            else:
                resAddType = "ACDEFGHIKLMNPQRSTVWY"
            self.resfile[self.currSelection][6] = resAddType
            self.grdResfile.SetCellValue(self.currSelection, 0, self.resfile[self.currSelection][6])
    
    def addRes(self, event):
        # Add the chosen AA to the list of residue options at the selected position
        if (self.grdCol == 0 and self.grdRow >= 0 and self.resEntityMenu.GetStringSelection() != ""):
            AA3 = self.resEntityMenu.GetStringSelection()
            logInfo("Added " + AA3 + " to the selected position")
            AAindx = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR ".find(AA3)
            AA1 = "ACDEFGHIKLMNPQRSTVWY"[AAindx/4]
            reslist = self.stateoptions[self.grdRow]
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
            self.stateoptions[self.grdRow] = reslist
            self.grdEntity.SetCellValue(self.grdRow, 0, self.stateoptions[self.grdRow])
    
    def removeRes(self, event):
        # Remove the chosen AA to the list of residue options at the selected position
        if (self.grdCol == 0 and self.grdRow >= 0 and self.resEntityMenu.GetStringSelection() != ""):
            AA3 = self.resEntityMenu.GetStringSelection()
            logInfo("Removed " + AA3 + " from the selected position")
            AAindx = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR ".find(AA3)
            AA1 = "ACDEFGHIKLMNPQRSTVWY"[AAindx/4]
            reslist = self.stateoptions[self.grdRow]
            todeleteindx = reslist.find(AA1)
            if (todeleteindx >= 0 and len(reslist) > 1):
                reslist = reslist[0:todeleteindx] + reslist[todeleteindx+1:]
                self.stateoptions[self.grdRow] = reslist
                self.grdEntity.SetCellValue(self.grdRow, 0, self.stateoptions[self.grdRow])
            elif (todeleteindx >= 0 and len(reslist) == 1):
                # This was the only residue at this position so remove this entity resfile entry
                self.entityMenu.SetSelection(self.grdRow)
                self.removeEntityPosition(None)
    
    def MSDClick(self, event):
        # This is also the "Finalize!" button
        logInfo("MSD button clicked")
        goToSandbox()
        # Make sure all of the data has been entered
        # First check to make sure the entity resfile has some data in it
        entityOkay = False
        nEntities = 0
        if (len(self.states) == 0):
            entityOkay = False
        else:
            for i in range(1, len(self.states[0])):
                for j in range(0, len(self.states)):
                    if (len(self.states[j][i].strip()) > 0):
                        entityOkay = True
                        nEntities = nEntities + 1
                        break
        if (not(entityOkay)):
            dlg = wx.MessageDialog(self, "You have not entered any information into the entity resfile!  Please specify some correspondences.", "Entity Resfile Empty", wx.OK|wx.ICON_EXCLAMATION)
            dlg.ShowModal()
            dlg.Destroy()
            return
        # Now make sure the fitness function is validly defined
        # There should be no ? in it
        if ("?" in self.fitness):
            dlg = wx.MessageDialog(self, "The specified fitness function is invalid.  Please fill in or remove any remaining \"?\".", "Fitness Function Undefined", wx.OK|wx.ICON_EXCLAMATION)
            dlg.ShowModal()
            dlg.Destroy()
            return
        # Now make sure that no unused state variables ended up in the fitness function
        for state in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[0:len(self.states)]:
            if (not(state in self.stateAssignments) and state in self.fitness):
                dlg = wx.MessageDialog(self, "The fitness function contains a reference to " + state + " but no structures are in this state!  Remove this state from the function or put a structure in this state.", "Fitness Function Invalid", wx.OK|wx.ICON_EXCLAMATION)
                dlg.ShowModal()
                dlg.Destroy()
                return
        # All the rest of the data is either optional or forced to be specified and valid
        # Now we have to generate a ton of input files and send them to a remote server
        # This is all going to go into a large text file that will be parsed by the server
        f = open("msdinputtemp", "w")
        # First write all the PDB data
        # Find all the PDBs that have data specified in the entity resfile, otherwise leave them out
        validmodels = []
        validmodelnames = []
        for i in range(0, len(self.states)):
            modelname = self.states[i][0][0:len(self.states[i][0])-2]
            for j in range(1, len(self.states[i])):
                if (len(self.states[i][j].strip()) > 0 and not(modelname in validmodelnames)):
                    validmodels.append(i)
                    validmodelnames.append(modelname)
                    break
        # Write the PDBs
        for modelname in validmodelnames:
            #modelname = self.grdEntity.GetColLabelValue(validmodel+1) # Remember col 0 is the options list
            # Dump the PDB from PyMOL first in case the coordinates were altered by the user
            self.cmd.save(modelname + ".pdb", "model " + modelname)
            fixPyMOLSave(modelname + ".pdb")
            f2 = open(modelname + ".pdb", "r")
            if (modelname.startswith("msd_output")):
                # The server daemon expects the output files to start with this string, so change it if
                # the same string is in the inputs
                f.write("BEGIN PDB 1" + modelname + ".pdb\n")
            else:
                f.write("BEGIN PDB " + modelname + ".pdb\n")
            for aline in f2:
                f.write(aline.strip() + "\n")
            f.write("END PDB " + modelname + ".pdb\n")
            f2.close()
        # Now write out the entity resfile data
        f.write("BEGIN ENTITY RESFILE entity.resfile\n")
        f.write(str(nEntities) + "\n")
        # Default to ALLAAxc for now
        f.write("ALLAAxc EX 1 EX ARO 2\n")
        f.write("start\n")
        # Write out the specified options for each position
        ires = 1
        validentities = []
        for i in range(1, len(self.states[0])):
            # Is this line empty?
            empty = True
            for j in range(0, len(self.states)):
                if (len(self.states[j][i].strip()) > 0):
                    empty = False
                    break
            if (not(empty)):
                validentities.append((i, ires))
                f.write(str(ires) + " A PIKAA " + self.stateoptions[i-1] + " EX 1 EX ARO 2\n")
                ires = ires + 1
        f.write("END ENTITY RESFILE entity.resfile\n")
        # Now write out the correspondence files that tells MSD what the entities refers to in each of the structures
        for modelname in validmodelnames:
            #modelname = self.grdEntity.GetColLabelValue(validmodel+1)
            f.write("BEGIN CORRESPONDENCE " + modelname + ".corr\n")
            for (i, ires) in validentities:
                # Is there a correspondence here? Otherwise, skip
                # There may be multiple positions for this model on multiple chains, look for that
                for j in range(0, len(self.states)):
                    if (self.states[j][0][0:len(self.states[j][0])-2] == modelname and len(self.states[j][i].strip()) > 0):
                        f.write(str(ires) + " " + self.states[j][i][3:].strip() + " " + self.states[j][i][0] + "\n")
                #if (len(self.states[validmodel][i].strip()) > 0):
                    #f.write(str(ires) + " " + self.states[validmodel][i][3:].strip() + " " + self.states[validmodel][i][0] + "\n")
            f.write("END CORRESPONDENCE " + modelname + ".corr\n")
        # Now write out the secondary resfiles
        for validmodel in validmodels:
            ID = self.grdEntity.GetColLabelValue(validmodel+1)
            modelname = ID[0:len(ID)-2]
            f.write("BEGIN 2RESFILE " + modelname + ".2resfile\n")
            # Default to NATRO unless otherwise specified to NATAA
            f.write("USE_INPUT_SC\n")
            f.write("NATRO\n")
            f.write("start\n")
            for [a, b, seqpos, c, chainID, d, e] in self.resfile[validmodel]:
                f.write(str(seqpos) + " " + chainID + " NATAA\n")
            f.write("END 2RESFILE " + modelname + ".2resfile\n")
        # Now write out the state files that declare the grouping of PDB to correspondence files and 2resfiles
        for state in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[0:len(self.states)]:
            if (not(state in self.stateAssignments)):
                continue
            f.write("BEGIN STATES " + state + ".states\n")
            for i in range(0, self.grdStates.NumberRows):
                if (self.grdStates.GetCellValue(i, 0) != state):
                    continue
                modelname = self.grdStates.GetRowLabelValue(i)
                #skip = True
                #for validmodel in validmodels:
                    #if (self.states[validmodel][0][0: == modelname):
                        #skip = False
                        #break
                if (modelname in validmodelnames):
                    f.write(modelname + ".pdb " + modelname + ".corr " + modelname + ".2resfile\n")
            f.write("END STATES " + state + ".states\n")
        # Now write out the fitness file that will define the vector operations on the states and the fitness function
        f.write("BEGIN FITNESS fitness.daf\n")
        for state in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[0:len(self.states)]:
            if (not(state in self.stateAssignments)):
                continue
            f.write("STATE_VECTOR vec" + state + " " + state + ".states\n")
        f.write("\n")
        for i in range(0, len(self.vectorFunctions)):
            state = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[i]
            if (not(state in self.stateAssignments)):
                continue
            f.write("SCALAR_EXPRESSION " + state + " = ")
            if (self.vectorFunctions[i] == "MAX"):
                f.write("vmax( vec" + state + " )\n")
            else:
                f.write("vmin( vec" + state + " )\n")
        f.write("\n")
        f.write("FITNESS " + self.fitness + "\n")
        f.write("END FITNESS fitness.daf\n")
        # Finally, write out the flags file that will be used for submission on the server
        f.write("BEGIN FLAGS flags\n")
        f.write("-entity_resfile entity.resfile\n")
        f.write("-fitness_file fitness.daf\n")
        f.write("-ms::pop_size 100\n")
        f.write("-ms::generations " + str(15*nEntities) + "\n")
        f.write("-ms::numresults 10\n")
        f.write("-no_his_his_pairE\n")
        f.write("-ms::fraction_by_recombination .02\n")
        f.write("-msd::double_lazy_ig_mem_limit 800\n")
        f.write("END FLAGS flags\n")
        f.close()
        appendScorefxnParamsInfoToFile("msdinputtemp", self.selectWin.weightsfile)
        # Send it to the server
        try: 
            self.ID = sendToServer("msdinput")
            dlg = wx.TextEntryDialog(None, "Enter a description for this submission:", "Job Description", "")
            if (dlg.ShowModal() == wx.ID_OK):
                desc = dlg.GetValue()
                desc = desc.replace("\t", " ").replace("\n", " ").strip()
            else:
                desc = self.ID
            # Now write this to a text file
            # A timer on the SequenceWindow will use this information to look for the files to appear on the server
            goToSandbox()
            # First make sure this isn't a duplicate
            alreadythere = False
            try:
                f = open("downloadwatch", "r")
                for aline in f:
                    if (len(aline.split("\t")) >= 2 and aline.split("\t")[0] == "MSD" and aline.split("\t")[1] == self.ID.strip()):
                        alreadythere = True
                        break
                f.close()
            except:
                pass
            if (not(alreadythere)):
                f = open("downloadwatch", "a")
                f.write("MSD\t" + self.ID.strip() + "\t" + str(datetime.datetime.now().strftime("%A, %B %d - %I:%M:%S %p")) + "\t" + getServerName() + "\t" + desc + "\n")
                f.close()
            dlg = wx.MessageDialog(self, "InteractiveROSETTA is now watching the server for job ID " + desc.strip() + ".  You will be notified when the package is available for download.", "Listening for Download", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
            dlg.ShowModal()
            dlg.Destroy()
            logInfo("MSD design input sent to server daemon with ID " + self.ID)
        except Exception as e:
            logInfo("Server daemon not available")
            f = open("errreport", "w")
            f.write(e.message.strip())
            f.close()
            self.recoverFromError()
    
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
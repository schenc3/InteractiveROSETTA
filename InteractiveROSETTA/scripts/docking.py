import wx
import wx.grid
import wx.lib.scrolledpanel
import os
import os.path
import time
import platform
import math
import multiprocessing
import Bio.PDB
import webbrowser
import datetime
import gzip
import numpy
from threading import Thread
from tools import *

class SettingsDialog(wx.Dialog):
    def __init__(self, parent):
        if (platform.system() != "Linux"):
            wx.Dialog.__init__(self, parent, -1, "Constraint Default Settings", size=(330, 270))
        else:
            wx.Dialog.__init__(self, parent, -1, "Constraint Default Settings", size=(330, 270))
        
        # Structure label
        self.lblTitle = wx.StaticText(self, -1, "Constraint Default Settings", (0, 10), (320, 30))
        self.lblTitle.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
        resizeTextControlForUNIX(self.lblTitle, 0, 320)
        
        ypos = 35
        if (platform.system() == "Windows"):
            self.lblFunction = wx.StaticText(self, -1, "Function Type:", (0, ypos+23), (160, 20), wx.ALIGN_CENTRE)
            self.lblFunction.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblFunction = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/lblFunction_Dialog.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, ypos+23), size=(160, 20))
        else:
            self.lblFunction = wx.StaticText(self, -1, "Function Type:", (0, ypos+23), style=wx.ALIGN_CENTRE)
            self.lblFunction.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblFunction, 0, 160)
        self.menuFunctions = wx.ComboBox(self, pos=(160, ypos+20), size=(160, 25), choices=["Harmonic", "Gaussian", "Bounded", "Sigmoid"], style=wx.CB_READONLY)
        self.menuFunctions.SetToolTipString("List of supported constraint function types")
        self.menuFunctions.SetStringSelection(parent.dFunction)
        if (platform.system() == "Windows"):
            self.lblIdealD = wx.StaticText(self, -1, "Ideal Distance (A):", (0, ypos+53), (200, 20), wx.ALIGN_CENTRE)
            self.lblIdealD.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblIdealD = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/lblIdealD_Dialog.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, ypos+53), size=(200, 20))
        else:
            self.lblIdealD = wx.StaticText(self, -1, "Ideal Distance (A):", (0, ypos+53), style=wx.ALIGN_CENTRE)
            self.lblIdealD.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblIdealD, 0, 200)
        self.txtIdealD = wx.TextCtrl(self, -1, pos=(200, ypos+50), size=(120, 25))
        self.txtIdealD.SetValue(str(parent.dIdeal))
        self.txtIdealD.SetToolTipString("Ideal distance between the residue and its partner residue/chain")
        if (platform.system() == "Windows"):
            self.lblMaxD = wx.StaticText(self, -1, "Maximum Distance (A):", (0, ypos+83), (200, 20), wx.ALIGN_CENTRE)
            self.lblMaxD.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblMaxD = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/lblMaxD_Dialog.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, ypos+83), size=(200, 20))
        else:
            self.lblMaxD = wx.StaticText(self, -1, "Maximum Distance (A):", (0, ypos+83), style=wx.ALIGN_CENTRE)
            self.lblMaxD.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblMaxD, 0, 200)
        self.txtMaxD = wx.TextCtrl(self, -1, pos=(200, ypos+80), size=(120, 25))
        self.txtMaxD.SetValue(str(parent.dMax))
        self.txtMaxD.SetToolTipString("Maximum distance between the residue and its partner residue/chain")
        if (platform.system() == "Windows"):
            self.lblMinD = wx.StaticText(self, -1, "Minimum Distance (A):", (0, ypos+113), (200, 20), wx.ALIGN_CENTRE)
            self.lblMinD.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblMinD = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/lblMinD_Dialog.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, ypos+113), size=(200, 20))
        else:
            self.lblMinD = wx.StaticText(self, -1, "Minimum Distance (A):", (0, ypos+113), style=wx.ALIGN_CENTRE)
            self.lblMinD.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblMinD, 0, 200)
        self.txtMinD = wx.TextCtrl(self, -1, pos=(200, ypos+110), size=(120, 25))
        self.txtMinD.SetValue(str(parent.dMin))
        self.txtMinD.SetToolTipString("Minimum distance between the residue and its partner residue/chain")
        if (platform.system() == "Windows"):
            self.lblWeight = wx.StaticText(self, -1, "Weight:", (0, ypos+143), (200, 20), wx.ALIGN_CENTRE)
            self.lblWeight.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblWeight = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/lblWeight_Dialog.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, ypos+143), size=(200, 20))
        else:
            self.lblWeight = wx.StaticText(self, -1, "Weight:", (0, ypos+143), style=wx.ALIGN_CENTRE)
            self.lblWeight.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblWeight, 0, 200)
        self.txtWeight = wx.TextCtrl(self, -1, pos=(200, ypos+140), size=(120, 25))
        self.txtWeight.SetValue(str(parent.dWeight))
        self.txtWeight.SetToolTipString("Relative weight of this constraint (constraints with higher weights are favored over those with lower weights)")
        
        # OK and Cancel buttons
        self.btnOK = wx.Button(self, id=-1, label="OK", pos=(40, ypos+180), size=(100, 30))
        self.btnOK.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnOK.Bind(wx.EVT_BUTTON, self.okDialog)
        self.btnOK.SetToolTipString("Confirm current indicated default settings")
        self.btnCancel = wx.Button(self, id=-1, label="Cancel", pos=(180, ypos+180), size=(100, 30))
        self.btnCancel.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnCancel.Bind(wx.EVT_BUTTON, self.cancelDialog)
        self.btnCancel.SetToolTipString("Cancel this operation")
        
        # Center the dialog in the middle of the screen
        self.SetPosition((wx.GetDisplaySize()[0]/2-150, wx.GetDisplaySize()[1]/2-150))
            
    # Return codes, the main script knows how to interpret these
    def okDialog(self, event):
        self.SetReturnCode(wx.OK)
        self.EndModal(wx.OK)
        
    def cancelDialog(self, event):
        self.SetReturnCode(wx.CANCEL)
        self.EndModal(wx.CANCEL)

class DockingPanel(wx.lib.scrolledpanel.ScrolledPanel):
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
            self.lblProt = wx.StaticText(self, -1, "Docking", (25, 15), (270, 25), wx.ALIGN_CENTRE)
            self.lblProt.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblProt = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/lblDocking.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 15), size=(320, 25))
        else:
            self.lblProt = wx.StaticText(self, -1, "Docking", (70, 15), style=wx.ALIGN_CENTRE)
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
            self.lblInst = wx.StaticText(self, -1, "Dock proteins and/or ligands to receptor proteins", (0, 45), (320, 25), wx.ALIGN_CENTRE)
            self.lblInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
        elif (platform.system() == "Darwin"):
            self.lblInst = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/lblInstDocking.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 45), size=(320, 25))
        else:
            self.lblInst = wx.StaticText(self, -1, "Dock proteins and/or ligands to receptor proteins", (5, 45), style=wx.ALIGN_CENTRE)
            self.lblInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
            resizeTextControlForUNIX(self.lblInst, 0, self.GetSize()[0]-20)
        self.lblInst.SetForegroundColour("#FFFFFF")
        
        if (platform.system() == "Windows"):
            self.lblStatic = wx.StaticText(self, -1, "Static Chains", (0, 90), (155, 20), wx.ALIGN_CENTRE)
            self.lblStatic.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblStatic = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/lblStatic.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 90), size=(155, 20))
        else:
            self.lblStatic = wx.StaticText(self, -1, "Static Chains", (0, 90), style=wx.ALIGN_CENTRE)
            self.lblStatic.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblStatic, 0, 155)
        self.lblStatic.SetForegroundColour("#FFFFFF")
        self.staticMenu = wx.ComboBox(self, pos=(0, 110), size=(155, 25), choices=[], style=wx.CB_READONLY)
        self.staticMenu.Bind(wx.EVT_COMBOBOX, self.staticMenuSelect)
        self.staticMenu.SetToolTipString("Select chains that will be fixed receptors")
        if (platform.system() == "Windows"):
            self.lblMoving = wx.StaticText(self, -1, "Movable Chains", (165, 90), (155, 20), wx.ALIGN_CENTRE)
            self.lblMoving.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblMoving = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/lblMoving.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(165, 90), size=(155, 20))
        else:
            self.lblMoving = wx.StaticText(self, -1, "Movable Chains", (165, 90), style=wx.ALIGN_CENTRE)
            self.lblMoving.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblMoving, 165, 155)
        self.lblMoving.SetForegroundColour("#FFFFFF")
        self.movingMenu = wx.ComboBox(self, pos=(165, 110), size=(155, 25), choices=[], style=wx.CB_READONLY)
        self.movingMenu.Bind(wx.EVT_COMBOBOX, self.movingMenuSelect)
        self.movingMenu.SetToolTipString("Select chains that will be docked to the static chains")
        
        if (platform.system() == "Darwin"):
            self.btnAddStatic = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/btnAddChain.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 140), size=(77, 25))
        else:
            self.btnAddStatic = wx.Button(self, id=-1, label="Add", pos=(0, 140), size=(77, 25))
            self.btnAddStatic.SetForegroundColour("#000000")
            self.btnAddStatic.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAddStatic.Bind(wx.EVT_BUTTON, self.addStatic)
        self.btnAddStatic.SetToolTipString("Add selected chain to the list of static receptor chains")
        if (platform.system() == "Darwin"):
            self.btnRemoveStatic = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/btnRemoveChain.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(77, 140), size=(78, 25))
        else:
            self.btnRemoveStatic = wx.Button(self, id=-1, label="Remove", pos=(77, 140), size=(78, 25))
            self.btnRemoveStatic.SetForegroundColour("#000000")
            self.btnRemoveStatic.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnRemoveStatic.Bind(wx.EVT_BUTTON, self.removeStatic)
        self.btnRemoveStatic.SetToolTipString("Remove selected chain from the list of movable docked chains")
        self.staticChains = []
        
        if (platform.system() == "Darwin"):
            self.btnAddMoving = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/btnAddChain.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(165, 140), size=(77, 25))
        else:
            self.btnAddMoving = wx.Button(self, id=-1, label="Add", pos=(165, 140), size=(77, 25))
            self.btnAddMoving.SetForegroundColour("#000000")
            self.btnAddMoving.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAddMoving.Bind(wx.EVT_BUTTON, self.addMoving)
        self.btnAddMoving.SetToolTipString("Add selected chain to the list of static receptor chains")
        if (platform.system() == "Darwin"):
            self.btnRemoveMoving = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/btnRemoveChain.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(242, 140), size=(78, 25))
        else:
            self.btnRemoveMoving = wx.Button(self, id=-1, label="Remove", pos=(242, 140), size=(78, 25))
            self.btnRemoveMoving.SetForegroundColour("#000000")
            self.btnRemoveMoving.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnRemoveMoving.Bind(wx.EVT_BUTTON, self.removeMoving)
        self.btnRemoveMoving.SetToolTipString("Remove selected chain from the list of movable docked chains")
        self.movingChains = []
        
        self.grdDocking = wx.grid.Grid(self)
        self.grdDocking.CreateGrid(0, 1)
        self.grdDocking.SetSize((320, 100))
        self.grdDocking.SetPosition((0, 175))
        self.grdDocking.SetLabelFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.grdDocking.DisableDragColSize()
        self.grdDocking.DisableDragRowSize()
        self.grdDocking.SetColLabelValue(0, "Movable Chains")
        self.grdDocking.SetRowLabelSize(160)
        self.grdDocking.SetColSize(0, 160)
        
        ypos = self.grdDocking.GetPosition()[1] + self.grdDocking.GetSize()[1] + 10
        if (platform.system() == "Windows"):
            self.lblProt2 = wx.StaticText(self, -1, "Constraints (Recommended)", (0, ypos), (320, 25), wx.ALIGN_CENTRE)
            self.lblProt2.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblProt2 = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/lblConstraints.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(25, ypos), size=(270, 25))
        else:
            self.lblProt2 = wx.StaticText(self, -1, "Constraints (Recommended)", (70, ypos), style=wx.ALIGN_CENTRE)
            self.lblProt2.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
            resizeTextControlForUNIX(self.lblProt2, 0, self.GetSize()[0]-20)
        self.lblProt2.SetForegroundColour("#FFFFFF")
        
        if (platform.system() == "Windows"):
            self.lblInst2 = wx.StaticText(self, -1, "Specify which elements of the proteins interact.\n\nAtom Pair Constraint - Specify residue contacts\nSite Constraint - Residue to chain interactions\nUse the default settings unless you understand\nwhat you are doing.", (0, ypos+30), (320, 25), wx.ALIGN_CENTRE)
            self.lblInst2.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
        elif (platform.system() == "Darwin"):
            self.lblInst2 = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/lblInstConstraints.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, ypos+30), size=(320, 95))
        else:
            self.lblInst2 = wx.StaticText(self, -1, "Specify which elements of the proteins interact.\n\nAtom Pair Constraint - Specify residue contacts\nSite Constraint - Residue to chain interactions\nUse the default settings unless you understand\nwhat you are doing.", (5, ypos+30), style=wx.ALIGN_CENTRE)
            self.lblInst2.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
            resizeTextControlForUNIX(self.lblInst2, 0, self.GetSize()[0]-20)
        self.lblInst2.SetForegroundColour("#FFFFFF")
        
        if (platform.system() == "Darwin"):
            self.btnAdd = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/btnAddConstraint.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, ypos+130), size=(100, 25))
        else:
            self.btnAdd = wx.Button(self, id=-1, label="Add", pos=(0, ypos+130), size=(100, 25))
            self.btnAdd.SetForegroundColour("#000000")
            self.btnAdd.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAdd.Bind(wx.EVT_BUTTON, self.add)
        self.btnAdd.SetToolTipString("Add the selected residues to the list of constraints")
        if (platform.system() == "Darwin"):
            self.btnRemove = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/btnRemoveConstraint.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(110, ypos+130), size=(100, 25))
        else:
            self.btnRemove = wx.Button(self, id=-1, label="Remove", pos=(110, ypos+130), size=(100, 25))
            self.btnRemove.SetForegroundColour("#000000")
            self.btnRemove.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnRemove.Bind(wx.EVT_BUTTON, self.remove)
        self.btnRemove.SetToolTipString("Remove the selected residues from the list of constraints")
        if (platform.system() == "Darwin"):
            self.btnClear = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/btnClearConstraints.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(220, ypos+130), size=(100, 25))
        else:
            self.btnClear = wx.Button(self, id=-1, label="Clear", pos=(220, ypos+130), size=(100, 25))
            self.btnClear.SetForegroundColour("#000000")
            self.btnClear.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnClear.Bind(wx.EVT_BUTTON, self.clear)
        self.btnClear.SetToolTipString("Clear the list of constraints")
        
        if (platform.system() == "Darwin"):
            self.btnAtomPairConstraint = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/btnAtomPairConstraint_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, ypos+160), size=(150, 25))
        else:
            self.btnAtomPairConstraint = wx.Button(self, id=-1, label="Atom Pair", pos=(5, ypos+160), size=(150, 25))
            self.btnAtomPairConstraint.SetForegroundColour("#FF0000")
            self.btnAtomPairConstraint.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAtomPairConstraint.Bind(wx.EVT_BUTTON, self.toggleAtomPair)
        self.btnAtomPairConstraint.SetToolTipString("Constrain a contact between two residues")
        if (platform.system() == "Darwin"):
            self.btnSiteConstraint = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/btnSiteConstraint.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(165, ypos+160), size=(150, 25))
        else:
            self.btnSiteConstraint = wx.Button(self, id=-1, label="Site", pos=(165, ypos+160), size=(150, 25))
            self.btnSiteConstraint.SetForegroundColour("#000000")
            self.btnSiteConstraint.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnSiteConstraint.Bind(wx.EVT_BUTTON, self.toggleSite)
        self.btnSiteConstraint.SetToolTipString("Constrain a contact between a residue and a chain")
        self.constraintType = "AtomPair"
        self.menuPartners = wx.ComboBox(self, pos=(0, ypos+190), size=(320, 25), choices=[], style=wx.CB_READONLY)
        self.menuPartners.Bind(wx.EVT_COMBOBOX, self.partnerMenuSelect)
        self.menuPartners.SetToolTipString("List of potential partners for the selected residue in the constraints list")
        
        self.grdConstraints = wx.grid.Grid(self)
        self.grdConstraints.CreateGrid(0, 7)
        self.grdConstraints.SetSize((320, 200))
        self.grdConstraints.SetPosition((0, ypos+220))
        self.grdConstraints.SetLabelFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.grdConstraints.DisableDragColSize()
        self.grdConstraints.DisableDragRowSize()
        self.grdConstraints.SetColLabelValue(0, "Group")
        self.grdConstraints.SetColLabelValue(1, "Residue")
        self.grdConstraints.SetColLabelValue(2, "Model")
        self.grdConstraints.SetColLabelValue(3, "Partner")
        self.grdConstraints.SetColLabelValue(4, "Model")
        self.grdConstraints.SetColLabelValue(5, "Function")
        self.grdConstraints.SetColLabelValue(6, "Arguments")
        self.grdConstraints.SetRowLabelSize(80)
        self.grdConstraints.SetColSize(0, 50)
        self.grdConstraints.SetColSize(1, 100)
        self.grdConstraints.SetColSize(2, 100)
        self.grdConstraints.SetColSize(3, 100)
        self.grdConstraints.SetColSize(4, 100)
        self.grdConstraints.SetColSize(5, 100)
        self.grdConstraints.SetColSize(6, 200)
        self.grdConstraints.Bind(wx.grid.EVT_GRID_CELL_CHANGE, self.gridChange)
        self.grdConstraints.Bind(wx.grid.EVT_GRID_CELL_LEFT_CLICK, self.gridClick)
        self.grdConstraints.Bind(wx.grid.EVT_GRID_CELL_RIGHT_CLICK, self.gridRClick)
        self.constraints = []
        self.selectedr = -1
        
        ypos = self.grdConstraints.GetPosition()[1] + self.grdConstraints.GetSize()[1] + 10
        if (platform.system() == "Darwin"):
            self.btnLoadCST = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/btnLoadCST.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(20, ypos), size=(120, 25))
        else:
            self.btnLoadCST = wx.Button(self, id=-1, label="Load CST", pos=(20, ypos), size=(120, 25))
            self.btnLoadCST.SetForegroundColour("#000000")
            self.btnLoadCST.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.btnLoadCST.Bind(wx.EVT_BUTTON, self.loadCST)
        self.btnLoadCST.SetToolTipString("Load the data in a premade constraints file")
        if (platform.system() == "Darwin"):
            self.btnSaveCST = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/btnSaveCST.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(175, ypos), size=(120, 25))
        else:
            self.btnSaveCST = wx.Button(self, id=-1, label="Save CST", pos=(175, ypos), size=(120, 25))
            self.btnSaveCST.SetForegroundColour("#000000")
            self.btnSaveCST.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.btnSaveCST.Bind(wx.EVT_BUTTON, self.saveCST)
        self.btnSaveCST.SetToolTipString("Save the current constraints data to a real Rosetta constraints file")
        ypos = self.btnSaveCST.GetPosition()[1] + self.btnSaveCST.GetSize()[1] + 10
        
        self.btnDefaults = wx.BitmapButton(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/wrench.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, ypos-5), size=(25, 25))
        self.btnDefaults.Bind(wx.EVT_BUTTON, self.configureDefaults)
        self.btnDefaults.SetToolTipString("Configure the default settings for docking constraints")
        self.dFunction = "Bounded"
        self.dMax = 12.0
        self.dMin = 6.0
        self.dIdeal = 9.0
        self.dWeight = 1.0
        # Try to read the default settings from the cfg file
        goToSandbox()
        try:
            fin = open("docking.cfg", "r")
            for aline in fin:
                if (len(aline.strip()) == 0):
                    continue
                if (aline.startswith("[FUNCTION]")):
                    self.dFunction = aline.split("\t")[1].strip()
                elif (aline.startswith("[IDEAL]")):
                    self.dIdeal = float(aline.split("\t")[1].strip())
                elif (aline.startswith("[MAX]")):
                    self.dMax = float(aline.split("\t")[1].strip())
                elif (aline.startswith("[MIN]")):
                    self.dMin = float(aline.split("\t")[1].strip())
                elif (aline.startswith("[WEIGHT]")):
                    self.dWeight = float(aline.split("\t")[1].strip())
            fin.close()
        except:
            pass
        if (platform.system() == "Windows"):
            self.lblAdvanced = wx.StaticText(self, -1, "Advanced Options", (0, ypos), (320, 20), wx.ALIGN_CENTRE)
            self.lblAdvanced.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblAdvanced = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/lblAdvanced.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, ypos), size=(320, 20))
        else:
            self.lblAdvanced = wx.StaticText(self, -1, "Advanced Options", (0, ypos), style=wx.ALIGN_CENTRE)
            self.lblAdvanced.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblAdvanced, 0, 320)
        self.lblAdvanced.SetForegroundColour("#FFFFFF")
        if (platform.system() == "Windows"):
            self.lblFunction = wx.StaticText(self, -1, "Function Type:", (0, ypos+23), (160, 20), wx.ALIGN_CENTRE)
            self.lblFunction.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblFunction = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/lblFunction.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, ypos+23), size=(160, 20))
        else:
            self.lblFunction = wx.StaticText(self, -1, "Function Type:", (0, ypos+23), style=wx.ALIGN_CENTRE)
            self.lblFunction.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblFunction, 0, 160)
        self.lblFunction.SetForegroundColour("#FFFFFF")
        self.menuFunctions = wx.ComboBox(self, pos=(160, ypos+20), size=(160, 25), choices=["Harmonic", "Gaussian", "Bounded", "Sigmoid"], style=wx.CB_READONLY)
        self.menuFunctions.Bind(wx.EVT_COMBOBOX, self.functionMenuSelect)
        self.menuFunctions.SetToolTipString("List of supported constraint function types")
        self.menuFunctions.Disable()
        if (platform.system() == "Windows"):
            self.lblIdealD = wx.StaticText(self, -1, "Ideal Distance (A):", (0, ypos+53), (200, 20), wx.ALIGN_CENTRE)
            self.lblIdealD.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblIdealD = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/lblIdealD.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, ypos+53), size=(200, 20))
        else:
            self.lblIdealD = wx.StaticText(self, -1, "Ideal Distance (A):", (0, ypos+53), style=wx.ALIGN_CENTRE)
            self.lblIdealD.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblIdealD, 0, 200)
        self.lblIdealD.SetForegroundColour("#FFFFFF")
        self.txtIdealD = wx.TextCtrl(self, -1, pos=(200, ypos+50), size=(120, 25))
        self.txtIdealD.SetValue("")
        self.txtIdealD.Disable()
        self.txtIdealD.Bind(wx.EVT_TEXT, self.advancedTextUpdate)
        self.txtIdealD.SetToolTipString("Ideal distance between the residue and its partner residue/chain")
        if (platform.system() == "Windows"):
            self.lblMaxD = wx.StaticText(self, -1, "Maximum Distance (A):", (0, ypos+83), (200, 20), wx.ALIGN_CENTRE)
            self.lblMaxD.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblMaxD = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/lblMaxD.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, ypos+83), size=(200, 20))
        else:
            self.lblMaxD = wx.StaticText(self, -1, "Maximum Distance (A):", (0, ypos+83), style=wx.ALIGN_CENTRE)
            self.lblMaxD.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblMaxD, 0, 200)
        self.lblMaxD.SetForegroundColour("#FFFFFF")
        self.txtMaxD = wx.TextCtrl(self, -1, pos=(200, ypos+80), size=(120, 25))
        self.txtMaxD.SetValue("")
        self.txtMaxD.Disable()
        self.txtMaxD.Bind(wx.EVT_TEXT, self.advancedTextUpdate)
        self.txtMaxD.SetToolTipString("Maximum distance between the residue and its partner residue/chain")
        if (platform.system() == "Windows"):
            self.lblMinD = wx.StaticText(self, -1, "Minimum Distance (A):", (0, ypos+113), (200, 20), wx.ALIGN_CENTRE)
            self.lblMinD.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblMinD = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/lblMinD.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, ypos+113), size=(200, 20))
        else:
            self.lblMinD = wx.StaticText(self, -1, "Minimum Distance (A):", (0, ypos+113), style=wx.ALIGN_CENTRE)
            self.lblMinD.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblMinD, 0, 200)
        self.lblMinD.SetForegroundColour("#FFFFFF")
        self.txtMinD = wx.TextCtrl(self, -1, pos=(200, ypos+110), size=(120, 25))
        self.txtMinD.SetValue("")
        self.txtMinD.Disable()
        self.txtMinD.Bind(wx.EVT_TEXT, self.advancedTextUpdate)
        self.txtMinD.SetToolTipString("Minimum distance between the residue and its partner residue/chain")
        if (platform.system() == "Windows"):
            self.lblWeight = wx.StaticText(self, -1, "Weight:", (0, ypos+143), (200, 20), wx.ALIGN_CENTRE)
            self.lblWeight.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblWeight = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/lblWeight.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, ypos+143), size=(200, 20))
        else:
            self.lblWeight = wx.StaticText(self, -1, "Weight:", (0, ypos+143), style=wx.ALIGN_CENTRE)
            self.lblWeight.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblWeight, 0, 200)
        self.lblWeight.SetForegroundColour("#FFFFFF")
        self.txtWeight = wx.TextCtrl(self, -1, pos=(200, ypos+140), size=(120, 25))
        self.txtWeight.SetValue("")
        self.txtWeight.Disable()
        self.txtWeight.Bind(wx.EVT_TEXT, self.advancedTextUpdate)
        self.txtWeight.SetToolTipString("Relative weight of this constraint (constraints with higher weights are favored over those with lower weights)")
        if (platform.system() == "Windows"):
            self.lblResidueAtom = wx.StaticText(self, -1, "Residue Atom", (0, ypos+170), (160, 20), wx.ALIGN_CENTRE)
            self.lblResidueAtom.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblResidueAtom = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/lblResidueAtom.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, ypos+170), size=(160, 20))
        else:
            self.lblResidueAtom = wx.StaticText(self, -1, "Residue Atom", (0, ypos+170), style=wx.ALIGN_CENTRE)
            self.lblResidueAtom.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblResidueAtom, 0, 160)
        self.lblResidueAtom.SetForegroundColour("#FFFFFF")
        if (platform.system() == "Windows"):
            self.lblPartnerAtom = wx.StaticText(self, -1, "Partner Atom", (160, ypos+170), (160, 20), wx.ALIGN_CENTRE)
            self.lblPartnerAtom.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblPartnerAtom = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/lblPartnerAtom.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(160, ypos+170), size=(160, 20))
        else:
            self.lblPartnerAtom = wx.StaticText(self, -1, "Partner Atom", (160, ypos+170), style=wx.ALIGN_CENTRE)
            self.lblPartnerAtom.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblPartnerAtom, 160, 160)
        self.lblPartnerAtom.SetForegroundColour("#FFFFFF")
        self.menuResidueAtoms = wx.ComboBox(self, pos=(5, ypos+190), size=(150, 25), choices=[], style=wx.CB_READONLY)
        self.menuResidueAtoms.Bind(wx.EVT_COMBOBOX, self.residueAtomMenuSelect)
        self.menuResidueAtoms.SetToolTipString("List of atoms for the selected residue")
        self.menuResidueAtoms.Disable()
        self.menuPartnerAtoms = wx.ComboBox(self, pos=(165, ypos+190), size=(150, 25), choices=[], style=wx.CB_READONLY)
        self.menuPartnerAtoms.Bind(wx.EVT_COMBOBOX, self.partnerAtomMenuSelect)
        self.menuPartnerAtoms.SetToolTipString("List of atoms for the partner residue")
        self.menuPartnerAtoms.Disable()
        
        if (platform.system() == "Windows"):
            self.lblLine = wx.StaticText(self, -1, "==========================", (0, ypos+220), (320, 20), wx.ALIGN_CENTRE)
            self.lblLine.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.NORMAL))
        elif (platform.system() == "Darwin"):
            self.lblLine = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/lblLine.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, ypos+220), size=(320, 20))
        else:
            self.lblLine = wx.StaticText(self, -1, "==========================", (0, ypos+220), style=wx.ALIGN_CENTRE)
            self.lblLine.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.NORMAL))
            resizeTextControlForUNIX(self.lblLine, 20, 120)
        self.lblLine.SetForegroundColour("#FFFFFF")
        
        ypos = self.lblLine.GetPosition()[1] + self.lblLine.GetSize()[1] + 10
        if (platform.system() == "Windows"):
            self.lblProt3 = wx.StaticText(self, -1, "Ensemble Docking (Optional)", (0, ypos), (320, 25), wx.ALIGN_CENTRE)
            self.lblProt3.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblProt3 = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/lblEnsembleDocking.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, ypos), size=(320, 25))
        else:
            self.lblProt3 = wx.StaticText(self, -1, "Ensemble Docking (Optional)", (0, ypos), style=wx.ALIGN_CENTRE)
            self.lblProt3.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
            resizeTextControlForUNIX(self.lblProt3, 0, self.GetSize()[0]-20)
        self.lblProt3.SetForegroundColour("#FFFFFF")
        
        if (platform.system() == "Windows"):
            self.lblInst3 = wx.StaticText(self, -1, "Rosetta can represent the static and/or moving\nmodels as ensembles rather than single models.\nProvide an ensemble archive (.ensb) as input for\neither structure to activate ensemble docking.", (0, ypos+30), (320, 25), wx.ALIGN_CENTRE)
            self.lblInst3.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
        elif (platform.system() == "Darwin"):
            self.lblInst3 = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/lblInstEnsembleDocking.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, ypos+30), size=(320, 75))
        else:
            self.lblInst3 = wx.StaticText(self, -1, "Rosetta can represent the static and/or moving\nmodels as ensembles rather than single models.\nProvide an ensemble archive (.ensb) as input for\neither structure to activate ensemble docking.", (5, ypos+30), style=wx.ALIGN_CENTRE)
            self.lblInst3.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
            resizeTextControlForUNIX(self.lblInst3, 0, self.GetSize()[0]-20)
        self.lblInst3.SetForegroundColour("#FFFFFF")
        
        if (platform.system() == "Windows"):
            self.lblStaticEnsb = wx.StaticText(self, -1, "Static Ensemble:", (0, ypos+103), (160, 20), wx.ALIGN_CENTRE)
            self.lblStaticEnsb.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblStaticEnsb = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/lblStaticEnsb.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, ypos+103), size=(160, 20))
        else:
            self.lblStaticEnsb = wx.StaticText(self, -1, "Static Ensemble:", (0, ypos+103), style=wx.ALIGN_CENTRE)
            self.lblStaticEnsb.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblStaticEnsb, 0, 160)
        self.lblStaticEnsb.SetForegroundColour("#FFFFFF")
        if (platform.system() == "Darwin"):
            self.btnLoadStaticEnsb = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/btnLoadStaticEnsb.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(170, ypos+100), size=(70, 25))
        else:
            self.btnLoadStaticEnsb = wx.Button(self, id=-1, label="Load", pos=(170, ypos+100), size=(70, 25))
            self.btnLoadStaticEnsb.SetForegroundColour("#000000")
            self.btnLoadStaticEnsb.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnLoadStaticEnsb.Bind(wx.EVT_BUTTON, self.loadStaticEnsb)
        self.btnLoadStaticEnsb.SetToolTipString("Load an ensemble archive for the static chains")
        if (platform.system() == "Darwin"):
            self.btnDeleteStaticEnsb = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/btnDeleteStaticEnsb.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(250, ypos+100), size=(70, 25))
        else:
            self.btnDeleteStaticEnsb = wx.Button(self, id=-1, label="Delete", pos=(250, ypos+100), size=(70, 25))
            self.btnDeleteStaticEnsb.SetForegroundColour("#000000")
            self.btnDeleteStaticEnsb.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnDeleteStaticEnsb.Bind(wx.EVT_BUTTON, self.deleteStaticEnsb)
        self.btnDeleteStaticEnsb.SetToolTipString("Delete the loaded ensemble archive for the static chains")
        if (platform.system() == "Windows"):
            self.lblSelStaticEnsb = wx.StaticText(self, -1, "No Ensemble Specified", (0, ypos+133), (320, 20), wx.ALIGN_CENTRE)
            self.lblSelStaticEnsb.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        else:
            self.lblSelStaticEnsb = wx.StaticText(self, -1, "No Ensemble Specified", (0, ypos+133), style=wx.ALIGN_CENTRE)
            self.lblSelStaticEnsb.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblSelStaticEnsb, 0, 320)
        self.lblSelStaticEnsb.SetForegroundColour("#FFFFFF")
        
        if (platform.system() == "Windows"):
            self.lblMovingEnsb = wx.StaticText(self, -1, "Moving Ensemble:", (0, ypos+163), (160, 20), wx.ALIGN_CENTRE)
            self.lblMovingEnsb.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblMovingEnsb = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/lblMovingEnsb.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, ypos+163), size=(160, 20))
        else:
            self.lblMovingEnsb = wx.StaticText(self, -1, "Moving Ensemble:", (0, ypos+163), style=wx.ALIGN_CENTRE)
            self.lblMovingEnsb.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblMovingEnsb, 0, 160)
        self.lblMovingEnsb.SetForegroundColour("#FFFFFF")
        if (platform.system() == "Darwin"):
            self.btnLoadMovingEnsb = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/btnLoadStaticEnsb.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(170, ypos+160), size=(70, 25))
        else:
            self.btnLoadMovingEnsb = wx.Button(self, id=-1, label="Load", pos=(170, ypos+160), size=(70, 25))
            self.btnLoadMovingEnsb.SetForegroundColour("#000000")
            self.btnLoadMovingEnsb.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnLoadMovingEnsb.Bind(wx.EVT_BUTTON, self.loadMovingEnsb)
        self.btnLoadMovingEnsb.SetToolTipString("Load an ensemble archive for the moving chains")
        if (platform.system() == "Darwin"):
            self.btnDeleteMovingEnsb = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/btnDeleteStaticEnsb.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(250, ypos+160), size=(70, 25))
        else:
            self.btnDeleteMovingEnsb = wx.Button(self, id=-1, label="Delete", pos=(250, ypos+160), size=(70, 25))
            self.btnDeleteMovingEnsb.SetForegroundColour("#000000")
            self.btnDeleteMovingEnsb.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnDeleteMovingEnsb.Bind(wx.EVT_BUTTON, self.deleteMovingEnsb)
        self.btnDeleteMovingEnsb.SetToolTipString("Delete the loaded ensemble archive for the moving chains")
        if (platform.system() == "Windows"):
            self.lblSelMovingEnsb = wx.StaticText(self, -1, "No Ensemble Specified", (0, ypos+193), (320, 20), wx.ALIGN_CENTRE)
            self.lblSelMovingEnsb.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        else:
            self.lblSelMovingEnsb = wx.StaticText(self, -1, "No Ensemble Specified", (0, ypos+193), style=wx.ALIGN_CENTRE)
            self.lblSelMovingEnsb.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblSelMovingEnsb, 0, 320)
        self.lblSelMovingEnsb.SetForegroundColour("#FFFFFF")
        self.ensemble1 = None
        self.ensemble2 = None
        
        if (platform.system() == "Windows"):
            self.lblLine2 = wx.StaticText(self, -1, "==========================", (0, ypos+220), (320, 20), wx.ALIGN_CENTRE)
            self.lblLine2.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.NORMAL))
        elif (platform.system() == "Darwin"):
            self.lblLine2 = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/lblLine.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, ypos+220), size=(320, 20))
        else:
            self.lblLine2 = wx.StaticText(self, -1, "==========================", (0, ypos+220), style=wx.ALIGN_CENTRE)
            self.lblLine2.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.NORMAL))
            resizeTextControlForUNIX(self.lblLine2, 20, 120)
        self.lblLine2.SetForegroundColour("#FFFFFF")
        
        if (platform.system() == "Windows"):
            self.lblInst4 = wx.StaticText(self, -1, "If you know the interface on the static chains,\nselect it and reorient the interface to point\ntowards the ligand.", (0, ypos+240), (320, 25), wx.ALIGN_CENTRE)
            self.lblInst4.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
        elif (platform.system() == "Darwin"):
            self.lblInst4 = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/lblInstInterface.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, ypos+240), size=(320, 60))
        else:
            self.lblInst4 = wx.StaticText(self, -1, "If you know the interface on the static chains,\nselect it and reorient the interface to point\ntowards the ligand.", (5, ypos+240), style=wx.ALIGN_CENTRE)
            self.lblInst4.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
            resizeTextControlForUNIX(self.lblInst4, 0, self.GetSize()[0]-20)
        self.lblInst4.SetForegroundColour("#FFFFFF")
        
        if (platform.system() == "Darwin"):
            self.btnReorient = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/btnReorient.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(60, ypos+290), size=(200, 25))
        else:
            self.btnReorient = wx.Button(self, id=-1, label="Re-orient Partners", pos=(60, ypos+290), size=(200, 25))
            self.btnReorient.SetForegroundColour("#000000")
            self.btnReorient.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnReorient.Bind(wx.EVT_BUTTON, self.reorient)
        self.btnReorient.SetToolTipString("Reorient the static chain selection to point towards the ligand.  Then use Fix Stat for the docking mode.")
        
        if (platform.system() == "Windows"):
            self.lblCoarse = wx.StaticText(self, -1, "Coarse Models", (0, ypos+320), (155, 20), wx.ALIGN_CENTRE)
            self.lblCoarse.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblCoarse = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/lblCoarse.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, ypos+320), size=(155, 20))
        else:
            self.lblCoarse = wx.StaticText(self, -1, "Coarse Models", (0, ypos+320), style=wx.ALIGN_CENTRE)
            self.lblCoarse.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblCoarse, 0, 155)
        self.lblCoarse.SetForegroundColour("#FFFFFF")
        self.txtCoarse = wx.TextCtrl(self, -1, pos=(0, ypos+340), size=(155, 25))
        self.txtCoarse.SetValue("1000")
        self.txtCoarse.SetToolTipString("Number of decoys to generate in the coarse docking simulation")
        if (platform.system() == "Windows"):
            self.lblRefined = wx.StaticText(self, -1, "Refined Models", (165, ypos+320), (155, 20), wx.ALIGN_CENTRE)
            self.lblRefined.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblRefined = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/lblRefined.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(165, ypos+320), size=(155, 20))
        else:
            self.lblRefined = wx.StaticText(self, -1, "Refined Models", (165, ypos+320), style=wx.ALIGN_CENTRE)
            self.lblRefined.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblRefined, 165, 155)
        self.lblRefined.SetForegroundColour("#FFFFFF")
        self.txtRefined = wx.TextCtrl(self, -1, pos=(165, ypos+340), size=(155, 25))
        self.txtRefined.SetValue("10")
        self.txtRefined.SetToolTipString("Number of decoys to generate in the refined docking simulation that you will be able to view")
        
        if (platform.system() == "Windows"):
            self.lblPostDock = wx.StaticText(self, -1, "Post-Docking", (0, ypos+370), (320, 20), wx.ALIGN_CENTRE)
            self.lblPostDock.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblPostDock = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/lblPostDock.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, ypos+370), size=(320, 20))
        else:
            self.lblPostDock = wx.StaticText(self, -1, "Post-Docking", (0, ypos+370), style=wx.ALIGN_CENTRE)
            self.lblPostDock.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
            resizeTextControlForUNIX(self.lblPostDock, 0, self.GetSize()[0]-20)
        self.lblPostDock.SetForegroundColour("#FFFFFF")
        
        if (platform.system() == "Windows"):
            self.lblModelView = wx.StaticText(self, -1, "View Structures:", (20, ypos+398), (120, 20), wx.ALIGN_CENTRE)
            self.lblModelView.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblModelView = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/lblModelView.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(20, ypos+398), size=(120, 20))
        else:
            self.lblModelView = wx.StaticText(self, -1, "View Structures:", (20, ypos+398), style=wx.ALIGN_CENTRE)
            self.lblModelView.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblModelView, 20, 120)
        self.lblModelView.SetForegroundColour("#FFFFFF")
        self.viewMenu = wx.ComboBox(self, pos=(175, ypos+395), size=(120, 25), choices=[], style=wx.CB_READONLY)
        self.viewMenu.Bind(wx.EVT_COMBOBOX, self.viewMenuSelect)
        self.viewMenu.Disable()
        self.viewMenu.SetToolTipString("Select docked positions to view in PyMOL")
        
        if (platform.system() == "Darwin"):
            self.btnServerToggle = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/btnServer_Off.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, ypos+430), size=(100, 25))
        else:
            self.btnServerToggle = wx.Button(self, id=-1, label="Server Off", pos=(0, ypos+430), size=(100, 25))
            self.btnServerToggle.SetForegroundColour("#000000")
            self.btnServerToggle.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.btnServerToggle.Bind(wx.EVT_BUTTON, self.serverToggle)
        self.btnServerToggle.SetToolTipString("Perform docking simulations locally")
        self.serverOn = False
        if (platform.system() == "Darwin"):
            self.btnStarting = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/btnStarting_Global.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(110, ypos+430), size=(100, 25))
        else:
            self.btnStarting = wx.Button(self, id=-1, label="Global", pos=(110, ypos+430), size=(100, 25))
            self.btnStarting.SetForegroundColour("#000000")
            self.btnStarting.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.btnStarting.Bind(wx.EVT_BUTTON, self.toggleStarting)
        self.btnStarting.SetToolTipString("Perform a global dock, where the orientations of both partners is unknown")
        self.startingType = "Global"
        if (platform.system() == "Darwin"):
            self.btnDock = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/btnDock.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(220, ypos+430), size=(100, 25))
        else:
            self.btnDock = wx.Button(self, id=-1, label="Dock!", pos=(220, ypos+430), size=(100, 25))
            self.btnDock.SetForegroundColour("#000000")
            self.btnDock.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.btnDock.Bind(wx.EVT_BUTTON, self.dockClick)
        self.btnDock.SetToolTipString("Begin docking simulation with selected parameters")
        self.buttonState = "Dock!"
        
        self.scrollh = self.btnDock.GetPosition()[1] + self.btnDock.GetSize()[1] + 5
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
            browser.open(self.parent.parent.scriptdir + "/help/docking.html")
        else:
            webbrowser.open(self.parent.parent.scriptdir + "/help/docking.html")
    
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
        
    def scrolled(self, event):
        self.winscrollpos = self.GetScrollPos(wx.VERTICAL)
        event.Skip()
        
    def activate(self):
        # Get the list of all the chains in the sequence viewer
        chainList = []
        for r in range(0, self.seqWin.SeqViewer.NumberRows):
            chainList.append(self.seqWin.IDs[r])
        # Update the combobox list if the list has changed
        if (chainList != self.staticMenu.GetItems()):
            self.staticMenu.Clear()
            self.staticMenu.AppendItems(chainList)
            # Take out invalid chains if the list of models changed
            for i in range(len(self.staticChains)-1, -1, -1):
                if (self.staticChains[i] not in chainList):
                    self.staticChains.pop(i)
            self.updateGrid()
        if (chainList != self.movingMenu.GetItems()):
            self.movingMenu.Clear()
            self.movingMenu.AppendItems(chainList)
            # Take out invalid chains if the list of models changed
            for i in range(len(self.movingChains)-1, -1, -1):
                if (self.movingChains[i] not in chainList):
                    self.movingChains.pop(i)
            self.updateGrid()
        # Grab the current selection of residues for adding constraints
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
        self.pruneConstraints()
        self.Scroll(0, self.winscrollpos)
    
    def configureDefaults(self, event):
        dlg = SettingsDialog(self)
        if (dlg.ShowModal() == wx.OK):
            # Update the default settings, ignoring invalid values
            self.dFunction = dlg.menuFunctions.GetStringSelection()
            try:
                val = float(dlg.txtIdealD.GetValue())
                if (val <= 0):
                    raise Exception()
                self.dIdeal = val
            except:
                pass
            try:
                val = float(dlg.txtMaxD.GetValue())
                if (val <= 0):
                    raise Exception()
                self.dMax = val
            except:
                pass
            try:
                val = float(dlg.txtMinD.GetValue())
                if (val <= 0):
                    raise Exception()
                self.dMin = val
            except:
                pass
            try:
                val = float(dlg.txtWeight.GetValue())
                if (val <= 0):
                    raise Exception()
                self.dWeight = val
            except:
                pass
            # Save these settings
            goToSandbox()
            fout = open("docking.cfg", "w")
            fout.write("[FUNCTION]\t" + self.dFunction + "\n")
            fout.write("[IDEAL]\t" + str(self.dIdeal) + "\n")
            fout.write("[MAX]\t" + str(self.dMax) + "\n")
            fout.write("[MIN]\t" + str(self.dMin) + "\n")
            fout.write("[WEIGHT]\t" + str(self.dWeight) + "\n")
            fout.close()
        dlg.Destroy()
    
    def showChain(self, ID):
        model = ID[0:len(ID)-2]
        chain = ID[len(ID)-1]
        defaultPyMOLView(self.cmd)
        if (chain != "_"):
            self.cmd.select("chainsele", "model " + model + " and chain " + chain)
        else:
            self.cmd.select("chainsele", "model " + model)
        self.cmd.zoom("chainsele")
        self.cmd.hide("everything", "not chainsele")
        self.cmd.delete("chainsele")
        self.cmd.enable("seqsele")
        self.seqWin.selectUpdate(False)
    
    def staticMenuSelect(self, event):
        logInfo("Selected static model " + self.staticMenu.GetStringSelection())
        #self.showChain(self.staticMenu.GetStringSelection())
        
    def movingMenuSelect(self, event):
        logInfo("Selected movable model " + self.movingMenu.GetStringSelection())
        #self.showChain(self.movingMenu.GetStringSelection())
    
    def updateGrid(self):
        if (self.grdDocking.NumberRows > 0):
            self.grdDocking.DeleteRows(0, self.grdDocking.NumberRows)
        for i in range(0, max(len(self.staticChains), len(self.movingChains))):
            self.grdDocking.AppendRows(1)
            readOnly = wx.grid.GridCellAttr()
            readOnly.SetReadOnly(True)
            self.grdDocking.SetRowAttr(i, readOnly)
            self.grdDocking.SetCellAlignment(i, 0, wx.ALIGN_CENTRE, wx.ALIGN_CENTRE)
            self.grdDocking.SetRowLabelValue(i, "")
        for i in range(0, len(self.staticChains)):
            self.grdDocking.SetRowLabelValue(i, self.staticChains[i])
        for i in range(0, len(self.movingChains)):
            self.grdDocking.SetCellValue(i, 0, self.movingChains[i])
    
    def addStatic(self, event):
        # Add this to the list of chains that will be fixed
        # If this chain is already in the movable list then we need to take it out of that list and put it as static
        chain = str(self.staticMenu.GetStringSelection())
        if (len(chain.strip()) == 0):
            return
        if (not(chain in self.staticChains)):
            self.staticChains.append(chain)
        self.staticChains.sort()
        if (chain in self.movingChains):
            indx = self.movingChains.index(chain)
            self.movingChains.pop(indx)
        # Update the grid with this new information
        self.updateGrid()
        self.pruneConstraints()
        self.deleteStaticEnsb(None) # Ensemble data is not valid anymore, remove it
        logInfo("Added " + chain + " to the list of static chains")
        
    def removeStatic(self, event):
        # Just pop out the indicated chain
        chain = self.staticMenu.GetStringSelection()
        if (chain in self.staticChains):
            indx = self.staticChains.index(chain)
            self.staticChains.pop(indx)
        # Update the grid with this new information
        self.updateGrid()
        self.pruneConstraints()
        self.deleteStaticEnsb(None) # Ensemble data is not valid anymore, remove it
        logInfo("Removed " + chain + " from the list of static chains")
            
    def addMoving(self, event):
        # Add this to the list of chains that will be docked
        # If this chain is already in the static list then we need to take it out of that list and put it as movable
        chain = str(self.movingMenu.GetStringSelection())
        if (len(chain.strip()) == 0):
            return
        if (not(chain in self.movingChains)):
            self.movingChains.append(chain)
        self.movingChains.sort()
        if (chain in self.staticChains):
            indx = self.staticChains.index(chain)
            self.staticChains.pop(indx)
        # Update the grid with this new information
        self.updateGrid()
        self.pruneConstraints()
        self.deleteMovingEnsb(None) # Ensemble data is not valid anymore, remove it
        logInfo("Added " + chain + " to the list of moving chains")
    
    def removeMoving(self, event):
        # Just pop out the indicated chain
        chain = self.movingMenu.GetStringSelection()
        if (chain in self.movingChains):
            indx = self.movingChains.index(chain)
            self.movingChains.pop(indx)
        # Update the grid with this new information
        self.updateGrid()
        self.pruneConstraints()
        self.deleteMovingEnsb(None) # Ensemble data is not valid anymore, remove it
        logInfo("Removed " + chain + " from the list of static chains")
        
    def toggleAtomPair(self, event):
        if (platform.system() == "Darwin"):
            self.btnAtomPairConstraint.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/btnAtomPairConstraint_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            self.btnSiteConstraint.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/btnSiteConstraint.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
        else:
            self.btnAtomPairConstraint.SetForegroundColour("#FF0000")
            self.btnSiteConstraint.SetForegroundColour("#000000")
        self.constraintType = "AtomPair"
        logInfo("Set constraint type to AtomPair")
        
    def toggleSite(self, event):
        if (platform.system() == "Darwin"):
            self.btnAtomPairConstraint.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/btnAtomPairConstraint.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            self.btnSiteConstraint.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/btnSiteConstraint_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
        else:
            self.btnSiteConstraint.SetForegroundColour("#FF0000")
            self.btnAtomPairConstraint.SetForegroundColour("#000000")
        self.constraintType = "Site"
        logInfo("Set constraint type to Site")
        
    def partnerMenuSelect(self, event):
        # Set this selection as the partner for the selected row in the constraints grid
        if (self.grdConstraints.GetRowLabelValue(self.selectedr) == "AtomPair"):
            ID = self.menuPartners.GetStringSelection().split(":")[1]
            model = ID[2:]
            chain = ID[0]
            ID = self.menuPartners.GetStringSelection().split(":")[0]
            seqpos = ID[1:]
            resn = ID[0]
            if (resn == "G"):
                atom = "CA"
            else:
                atom = "CB"
            self.constraints[self.selectedr][2] = model + "|" + chain + ":" + seqpos + resn + ":" + atom
            self.constraintView(seqpos+resn, model+"|"+chain)
            self.populatePartnerAtoms()
        else:
            ID = self.menuPartners.GetStringSelection()
            self.constraints[self.selectedr][2] = ID
            self.showChain(ID)
        self.updateConstraints()

    def pruneConstraints(self):
        # This function is called whenever a change is made to the docking setup
        # It will take out entires that are no longer valid because chains have moved around
        updateNeeded = False
        for i in range(len(self.constraints)-1, -1, -1):
            if (self.grdConstraints.GetRowLabelValue(i) == "AtomPair"):
                # First make sure the chains are still in the right categories
                if (not(self.grdConstraints.GetCellValue(i, 2) in self.staticChains)):
                    self.constraints.pop(i)
                    continue
                elif (len(self.grdConstraints.GetCellValue(i, 4).strip()) > 0 and not(self.grdConstraints.GetCellValue(i, 4) in self.movingChains)):
                    self.constraints.pop(i)
                    continue
                # Now make sure the residues still exist
                model = self.grdConstraints.GetCellValue(i, 2)
                chain = model[len(model)-1]
                model = model[0:len(model)-2]
                seqpos = self.grdConstraints.GetCellValue(i, 1).split(":")[0]
                seqpos = int(seqpos[0:len(seqpos)-1])
                if (self.seqWin.doesResidueExist(model, chain, seqpos)):
                    pass
                else:
                    self.constraints.pop(i)
                    updateNeeded = True
                    continue
                if (len(self.grdConstraints.GetCellValue(i, 3).strip()) == 0):
                    continue
                model = self.grdConstraints.GetCellValue(i, 4)
                chain = model[len(model)-1]
                model = model[0:len(model)-2]
                seqpos = self.grdConstraints.GetCellValue(i, 3).split(":")[0]
                seqpos = int(seqpos[0:len(seqpos)-1])
                if (self.seqWin.doesResidueExist(model, chain, seqpos)):
                    pass
                else:
                    updateNeeded = True
                    self.constraints.pop(i)
                    continue
            else:
                # First make sure the chains are still in the right categories
                residueStatic = False
                residueMoving = False
                if (self.grdConstraints.GetCellValue(i, 2) in self.staticChains):
                    residueStatic = True
                if (self.grdConstraints.GetCellValue(i, 2) in self.movingChains):
                    residueMoving = True
                if ((residueStatic and residueMoving) or (not(residueStatic) and not(residueMoving))):
                    self.constraints.pop(i)
                    continue
                elif (residueStatic and len(self.grdConstraints.GetCellValue(i, 4)) > 0 and not(self.grdConstraints.GetCellValue(i, 4) in self.movingChains)):
                    self.constraints.pop(i)
                    continue
                elif (residueMoving and len(self.grdConstraints.GetCellValue(i, 4)) > 0 and not(self.grdConstraints.GetCellValue(i, 4) in self.staticChains)):
                    self.constraints.pop(i)
                    continue
                # Now make sure the residues still exist
                model = self.grdConstraints.GetCellValue(i, 2)
                chain = model[len(model)-1]
                model = model[0:len(model)-2]
                seqpos = self.grdConstraints.GetCellValue(i, 1).split(":")[0]
                seqpos = int(seqpos[0:len(seqpos)-1])
                if (self.seqWin.doesResidueExist(model, chain, seqpos)):
                    pass
                else:
                    self.constraints.pop(i)
                    updateNeeded = True
                    continue
                if (len(self.grdConstraints.GetCellValue(i, 3).strip()) == 0):
                    continue
                # Now check that this chain exists
                model = self.grdConstraints.GetCellValue(i, 4)
                chain = model[len(model)-1]
                if (chain == "_"):
                    chain = " "
                model = model[0:len(model)-2]
                try:
                    poseindx = self.seqWin.getPoseIndexForModel(model)
                    self.seqWin.poses[poseindx][0][chain]
                except:
                    self.constraints.pop(i)
                    updateNeeded = True
                    continue
        if (updateNeeded):
            self.updateConstraints()
            self.updatePartnerList()

    def updateConstraints(self):
        scrollpos = self.grdConstraints.GetScrollPos(wx.VERTICAL)
        if (len(self.constraints) != self.grdConstraints.NumberRows):
            if (self.grdConstraints.NumberRows > 0):
                self.grdConstraints.DeleteRows(0, self.grdConstraints.NumberRows)
            self.selectedr = -1
            self.refreshAdvanced()
            addRows = True
        else:
            addRows = False
        row = 0
        for [constraintType, residue, partner, functionType, functionArgs, group] in self.constraints:
            if (addRows):
                self.grdConstraints.AppendRows(1)
            resID = residue.split(":")[0]
            resmodel = resID #[0:len(resID)-2]
            resatom = residue.split(":")[1] + ":" + residue.split(":")[2]
            if (len(partner.strip()) > 0):
                if (constraintType == "AtomPair"):
                    partID = partner.split(":")[0]
                    partmodel = partID #[0:len(partID)-2]
                    partatom = partner.split(":")[1] + ":" + partner.split(":")[2]
                else:
                    partmodel = partner
                    partatom = partner
            else:
                partmodel = ""
                partatom = ""
            self.grdConstraints.SetRowLabelValue(row, constraintType)
            self.grdConstraints.SetCellValue(row, 0, group)
            self.grdConstraints.SetCellValue(row, 1, resatom)
            self.grdConstraints.SetCellValue(row, 2, resmodel)
            self.grdConstraints.SetCellValue(row, 3, partatom)
            self.grdConstraints.SetCellValue(row, 4, partmodel)
            self.grdConstraints.SetCellValue(row, 5, functionType)
            self.grdConstraints.SetCellValue(row, 6, functionArgs)
            # Very important note: you actually need to create a new GridCellAttr for each new row
            # You cannot just declare it outside of the loop and use the same one for each row otherwise you
            # get some pretty nasty crashes when you delete rows
            readOnly = wx.grid.GridCellAttr()
            readOnly.SetReadOnly(True)
            readOnly.SetAlignment(wx.ALIGN_CENTRE, wx.ALIGN_CENTRE)
            readOnly.SetBackgroundColour("#FFFFFF")
            self.grdConstraints.SetRowAttr(row, readOnly)
            # Now update the drop down menu so the user can tweak the settings of individual residues in the minmap
            row = row + 1
        # Resize columns if necessary
        fitGridColumn(self.grdConstraints, 0, 50)
        fitGridColumn(self.grdConstraints, 1, 100)
        fitGridColumn(self.grdConstraints, 2, 100)
        fitGridColumn(self.grdConstraints, 3, 100)
        fitGridColumn(self.grdConstraints, 4, 100)
        fitGridColumn(self.grdConstraints, 5, 100)
        fitGridColumn(self.grdConstraints, 6, 200)
        # Make the group column writable
        readOnly = wx.grid.GridCellAttr()
        readOnly.SetReadOnly(False)
        self.grdConstraints.SetColAttr(0, readOnly)
        # If the vertical scroll position is given, scroll to that location
        self.grdConstraints.Scroll(0, scrollpos)
        
    def IDtoInt(self, ID):
        # This function converts a residue atom ID to an integer for sorting purposes
        val = 0
        modelchain = ID.split(":")[0]
        for i in range(0, len(self.seqWin.IDs)):
            if (modelchain == self.seqWin.IDs[i]):
                val = val + 10000 * i
                break
        seqpos = ID.split(":")[1]
        val = val + int(seqpos[0:len(seqpos)-1])
        return val
        
    def insertConstraint(self, constraintType, resID, partID, functionType, functionArgs, group):
        # Useful function for inserting constraints into the constraint list in the proper ordering
        if (len(self.constraints) == 0):
            # List empty, add new element
            self.constraints.append([constraintType, resID, partID, functionType, functionArgs, group])
        else:
            notInYet = True
            for i in range(0, len(self.constraints)):
                # All the AtomPair constraints come first
                if (constraintType == "AtomPair" and self.constraints[i][0] == "Site"):
                    self.constraints.insert(i, [constraintType, resID, partID, functionType, functionArgs, group])
                    notInYet = False
                    break
                elif (constraintType == self.constraints[i][0] and self.IDtoInt(resID) <= self.IDtoInt(self.constraints[i][1])):
                    self.constraints.insert(i, [constraintType, resID, partID, functionType, functionArgs, group])
                    notInYet = False
                    break
            if (notInYet):
                # Belongs at the end
                self.constraints.append([constraintType, resID, partID, functionType, functionArgs, group])
        
    def add(self, event):
        #self.activate()
        logInfo("Add button clicked")
        if (len(self.selectedData) == 0):
            return
        # For each of the selected entries, first verify that this entry is not already in the resfile and if it
        # isn't then add it in
        for i in range(0, len(self.selectedData)):
            [indx, r, seqpos, poseindx, chainID, chainoffset] = self.selectedData[i]
            if (self.constraintType == "AtomPair"):
                # If this is not from the static chain, then ignore it
                # The first residue in an atom pair should always be from the static chain to prevent redundancies
                if (not(self.seqWin.IDs[r] in self.staticChains)):
                    continue
            else:
                # For site constraints, just make sure the chain has been added to either the static or moving chains
                if (not(self.seqWin.IDs[r] in self.staticChains) and not(self.seqWin.IDs[r] in self.movingChains)):
                    continue
            # Make sure this is a CAA
            if (not(self.seqWin.getIsCanonicalAA(r, indx))):
                continue
            # The constrained atom defaults to CB, unless the residue is a GLY, in which case it is CA
            if (self.seqWin.sequences[r][indx] == "G"):
                residueAtom = "CA"
            else:
                residueAtom = "CB"
            thisID = self.seqWin.IDs[r] + ":" + str(seqpos) + self.seqWin.sequences[r][indx] + ":" + residueAtom
            alreadyIn = False
            #for j in range(0, len(self.constraints)):
                #[constraintType, residue, partner, functionType, functionArgs, group] = self.constraints[j]
                #if (residue == thisID and self.constraintType == constraintType):
                #    alreadyIn = True
                #    break
            if (not(alreadyIn)):
                # Default parameters are added and the user can change these later
                functionType = self.dFunction
                if (functionType == "Bounded"):
                    functionArgs = "Max: " + str(self.dMax) + ", Min: " + str(self.dMin) + ", Weight: " + str(self.dWeight)
                else:
                    functionArgs = "Ideal: " + str(self.dIdeal) + ", Weight: " + str(self.dWeight)
                self.insertConstraint(self.constraintType, thisID, "", functionType, functionArgs, "")
        self.updateConstraints()
        
    def remove(self, event):
        # For this function, remove the selected constraint
        self.activate()
        logInfo("Remove button clicked")
        if (self.selectedr >= 0 and self.selectedr < len(self.constraints)):
            self.constraints.pop(self.selectedr)
            self.selectedr = -1
        #for i in range(0, len(self.selectedData)):
        #    [indx, r, seqpos, poseindx, chainID, chainoffset] = self.selectedData[i]
        #    ID = self.seqWin.IDs[r] + ":" + str(seqpos) + self.seqWin.sequences[r][indx]
        #    for j in range(len(self.constraints)-1, -1, -1):
        #        [constraintType, residue, partner, functionType, functionArgs, group] = self.constraints[j]
        #        resID = residue.split(":")[0] + ":" + residue.split(":")[1]
        #        if (len(partner.strip()) > 0):
        #            partID = partner.split(":")[0] + ":" + partner.split(":")[1]
        #        else:
        #            partID = ""
        #        if (ID == resID or ID == partID):
        #            self.constraints.pop(j)
        self.updateConstraints()
        
    def clear(self, event):
        logInfo("Clear button clicked")
        # Remove everything
        self.constraints = []
        self.updateConstraints()

    def refreshAdvanced(self):
        # Function to update the advanced controls with whatever the current selection in the grid is
        if (self.selectedr < 0):
            self.menuFunctions.Disable()
            self.txtIdealD.SetValue("0")
            self.txtIdealD.Disable()
            self.txtMaxD.SetValue("0")
            self.txtMaxD.Disable()
            self.txtMinD.SetValue("0")
            self.txtMinD.Disable()
            self.txtWeight.SetValue("0")
            self.txtWeight.Disable()
            self.menuResidueAtoms.Disable()
            self.menuPartnerAtoms.Disable()
        else:
            args = self.constraints[self.selectedr][4].split(",")
            self.menuFunctions.Enable()
            self.txtIdealD.SetValue("0")
            self.txtIdealD.Enable()
            self.txtMaxD.SetValue("0")
            self.txtMaxD.Enable()
            self.txtMinD.SetValue("0")
            self.txtMinD.Enable()
            self.txtWeight.SetValue("0")
            self.txtWeight.Enable()
            self.menuResidueAtoms.Enable()
            self.menuPartnerAtoms.Disable()
            self.menuFunctions.SetValue(self.constraints[self.selectedr][3])
            self.functionMenuTooltip()
            if (self.constraints[self.selectedr][3] == "Bounded"):
                self.txtIdealD.Disable()
            else:
                self.txtMaxD.Disable()
                self.txtMinD.Disable()
            for arg in args:
                if (arg.split(":")[0].strip() == "Ideal"):
                    self.txtIdealD.SetValue(arg.split(":")[1].strip())
                elif (arg.split(":")[0].strip() == "Max"):
                    self.txtMaxD.SetValue(arg.split(":")[1].strip())
                elif (arg.split(":")[0].strip() == "Min"):
                    self.txtMinD.SetValue(arg.split(":")[1].strip())
                elif (arg.split(":")[0].strip() == "Weight"):
                    self.txtWeight.SetValue(arg.split(":")[1].strip())
            # Get all the atoms in this residue
            atoms = []
            model = self.constraints[self.selectedr][1].split(":")[0]
            seqpos = self.constraints[self.selectedr][1].split(":")[1]
            chain = model[len(model)-1]
            model = model[0:len(model)-2]
            if (chain == "_"):
                chain = " "
            seqpos = int(seqpos[0:len(seqpos)-1]) # Trim off the AA
            poseindx = self.seqWin.getPoseIndexForModel(model)
            for a in self.seqWin.poses[poseindx][0][chain][seqpos]:
                atoms.append(a.id)
            self.menuResidueAtoms.Clear()
            self.menuResidueAtoms.AppendItems(atoms)
            self.menuResidueAtoms.SetValue(self.constraints[self.selectedr][1].split(":")[2].strip())
            if (self.constraints[self.selectedr][0] == "AtomPair"):
                if (len(self.constraints[self.selectedr][2].strip()) > 0):
                    self.populatePartnerAtoms()
                else:
                    self.menuPartnerAtoms.SetSelection(0)

    def populatePartnerAtoms(self):
        # Get all the atoms in this residue
        atoms = []
        model = self.constraints[self.selectedr][2].split(":")[0]
        seqpos = self.constraints[self.selectedr][2].split(":")[1]
        chain = model[len(model)-1]
        model = model[0:len(model)-2]
        if (chain == "_"):
            chain = " "
        seqpos = int(seqpos[0:len(seqpos)-1]) # Trim off the AA
        poseindx = self.seqWin.getPoseIndexForModel(model)
        for a in self.seqWin.poses[poseindx][0][chain][seqpos]:
            atoms.append(a.id)
        self.menuPartnerAtoms.Clear()
        self.menuPartnerAtoms.AppendItems(atoms)
        self.menuPartnerAtoms.SetValue(self.constraints[self.selectedr][2].split(":")[2].strip())
        self.menuPartnerAtoms.Enable()

    def updatePartnerList(self):
        # Update the partner menu with the appropriate choices
        if (self.selectedr < 0):
            return
        if (self.grdConstraints.GetRowLabelValue(self.selectedr) == "AtomPair"):
            # Get a list of all the atoms in the moving chains
            residueList = []
            for ID in self.movingChains:
                model = ID[0:len(ID)-2]
                chain = ID[len(ID)-1]
                if (chain == "_"):
                    chain = " "
                poseindx = self.seqWin.getPoseIndexForModel(model)
                for r in self.seqWin.poses[poseindx][0][chain]:
                    seqpos = str(r.id[1])
                    resn = AA3to1(r.resname)
                    if (not(resn in "ACDEFGHIKLMNPQRSTVWY")):
                        continue
                    residueList.append(resn + seqpos + ":" + chain + "|" + model)
            self.menuPartners.Clear()
            self.menuPartners.AppendItems(residueList)
        else:
            # This is easy, all we have to do is add chains from the opposite list of chains
            if (self.grdConstraints.GetCellValue(self.selectedr, 1) in self.staticChains):
                self.menuPartners.Clear()
                self.menuPartners.AppendItems(self.movingChains)
            else:
                self.menuPartners.Clear()
                self.menuPartners.AppendItems(self.staticChains)

    def gridChange(self, event):
        (r, c) = event.GetRow(), event.GetCol()
        # Update the group attribute
        self.constraints[r][5] = self.grdConstraints.GetCellValue(r, c)

    def gridClick(self, event):
        # Set the selected residue's row to blue so it is easy to see what the selection is
        self.selectedr = event.GetRow()
        if (self.selectedr >= self.grdConstraints.NumberRows):
            self.selectedr = -1
        self.refreshAdvanced()
        for r in range(0, self.grdConstraints.NumberRows):
            if (r == self.selectedr):
                for c in range(0, self.grdConstraints.NumberCols):
                    self.grdConstraints.SetCellBackgroundColour(r, c, "light blue")
            else:
                for c in range(0, self.grdConstraints.NumberCols):
                    self.grdConstraints.SetCellBackgroundColour(r, c, "white")
        self.grdConstraints.Refresh()
        if (self.selectedr >= 0):
            #self.constraintView(self.grdConstraints.GetCellValue(self.selectedr, 0).split(":")[0], self.grdConstraints.GetCellValue(self.selectedr, 1))
            self.updatePartnerList()
        event.Skip()
        
    def gridRClick(self, event):
        # This will be a handy function for quickly filling in data into the partners column of constraints using
        # currently selected residues
        indx = 0
        r = event.GetRow()
        while (self.selectedData):
            if (r >= self.grdConstraints.NumberRows):
                break
            if (self.grdConstraints.GetRowLabelValue(r) == "AtomPair"):
                while (True):
                    if (indx >= len(self.selectedData) or self.seqWin.IDs[self.selectedData[indx][1]] in self.movingChains):
                        break
                    indx = indx + 1
                if (indx >= len(self.selectedData)):
                    break
                if (not(self.seqWin.sequences[self.selectedData[indx][1]][self.selectedData[indx][0]] in "ACDEFGHIKLMNPQRSTVWY")):
                    indx = indx + 1
                    continue
                modelID = self.seqWin.IDs[self.selectedData[indx][1]]
                if (self.seqWin.sequences[self.selectedData[indx][1]][self.selectedData[indx][0]] == "G"):
                    atomID = self.selectedData[indx][2] + self.seqWin.sequences[self.selectedData[indx][1]][self.selectedData[indx][0]] + ":CA"
                else:
                    atomID = self.selectedData[indx][2] + self.seqWin.sequences[self.selectedData[indx][1]][self.selectedData[indx][0]] + ":CB"
                self.grdConstraints.SetCellValue(r, 3, atomID)
                self.grdConstraints.SetCellValue(r, 4, modelID)
                self.constraints[r][2] = modelID + ":" + atomID
            else:
                while (True):
                    if (indx >= len(self.selectedData)):
                        break
                    elif (self.grdConstraints.GetCellValue(r, 2) in self.staticChains and self.seqWin.IDs[self.selectedData[indx][1]] in self.movingChains):
                        break
                    elif (self.grdConstraints.GetCellValue(r, 2) in self.movingChains and self.seqWin.IDs[self.selectedData[indx][1]] in self.staticChains):
                        break
                    indx = indx + 1
                if (indx >= len(self.selectedData)):
                    break
                if (not(self.seqWin.sequences[self.selectedData[indx][1]][self.selectedData[indx][0]] in "ACDEFGHIKLMNPQRSTVWY")):
                    indx = indx + 1
                    continue
                modelID = self.seqWin.IDs[self.selectedData[indx][1]]
                self.grdConstraints.SetCellValue(r, 3, modelID)
                self.grdConstraints.SetCellValue(r, 4, modelID)
                self.constraints[r][2] = modelID
            r = r + 1
            indx = indx + 1
        # Resize columns if necessary
        fitGridColumn(self.grdConstraints, 0, 50)
        fitGridColumn(self.grdConstraints, 1, 100)
        fitGridColumn(self.grdConstraints, 2, 100)
        fitGridColumn(self.grdConstraints, 3, 100)
        fitGridColumn(self.grdConstraints, 4, 100)
        fitGridColumn(self.grdConstraints, 5, 100)
        fitGridColumn(self.grdConstraints, 6, 200)
        event.Skip()
        
    def updateFunctionArguments(self):
        # Function for updating the arguments for this constraint's function
        args = ""
        if (len(self.txtIdealD.GetValue().strip()) > 0 and float(self.txtIdealD.GetValue()) > 0):
            args = args + "Ideal: " + self.txtIdealD.GetValue().strip() + ", "
        if (len(self.txtMaxD.GetValue().strip()) > 0 and float(self.txtMaxD.GetValue()) > 0):
            args = args + "Max: " + self.txtMaxD.GetValue().strip() + ", "
        if (len(self.txtMinD.GetValue().strip()) > 0 and float(self.txtMinD.GetValue()) > 0):
            args = args + "Min: " + self.txtMinD.GetValue().strip() + ", "
        if (len(self.txtWeight.GetValue().strip()) > 0 and float(self.txtWeight.GetValue()) > 0):
            args = args + "Weight: " + self.txtWeight.GetValue().strip() + ", "
        if (len(args) > 0):
            args = args[0:len(args)-2] # Trim off the last ,
        self.constraints[self.selectedr][4] = args
        self.grdConstraints.SetCellValue(self.selectedr, 6, args)
        
    def functionMenuTooltip(self):
        if (self.menuFunctions.GetStringSelection() == "Bounded"):
            self.menuFunctions.SetToolTipString("Constraint will obey a bounded function, such that within the minimum/maximum values there is no penalty, but outside the bounds there is a parabolic penalty")
        elif (self.menuFunctions.GetStringSelection() == "Harmonic"):
            self.menuFunctions.SetToolTipString("Constraint will obey a harmonic function, such that at the ideal value there is no penalty, otherwise there is a parabolic penalty")
        elif (self.menuFunctions.GetStringSelection() == "Gaussian"):
            self.menuFunctions.SetToolTipString("Constraint will obey a gaussian function, such that at the ideal value there is no penalty, otherwise there is a penalty that obeys a gaussian curve")
        else:
            self.menuFunctions.SetToolTipString("Constraint will obey a sigmoid function, such that there is no penalty at the ideal value, a slight bonus for being less than the ideal value, and a slight penalty for being greater (you can reverse this behavior by setting a negative weight)")
        
    def functionMenuSelect(self, event):
        # Certain advanced options are disabled/enabled depending on what the function type is
        if (self.menuFunctions.GetStringSelection() == "Bounded"):
            self.txtIdealD.SetValue("0")
            self.txtIdealD.Disable()
            self.txtMaxD.SetValue(str(self.dMax))
            self.txtMaxD.Enable()
            self.txtMinD.SetValue(str(self.dMin))
            self.txtMinD.Enable()
        else:
            self.txtIdealD.SetValue(str(self.dIdeal))
            self.txtIdealD.Enable()
            self.txtMaxD.SetValue("0")
            self.txtMaxD.Disable()
            self.txtMinD.SetValue("0")
            self.txtMinD.Disable()
        self.constraints[self.selectedr][3] = self.menuFunctions.GetStringSelection()
        self.grdConstraints.SetCellValue(self.selectedr, 4, self.menuFunctions.GetStringSelection())
        self.updateFunctionArguments()
        self.functionMenuTooltip()

    def advancedTextUpdate(self, event):
        # Don't let the user enter invalid data
        if (self.selectedr >= 0):
            try:
                if (float(event.GetEventObject().GetValue()) < 0):
                    raise Exception()
                self.updateFunctionArguments()
            except:
                event.GetEventObject().SetValue("0")
        event.Skip()
        
    def residueAtomMenuSelect(self, event):
        # Update the data with the new selection
        resID = self.constraints[self.selectedr][1]
        self.constraints[self.selectedr][1] = resID.split(":")[0] + ":" + resID.split(":")[1] + ":" + self.menuResidueAtoms.GetStringSelection()
        self.grdConstraints.SetCellValue(self.selectedr, 1, resID.split(":")[1] + ":" + self.menuResidueAtoms.GetStringSelection())
        self.constraintView(self.grdConstraints.GetCellValue(self.selectedr, 0).split(":")[0], self.grdConstraints.GetCellValue(self.selectedr, 1), self.grdConstraints.GetCellValue(self.selectedr, 0).split(":")[1])

    def partnerAtomMenuSelect(self, event):
        # Update the data with the new selection
        resID = self.constraints[self.selectedr][2]
        self.constraints[self.selectedr][2] = resID.split(":")[0] + ":" + resID.split(":")[1] + ":" + self.menuPartnerAtoms.GetStringSelection()
        self.grdConstraints.SetCellValue(self.selectedr, 3, resID.split(":")[1] + ":" + self.menuPartnerAtoms.GetStringSelection())
        self.constraintView(self.grdConstraints.GetCellValue(self.selectedr, 2).split(":")[0], self.grdConstraints.GetCellValue(self.selectedr, 3), self.grdConstraints.GetCellValue(self.selectedr, 2).split(":")[1])

    def constraintView(self, posID, origmodel, atom=None):
        self.cmd.set("sphere_scale", 1, "all")
        chain = origmodel[len(origmodel)-1]
        firstmodel = origmodel[0:len(origmodel)-2]
        seqpos = posID[0:len(posID)-1]
        self.cmd.hide("all")
        if (chain == " " or chain == "_"):
            self.cmd.select("viewsele", "resi " + seqpos + " and model " + firstmodel)
        else:
            self.cmd.select("viewsele", "resi " + seqpos + " and model " + firstmodel + " and chain " + chain)
        self.cmd.select("exviewsele", "all within 12 of viewsele")
        self.cmd.show("cartoon", "exviewsele")
        self.cmd.hide("ribbon", "exviewsele")
        self.cmd.show("sticks", "exviewsele")
        self.cmd.set_bond("stick_radius", 0.1, "exviewsele")
        self.cmd.zoom("exviewsele")
        self.cmd.show("sticks", "viewsele")
        self.cmd.set_bond("stick_radius", 0.25, "viewsele")
        self.cmd.select("sele", "viewsele")
        if (atom):
            # Enlarge the selected atom
            if (chain == " " or chain == "_"):
                self.cmd.select("viewsele", "resi " + seqpos + " and model " + firstmodel + " and name " + atom)
            else:
                self.cmd.select("viewsele", "resi " + seqpos + " and model " + firstmodel + " and chain " + chain + " and name " + atom)
            self.cmd.set("sphere_scale", 0.3, "viewsele")
            self.cmd.show("spheres", "viewsele")
        # Highlight this residue in PyMOL
        self.cmd.enable("sele")
        self.cmd.delete("viewsele")
        self.cmd.select("exviewsele", "solvent")
        self.cmd.hide("everything", "exviewsele")
        self.cmd.delete("exviewsele")
        self.seqWin.selectUpdate(False)

    def viewMenuSelect(self, event):
        del self.dockView
        self.cmd.remove("dock_view")
        self.cmd.delete("dock_view")
        self.cmd.load(self.viewMenu.GetStringSelection(), "dock_view")
        self.selectedModel = self.viewMenu.GetStringSelection()
        self.dockView = self.seqWin.pdbreader.get_structure("dock_view", str(self.selectedModel))
        self.focusView(self.viewMenu.GetStringSelection(), "dock_view")
        logInfo("Viewing " + self.viewMenu.GetStringSelection())
    
    def focusView(self, posID, origmodel, newmodel=None):
        model = origmodel
        chainstr = "model " + origmodel + " and (chain " + self.jumpconfig.split("_")[0][0]
        for chain in self.jumpconfig.split("_")[0][1:]:
            chainstr = chainstr + " or chain " + chain
        chainstr = chainstr + ")"
        self.cmd.color("cyan", chainstr)
        chainstr = "model " + origmodel + " and (chain " + self.jumpconfig.split("_")[1][0]
        for chain in self.jumpconfig.split("_")[1][1:]:
            chainstr = chainstr + " or chain " + chain
        chainstr = chainstr + ")"
        self.cmd.color("purple", chainstr)
        if (posID != "Whole View" and not(".pdb" in posID)):
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
            relabelEnergies(self.dockView, self.residue_E, newmodel, self.scoretypeMenu.GetStringSelection(), self.cmd, seqpos)
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
        recolorEnergies(self.dockView, self.residue_E, "dock_view", self.scoretypeMenu.GetStringSelection(), self.cmd)
        self.viewMenuSelect(event) # To update all the labels
    
    def saveCST(self, event):
        if (len(self.constraints) == 0):
            wx.MessageBox("There's nothing to save to a constraints file!", "No Constraints", wx.OK|wx.ICON_EXCLAMATION)
            return
        logInfo("Clicked the Save CST button")
        # Save the data in the resfile graph as an actual resfile for the user's convenience
        dlg = wx.FileDialog(
            self, message="Save a Constraints File",
            defaultDir=self.seqWin.cwd,
            defaultFile="",
            wildcard="Constraints (*.cst)|*.cst",
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
            filename = str(paths[0]).split(".cst")[0] + ".cst"
            # Does it exist already?  If so, ask if the user really wants to overwrite it
            if (os.path.isfile(filename)):
                dlg2 = wx.MessageDialog(self, "The file " + filename + " already exists.  Overwrite it?", "Filename Already Exists", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
                if (dlg2.ShowModal() == wx.ID_NO):
                    dlg2.Destroy()
                    logInfo("Canceled save operation due to filename already existing")
            self.saveConstraints(filename)
        else:
            logInfo("Cancelled save constraints operation")
    
    def loadCST(self, event):
        logInfo("Load CST button clicked")
        # Load data from an existing resfile into the resfile window
        # If there is already data in the graph, notify the user that this data will be erased if they
        # proceed further
        if (len(self.constraints) > 0):
            dlg = wx.MessageDialog(self, "The data in your current workflow will be lost and replaced with a loaded constraints file.  Are you sure you want to proceed?", "Current Constraints Data Will Be Lost", wx.YES_NO | wx.ICON_EXCLAMATION | wx.CENTRE)
            if (dlg.ShowModal() == wx.ID_YES):
                # Don't clear it out until the user actually selects a resfile (instead of cancelling
                # at the file selection dialog
                pass
            else:
                logInfo("Load CST operation cancelled due to data already being in the resfile")
                return
        dlg = wx.FileDialog(
            self, message="Choose a File",
            defaultDir=self.seqWin.cwd,
            defaultFile="",
            wildcard="Constraints (*.cst)|*.cst",
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
            logInfo("Loaded data from a constraints file", filename)
            self.constraints = []
            # Read the data and only add constraints for currently valid information
            # If models are missing for constraints, do not load those constraints
            group = ""
            groupid = 1
            for aline in f:
                if (len(aline.strip()) == 0):
                    continue
                if (aline.strip().lower().startswith("ambiguousconstraint")):
                    group = str(groupid)
                elif (aline.strip().lower().startswith("end")):
                    group = ""
                    groupid += 1
                try:
                    constraintType = aline.split()[0]
                    if (not(constraintType in ["AtomPair", "SiteConstraint"])):
                        raise Exception()
                    if (constraintType == "SiteConstraint"):
                        constraintType = "Site"
                    resatom = aline.split()[1]
                    resloc = aline.split()[2]
                    if (resloc[len(resloc)-1] in "0123456789"):
                        # Undefined chain ID
                        reschain = "_"
                        resseqpos = resloc
                    else:
                        reschain = resloc[len(resloc)-1]
                        resseqpos = resloc[0:len(resloc)-1]
                    if (constraintType == "Site"):
                        partchain = aline.split()[3]
                        funcstart = 4
                    else:
                        partatom = aline.split()[3]
                        partloc = aline.split()[4]
                        if (partloc[len(partloc)-1] in "0123456789"):
                            # Undefined chain ID
                            partchain = "_"
                            partseqpos = partloc
                        else:
                            partchain = partloc[len(partloc)-1]
                            partseqpos = partloc[0:len(partloc)-1]
                        funcstart = 5
                    functype = aline.split()[funcstart]
                    functype = functype[0] + functype[1:].lower() # To get the casing right
                    if (not(functype in ["Bounded", "Sigmoid", "Gaussian", "Harmonic"])):
                        raise Exception()
                    funcargs = None
                    for i in range(funcstart+1, len(aline.split())):
                        if (aline.split()[i][0] == "#"):
                            funcargs = aline.split()[funcstart+1:i]
                            break
                    if (not(funcargs)):
                        raise Exception()
                    commentstart = aline.index("#")
                    resmodel = aline[commentstart:].split("\t")[1].strip()
                    partmodel = aline[commentstart:].split("\t")[3].strip()
                    rescoord = self.seqWin.doesAtomExist(resmodel, reschain, resseqpos, resatom)
                    if (constraintType == "Site"):
                        partcoord = True
                    else:
                        partcoord = self.seqWin.doesAtomExist(partmodel, partchain, partseqpos, partatom)
                    if (not(rescoord) or not(partcoord)):
                        raise Exception()
                    if (constraintType == "Site"):
                        allchains = []
                        allchains.extend(self.staticChains)
                        allchains.extend(self.movingChains)
                        if (not(resmodel + "|" + reschain in allchains) or not(partmodel + "|" + partchain in allchains)):
                            raise Exception()
                    else:
                        if (not(resmodel + "|" + reschain in self.staticChains) or not(partmodel + "|" + partchain in self.movingChains)):
                            raise Exception()
                    # If we've survived this far, all of the data should be valid
                    resID = resmodel + "|" + reschain + ":" + resseqpos + self.seqWin.sequences[rescoord[0]][rescoord[1]] + ":" + resatom
                    if (constraintType == "Site"):
                        partID = partmodel + "|" + partchain
                    else:
                        partID = partmodel + "|" + partchain + ":" + partseqpos + self.seqWin.sequences[partcoord[0]][partcoord[1]] + ":" + partatom
                    # Back calculate what the weight should be (because in the constraints file it was translated to an SD value
                    if (functype == "Sigmoid"):
                        weight = funcargs[1]
                    elif (functype == "Bounded"):
                        weight = str(1.0 / (float(funcargs[2])**2))
                    else:
                        weight = str(1.0 / (float(funcargs[1])**2))
                    if (functype == "Bounded"):
                        funcstr = "Max: " + funcargs[0] + ", Min: " + funcargs[1] + ", Weight: " + weight
                    else:
                        funcstr = "Ideal: " + funcargs[0] + ", Weight: " + weight
                    self.insertConstraint(constraintType, resID, partID, functype, funcstr, group)
                except:
                    continue
            f.close()
            self.updateConstraints()
    
    def writeConstraint(self, f, constraintType, residue, partner, functionType, functionArgs, group, convertChains):
        # Writes out the constraint
        trailingcomment = ""
        residueAtom = residue.split(":")[2].strip()
        if (convertChains):
            residueChain = self.realChainConversion[residue.split(":")[0].strip()]
        else:
            residueChain = residue.split(":")[0].strip()
            trailingcomment = "# MODEL1\t" + residueChain[0:len(residueChain)-2]
            residueChain = residueChain[len(residueChain)-1]
        residueSeqpos = residue.split(":")[1]
        residueSeqpos = residueSeqpos[0:len(residueSeqpos)-1]
        if (convertChains):
            partnerChain = self.realChainConversion[partner.split(":")[0].strip()]
        else:
            partnerChain = partner.split(":")[0].strip()
            trailingcomment += "\tMODEL2\t" + partnerChain[0:len(partnerChain)-2]
            partnerChain = partnerChain[len(partnerChain)-1]
        if (constraintType == "AtomPair"):
            partnerAtom = partner.split(":")[2].strip()
            partnerSeqpos = partner.split(":")[1]
            partnerSeqpos = partnerSeqpos[0:len(partnerSeqpos)-1]
            f.write("AtomPair " + residueAtom + " " + residueSeqpos + residueChain + " " + partnerAtom + " " + partnerSeqpos + partnerChain + " ")
        else:
            f.write("SiteConstraint " + residueAtom + " " + residueSeqpos + residueChain + " " + partnerChain + " ")
        # Now make sure the arguments are all validly defined
        if (functionType == "Gaussian"):
            f.write(functionType.upper() + "FUNC ")
        else:
            f.write(functionType.upper() + " ")
        args = functionArgs.split(",")
        if (functionType.upper() == "BOUNDED"):
            maxD = 12
            minD = 6
            weight = 1
            for arg in args:
                field = arg.split(":")[0].strip()
                val = float(arg.split(":")[1].strip())
                if (field == "Max"):
                    maxD = val
                elif (field == "Min"):
                    minD = val
                elif (field == "Weight"):
                    weight = val
            # Fix any issues
            if (weight <= 0):
                weight = 1
            if (minD < 0):
                minD = 0
            if (maxD < 0):
                maxD = 0
            if (maxD < minD):
                temp = maxD
                maxD = minD
                minD = temp
            sd = 1.0 / math.sqrt(weight)
            #f.write(str(minD) + " " + str(maxD) + " " + str(sd) + " constraint" + str(tagno))
            f.write(str(minD) + " " + str(maxD) + " " + str(sd))
            #tagno = tagno + 1
        else:
            idealD = 9
            weight = 1
            for arg in args:
                field = arg.split(":")[0].strip()
                val = float(arg.split(":")[1].strip())
                if (field == "Ideal"):
                    idealD = val
                elif (field == "Weight"):
                    weight = val
            # Fix any issues
            if (weight <= 0):
                weight = 1
            if (idealD < 0):
                idealD = 0
            sd = 1.0 / math.sqrt(weight)
            if (functionType.upper() == "SIGMOID"):
                f.write(str(idealD) + str(weight))
            elif (functionType.upper() == "HARMONIC"):
                f.write(str(idealD) + str(sd))
            else:
                #f.write(str(idealD) + str(sd) + " constraint" + str(tagno))
                f.write(str(idealD) + str(sd))
                #tagno = tagno + 1
        if (len(trailingcomment) > 0):
            f.write(" " + trailingcomment + "\n")
        else:
            f.write("\n")
    
    def saveConstraints(self, filename, convertChains=False):
        f = open(filename, "w")
        ambig_groups = {}
        tagno = 1
        # Ambiguous constraints have to come first before non-ambiguous ones, otherwise it
        # crashes for no reason
        for [constraintType, residue, partner, functionType, functionArgs, group] in self.constraints:
            if (len(partner.strip()) == 0):
                # Ignore incomplete constraints when saving
                continue
            # Let's skip the grouped constraints, which will end up in ambiguous constraints eventually
            if (len(group.strip()) > 0):
                try:
                    ambig_groups[group].append([constraintType, residue, partner, functionType, functionArgs, group])
                except:
                    ambig_groups[group] = []
                    ambig_groups[group].append([constraintType, residue, partner, functionType, functionArgs, group])
                continue
            # Don't add entries for positions that do not have partners specified
            if (len(partner.strip()) == 0):
                continue
            #self.writeConstraint(f, constraintType, residue, partner, functionType, functionArgs, group, convertChains)
        # Now write out the ambiguous constraints for grouped constraints
        for group in ambig_groups:
            # "group" is a key
            if (len(ambig_groups[group]) > 1):
                f.write("AmbiguousConstraint\n")
            for [constraintType, residue, partner, functionType, functionArgs, group] in ambig_groups[group]:
                self.writeConstraint(f, constraintType, residue, partner, functionType, functionArgs, group, convertChains)
            if (len(ambig_groups[group]) > 1):
                f.write("END\n")
        for [constraintType, residue, partner, functionType, functionArgs, group] in self.constraints:
            if (len(partner.strip()) == 0):
                # Ignore incomplete constraints when saving
                continue
            # Let's skip the grouped constraints
            if (len(group.strip()) > 0):
                continue
            # Don't add entries for positions that do not have partners specified
            if (len(partner.strip()) == 0):
                continue
            self.writeConstraint(f, constraintType, residue, partner, functionType, functionArgs, group, convertChains)
        f.close()
    
    def writeSinglePDB(self):
        # We need to get everything that will be docked into one PDB file, and make sure all the chains are distinct
        # The docker seems to work on a single pose
        # First reload all the structures from PyMOL since the user might have been rotating structures around
        for i in range(0, len(self.seqWin.poses)):
            if (self.seqWin.poses[i]):
                self.seqWin.reloadFromPyMOL(i)
        chains = self.staticChains[:]
        chains.extend(self.movingChains[:])
        self.jumpconfig = ""
        self.realChainConversion = {}
        usedChains = ""
        for i in range(0, len(chains)):
            if (i == len(self.staticChains)):
                self.jumpconfig = self.jumpconfig + "_"
            # Try to preserve chain IDs, but rename them if there are conflicts with two chains from different
            # models having the same chain ID
            chainID = chains[i][len(chains[i])-1]
            if (chainID in usedChains or chainID == "_"):
                # Find a new one
                for char in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
                    if (not(char in usedChains)):
                        chainID = char
                        break
            usedChains += chainID
            self.jumpconfig = self.jumpconfig + chainID # This is a string that helps Rosetta figure out what is static/movable
            self.realChainConversion[chains[i]] = chainID # For the constraints file to be able to convert pre-dock chains to the actual chainIDs in the docked PDB
            posechain = self.seqWin.getChainPose(chains[i])
            posechain.id = chainID
            # Rename the chainIDs
            if (i == 0):
                m = Bio.PDB.Model.Model(0)
                m.add(posechain)
                pose = Bio.PDB.Structure.Structure("dock")
                pose.add(m)
            else:
                pose[0].add(posechain)
        self.dockpdbname = str(chains[0][0:len(chains[0])-2]) + "_dock.pdb"
        self.seqWin.pdbwriter.set_structure(pose)
        self.seqWin.pdbwriter.save(self.dockpdbname)
        # Now let's write the constraints file, if data exists for one
        try:
            os.remove("constraints.cst")
        except:
            pass
        if (len(self.constraints) > 0):
            self.saveConstraints("constraints.cst", convertChains=True)
    
    def toggleStarting(self, event):
        # Toggle the randomize/use current orientation option
        if (self.startingType == "Global"):
            self.startingType = "Fix Stat"
            if (platform.system() == "Darwin"):
                self.btnStarting.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/btnStarting_FixStat.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnStarting.SetLabel(self.startingType)
            self.btnStarting.SetToolTipString("Use the current orientation of the static chains, but the orientation of the moving chains is unknown")
        elif (self.startingType == "Fix Stat"):
            self.startingType = "Fix Mov"
            if (platform.system() == "Darwin"):
                self.btnStarting.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/btnStarting_FixMov.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnStarting.SetLabel(self.startingType)
            self.btnStarting.SetToolTipString("Use the current orientation of the moving chains, but the orientation of the static chains is unknown")
        elif (self.startingType == "Fix Mov"):
            self.startingType = "Fix Both"
            if (platform.system() == "Darwin"):
                self.btnStarting.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/btnStarting_FixBoth.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnStarting.SetLabel(self.startingType)
            self.btnStarting.SetToolTipString("Use the current orientations of both partners")
        else:
            self.startingType = "Global"
            if (platform.system() == "Darwin"):
                self.btnStarting.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/btnStarting_Global.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnStarting.SetLabel(self.startingType)
            self.btnStarting.SetToolTipString("Perform a global dock, where the orientations of both partners is unknown")
        logInfo("Set orientation mode to " + self.startingType)
    
    def serverToggle(self, event):
        if (self.serverOn):
            self.serverOn = False
            if (platform.system() == "Darwin"):
                self.btnServerToggle.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/btnServer_Off.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnServerToggle.SetLabel("Server Off")
            self.btnServerToggle.SetToolTipString("Perform docking simulations locally")
            logInfo("Turned off docking server usage")
        else:
            self.serverOn = True
            if (platform.system() == "Darwin"):
                self.btnServerToggle.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/btnServer_On.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnServerToggle.SetLabel("Server On")
            self.btnServerToggle.SetToolTipString("Perform docking simulations on a remote server")
            logInfo("Turned on docking server usage")
    
    def saveClick(self, event):
        # Save the currently loaded docking structure
        # The Finalize button will only save the one that is currently loaded and the others will be lost
        # This button lets the user save some of the other structures in case they like more than one model
        logInfo("Clicked the Save button")
        # Save the data in the resfile graph as an actual resfile for the user's convenience
        dlg = wx.FileDialog(
            self, message="Save a Docked Structure",
            defaultDir=self.seqWin.cwd,
            defaultFile="",
            wildcard="PDBs (*.pdb)|*.pdb",
            style=wx.SAVE | wx.CHANGE_DIR)
        if (dlg.ShowModal() == wx.ID_OK):
            paths = dlg.GetPaths()
            # Change cwd to the last opened file
            lastDirIndx = paths[len(paths)-1].rfind("/")
            self.seqWin.cwd = str(paths[len(paths)-1][0:lastDirIndx])
            self.seqWin.saveWindowData(None)
            filename = str(paths[0]).split(".pdb")[0] + ".pdb"
            if (os.path.isfile(filename)):
                dlg2 = wx.MessageDialog(self, "The file " + filename + " already exists.  Overwrite it?", "Filename Already Exists", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
                if (dlg2.ShowModal() == wx.ID_NO):
                    dlg2.Destroy()
                    logInfo("Cancelled save operation due to filename already existing")
                    return
                dlg2.Destroy() 
            goToSandbox()
            f = open(str(self.viewMenu.GetStringSelection().strip()), "r")
            f2 = open(filename, "w")
            for aline in f:
                f2.write(aline.strip() + "\n")
            f.close()
            f2.close()
            logInfo("Saved model " + self.viewMenu.GetStringSelection() + " to " + filename)
        else:
            logInfo("Cancelled save operation")
    
    def loadStaticEnsb(self, event):
        logInfo("Load Static Ensemble button clicked")
        # If no static chain has been specified, abort because we cannot check to make sure the ensemble
        # and the static chains are compatible
        if (len(self.staticChains) == 0):
            dlg = wx.MessageDialog(self, "Please specify which chains are in the static structure before specifying an ensemble.", "No Static Structure", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
            dlg.ShowModal()
            dlg.Destroy()
            return
        dlg = wx.FileDialog(
            self, message="Choose a File",
            defaultDir=self.seqWin.cwd,
            defaultFile="",
            wildcard="Ensemble Archives (*.ensb)|*.ensb",
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
            filelabel = filename[lastDirIndx+1:]
            # Does it open?  If yes, then erase the resfile data and continue
            try:
                fin = gzip.open(filename, "r")
            except:
                wx.MessageBox("The file " + filename.strip() + " cannot be opened!", "File Cannot Be Read", wx.OK|wx.ICON_EXCLAMATION)
                return
            logInfo("Loaded data from an ensemble archive", filename)
            # Now we have to read the first model in the archive, find its chains and sequences, and make
            # sure they all match what is already loaded in the static chains
            readingData = False
            modeldata = {}
            lastres = " 0000"
            for aline in fin:
                if (aline.startswith("BEGIN PDB")):
                    readingData = True
                elif (aline.startswith("END PDB")):
                    break
                elif (readingData and (aline.startswith("ATOM") or aline.startswith("HETATM"))):
                    chain = aline[21]
                    if (len(chain.strip()) == 0):
                        chain = "_"
                    if (chain + aline[22:27] != lastres):
                        lastres = chain + aline[22:27]
                        try:
                            modeldata[chain] += AA3to1(aline[17:20])
                        except:
                            modeldata[chain] = AA3to1(aline[17:20])
            fin.close()
            # Now check against the static chains loaded
            for modelchain in self.staticChains:
                chain = modelchain[len(modelchain)-1]
                for i in range(0, len(self.seqWin.IDs)):
                    if (modelchain == self.seqWin.IDs[i]):
                        try:
                            if (modeldata[chain] != self.seqWin.sequences[i]):
                                dlg = wx.MessageDialog(self, "The sequences for chain " + chain + " do not match between the specified static chains and the selected ensemble.", "Ensemble Data Mismatch", wx.OK | wx.ICON_ERROR | wx.CENTRE)
                                dlg.ShowModal()
                                dlg.Destroy()
                                return
                        except:
                            dlg = wx.MessageDialog(self, "Chain " + chain + " in the static chains is not in the ensemble", "Ensemble Data Mismatch", wx.OK | wx.ICON_ERROR | wx.CENTRE)
                            dlg.ShowModal()
                            dlg.Destroy()
                            return
                        modeldata.pop(chain)
            # Are there extra chains left over from the ensemble?
            if (len(modeldata.keys()) > 0):
                dlg = wx.MessageDialog(self, "Chain " + modeldata.keys()[0] + " does not exist in the specified static chains.", "Ensemble Data Mismatch", wx.OK | wx.ICON_ERROR | wx.CENTRE)
                dlg.ShowModal()
                dlg.Destroy()
                return
            # If we got this far, then the ensemble should be okay
            self.ensemble1 = filename
            self.lblSelStaticEnsb.SetLabel(filelabel)
            if (platform.system() == "Linux"):
                resizeTextControlForUNIX(self.lblSelStaticEnsb, 0, 320)
        # Turn the server on by default if it is not already on
        if (not(self.serverOn)):
            self.serverToggle(event)
    
    def deleteStaticEnsb(self, event):
        self.ensemble1 = None
        self.lblSelStaticEnsb.SetLabel("No Ensemble Specified")
        if (platform.system() == "Linux"):
            resizeTextControlForUNIX(self.lblSelStaticEnsb, 0, 320)
    
    def loadMovingEnsb(self, event):
        logInfo("Load Moving Ensemble button clicked")
        # If no static chain has been specified, abort because we cannot check to make sure the ensemble
        # and the static chains are compatible
        if (len(self.movingChains) == 0):
            dlg = wx.MessageDialog(self, "Please specify which chains are in the moving structure before specifying an ensemble.", "No Moving Structure", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
            dlg.ShowModal()
            dlg.Destroy()
            return
        dlg = wx.FileDialog(
            self, message="Choose a File",
            defaultDir=self.seqWin.cwd,
            defaultFile="",
            wildcard="Ensemble Archives (*.ensb)|*.ensb",
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
            filelabel = filename[lastDirIndx+1:]
            # Does it open?  If yes, then erase the resfile data and continue
            try:
                fin = gzip.open(filename, "r")
            except:
                wx.MessageBox("The file " + filename.strip() + " cannot be opened!", "File Cannot Be Read", wx.OK|wx.ICON_EXCLAMATION)
                return
            logInfo("Loaded data from an ensemble archive", filename)
            # Now we have to read the first model in the archive, find its chains and sequences, and make
            # sure they all match what is already loaded in the static chains
            readingData = False
            modeldata = {}
            lastres = " 0000"
            for aline in fin:
                if (aline.startswith("BEGIN PDB")):
                    readingData = True
                elif (aline.startswith("END PDB")):
                    break
                elif (readingData and (aline.startswith("ATOM") or aline.startswith("HETATM"))):
                    chain = aline[21]
                    if (len(chain.strip()) == 0):
                        chain = "_"
                    if (chain + aline[22:27] != lastres):
                        lastres = chain + aline[22:27]
                        try:
                            modeldata[chain] += AA3to1(aline[17:20])
                        except:
                            modeldata[chain] = AA3to1(aline[17:20])
            fin.close()
            # Now check against the moving chains loaded
            for modelchain in self.movingChains:
                chain = modelchain[len(modelchain)-1]
                for i in range(0, len(self.seqWin.IDs)):
                    if (modelchain == self.seqWin.IDs[i]):
                        try:
                            if (modeldata[chain] != self.seqWin.sequences[i]):
                                dlg = wx.MessageDialog(self, "The sequences for chain " + chain + " do not match between the specified moving chains and the selected ensemble.", "Ensemble Data Mismatch", wx.OK | wx.ICON_ERROR | wx.CENTRE)
                                dlg.ShowModal()
                                dlg.Destroy()
                                return
                        except:
                            dlg = wx.MessageDialog(self, "Chain " + chain + " in the moving chains is not in the ensemble", "Ensemble Data Mismatch", wx.OK | wx.ICON_ERROR | wx.CENTRE)
                            dlg.ShowModal()
                            dlg.Destroy()
                            return
                        modeldata.pop(chain)
            # Are there extra chains left over from the ensemble?
            if (len(modeldata.keys()) > 0):
                dlg = wx.MessageDialog(self, "Chain " + modeldata.keys()[0] + " does not exist in the specified moving chains.", "Ensemble Data Mismatch", wx.OK | wx.ICON_ERROR | wx.CENTRE)
                dlg.ShowModal()
                dlg.Destroy()
                return
            # If we got this far, then the ensemble should be okay
            self.ensemble2 = filename
            self.lblSelMovingEnsb.SetLabel(filelabel)
            if (platform.system() == "Linux"):
                resizeTextControlForUNIX(self.lblSelMovingEnsb, 0, 320)
        # Turn the server on by default if it is not already on
        if (not(self.serverOn)):
            self.serverToggle(event)
    
    def deleteMovingEnsb(self, event):
        self.ensemble2 = None
        self.lblSelMovingEnsb.SetLabel("No Ensemble Specified")
        if (platform.system() == "Linux"):
            resizeTextControlForUNIX(self.lblSelMovingEnsb, 0, 320)
    
    def reorient(self, event):
        #from pymol.cgo import *
        # First find the interface in the current selection that are on the static chains
        staticchainstr = "("
        for modelchain in self.staticChains:
            model = modelchain[0:len(modelchain)-2]
            chainID = modelchain[len(modelchain)-1]
            if (chainID == "_"):
                staticchainstr += "(model " + model + ") or "
            else:
                staticchainstr += "(model " + model + " and chain " + chainID + ") or "
        staticchainstr = staticchainstr[0:len(staticchainstr)-4] + ")"
        # And now the moving chains
        movingchainstr = "("
        for modelchain in self.movingChains:
            model = modelchain[0:len(modelchain)-2]
            chainID = modelchain[len(modelchain)-1]
            if (chainID == "_"):
                movingchainstr += "(model " + model + ") or "
            else:
                movingchainstr += "(model " + model + " and chain " + chainID + ") or "
        movingchainstr = movingchainstr[0:len(movingchainstr)-4] + ")"
        try:
            # Now find the average center point of all atoms in this chainset
            self.stored.sum_x = 0.0
            self.stored.sum_y = 0.0
            self.stored.sum_z = 0.0
            self.stored.n = 0
            self.cmd.iterate_state(1, staticchainstr + " and name ca", "stored.sum_x += x; stored.sum_y +=y; stored.sum_z += z; stored.n += 1")
            center_x = self.stored.sum_x / float(self.stored.n)
            center_y = self.stored.sum_y / float(self.stored.n)
            center_z = self.stored.sum_z / float(self.stored.n)
            vcenter = numpy.array([center_x, center_y, center_z])            
            # Now get the center of the interface only
            self.stored.sum_x = 0.0
            self.stored.sum_y = 0.0
            self.stored.sum_z = 0.0
            self.stored.n = 0
            self.cmd.iterate_state(1, staticchainstr + " and seqsele and name ca", "stored.sum_x += x; stored.sum_y +=y; stored.sum_z += z; stored.n += 1")
            interface_x = self.stored.sum_x / float(self.stored.n)
            interface_y = self.stored.sum_y / float(self.stored.n)
            interface_z = self.stored.sum_z / float(self.stored.n)
            vinterface = numpy.array([interface_x, interface_y, interface_z])
            # Now calculate the destination of the opposing chain
            # It will be along the line between these previous two points, four times the distance between them
            vtranslate = vinterface - vcenter
            vdestination = vinterface + 4*vtranslate
            # Now calculate the center of mass of the other chainset
            self.stored.sum_x = 0.0
            self.stored.sum_y = 0.0
            self.stored.sum_z = 0.0
            self.stored.n = 0
            self.cmd.iterate_state(1, movingchainstr + " and name ca", "stored.sum_x += x; stored.sum_y +=y; stored.sum_z += z; stored.n += 1")
            center2_x = self.stored.sum_x / float(self.stored.n)
            center2_y = self.stored.sum_y / float(self.stored.n)
            center2_z = self.stored.sum_z / float(self.stored.n)
            vcenter2 = numpy.array([center2_x, center2_y, center2_z])
            # Now calculate the translation vector
            vtranslate = vdestination - vcenter2
            # Now translate
            self.cmd.translate(list(vtranslate), movingchainstr, camera=0)
            self.cmd.center(staticchainstr + " or " + movingchainstr)
            self.cmd.zoom(staticchainstr + " or " + movingchainstr)
            #self.stored.sum_x = 0.0
            #self.stored.sum_y = 0.0
            #self.stored.sum_z = 0.0
            #self.stored.n = 0
            #self.cmd.iterate_state(1, partner2 + " and name ca", "stored.sum_x += x; stored.sum_y +=y; stored.sum_z += z; stored.n += 1")
            #center3_x = self.stored.sum_x / float(self.stored.n)
            #center3_y = self.stored.sum_y / float(self.stored.n)
            #center3_z = self.stored.sum_z / float(self.stored.n)
            #spherelist = [
            #    COLOR, 1.0, 0.0, 0.0,
            #    SPHERE, center_x, center_y, center_z, 3,
            #    COLOR, 1.0, 0.0, 1.0,
            #    SPHERE, center2_x, center2_y, center2_z, 3,
            #    COLOR, 1.0, 1.0, 0.0,
            #    SPHERE, center3_x, center3_y, center3_z, 3,
            #    COLOR, 0.0, 0.0, 1.0,
            #    SPHERE, vdestination[0], vdestination[1], vdestination[2], 3,
            #    COLOR, 0.0, 1.0, 0.0,
            #    SPHERE, interface_x, interface_y, interface_z, 3]
            #self.cmd.load_cgo(spherelist, "segment", 1)
        except:
            pass
    
    def cancelDock(self):
        logInfo("Canceled docking operation")
        try:
            self.progress.Destroy()
        except:
            pass
        try:
            os.remove("coarsedockinput")
        except:
            pass
        try:
            os.remove("coarsedockinputtemp")
        except:
            pass
        try:
            os.remove("repackmedock_0.pdb")
        except:
            pass
        try:
            os.remove("finedockinput")
        except:
            pass
        self.tmrDock.Stop()
        self.seqWin.cannotDelete = False
        self.parent.GoBtn.Enable()
        self.btnAddStatic.Enable()
        self.btnRemoveStatic.Enable()
        self.btnAddMoving.Enable()
        self.btnRemoveMoving.Enable()
        self.btnStarting.Enable()
        if (platform.system() == "Darwin"):
            self.btnDock.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/btnDock.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
        else:
            self.btnDock.SetLabel("Dock!")
        self.buttonState = "Dock!"
        self.btnDock.SetToolTipString("Perform docking simulation with selected parameters")
        deleteInputFiles()
        self.parent.parent.restartDaemon()
        # Get rid of the messages
        for i in range(0, len(self.seqWin.msgQueue)):
            if (self.seqWin.msgQueue[i].find("Performing protein docking") >= 0):
                self.seqWin.msgQueue.pop(i)
                break
        for i in range(0, len(self.seqWin.msgQueue)):
            if (self.seqWin.msgQueue[i].find("Performing rotamer repacking") >= 0):
                self.seqWin.msgQueue.pop(i)
                break
        for i in range(0, len(self.seqWin.msgQueue)):
            if (self.seqWin.msgQueue[i].find("Performing refined protein docking") >= 0):
                self.seqWin.msgQueue.pop(i)
                break
        if (len(self.seqWin.msgQueue) > 0):
            self.seqWin.labelMsg.SetLabel(self.seqWin.msgQueue[len(self.seqWin.msgQueue)-1])
        else:
            self.seqWin.labelMsg.SetLabel("")
        self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")
    
    def dockClick(self, event):
        # This is also the "Finalize!" button
        if (self.buttonState == "Dock!"):
            # Make sure that there is at least one chain in each of the lists and that the fixed chain list includes
            # a polypeptide chain
            if (len(self.staticChains) == 0):
                wx.MessageBox("Please indicate a chain that will serve as the receptor!", "Static Chain Required", wx.OK|wx.ICON_EXCLAMATION)
                return
            if (len(self.movingChains) == 0):
                wx.MessageBox("Please indicate a chain that will serve as the docked structure!", "Movable Chain Required", wx.OK|wx.ICON_EXCLAMATION)
                return
            ispolypeptide = False
            for ID in self.staticChains:
                ch = self.seqWin.getChainPose(ID)
                for residue in ch:
                    if (residue.resname in "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR"):
                        ispolypeptide = True
                        break
                if (ispolypeptide):
                    break
            if (not(ispolypeptide)):
                wx.MessageBox("The fixed structure is not a protein!  Please add a polypeptide chain to the static model.", "Static Structure Is Not a Protein", wx.OK|wx.ICON_EXCLAMATION)
                return
            # Are the decoy numbers valid?
            try:
                n = int(self.txtCoarse.GetValue())
                if (n <= 0):
                    raise Exception("Coarse decoy number not a positive integer")
            except:
                wx.MessageBox("Please enter a positive integer for the number of coarse decoys.", "Invalid Number of Coarse Decoys", wx.OK|wx.ICON_EXCLAMATION)
                return
            try:
                n = int(self.txtRefined.GetValue())
                if (n <= 0):
                    raise Exception("Refined decoy number not a positive integer")
            except:
                wx.MessageBox("Please enter a positive integer for the number of refined decoys.", "Invalid Number of Refined Decoys", wx.OK|wx.ICON_EXCLAMATION)
                return
            # Ensemble docking is only available on the server
            if ((self.ensemble1 or self.ensemble2) and not(self.serverOn)):
                dlg = wx.MessageDialog(self, "Ensemble docking is only available through the Rosetta server.  Would you like to enable server mode?", "Server Needed for EnsembleDock", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
                if (dlg.ShowModal() == wx.ID_YES):
                    self.serverToggle(None)
                else:
                    return
                dlg.Destroy()
            self.seqWin.labelMsg.SetLabel("Performing protein docking, please be patient...")
            self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
            self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")
            self.seqWin.msgQueue.append("Performing protein docking, please be patient...")
            self.seqWin.cannotDelete = True
            self.parent.GoBtn.Disable()
            self.btnAddStatic.Disable()
            self.btnRemoveStatic.Disable()
            self.btnAddMoving.Disable()
            self.btnRemoveMoving.Disable()
            self.btnStarting.Disable()
            if (platform.system() == "Darwin"):
                self.btnDock.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/btnDock_Cancel.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnDock.SetLabel("Cancel!")
            self.buttonState = "Cancel!"
            self.btnDock.SetToolTipString("Cancel the docking simulation")
            self.stage = 1
            #thrKIC = Thread(target=self.threadKIC, args=())
            #thrKIC.start()
            logInfo("Clicked the Dock button")
            goToSandbox()
            self.tmrDock = wx.Timer(self)
            self.Bind(wx.EVT_TIMER, self.threadDock, self.tmrDock)
            self.writeSinglePDB()
            self.parent.parent.restartDaemon()
            self.tmrDock.Start(1000)
        elif (self.buttonState == "Cancel!"):
            dlg = wx.MessageDialog(self, "Are you sure you want to cancel the docking simulation?  All progress will be lost.", "Cancel Docking Simulation", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
            result = dlg.ShowModal()
            if (result == wx.ID_YES):
                self.cancelDock()
            dlg.Destroy()
        else:
            # Finalize button, ask whether the changes will be accepted or rejected
            dlg = wx.MessageDialog(self, "Do you want to accept the results of this docking session?", "Accept/Reject Model", wx.YES_NO | wx.CANCEL | wx.ICON_QUESTION | wx.CENTRE)
            result = dlg.ShowModal()
            if (result == wx.ID_YES):
                logInfo("Accepted docked model")
                accept = True
                # Get a filename prefix for these models
                while (True):
                    dlg3 = wx.FileDialog(
                        self, message="Save a PDB File",
                        defaultDir=self.seqWin.cwd,
                        defaultFile=self.staticChains[0][0:len(self.staticChains[0])-2],
                        wildcard="PDB Files (*.pdb)|*.pdb",
                        style=wx.SAVE | wx.CHANGE_DIR)
                    if (dlg3.ShowModal() == wx.ID_OK):
                        path = dlg3.GetPath()
                        # Change cwd to the last opened file
                        if (platform.system() == "Windows"):
                            lastDirIndx = path.rfind("\\")
                        else:
                            lastDirIndx = path.rfind("/")
                        self.seqWin.cwd = str(path[0:lastDirIndx])
                        self.seqWin.saveWindowData(None)
                        # Load the PDBs into PyMOL
                        filename = str(path).split(".pdb")[0] + "_0001.pdb"
                        # Does it exist already?  If so, ask if the user really wants to overwrite it
                        if (os.path.isfile(filename)):
                            dlg2 = wx.MessageDialog(self, "The file " + filename + " already exists.  Overwrite it?", "Filename Already Exists", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
                            if (dlg2.ShowModal() == wx.ID_NO):
                                dlg2.Destroy()
                                continue
                            dlg2.Destroy()
                    else:
                        # Default to cancel behavior
                        dlg.Destroy()
                        return
                    filename = filename.split("_0001.pdb")[0] + ".pdb"
                    break
            elif (result == wx.ID_NO):
                logInfo("Rejected docked model")
                accept = False
            else:
                logInfo("Canceled Finalize operation")
                dlg.Destroy()
                return
            dlg.Destroy()
            self.viewMenu.Disable()
            self.parent.GoBtn.Enable()
            self.btnAddStatic.Enable()
            self.btnRemoveStatic.Enable()
            self.btnAddMoving.Enable()
            self.btnRemoveMoving.Enable()
            self.btnStarting.Enable()
            self.grdDocking.ClearGrid()
            if (platform.system() == "Darwin"):
                self.btnDock.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/btnDock.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnDock.SetLabel("Dock!")
            self.buttonState = "Dock!"
            self.btnDock.SetToolTipString("Perform docking simulation with selected parameters")
            #self.btnSave.Disable()
            self.cmd.label("all", "")
            self.seqWin.cannotDelete = False
            if (not(accept)):
                self.cmd.remove("dock_view")
                self.cmd.delete("dock_view")
                return
            # Get rid of the original chains, save the docked pose, and reload the structure in PyMOL
            for ID in self.staticChains:
                row = self.seqWin.IDs.index(ID)
                self.seqWin.deleteChain(row)
            for ID in self.movingChains:
                row = self.seqWin.IDs.index(ID)
                self.seqWin.deleteChain(row)
            newname = filename.split(".pdb")[0]
            if (platform.system() == "Windows"):
                newID = newname[newname.rfind("\\")+1:]
            else:
                newID = newname[newname.rfind("/")+1:]
            self.movingChains = []
            self.staticChains = []
            #if (platform.system() == "Windows"):
                #newname = os.path.expanduser("~") + "\\InteractiveROSETTA\\" + newID
            #else:
                #newname = os.path.expanduser("~") + "/.InteractiveROSETTA/" + newID
            for i in range(0, len(self.viewMenu.GetItems())):
                if (self.selectedModel == self.viewMenu.GetItems()[i]):
                    realID = newID + "_%4.4i" % (i+1)
                    realname = newname + "_%4.4i" % (i+1) + ".pdb"
                try:
                    os.rename(self.viewMenu.GetItems()[i], newname + "_%4.4i" % (i+1) + ".pdb")
                except:
                    # Windows may complain if the file already exists
                    try:
                        os.remove(newname + "_%4.4i" % (i+1) + ".pdb")
                        os.rename(self.viewMenu.GetItems()[i], newname + "_%4.4i" % (i+1) + ".pdb")
                    except:
                        print "ERROR: Could not load the docked structure into the Sequence Window"
                        if (platform.system() == "Windows"):
                            print "The model is currently at " + self.selectedModel + " in " + os.path.expanduser("~") + "\\InteractiveROSETTA"
                        else:
                            print "The model is currently at " + self.selectedModel + " in " + os.path.expanduser("~") + "/.InteractiveROSETTA"
            try:
                self.cmd.remove("dock_view")
                self.cmd.delete("dock_view")
                self.cmd.load(realname, realID)
                self.seqWin.PyMOLPDBLoad(0, realname)
                defaultPyMOLView(self.cmd, realID)
                del self.dockView
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
            sessioninfo = os.path.expanduser("~") + "/.InteractiveRosetta/sessionlog"
        errmsg = errmsg + "\n\nIf you don't know what caused this, send the file " + sessioninfo + " to a developer along with an explanation of what you did."
        # You have to use a MessageDialog because the MessageBox doesn't always work for some reason
        dlg = wx.MessageDialog(self, errmsg, "Error Encountered", wx.OK|wx.ICON_EXCLAMATION)
        dlg.ShowModal()
        dlg.Destroy()
        self.seqWin.cannotDelete = False
        self.parent.GoBtn.Enable()
        self.btnAddStatic.Enable()
        self.btnRemoveStatic.Enable()
        self.btnAddMoving.Enable()
        self.btnRemoveMoving.Enable()
        self.btnDock.Enable()
        if (platform.system() == "Darwin"):
            self.btnDock.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/btnDock.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
        else:
            self.btnDock.SetLabel("Dock!")
        self.buttonState = "Dock!"
        self.btnDock.SetToolTipString("Perform docking simulation with selected parameters")
        # Get rid of the messages
        for i in range(0, len(self.seqWin.msgQueue)):
            if (self.seqWin.msgQueue[i].find("Performing protein docking") >= 0):
                self.seqWin.msgQueue.pop(i)
                break
        for i in range(0, len(self.seqWin.msgQueue)):
            if (self.seqWin.msgQueue[i].find("Performing rotamer repacking") >= 0):
                self.seqWin.msgQueue.pop(i)
                break
        for i in range(0, len(self.seqWin.msgQueue)):
            if (self.seqWin.msgQueue[i].find("Performing refined protein docking") >= 0):
                self.seqWin.msgQueue.pop(i)
                break
        if (len(self.seqWin.msgQueue) > 0):
            self.seqWin.labelMsg.SetLabel(self.seqWin.msgQueue[len(self.seqWin.msgQueue)-1])
        else:
            self.seqWin.labelMsg.SetLabel("")
        self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")
    
    def threadDock(self, event):
        # Dump a file with the loop modeling parameters for the daemon to pick up
        goToSandbox()
        if (self.stage == 1):
            self.tmrDock.Stop()
            self.timeoutCount = 0
            f = open("coarsedockinputtemp", "w")
            f.write("PDBFILE\t" + self.dockpdbname.strip() + "\n")
            try:
                f2 = open(self.dockpdbname.strip(), "r")
            except:
                # Might not have finished writing it yet
                return
            f.write("BEGIN PDB DATA\n")
            for aline in f2:
                f.write(aline.strip() + "\n")
            f.write("END PDB DATA\n")
            f2.close()
            try:
                f2 = open("constraints.cst", "r")
                f.write("BEGIN CST DATA\n")
                for aline in f2:
                    f.write(aline.strip() + "\n")
                f.write("END CST DATA\n")
                f2.close()
            except:
                pass
            f.write("JUMPCONFIG\t" + self.jumpconfig + "\n")
            f.write("ORIENT\t" + self.startingType + "\n")
            f.write("COARSEDECOYS\t" + self.txtCoarse.GetValue() + "\n")
            f.write("REFINEDDECOYS\t" + self.txtRefined.GetValue() + "\n")
            if (self.ensemble1):
                fin = gzip.open(self.ensemble1, "r")
                f.write("BEGIN ENSEMBLE1 DATA\n")
                for aline in fin:
                    f.write(aline)
                f.write("END ENSEMBLE1 DATA\n")
                fin.close()
            if (self.ensemble2):
                fin = gzip.open(self.ensemble2, "r")
                f.write("BEGIN ENSEMBLE2 DATA\n")
                for aline in fin:
                    f.write(aline)
                f.write("END ENSEMBLE2 DATA\n")
                fin.close()
            f.close()
            appendScorefxnParamsInfoToFile("coarsedockinputtemp", self.selectWin.weightsfile)
            if (self.serverOn):
                try: 
                    self.ID = sendToServer("coarsedockinput")
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
                            if (len(aline.split("\t")) >= 2 and aline.split("\t")[0] == "DOCK" and aline.split("\t")[1] == self.ID.strip()):
                                alreadythere = True
                                break
                        f.close()
                    except:
                        pass
                    if (not(alreadythere)):
                        f = open("downloadwatch", "a")
                        f.write("DOCK\t" + self.ID.strip() + "\t" + str(datetime.datetime.now().strftime("%A, %B %d - %I:%M:%S %p")) + "\t" + getServerName() + "\t" + desc + "\n")
                        f.close()
                    dlg = wx.MessageDialog(self, "InteractiveROSETTA is now watching the server for job ID " + desc + ".  You will be notified when the package is available for download.", "Listening for Download", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
                    dlg.ShowModal()
                    dlg.Destroy()
                    # Re-enable everything since we're not waiting for the local daemon to do anything
                    self.viewMenu.Disable()
                    self.parent.GoBtn.Enable()
                    self.btnAddStatic.Enable()
                    self.btnRemoveStatic.Enable()
                    self.btnAddMoving.Enable()
                    self.btnRemoveMoving.Enable()
                    self.btnStarting.Enable()
                    if (platform.system() == "Darwin"):
                        self.btnDock.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/btnDock.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                    else:
                        self.btnDock.SetLabel("Dock!")
                    self.buttonState = "Dock!"
                    self.btnDock.SetToolTipString("Perform docking simulation with selected parameters")
                    self.seqWin.cannotDelete = False
                    # Pop this message out of the queue
                    for i in range(0, len(self.seqWin.msgQueue)):
                        if (self.seqWin.msgQueue[i].find("Performing protein docking") >= 0):
                            self.seqWin.msgQueue.pop(i)
                            break
                    if (len(self.seqWin.msgQueue) > 0):
                        self.seqWin.labelMsg.SetLabel(self.seqWin.msgQueue[len(self.seqWin.msgQueue)-1])
                    else:
                        self.seqWin.labelMsg.SetLabel("")
                    logInfo("Coarse docking input sent to server daemon with ID " + self.ID)
                    return
                except Exception as e:
                    print e.message
                    dlg = wx.MessageDialog(self, "The server could not be reached!  Ensure that you have specified a valid server and that you have an network connection.", "Server Could Not Be Reached", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
                    dlg.ShowModal()
                    dlg.Destroy()
                    return
            else:
                os.rename("coarsedockinputtemp", "coarsedockinput")
                self.usingServer = False
                logInfo("Coarse docking input uploaded locally at coarsedockinput")
                self.stage = 2
            self.looptimecount = 0
            self.timeout = 180
            self.tmrDock.Start(1000)
            try:
                os.remove("dock_progress")
            except:
                pass
            self.progress = wx.ProgressDialog("Coarse Docking", "Performing coarse protein docking...", int(self.txtCoarse.GetValue()), style=wx.PD_CAN_ABORT | wx.PD_APP_MODAL | wx.PD_ELAPSED_TIME | wx.PD_REMAINING_TIME)
        elif (self.stage == 2):
            # This is really annoying, here's the ugly memory problem again
            # So first we have to do a coarse KIC job in the daemon
            # This involves using centroid residues, so those have to be repacked in another 
            # instance of the daemon process because the repacking step pushes the memory usage too
            # high, so first wait for the "repackmetemp.pdb" structure to show up, kill the daemon
            # and restart it to do the repacking step
            if (os.path.isfile("repackmedocktemp_" + str(int(self.txtRefined.GetValue())-1) + ".pdb")):
                self.tmrDock.Stop()
                try:
                    self.progress.Destroy()
                except:
                    pass
                # Pop this message out of the queue
                for i in range(0, len(self.seqWin.msgQueue)):
                    if (self.seqWin.msgQueue[i].find("Performing protein docking") >= 0):
                        self.seqWin.msgQueue.pop(i)
                        break
                self.seqWin.labelMsg.SetLabel("Performing rotamer repacking, please be patient...")
                self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
                self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")
                self.seqWin.msgQueue.append("Performing rotamer repacking, please be patient...")
                self.parent.parent.restartDaemon()
                for i in range(0, int(self.txtRefined.GetValue())):
                    os.rename("repackmedocktemp_" + str(i) + ".pdb", "repackmedock_" + str(i) + ".pdb") # So the new daemon sees it
                logInfo("repackmedocktemp.pdb sent to be rotamer repacked")
                self.stage = 3
                self.tmrDock.Start(1000)
            elif (os.path.isfile("errreport")):
                # Something went wrong, tell the user about it (loop sequence probably too short)
                self.tmrDock.Stop()
                try:
                    self.progress.Destroy()
                except:
                    pass
                self.parent.parent.restartDaemon() # Has to happen because coarse KIC is threaded
                self.recoverFromError()
            elif (os.path.isfile("dock_progress")):
                f = open("dock_progress", "r")
                for aline in f:
                    count = int(aline)
                f.close()
                if (count == int(self.txtCoarse.GetValue())):
                    try:
                        self.progress.Destroy()
                    except:
                        pass
                else:
                    (keepGoing, skip) = self.progress.Update(count)
                    if (not(keepGoing)):
                        # User clicked "Cancel" on the progress bar
                        self.cancelDock()
                        self.progress.Destroy()
        elif (self.stage == 3):
            # Now we have to wait for the output of the repacking step and restart the daemon again
            # so we can finish up with a fine-grained KIC step
            if (os.path.isfile("repackeddock_" + str(int(self.txtRefined.GetValue())-1) + ".pdb")):
                self.tmrDock.Stop()
                # Pop this message out of the queue
                for i in range(0, len(self.seqWin.msgQueue)):
                    if (self.seqWin.msgQueue[i].find("Performing rotamer repacking") >= 0):
                        self.seqWin.msgQueue.pop(i)
                        break
                self.seqWin.labelMsg.SetLabel("Performing refined protein docking, please be patient...")
                self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
                self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")
                self.seqWin.msgQueue.append("Performing refined protein docking, please be patient...")
                self.parent.parent.restartDaemon()
                os.rename("finedockinputtemp", "finedockinput") # So the new daemon sees it
                logInfo("Repacked coarse structure sent to fine grained docking")
                self.stage = 4
                self.tmrDock.Start(1000)
            elif (os.path.isfile("errreport")):
                # Something went wrong, tell the user about it
                self.tmrDock.Stop()
                self.recoverFromError()
        elif (self.stage == 4):
            if (self.usingServer):
                # See if the file has been uploaded yet and bring it here if so
                queryServerForResults("dockoutput-" + self.ID)
                queryServerForResults("coarsedockoutput-" + self.ID)
                self.timeoutCount = self.timeoutCount + 1
            if (self.timeoutCount >= serverTimeout):
                self.tmrDock.Stop()
                # If this is taking too long, maybe there's something wrong with the server
                # Ask the user if they want to continue waiting or use the local daemon instead
                dlg = wx.MessageDialog(self, "The server is taking a long time to respond.  Continue to wait?  Pressing No will run the calculations locally.", "Delayed Server Response", wx.YES_NO | wx.ICON_EXCLAMATION | wx.CENTRE)
                if (dlg.ShowModal() == wx.ID_YES):
                    # Reset the counter
                    self.timeoutCount = 0
                else:
                    self.usingServer = False
                    self.timeoutCount = 0
                    os.rename("coarsedockinputtemp", "coarsedockinput")
                    logInfo("Server took too long to respond so the local daemon was used")
                    self.stage = 2
                dlg.Destroy()
                self.tmrDock.Start(1000)
            # Read the output dumped by the child process (finally!)
            if (os.path.isfile("repackeddocktemp.pdb")):
                # Flip back so the timer sees repacked.pdb and runs the local daemon
                os.rename("coarsedockinputtemp", "finedockinputtemp")
                os.rename("repackeddocktemp.pdb", "repackeddock.pdb")
                # Pop this message out of the queue
                for i in range(0, len(self.seqWin.msgQueue)):
                    if (self.seqWin.msgQueue[i].find("Performing protein docking") >= 0):
                        self.seqWin.msgQueue.pop(i)
                        break
                self.usingServer = False
                self.timeoutCount = 0
                self.stage = 3
            elif (os.path.isfile("dockoutput")):
                self.tmrDock.Stop()
                self.residue_E = []
                self.dockmodels = []
                f = open("dockoutput", "r")
                for aline in f:
                    if (aline[0:6] == "OUTPUT"):
                        pdbfile = aline.split("\t")[1].strip()
                        self.dockmodels.append(pdbfile)
                        self.residue_E.append([])
                        #self.dockView = pose_from_pdb(pdbfile)
                    elif (aline[0:6] == "ENERGY"):
                        if (aline.split()[1] == "total_score"):
                            # This is the scoretype line, row 0 in residue_E
                            self.residue_E[len(self.residue_E)-1].append(aline.split()[1:])
                        else:
                            self.residue_E[len(self.residue_E)-1].append([])
                            indx = len(self.residue_E[len(self.residue_E)-1]) - 1
                            for E in aline.split()[1:]:
                                self.residue_E[len(self.residue_E)-1][indx].append(float(E))
                f.close()
                self.dockView = self.seqWin.pdbreader.get_structure("dock_view", self.dockmodels[0])
                self.selectedModel = self.dockmodels[0]
                logInfo("Found docking output at dockoutput")
                # Pop this message out of the queue
                for i in range(0, len(self.seqWin.msgQueue)):
                    if (self.seqWin.msgQueue[i].find("Performing refined protein docking") >= 0):
                        self.seqWin.msgQueue.pop(i)
                        break
                if (len(self.seqWin.msgQueue) > 0):
                    self.seqWin.labelMsg.SetLabel(self.seqWin.msgQueue[len(self.seqWin.msgQueue)-1])
                else:
                    self.seqWin.labelMsg.SetLabel("")
                self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
                self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")
                # Add these loop residues to the view menu so the user can look at the new loop
                self.viewMenu.Clear()
                self.viewMenu.AppendItems(self.dockmodels)
                self.viewMenu.SetSelection(0)
                self.viewMenu.Enable()
                self.parent.GoBtn.Enable()
                self.btnDock.Enable()
                #self.btnSave.Enable()
                #self.enableControls()
                #self.selectedModel = ""
                if (platform.system() == "Darwin"):
                    self.btnDock.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/docking/btnDock_Finalize.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                else:
                    self.btnDock.SetLabel("Finalize!")
                self.buttonState = "Finalize!"
                self.btnDock.SetToolTipString("Accept or reject protocol results")
                os.remove("dockoutput")
                # Load the docked pose as the "dock_view" model so the user can look at the results
                self.cmd.load(self.dockmodels[0], "dock_view")
                self.cmd.hide("everything", "model dock_view")
                #recolorEnergies(self.dockView, self.residue_E[0], "dock_view", self.scoretypeMenu.GetStringSelection(), self.cmd)
                self.focusView(self.viewMenu.GetStringSelection(), "dock_view")
            elif (os.path.isfile("errreport")):
                # Something went wrong, tell the user about it
                self.tmrDock.Stop()
                self.recoverFromError()
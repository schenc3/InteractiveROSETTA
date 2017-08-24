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

class MinimizationPanel(wx.lib.scrolledpanel.ScrolledPanel):
    def __init__(self, parent, W, H):
        #if (platform.system() == "Windows"):
        wx.lib.scrolledpanel.ScrolledPanel.__init__(self, parent, id=-1, pos=(10, 60), size=(340, H-330), name="ProtMinimization")
        winh = H-330
        #else:
            #wx.lib.scrolledpanel.ScrolledPanel.__init__(self, parent, id=-1, pos=(10, 60), size=(340, H-330), name="ProtMinimization")
            #winh = H-290
        self.SetBackgroundColour("#333333")
        self.parent = parent
        self.sizer = wx.GridBagSizer(1, 1)
        self.SetSizer(self.sizer)

        if (platform.system() == "Windows"):
            self.lblProt = wx.StaticText(self, -1, "Energy Minimization", (25, 15), (270, 25), wx.ALIGN_CENTRE)
            self.lblProt.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblProt = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/lblMinimization.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(25, 15), size=(270, 25))
        else:
            self.lblProt = wx.StaticText(self, -1, "Energy Minimization", (70, 15), style=wx.ALIGN_CENTRE)
            self.lblProt.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
            resizeTextControlForUNIX(self.lblProt, 0, self.GetSize()[0])
        self.lblProt.SetForegroundColour("#FFFFFF")
        self.sizer.Add(self.lblProt, (0, 0), span=(1, 2), flag = wx.ALIGN_CENTER | wx.EXPAND, border=5)

        if (platform.system() == "Darwin"):
            self.HelpBtn = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/HelpBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(295, 10), size=(25, 25))
        else:
            self.HelpBtn = wx.Button(self, id=-1, label="?", pos=(295, 10), size=(25, 25))
            self.HelpBtn.SetForegroundColour("#0000FF")
            self.HelpBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.HelpBtn.Bind(wx.EVT_BUTTON, self.showHelp)
        self.HelpBtn.SetToolTipString("Display the help file for this window")
        self.sizer.Add(self.HelpBtn,(0,2),(1,1))

        if (platform.system() == "Windows"):
            self.lblInst = wx.StaticText(self, -1, "Highlight residues to add/remove to minimization", (0, 45), (320, 25), wx.ALIGN_CENTRE)
            self.lblInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
        elif (platform.system() == "Darwin"):
            self.lblInst = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/lblInstMinimization.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 45), size=(320, 25))
        else:
            self.lblInst = wx.StaticText(self, -1, "Highlight residues to add/remove to minimization", (5, 45), style=wx.ALIGN_CENTRE)
            self.lblInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
            resizeTextControlForUNIX(self.lblInst, 0, self.GetSize()[0])
        self.lblInst.SetForegroundColour("#FFFFFF")
        self.sizer.Add(self.lblInst, (1, 0), span=(1, 5), flag = wx.ALIGN_CENTER | wx.EXPAND, border=5)

        if (platform.system() == "Darwin"):
            self.btnAddBB = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnAddBB.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(7, 70), size=(95, 25))
        else:
            self.btnAddBB = wx.Button(self, id=-1, label="BB", pos=(7, 70), size=(95, 25))
            self.btnAddBB.SetForegroundColour("#000000")
            self.btnAddBB.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAddBB.Bind(wx.EVT_BUTTON, self.selectBB)
        self.btnAddBB.SetToolTipString("Residues added will default to minimize the backbone only")
        if (platform.system() == "Darwin"):
            self.btnAddChi = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnAddChi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(112, 70), size=(95, 25))
        else:
            self.btnAddChi = wx.Button(self, id=-1, label="Chi", pos=(112, 70), size=(95, 25))
            self.btnAddChi.SetForegroundColour("#000000")
            self.btnAddChi.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAddChi.Bind(wx.EVT_BUTTON, self.selectChi)
        self.btnAddChi.SetToolTipString("Residues added will default to minimize the sidechains only")
        if (platform.system() == "Darwin"):
            self.btnAddBoth = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnAddBoth_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(217, 70), size=(95, 25))
        else:
            self.btnAddBoth = wx.Button(self, id=-1, label="Both", pos=(217, 70), size=(95, 25))
            self.btnAddBoth.SetForegroundColour("#FF0000")
            self.btnAddBoth.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAddBoth.Bind(wx.EVT_BUTTON, self.selectBoth)
        self.btnAddBoth.SetToolTipString("Residues added will default to minimize both the backbone and sidechains")
        self.sizer.Add(self.btnAddBB, (2, 0), span=(1,1), flag = wx.ALIGN_CENTER | wx.EXPAND, border=5)
        self.sizer.Add(self.btnAddChi, (2, 1), span=(1,1), flag = wx.ALIGN_CENTER | wx.EXPAND, border=5)
        self.sizer.Add(self.btnAddBoth, (2, 2), span=(1,1), flag = wx.ALIGN_CENTER | wx.EXPAND, border=5)
        self.addType = "BB+Chi"

        #Add constraint menu here????
        self.btnCst = wx.Button(self,id=-1,label="Constraints")
        self.btnCst.SetFont(wx.Font(10,wx.DEFAULT,wx.NORMAL,wx.BOLD))
        self.btnCst.SetForegroundColour("#000000")
        self.btnCst.Bind(wx.EVT_BUTTON,self.open_csts)
        self.sizer.Add(self.btnCst,(4,0),(1,1),flag = wx.ALIGN_CENTER | wx.EXPAND,border=5)
        self.Layout()
        self.ConstraintSet = []

        if (platform.system() == "Darwin"):
            self.btnAdd = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnAdd.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(7, 100), size=(57, 25))
        else:
            self.btnAdd = wx.Button(self, id=-1, label="Add", pos=(7, 100), size=(57, 25))
            self.btnAdd.SetForegroundColour("#000000")
            self.btnAdd.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAdd.Bind(wx.EVT_BUTTON, self.add)
        self.btnAdd.SetToolTipString("Add selected residues to the minimize map")
        if (platform.system() == "Darwin"):
            self.btnRemove = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnRemove.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(69, 100), size=(57, 25))
        else:
            self.btnRemove = wx.Button(self, id=-1, label="Remove", pos=(69, 100), size=(57, 25))
            self.btnRemove.SetForegroundColour("#000000")
            self.btnRemove.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnRemove.Bind(wx.EVT_BUTTON, self.remove)
        self.btnRemove.SetToolTipString("Remove selected residues from the minimize map")
        if (platform.system() == "Darwin"):
            self.btnRestrict = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnRestrict.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(131, 100), size=(57, 25))
        else:
            self.btnRestrict = wx.Button(self, id=-1, label="Restrict", pos=(131, 100), size=(57, 25))
            self.btnRestrict.SetForegroundColour("#000000")
            self.btnRestrict.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnRestrict.Bind(wx.EVT_BUTTON, self.restrict)
        self.btnRestrict.SetToolTipString("Set the minimize map to exclusively contain only the selected residues")
        if (platform.system() == "Darwin"):
            self.btnAll = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnAll.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(193, 100), size=(57, 25))
        else:
            self.btnAll = wx.Button(self, id=-1, label="All", pos=(193, 100), size=(57, 25))
            self.btnAll.SetForegroundColour("#000000")
            self.btnAll.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAll.Bind(wx.EVT_BUTTON, self.addAll)
        self.btnAll.SetToolTipString("Add all residues to the minimize map")
        if (platform.system() == "Darwin"):
            self.btnClear = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnClear.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(255, 100), size=(57, 25))
        else:
            self.btnClear = wx.Button(self, id=-1, label="Clear", pos=(255, 100), size=(57, 25))
            self.btnClear.SetForegroundColour("#000000")
            self.btnClear.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnClear.Bind(wx.EVT_BUTTON, self.clear)
        self.btnClear.SetToolTipString("Clear all entries from the minimize map")
        self.sizer.Add(self.btnAdd, (3, 0), span=(1, 1), flag = wx.ALIGN_CENTER | wx.EXPAND, border=5)
        self.sizer.Add(self.btnRemove, (3, 1), span=(1, 1), flag = wx.ALIGN_CENTER | wx.EXPAND, border=5)
        self.sizer.Add(self.btnRestrict, (3, 2), span=(1, 1), flag = wx.ALIGN_CENTER | wx.EXPAND, border=5)
        self.sizer.Add(self.btnAll, (4, 1), span=(1, 1), flag = wx.ALIGN_CENTER | wx.EXPAND, border=5)
        self.sizer.Add(self.btnClear, (4, 2), span=(1, 1), flag = wx.ALIGN_CENTER | wx.EXPAND, border=5)
        self.selectedData = []

        self.minMenu = wx.ComboBox(self, pos=(7, 130), size=(119, 25), choices=[], style=wx.CB_READONLY)
        self.minMenu.Bind(wx.EVT_COMBOBOX, self.minMenuSelect)
        self.minMenu.SetToolTipString("Select minimize map entries to edit")
        self.selectedModel = ""
        if (platform.system() == "Darwin"):
            self.btnBB = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnBB.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(131, 130), size=(57, 25))
        else:
            self.btnBB = wx.Button(self, id=-1, label="BB", pos=(131, 130), size=(57, 25))
            self.btnBB.SetForegroundColour("#000000")
            self.btnBB.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnBB.Bind(wx.EVT_BUTTON, self.changeBB)
        self.btnBB.SetToolTipString("Set current minimize map selection to minimize the backbone only")
        if (platform.system() == "Darwin"):
            self.btnChi = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnChi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(193, 130), size=(57, 25))
        else:
            self.btnChi = wx.Button(self, id=-1, label="Chi", pos=(193, 130), size=(57, 25))
            self.btnChi.SetForegroundColour("#000000")
            self.btnChi.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnChi.Bind(wx.EVT_BUTTON, self.changeChi)
        self.btnChi.SetToolTipString("Set current minimize map selection to minimize the sidechain only")
        if (platform.system() == "Darwin"):
            self.btnBoth = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnBoth.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(255, 130), size=(57, 25))
        else:
            self.btnBoth = wx.Button(self, id=-1, label="Both", pos=(255, 130), size=(57, 25))
            self.btnBoth.SetForegroundColour("#000000")
            self.btnBoth.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnBoth.Bind(wx.EVT_BUTTON, self.changeBoth)
        self.btnBoth.SetToolTipString("Set current minimize map selection to minimize both the backbone and sidechain")
        self.sizer.Add(self.minMenu, (5, 0), span=(1, 3), flag = wx.ALIGN_CENTER | wx.EXPAND, border=5)
        self.sizer.Add(self.btnBB, (6, 0), span=(1, 1), flag = wx.ALIGN_CENTER | wx.EXPAND, border=5)
        self.sizer.Add(self.btnChi, (6, 1), span=(1, 1), flag = wx.ALIGN_CENTER | wx.EXPAND, border=5)
        self.sizer.Add(self.btnBoth, (6, 2), span=(1, 1), flag = wx.ALIGN_CENTER | wx.EXPAND, border=5)



        if (platform.system() == "Darwin"):
            self.scoretypeMenu = wx.ComboBox(self, pos=(7, 160), size=(305, 25), choices=[], style=wx.CB_READONLY)
        else:
            self.scoretypeMenu = wx.ComboBox(self, pos=(7, 160), size=(305, 25), choices=[], style=wx.CB_READONLY | wx.CB_SORT)
        self.scoretypeMenu.Bind(wx.EVT_COMBOBOX, self.scoretypeMenuSelect)
        self.scoretypeMenu.Disable() # Is only enabled after a minimization and before accepting it
        self.scoretypeMenu.SetToolTipString("Set the scoretype by which PyMOL residues will be colored")
        self.sizer.Add(self.scoretypeMenu, (7, 0), span=(1, 3), flag = wx.ALIGN_CENTER | wx.EXPAND, border=5)

        self.grdMinMap = wx.grid.Grid(self)
        self.grdMinMap.CreateGrid(0, 2)
        if (winh-235 > 200):
            self.grdMinMap.SetSize((320, winh-235))
        else:
            self.grdMinMap.SetSize((320, 200))
        self.grdMinMap.SetPosition((0, 190))
        self.grdMinMap.SetLabelFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.grdMinMap.DisableDragColSize()
        self.grdMinMap.DisableDragRowSize()
        self.grdMinMap.SetColLabelValue(0, "Type")
        self.grdMinMap.SetColLabelValue(1, "Model")
        self.grdMinMap.SetRowLabelSize(80)
        self.grdMinMap.SetColSize(0, 150)
        self.grdMinMap.SetColSize(1, 90)
        self.grdMinMap.Bind(wx.grid.EVT_GRID_CELL_LEFT_CLICK, self.gridClick)
        self.selectedr = -1
        self.sizer.Add(self.grdMinMap, (8, 0), span=(1, 3), flag = wx.ALIGN_CENTER | wx.EXPAND, border=5)
        self.minmap = []

        ypos = self.grdMinMap.GetPosition()[1] + self.grdMinMap.GetSize()[1] + 10
        if (platform.system() == "Darwin"):
            self.btnMinType = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnMinType_Torsion.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(20, ypos), size=(120, 25))
        else:
            self.btnMinType = wx.Button(self, id=-1, label="Torsion", pos=(20, ypos), size=(120, 25))
            self.btnMinType.SetForegroundColour("#000000")
            self.btnMinType.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.btnMinType.Bind(wx.EVT_BUTTON, self.changeMinType)
        self.btnMinType.SetToolTipString("Minimize the models in torsion space (faster)")
        if (platform.system() == "Darwin"):
            self.btnMinimize = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnMinimize.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(175, ypos), size=(120, 25))
        else:
            self.btnMinimize = wx.Button(self, id=-1, label="Minimize!", pos=(175, ypos), size=(120, 25))
            self.btnMinimize.SetForegroundColour("#000000")
            self.btnMinimize.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.btnMinimize.Bind(wx.EVT_BUTTON, self.minimizeClick)
        self.btnMinimize.SetToolTipString("Perform energy minimization")
        self.buttonState = "Minimize!"
        self.sizer.Add(self.btnMinType, (9, 0), span=(1, 2), flag = wx.ALIGN_CENTER | wx.EXPAND, border=5)
        self.sizer.Add(self.btnMinimize, (9, 2), span=(1, 1), flag = wx.ALIGN_CENTER | wx.EXPAND, border=5)
        self.minType = "Torsion"

        self.sizer.AddGrowableRow(8)
        #self.SetSizerAndFit(self.sizer)
        #self.SetupScrolling()


        self.scrollh = self.btnMinimize.GetPosition()[1] + self.btnMinimize.GetSize()[1] + 5
        # self.SetScrollbars(1, 1, 320, self.scrollh)
        self.winscrollpos = 0
        self.Bind(wx.EVT_SCROLLWIN, self.scrolled)
        self.SetupScrolling()
        self.Layout()

    def open_csts(self,event):
        '''Creates the Constraints menu and allows constraints to be added.
        Each time a constraint is added, it is put in the local constraints set,
        so constraints are maintained as the menu is destroyed and recreated'''
        try:
         import constraints
         # print 'constraints imported'
         self.frame = wx.Frame(None,-1,title="Constraints Menu")
         # print 'frame generated'
         self.ConstraintPanel=constraints.ConstraintPanel(self.frame,self)
         self.frame.Fit()
         # print 'constraintpanel created'
         self.frame.Show()
         # print 'showing frame'

         self.ConstraintPanel.setSelectWin(self.selectWin)
         self.ConstraintPanel.setSeqWin(self.seqWin)
         self.ConstraintPanel.setPyMOL(self.pymol)
        except Exception as e:
         import traceback
         # print 'Error importing constraints',e.message
         traceback.print_tb(sys.exc_info()[2])
         pass

    def showHelp(self, event):
        # Open the help page
        if (platform.system() == "Darwin"):
            try:
                browser = webbrowser.get("Safari")
            except:
                print "Could not load Safari!  The help files are located at " + self.scriptdir + "/help"
                return
            browser.open(self.parent.parent.scriptdir + "/help/minimization.html")
        else:
            webbrowser.open(self.parent.parent.scriptdir + "/help/minimization.html")

    def setSeqWin(self, seqWin):
        self.seqWin = seqWin
        # So the sequence window knows about what model "minimized_view" really is
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
        self.minMenu.SetSelection(event.GetRow())
        self.minMenuSelect(event)
        event.Skip()

    def updateMinMap(self):
        # This function redraws the minmap grid to reflect changes to self.minmap
        scrollpos = self.grdMinMap.GetScrollPos(wx.VERTICAL)
        self.minMenu.Clear()
        self.selectedModel = ""
        if ("minimized_view" in self.cmd.get_names("objects")):
            self.cmd.remove("minimized_view")
            self.cmd.delete("minimized_view")
        if (self.grdMinMap.NumberRows > 0):
            self.grdMinMap.DeleteRows(0, self.grdMinMap.NumberRows)
        row = 0
        for [indx, r, seqpos, poseindx, chainoffset, mtype,r_indx] in self.minmap:
            self.grdMinMap.AppendRows(1)
            ID = self.seqWin.IDs[r]
            resn = self.seqWin.SeqViewer.GetCellValue(r, indx)
            label = str(row+1) + ": " + str(resn) + str(seqpos)
            self.grdMinMap.SetRowLabelValue(row, label)
            self.grdMinMap.SetCellValue(row, 0, mtype)
            self.grdMinMap.SetCellValue(row, 1, ID)
            # This needs to happen before the readOnly attr is set otherwise it doesn't apply the center on the
            # last cell for some bizarre reason
            self.grdMinMap.SetCellAlignment(row, 0, wx.ALIGN_CENTRE, wx.ALIGN_CENTRE)
            self.grdMinMap.SetCellAlignment(row, 1, wx.ALIGN_CENTRE, wx.ALIGN_CENTRE)
            # Very important note: you actually need to create a new GridCellAttr for each new row
            # You cannot just declare it outside of the loop and use the same one for each row otherwise you
            # get some pretty nasty crashes when you delete rows
            readOnly = wx.grid.GridCellAttr()
            readOnly.SetReadOnly(True)
            self.grdMinMap.SetRowAttr(row, readOnly)
            # Now update the drop down menu so the user can tweak the settings of individual residues in the minmap
            self.minMenu.AppendItems([label])
            row = row + 1
        # Resize the label width to fit long PDB filenames
        font = self.grdMinMap.GetFont()
        dc = wx.WindowDC(self.grdMinMap)
        dc.SetFont(font)
        maxwidth = 80
        for i in range(0, self.grdMinMap.NumberRows):
            (w, h) = dc.GetTextExtent(self.grdMinMap.GetRowLabelValue(i))
            if (w > maxwidth):
                maxwidth = w
        self.grdMinMap.SetRowLabelSize(maxwidth+10)
        # Resize columns if necessary
        fitGridColumn(self.grdMinMap, 0, 150)
        fitGridColumn(self.grdMinMap, 1, 90)
        self.grdMinMap.Scroll(0, scrollpos)
        # Recolor the grid if after a minimization
        if (self.buttonState != "Minimize!"):
            self.recolorGrid(self.minposes, self.residue_E, self.grdMinMap, self.scoretypeMenu.GetStringSelection())

    def scrolled(self, event):
        self.winscrollpos = self.GetScrollPos(wx.VERTICAL)
        event.Skip()

    def activate(self):
        # It's possible that the user could have deleted chains/residues that are currently in the minmap
        # Let's first make sure everything in the minmap still exists
        redrawNeeded = False
        for i in range(len(self.minmap)-1, -1, -1):
            [indx, r, seqpos, p, offset, mtype,r_indx] = self.minmap[i]
            ID = self.grdMinMap.GetCellValue(i, 1)
            #ID = ID.split(":")[0].strip()
            if (r >= len(self.seqWin.IDs) or ID != self.seqWin.IDs[r]):
                # Chain was deleted, pop this item
                self.minmap.pop(i)
                redrawNeeded = True
            elif (indx >= len(self.seqWin.indxToSeqPos[r]) or self.seqWin.indxToSeqPos[r][indx][1] != int(seqpos)):
                # Residue was deleted, pop this item
                self.minmap.pop(i)
                redrawNeeded = True
        if (redrawNeeded):
            self.updateMinMap()
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
                    chain = self.seqWin.IDs[r][len(self.seqWin.IDs[r])-1]
                    if (chain == "_"):
                        chain = " "
                    # Don't add any NCAAs or HETATMs for now
                    if ("ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".find(self.seqWin.poses[poseindx][0][chain][self.seqWin.indxToSeqPos[r][c]].resname) >= 0):
                        self.selectedData.append([indx, r, seqpos, poseindx, chainoffset])
        self.Scroll(0, self.winscrollpos)

    def selectBB(self, event):
        self.addType = "BB"
        if (platform.system() == "Darwin"):
            self.btnAddBB.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnAddBB_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            self.btnAddChi.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnAddChi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            self.btnAddBoth.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnAddBoth.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
        else:
            self.btnAddBB.SetForegroundColour("#FF0000")
            self.btnAddChi.SetForegroundColour("#000000")
            self.btnAddBoth.SetForegroundColour("#000000")
        logInfo("The add type was changed to BB")

    def selectChi(self, event):
        self.addType = "Chi"
        if (platform.system() == "Darwin"):
            self.btnAddBB.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnAddBB.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            self.btnAddChi.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnAddChi_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            self.btnAddBoth.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnAddBoth.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
        else:
            self.btnAddBB.SetForegroundColour("#000000")
            self.btnAddChi.SetForegroundColour("#FF0000")
            self.btnAddBoth.SetForegroundColour("#000000")
        logInfo("The add type was changed to Chi")

    def selectBoth(self, event):
        self.addType = "BB+Chi"
        if (platform.system() == "Darwin"):
            self.btnAddBB.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnAddBB.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            self.btnAddChi.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnAddChi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            self.btnAddBoth.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnAddBoth_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
        else:
            self.btnAddBB.SetForegroundColour("#000000")
            self.btnAddChi.SetForegroundColour("#000000")
            self.btnAddBoth.SetForegroundColour("#FF0000")
        logInfo("The add type was changed to Both")

    def getR_indx(self,selectedData):
      [indx,r,seqpos,poseindx,co] = selectedData
      model = self.seqWin.poses[poseindx][0]
      offset = 0
      chain = self.seqWin.IDs[r][len(self.seqWin.IDs[r])-1]
      if chain == '_':
        chain == ' '
    #   print self.seqWin.IDs
      for ch in model:
        chain2 = ch.get_id()
        if chain2 =='_':
          chain2 = ' '
        if chain2 == chain:
        #   print offset+indx+1,'R_indx'
          return offset + indx + 1
        else:
          offset += len(model[chain2])
      #if we got this far, there is an error
      return None

    def add(self, event, updateSelection=True):
        if (updateSelection):
            self.activate()
        # For each of the selected entries, first verify that this entry is not already in the minmap and if it
        # isn't then add it in
        logInfo("Add button clicked")
        for i in range(0, len(self.selectedData)):
            [indx, r, seqpos, poseindx, chainoffset] = self.selectedData[i]
            # print "selectedData:",indx, r, seqpos, poseindx, chainoffset
            r_indx = self.getR_indx(self.selectedData[i])
            # Make sure this is a CAA
            if (not(self.seqWin.getIsCanonicalAA(r, indx))):
                continue
            alreadyIn = False
            for j in range(0, len(self.minmap)):
                [mindx, mr, mseqpos, mposeindx, mchainoffset, mtype, mr_indx] = self.minmap[j]
                if (r == mr and indx == mindx):
                    # Just switch the BB/Chi/Both flag in case it is different in this selection
                    self.minmap[j][5] = self.addType
                    alreadyIn = True
                    break
            if (not(alreadyIn)):
                # Now figure out where this belongs so the list stays sorted
                if (len(self.minmap) == 0):
                    # List empty, add new element
                    self.minmap.append([indx, r, seqpos, poseindx, chainoffset, self.addType,r_indx])
                elif (r < self.minmap[0][1] or (r == self.minmap[0][1] and indx < self.minmap[0][0])):
                    # Belongs first
                    self.minmap.insert(0, [indx, r, seqpos, poseindx, chainoffset, self.addType,r_indx])
                else:
                    notInYet = True
                    # Maybe it belongs somewhere in the middle?
                    for i in range(0, len(self.minmap)-1):
                        [indx1, r1, seqpos1, poseindx1, chainoffset1, type1,r_indx1] = self.minmap[i]
                        [indx2, r2, seqpos2, poseindx2, chainoffset2, type2,r_indx2] = self.minmap[i+1]
                        if (r == r1 and r == r2 and indx > indx1 and indx < indx2):
                            notInYet = False
                            self.minmap.insert(i+1, [indx, r, seqpos, poseindx, chainoffset, self.addType,r_indx])
                        elif (r == r1 and r < r2 and indx > indx1):
                            notInYet = False
                            self.minmap.insert(i+1, [indx, r, seqpos, poseindx, chainoffset, self.addType,r_indx])
                    if (notInYet):
                        # Belongs at the end
                        self.minmap.append([indx, r, seqpos, poseindx, chainoffset, self.addType,r_indx])
        self.updateMinMap()

    def remove(self, event):
        self.activate()
        # For each of the selected entries, find out if it is in the minmap and remove it if it is
        logInfo("Remove button clicked")
        if (self.selectedr >= 0 and self.selectedr < len(self.minmap)):
            self.minmap.pop(self.selectedr)
            self.selectedr = -1
        #for i in range(0, len(self.selectedData)):
        #    [indx, r, seqpos, poseindx, chainoffset] = self.selectedData[i]
        #    for j in range(0, len(self.minmap)):
        #        [mindx, mr, mseqpos, mposeindx, mchainoffset, mtype] = self.minmap[j]
        #        if (r == mr and indx == mindx):
        #            self.minmap.pop(j)
        #            break
        self.updateMinMap()

    def restrict(self, event):
        self.activate()
        # Remove everything and add only the selected residues
        logInfo("Restrict button clicked")
        self.minmap = []
        self.add(event)

    def addAll(self, event):
        # Add everything that is in the sequence viewer
        logInfo("All button clicked")
        self.minmap = []
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
                chain = self.seqWin.IDs[r][len(self.seqWin.IDs[r])-1]
                if (chain == "_"):
                    chain = " "
                # Don't add any NCAAs or HETATMs for now
                if ("ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".find(self.seqWin.poses[poseindx][0][chain][self.seqWin.indxToSeqPos[r][c]].resname) >= 0):
                    allData.append([indx, r, seqpos, poseindx, chainoffset])
        saveit = self.selectedData
        self.selectedData = allData
        self.add(event, updateSelection=False)
        self.selectedData = saveit

    def clear(self, event):
        # Remove everything
        logInfo("Clear button clicked")
        self.minmap = []
        self.updateMinMap()
        pass

    def minMenuSelect(self, event):
        # Find this item in the minmap and highlight whatever the minimization type is
        selectedValue = self.minMenu.GetStringSelection()
        logInfo("The minmization menu was changed to " + selectedValue)
        row = 0
        for r in range(0, self.grdMinMap.NumberRows):
            if (self.grdMinMap.GetRowLabelValue(r) == selectedValue):
                row = r
                break
        self.currSelection = row
        # Set the selected residue's row to red so it is easy to see what the selection is
        for r in range(0, self.grdMinMap.NumberRows):
            if (r == row):
                for c in range(0, self.grdMinMap.NumberCols):
                    self.grdMinMap.SetCellBackgroundColour(r, c, "light blue")
            else:
                for c in range(0, self.grdMinMap.NumberCols):
                    self.grdMinMap.SetCellBackgroundColour(r, c, "white")
        self.grdMinMap.Refresh()
        r = row
        if (platform.system() == "Darwin"):
            if (self.minmap[r][5] == "BB"):
                self.btnBB.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnBB_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnChi.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnChi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnBoth.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnBoth.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            elif (self.minmap[r][5] == "Chi"):
                self.btnBB.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnBB.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnChi.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnChi_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnBoth.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnBoth.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnBB.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnBB.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnChi.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnChi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnBoth.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnBoth_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
        else:
            if (self.minmap[r][5] == "BB"):
                self.btnBB.SetForegroundColour("#FF0000")
                self.btnChi.SetForegroundColour("#000000")
                self.btnBoth.SetForegroundColour("#000000")
            elif (self.minmap[r][5] == "Chi"):
                self.btnBB.SetForegroundColour("#000000")
                self.btnChi.SetForegroundColour("#FF0000")
                self.btnBoth.SetForegroundColour("#000000")
            else:
                self.btnBB.SetForegroundColour("#000000")
                self.btnChi.SetForegroundColour("#000000")
                self.btnBoth.SetForegroundColour("#FF0000")
        # Do a neighborhood view on the selected position
        model = ""
        fields = self.grdMinMap.GetCellValue(r, 1).split("|")
        for field in fields[0:len(fields)-1]:
            model = model + field + "|"
        model = model[0:len(model)-1]
        chain = fields[len(fields)-1][0]
        if (model != self.selectedModel):
            self.selectedModel = model
            for i in range(0, len(self.seqWin.IDs)):
                if (self.seqWin.IDs[i].find(model) >= 0):
                    if (self.buttonState == "Minimize!"):
                        pass
                    else:
                        # After a minimization, this should display the minimized structure (which hasn't been
                        # committed yet) along with the original structure
                        # Find the corresponding minimized model
                        for j in range(0, len(self.minmodels)):
                            if (self.minmodels[j] == self.selectedModel):
                                self.seqWin.pdbwriter.set_structure(self.minposes[j])
                                self.seqWin.pdbwriter.save("minimized_view.pdb")
                                break
                        if ("minimized_view" in self.cmd.get_names("objects")):
                            self.cmd.remove("minimized_view")
                            self.cmd.delete("minimized_view")
                        self.cmd.load("minimized_view.pdb", "minimized_view")
                        # Grab the residue energies that were stored from the minimize job
                        indx = self.minmodels.index(str(self.selectedModel))
                        recolorEnergies(self.minposes[indx], self.residue_E[indx], "minimized_view", self.scoretypeMenu.GetStringSelection(), self.cmd)
                    break
        # Find the neighborhood view
        seqpos = selectedValue.split(":")[len(selectedValue.split(":"))-1]
        seqpos = seqpos.strip()
        seqpos = seqpos[1:]
        if (self.buttonState == "Minimize!"):
            firstmodel = model
        else:
            firstmodel = "minimized_view"
        self.cmd.hide("all")
        if (chain == " " or chain == "_"):
            self.cmd.select("minsele", "resi " + seqpos + " and model " + firstmodel)
        else:
            self.cmd.select("minsele", "resi " + seqpos + " and model " + firstmodel + " and chain " + chain)
        self.cmd.select("minsele", "model " + firstmodel + " within 12 of minsele")
        self.cmd.show("cartoon", "minsele")
        self.cmd.hide("ribbon", "minsele")
        self.cmd.show("sticks", "minsele")
        self.cmd.set_bond("stick_radius", 0.1, "minsele")
        # Display energy labels for minimized structures
        if (self.buttonState == "Finalize!"):
            # Grab the residue energies that were stored from the minimize job
            indx = self.minmodels.index(str(self.selectedModel))
            relabelEnergies(self.minposes[indx], self.residue_E[indx], "minimized_view", self.scoretypeMenu.GetStringSelection(), self.cmd, seqpos)
            self.cmd.label("not minsele", "")
        self.cmd.zoom("minsele")
        if (chain == " " or chain == "_"):
            self.cmd.select("minsele", "resi " + seqpos + " and model " + firstmodel)
        else:
            self.cmd.select("minsele", "resi " + seqpos + " and model " + firstmodel + " and chain " + chain)
        self.cmd.show("sticks", "minsele")
        self.cmd.set_bond("stick_radius", 0.25, "minsele")
        # Highlight this residue in PyMOL
        self.cmd.select("seqsele", "minsele")
        if (self.buttonState == "Finalize!"):
            # If this is after a minimization, also show the original structure in green for comparison
            self.cmd.select("minsele", "model " + self.selectedModel + " and symbol c")
            self.cmd.color("green", "minsele")
            self.cmd.set("cartoon_color", "green", "minsele")
            if (chain == " " or chain == "_"):
                self.cmd.select("minsele", "resi " + seqpos + " and model " + self.selectedModel)
            else:
                self.cmd.select("minsele", "resi " + seqpos + " and model " + self.selectedModel + " and chain " + chain)
            self.cmd.select("minsele", "model " + self.selectedModel + " within 12 of minsele")
            self.cmd.show("cartoon", "minsele")
            self.cmd.hide("ribbon", "minsele")
            self.cmd.show("sticks", "minsele")
            self.cmd.set_bond("stick_radius", 0.1, "minsele")
            self.cmd.zoom("minsele")
            if (chain == " " or chain == "_"):
                self.cmd.select("minsele", "resi " + seqpos + " and model " + self.selectedModel)
            else:
                self.cmd.select("minsele", "resi " + seqpos + " and model " + self.selectedModel + " and chain " + chain)
            self.cmd.show("sticks", "minsele")
            self.cmd.set_bond("stick_radius", 0.25, "minsele")
        self.cmd.enable("seqsele")
        self.cmd.delete("minsele")
        self.seqWin.selectUpdate(False)

    def changeBB(self, event):
        if (self.minMenu.GetStringSelection() != ""):
            # Change the selected position to this minimization type
            logInfo("The selected position was changed to BB")
            if (platform.system() == "Darwin"):
                self.btnBB.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnBB_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnChi.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnChi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnBoth.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnBoth.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnBB.SetForegroundColour("#FF0000")
                self.btnChi.SetForegroundColour("#000000")
                self.btnBoth.SetForegroundColour("#000000")
            self.minmap[self.currSelection][5] = "BB"
            self.grdMinMap.SetCellValue(self.currSelection, 0, self.minmap[self.currSelection][5])

    def changeChi(self, event):
        if (self.minMenu.GetStringSelection() != ""):
            # Change the selected position to this minimization type
            logInfo("The selected position was changed to Chi")
            if (platform.system() == "Darwin"):
                self.btnBB.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnBB.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnChi.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnChi_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnBoth.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnBoth.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnBB.SetForegroundColour("#000000")
                self.btnChi.SetForegroundColour("#FF0000")
                self.btnBoth.SetForegroundColour("#000000")
            self.minmap[self.currSelection][5] = "Chi"
            self.grdMinMap.SetCellValue(self.currSelection, 0, self.minmap[self.currSelection][5])

    def changeBoth(self, event):
        if (self.minMenu.GetStringSelection() != ""):
            # Change the selected position to this minimization type
            logInfo("The selected position was changed to Both")
            if (platform.system() == "Darwin"):
                self.btnBB.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnBB.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnChi.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnChi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                self.btnBoth.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnBoth_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnBB.SetForegroundColour("#000000")
                self.btnChi.SetForegroundColour("#000000")
                self.btnBoth.SetForegroundColour("#FF0000")
            self.minmap[self.currSelection][5] = "BB+Chi"
            self.grdMinMap.SetCellValue(self.currSelection, 0, self.minmap[self.currSelection][5])

    def recolorGrid(self, poses, allresidue_E, grid, selectedScoretype):
        # Useful function for recoloring the text color of a grid by energy so the user can see all the relative
        # energies just by looking at the colors in the grid
        for indx in range(0, len(poses)):
            model = self.minmodels[indx]
            found = False
            for (scoretypestr, name) in scoretypes.items():
                if (name == selectedScoretype or scoretypestr == selectedScoretype):
                    scoretype = scoretypestr
                    found = True
                    break
            if (found):
                sindx = allresidue_E[indx][0].index(scoretype)
            else:
                sindx = 0 # Default to total_score
            residue_E = []
            for i in range(1, len(allresidue_E[indx])):
                residue_E.append(allresidue_E[indx][i][sindx])
            residue_E = scale_list(residue_E)
            for row in range(0, grid.NumberRows):
                thismodel = grid.GetCellValue(row, 1)[0:len(grid.GetCellValue(row, 1))-2]
                chain = grid.GetCellValue(row, 1)[len(grid.GetCellValue(row, 1))-1]
                if (chain == "_"):
                    chain = " "
                if (thismodel != model):
                    continue
                # print grid.GetRowLabelValue(row).split()
                seqpos = grid.GetRowLabelValue(row).split()[1]
                seqpos = int(seqpos[1:len(seqpos)])
                # Find the rosetta index
                ires = 0
                found = False
                while not found:
                    for ch in poses[indx][0]:
                      if found:
                        break
                      for residue in ch:
                        if found:
                          break
                        # print chain,seqpos,"=",ch.id,residue.id[1],"?"
#                        ires = ires + 1
                        if (ch.id == chain and residue.id[1] == seqpos):
                        #   print 'Found!'
                          found = True
                        else:
                          ires += 1
#                                  break
                    break
                try:
                  r = residue_E[ires]
                  g = 0
                  b = 255 -r
                except:
                  r = 0
                  b = 0
                  g = 0
                for c in range(0, grid.NumberCols):
                    grid.SetCellTextColour(row, c, (r, g, b))
        grid.Refresh()

    def scoretypeMenuSelect(self, event):
        # Make sure there is even a PyMOL_Mover pose loaded
        if (self.selectedModel == ""):
            return
        logInfo("The scoretype view was changed to " + self.scoretypeMenu.GetStringSelection())
        if (self.buttonState != "Minimize!"):
            # Grab the residue energies that were stored from the minimize job
            indx = self.minmodels.index(str(self.selectedModel))
            recolorEnergies(self.minposes[indx], self.residue_E[indx], "minimized_view", self.scoretypeMenu.GetStringSelection(), self.cmd)
            self.recolorGrid(self.minposes, self.residue_E, self.grdMinMap, self.scoretypeMenu.GetStringSelection())
        self.minMenuSelect(event) # To update all the labels

    def changeMinType(self, event):
        if (self.minType == "Torsion"):
            self.minType = "Cartesian"
            if (platform.system() == "Darwin"):
                self.btnMinType.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnMinType_Cartesian.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnMinType.SetLabel(self.minType)
            self.btnMinType.SetToolTipString("Minimize the models in Cartesian space (slower)")
        else:
            self.minType = "Torsion"
            if (platform.system() == "Darwin"):
                self.btnMinType.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnMinType_Torsion.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnMinType.SetLabel(self.minType)
            self.btnMinType.SetToolTipString("Minimize the models in torsion space (faster)")
        logInfo("The minimization type was changed to " + self.minType)

    def cancelMinimization(self):
        logInfo("Canceled minimization operation")
        try:
            os.remove("minimizeinput")
        except:
            pass
        try:
            os.remove("minimizeinputtemp")
        except:
            pass
        self.tmrMinimize.Stop()
        self.seqWin.cannotDelete = False
        self.btnAddBB.Enable()
        self.btnAddChi.Enable()
        self.btnAddBoth.Enable()
        self.btnAdd.Enable()
        self.btnRemove.Enable()
        self.btnRestrict.Enable()
        self.btnAll.Enable()
        self.btnClear.Enable()
        self.btnBB.Enable()
        self.btnChi.Enable()
        self.btnBoth.Enable()
        self.btnMinType.Enable()
        self.btnCst.Enable()
        if (platform.system() == "Darwin"):
            self.btnMinimize.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnMinimize.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
        else:
            self.btnMinimize.SetLabel("Minimize!")
        self.buttonState = "Minimize!"
        self.btnMinimize.SetToolTipString("Perform energy minimization")
        deleteInputFiles()
        self.parent.parent.restartDaemon()
        self.parent.GoBtn.Enable()
        # Get rid of the messages
        for i in range(0, len(self.seqWin.msgQueue)):
            if (self.seqWin.msgQueue[i].find("Performing energy minimization") >= 0):
                self.seqWin.msgQueue.pop(i)
                break
        if (len(self.seqWin.msgQueue) > 0):
            self.seqWin.labelMsg.SetLabel(self.seqWin.msgQueue[len(self.seqWin.msgQueue)-1])
        else:
            self.seqWin.labelMsg.SetLabel("")
        self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")

    def minimizeClick(self, event):
        # This is also the "Finalize!" button
        logInfo("Minimize button clicked")
        if (self.buttonState == "Minimize!"):
            if (len(self.minmap) > 0):
                self.seqWin.labelMsg.SetLabel("Performing energy minimization, please be patient...")
                self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
                self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")
                self.seqWin.msgQueue.append("Performing energy minimization, please be patient...")
                self.parent.GoBtn.Disable()
                self.btnAddBB.Disable()
                self.btnAddChi.Disable()
                self.btnAddBoth.Disable()
                self.btnAdd.Disable()
                self.btnRemove.Disable()
                self.btnRestrict.Disable()
                self.btnAll.Disable()
                self.btnClear.Disable()
                self.btnBB.Disable()
                self.btnChi.Disable()
                self.btnBoth.Disable()
                self.btnMinType.Disable()
                self.btnCst.Disable()
                self.seqWin.cannotDelete = True
                #thrMinimize = Thread(target=self.threadMinimization, args=())
                #thrMinimize.start()
                self.stage = 1
                self.save_constraints()
                if (platform.system() == "Darwin"):
                    self.btnMinimize.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnMinimize_Cancel.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                else:
                    self.btnMinimize.SetLabel("Cancel!")
                self.buttonState = "Cancel!"
                self.btnMinimize.SetToolTipString("Cancel the energy minimization")
                self.tmrMinimize = wx.Timer(self)
                self.Bind(wx.EVT_TIMER, self.threadMinimization, self.tmrMinimize)
                self.tmrMinimize.Start(1000)
            else:
                wx.MessageBox("There's nothing to minimize!", "Nothing to Minimize", wx.OK|wx.ICON_EXCLAMATION)
        elif (self.buttonState == "Cancel!"):
            dlg = wx.MessageDialog(self, "Are you sure you want to cancel the energy minimization?  All progress will be lost.", "Cancel Energy Minimization", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
            result = dlg.ShowModal()
            if (result == wx.ID_YES):
                self.cancelMinimization()
            dlg.Destroy()
        else:
            # Finalize button, ask whether the changes will be accepted or rejected
            dlg = wx.MessageDialog(self, "Do you want to accept the results of this minimization?", "Accept/Reject Minimization", wx.YES_NO | wx.CANCEL | wx.ICON_QUESTION | wx.CENTRE)
            result = dlg.ShowModal()
            if (result == wx.ID_YES):
                logInfo("Minimization accepted")
                accept = True
            elif (result == wx.ID_NO):
                logInfo("Minimization rejected")
                accept = False
            else:
                logInfo("Finalization operation cancelled")
                dlg.Destroy()
                return
            dlg.Destroy()
            # Set the text colors back to black
            for r in range(0, self.grdMinMap.NumberRows):
                for c in range(0, self.grdMinMap.NumberCols):
                    self.grdMinMap.SetCellTextColour(r, c, "black")
            self.scoretypeMenu.Disable()
            if (platform.system() == "Darwin"):
                self.btnMinimize.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnMinimize.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnMinimize.SetLabel("Minimize!")
            self.buttonState = "Minimize!"
            self.btnMinimize.SetToolTipString("Perform energy minimization")
            self.cmd.label("all", "")
            self.seqWin.cannotDelete = False
            if (not(accept)):
                # try:
                #     self.cmd.remove("minimized_view")
                #     self.cmd.delete("minimized_view")
                # except: 
                #     pass
                try:
                    self.cmd.remove("protocol_view")
                    self.cmd.delete("protocol_view")
                except:
                    pass
                defaultPyMOLView(self.cmd)
                return
            # Get rid of the original poses, save the minimized poses, and reload their structures in PyMOL
            for i in range(0, len(self.minposes)):
                minmodel = str(self.minmodels[i])
                poseindx = -1
                for r in range(0, len(self.seqWin.IDs)):
                    if (self.seqWin.IDs[r].find(minmodel) >= 0):
                        poseindx = r
                        break
                try:
                    self.seqWin.poses[poseindx] = self.seqWin.pdbreader.get_structure(self.minmodels[i], self.minmodels[i] + "_M.pdb")
                    self.cmd.remove(minmodel)
                    self.cmd.delete(minmodel)
                    # try:
                    #     self.cmd.remove("minimized_view")
                    #     self.cmd.delete("minimized_view")
                    # except:
                    #     pass
                    self.cmd.load(minmodel + "_M.pdb", minmodel)
                    # IMPORTANT: You have to replace the model in the sandbox with the new minimized model
                    os.remove(minmodel + ".pdb")
                    os.rename(minmodel + "_M.pdb", minmodel + ".pdb")
                    self.seqWin.rerunDSSPForModel(minmodel)
                    defaultPyMOLView(self.cmd, minmodel)
                except:
                    # Some weird error happened, do nothing instead of crashing
                    print "Bug at accept button click"
                    pass
            self.seqWin.recolorResidues()
            try:
                self.cmd.remove("protocol_view")
                self.cmd.delete("protocol_view")
            except:
                print "Couldn't close protocol_view"
                pass
                

    def recoverFromError(self):
        # This function tells the user what the error was and tries to revert the protocol
        # back to the pre-daemon state so the main GUI can continue to be used
        f = open("errreport", "r")
        errmsg = "An error was encountered during the protocol:\n\n"
        for aline in f:
            errmsg = errmsg + str(aline)
        f.close()
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
        self.seqWin.cannotDelete = True
        self.parent.GoBtn.Enable()
        self.btnMinimize.Enable()
        self.btnAddBB.Enable()
        self.btnAddChi.Enable()
        self.btnAddBoth.Enable()
        self.btnAdd.Enable()
        self.btnRemove.Enable()
        self.btnRestrict.Enable()
        self.btnAll.Enable()
        self.btnClear.Enable()
        self.btnBB.Enable()
        self.btnChi.Enable()
        self.btnBoth.Enable()
        self.btnMinType.Enable()
        self.btnCst.Enable()
        if (platform.system() == "Darwin"):
            self.btnMinimize.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnMinimize.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
        else:
            self.btnMinimize.SetLabel("Minimize!")
        self.buttonState = "Minimize!"
        # Get rid of the messages
        for i in range(0, len(self.seqWin.msgQueue)):
            if (self.seqWin.msgQueue[i].find("Performing energy minimization") >= 0):
                self.seqWin.msgQueue.pop(i)
                break
        if (len(self.seqWin.msgQueue) > 0):
            self.seqWin.labelMsg.SetLabel(self.seqWin.msgQueue[len(self.seqWin.msgQueue)-1])
        else:
            self.seqWin.labelMsg.SetLabel("")
        self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")

    def threadMinimization(self, event):
        # Why am I using a Timer?  See the explanation in kic.py
        if (self.stage == 1):
            self.tmrMinimize.Stop()
            self.timeoutCount = 0
            # First find the number of distinct minimization jobs (since the user is allowed to have multiple
            # models in the minmap
            jobs = []
            currpose = self.minmap[0][3]
            laststart = 0
            for i in range(1, len(self.minmap)):
                if (self.minmap[i][3] != currpose):
                    jobs.append([currpose, laststart, i])
                    laststart = i
                    currpose = self.minmap[i][3]
            jobs.append([currpose, laststart, len(self.minmap)])
            f = open("minimizeinputtemp", "w")
            for [currpose, start, end] in jobs:
                constraintFile = ''
                constraints = self.ConstraintSet
                for [pdb,pose,constraint] in constraints:
                  if str(pose) == str(currpose):
                    constraintFile = str(pose)+".cst"
                    break
                fields = self.seqWin.IDs[currpose].split("|")
                pdbfile = ""
                for field in fields[0:len(fields)-1]:
                    pdbfile = pdbfile + field + "|"
                pdbfile = pdbfile[0:len(pdbfile)-1] + ".pdb"
                # Dump the PDB from PyMOL first in case the coordinates were altered by the user
                dumpmodel = pdbfile.split(".pdb")[0]
                self.cmd.save(pdbfile.strip(), "model " + dumpmodel)
                fixPyMOLSave(pdbfile.strip())
                f.write("JOB\t" + pdbfile + "\t" + str(start) + "\t" + str(end) + "\t"+ constraintFile+"\n")
                f2 = open(pdbfile, "r")
                f.write("BEGIN PDB DATA\n")
                for aline in f2:
                    f.write(aline.strip() + "\n")
                f.write("END PDB DATA\n")
                f2.close()
            for [indx, r, seqpos, p, co, mtype,r_indx] in self.minmap:
                f.write("MINMAP\t" + str(indx) + "\t" + str(r) + "\t" + str(seqpos) + "\t" + str(p) + "\t" + str(co) + "\t" + str(mtype) + "\t"+str(r_indx)+"\n")
            f.write("MINTYPE\t" + self.minType + "\n")
            f.write("SCOREFXN\t" + self.selectWin.weightsfile + "\n")
            f.close()
            # Get PyMOL ready for minimization viewing
            self.cmd.copy("protocol_view", pdbfile.split(".pdb")[0])
            self.cmd.hide("everything", "not protocol_view")
            appendScorefxnParamsInfoToFile("minimizeinputtemp", self.selectWin.weightsfile)
            if (useServer):
                try:
                    self.ID = sendToServer("minimizeinput")
                    self.usingServer = True
                    logInfo("Minimization input sent to server daemon with ID " + self.ID)
                except:
                    # Something failed, default to the local daemon
                    os.rename("minimizeinputtemp", "minimizeinput")
                    self.usingServer = False
                    logInfo("Server daemon not available, minimization input uploaded at minimizeinput")
            else:
                os.rename("minimizeinputtemp", "minimizeinput")
                self.usingServer = False
                logInfo("Minimization input uploaded locally at minimizeinput")
            self.stage = 2
            self.seqWin.protocol_view_active = True
            self.tmrMinimize.Start(1000)
        else:
            if (self.usingServer):
                # See if the file has been uploaded yet and bring it here if so
                queryServerForResults("minimizeoutput-" + self.ID)
                self.timeoutCount = self.timeoutCount + 1
            if (self.timeoutCount >= serverTimeout):
                self.tmrMinimize.Stop()
                # If this is taking too long, maybe there's something wrong with the server
                # Ask the user if they want to continue waiting or use the local daemon instead
                dlg = wx.MessageDialog(self, "The server is taking a long time to respond.  Continue to wait?  Pressing No will run the calculations locally.", "Delayed Server Response", wx.YES_NO | wx.ICON_EXCLAMATION | wx.CENTRE)
                if (dlg.ShowModal() == wx.ID_YES):
                    # Reset the counter
                    self.timeoutCount = 0
                else:
                    self.usingServer = False
                    self.timeoutCount = 0
                    os.rename("minimizeinputtemp", "minimizeinput")
                    logInfo("Server took too long to respond so the local daemon was used")
                self.tmrMinimize.Start(1000)
            if (os.path.isfile("minimizeoutput")):
                self.tmrMinimize.Stop()
                # Read the output dumped by the child process
                self.minposes = []
                self.minmodels = []
                self.residue_E = []
                pdbreader = Bio.PDB.PDBParser()
                f = open("minimizeoutput", "r")
                for aline in f:
                    if (aline[0:6] == "OUTPUT"):
                        pdbfile = aline.split("\t")[1].strip()
                        self.minposes.append(pdbreader.get_structure(pdbfile.split("_M.pdb")[0], pdbfile))
                        self.minmodels.append(pdbfile.split("_M.pdb")[0])
                        self.residue_E.append([])
                    elif (aline[0:6] == "ENERGY"):
                        indx = len(self.residue_E) - 1
                        if (aline.split()[1] == "total_score"):
                            # This is the scoretype line, row 0 in residue_E[indx]
                            self.residue_E[indx].append(aline.split()[1:])
                        else:
                            self.residue_E[indx].append([])
                            indx2 = len(self.residue_E[indx]) - 1
                            for E in aline.split()[1:]:
                                self.residue_E[indx][indx2].append(float(E))
                f.close()
                for i in range(0, len(self.minposes)):
                    # To get the energy values in the B-factors
                    recolorEnergies(self.minposes[i], self.residue_E[i], "nomodel", "Total Energy", self.cmd)
                    self.seqWin.pdbwriter.set_structure(self.minposes[i])
                    self.seqWin.pdbwriter.save(self.minmodels[i] + "_M.pdb")
                logInfo("Found minimization output at minimizeoutput")
                # Add the nonzero scoretypes to the energy viewing list from the current score function
                self.scoretypeMenu.Clear()
                #self.scoretypeMenu.Append("Total Energy")
                for scoretype in self.residue_E[0][0]:
                    try:
                        toAdd = scoretypes[str(scoretype)]
                    except:
                        toAdd = str(scoretype)
                    self.scoretypeMenu.Append(toAdd)
                self.scoretypeMenu.Enable()
                # Pop this message out of the queue
                for i in range(0, len(self.seqWin.msgQueue)):
                    if (self.seqWin.msgQueue[i].find("Performing energy minimization") >= 0):
                        self.seqWin.msgQueue.pop(i)
                        break
                if (len(self.seqWin.msgQueue) > 0):
                    self.seqWin.labelMsg.SetLabel(self.seqWin.msgQueue[len(self.seqWin.msgQueue)-1])
                else:
                    self.seqWin.labelMsg.SetLabel("")
                self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
                self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")
                self.parent.GoBtn.Enable()
                self.btnMinimize.Enable()
                self.btnAddBB.Enable()
                self.btnAddChi.Enable()
                self.btnAddBoth.Enable()
                self.btnAdd.Enable()
                self.btnRemove.Enable()
                self.btnRestrict.Enable()
                self.btnAll.Enable()
                self.btnClear.Enable()
                self.btnBB.Enable()
                self.btnChi.Enable()
                self.btnBoth.Enable()
                self.btnMinType.Enable()
                self.btnCst.Enable()
                self.selectedModel = ""
                if (platform.system() == "Darwin"):
                    self.btnMinimize.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/minimization/btnMinimize_Finalize.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                else:
                    self.btnMinimize.SetLabel("Finalize!")
                self.buttonState = "Finalize!"
                self.btnMinimize.SetToolTipString("Accept or reject the results of this protocol")
                os.remove("minimizeoutput")
                self.recolorGrid(self.minposes, self.residue_E, self.grdMinMap, self.scoretypeMenu.GetStringSelection())
                self.seqWin.protocol_view_active = False
            elif (os.path.isfile("errreport")):
                self.tmrMinimize.Stop()
                self.recoverFromError()
                self.seqWin.protocol_view_active = False
                try:
                    self.cmd.remove("protocol_view")
                except:
                    pass

    def save_constraints(self):
      constraints = self.ConstraintSet
      goToSandbox()
      outputs = {}
      for [pdb,poseindx,constraint] in constraints:
        if pdb not in outputs:
          outputs[pdb] = open('%s.cst'%(str(poseindx)),'w+')
          outputs[pdb].write('#%s\n'%(pdb))
        outputs[pdb].write('%s\n'%(constraint))
      for pdb in outputs:
        outputs[pdb].close()

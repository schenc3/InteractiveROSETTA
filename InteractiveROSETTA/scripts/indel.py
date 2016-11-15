### AUTHOR: William F. Hooper
### Affiliation: Rensselaer Polytechnic Institute
### Based off of kic.py, included with InteractiveROSETTA
### Additional dependencies: LoopHash, available at github.com/willhooper/LoopHash (might be packaged with this install already)
### Database files iRosetta_Lookup.exe, pdbselect.dat, looplist.dat, and grid.dat should be in sandbox

import wx
import wx.grid
import wx.lib.scrolledpanel
import wx.lib.intctrl
import os
import os.path
import time
import platform
import multiprocessing
import webbrowser
import datetime
from threading import Thread
from tools import *

class INDELmodelPanel(wx.lib.scrolledpanel.ScrolledPanel):
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
            self.lblProt = wx.StaticText(self, -1, "INDEL Loop Design", (25, 15), (270, 25), wx.ALIGN_CENTRE)
            self.lblProt.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        # elif (platform.system() == "Darwin"):
        #     self.lblProt = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/indel/label_INDEL.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(25, 15), size=(270, 25))
        else:
            self.lblProt = wx.StaticText(self, -1, "INDEL Loop Design", (70, 15), style=wx.ALIGN_CENTRE)
            self.lblProt.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
            resizeTextControlForUNIX(self.lblProt, 0, self.GetSize()[0]-20)
        self.lblProt.SetForegroundColour("#FFFFFF")

        # if (platform.system() == "Darwin"):
        #     self.HelpBtn = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/HelpBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(295, 10), size=(25, 25))
        # else:
        self.HelpBtn = wx.Button(self, id=-1, label="?", pos=(295, 10), size=(25, 25))
        self.HelpBtn.SetForegroundColour("#0000FF")
        self.HelpBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.HelpBtn.Bind(wx.EVT_BUTTON, self.showHelp)
        self.HelpBtn.SetToolTipString("Display the help file for this window")

        if (platform.system() == "Windows"):
            self.lblInst = wx.StaticText(self, -1, "Remodels loops via a \n  fragment database search", (0, 45), (320, 25), wx.ALIGN_CENTRE)
            self.lblInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
        # elif (platform.system() == "Darwin"):
        #     self.lblInst = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/indel/lbl_description_INDEL.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 45), size=(320, 25))
        else:
            self.lblInst = wx.StaticText(self, -1, "Remodels loops via a \n  fragment database search", (5, 45), style=wx.ALIGN_CENTRE)
            self.lblInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
            resizeTextControlForUNIX(self.lblInst, 0, self.GetSize()[0]-20)
        self.lblInst.SetForegroundColour("#FFFFFF")

        # Model selection

        if (platform.system() == "Windows"):
            self.lblModel = wx.StaticText(self, -1, "Model", (10, 90), (140, 20), wx.ALIGN_CENTRE)
            self.lblModel.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        # elif (platform.system() == "Darwin"):
        #     self.lblModel = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/kic/lblModelKIC.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(10, 90), size=(140, 20))
        else:
            self.lblModel = wx.StaticText(self, -1, "Model", (10, 90), style=wx.ALIGN_CENTRE)
            self.lblModel.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblModel, 10, 140)
        self.lblModel.SetForegroundColour("#FFFFFF")
        self.modelMenu = wx.ComboBox(self, pos=(10, 110), size=(140, 25), choices=[], style=wx.CB_READONLY)
        self.modelMenu.Bind(wx.EVT_COMBOBOX, self.modelMenuSelect)
        self.modelMenu.SetToolTipString("Model on which to perform loop modeling")
        self.selectedModel = ""

        # N-term anchor selection

        if (platform.system() == "Windows"):
            self.lblBegin = wx.StaticText(self, -1, "Loop Begin", (10, 140), (120, 20), wx.ALIGN_CENTRE)
            self.lblBegin.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        # elif (platform.system() == "Darwin"):
        #     self.lblBegin = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/kic/lblBegin.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(10, 140), size=(140, 20))
        else:
            self.lblBegin = wx.StaticText(self, -1, "Loop Begin", (10, 140), style=wx.ALIGN_CENTRE)
            self.lblBegin.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblBegin, 10, 140)
        self.lblBegin.SetForegroundColour("#FFFFFF")
        self.beginMenu = wx.ComboBox(self, pos=(10, 160), size=(140, 25), choices=[], style=wx.CB_READONLY)
        self.beginMenu.Bind(wx.EVT_COMBOBOX, self.beginMenuSelect)
        self.beginMenu.Bind(wx.EVT_RIGHT_DOWN, self.rightClick)
        self.beginMenu.SetToolTipString("Loop N-terminus")
        self.loopBegin = -1

        # C-term anchor selection

        if (platform.system() == "Windows"):
            self.lblEnd = wx.StaticText(self, -1, "Loop End", (170, 140), (140, 20), wx.ALIGN_CENTRE)
            self.lblEnd.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        # elif (platform.system() == "Darwin"):
        #     self.lblEnd = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/kic/lblEnd.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(170, 140), size=(140, 20))
        else:
            self.lblEnd = wx.StaticText(self, -1, "Loop End", (170, 140), style=wx.ALIGN_CENTRE)
            self.lblEnd.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblEnd, 170, 140)
        self.lblEnd.SetForegroundColour("#FFFFFF")
        self.endMenu = wx.ComboBox(self, pos=(170, 160), size=(140, 25), choices=[], style=wx.CB_READONLY)
        self.endMenu.Bind(wx.EVT_COMBOBOX, self.endMenuSelect)
        self.endMenu.Bind(wx.EVT_RIGHT_DOWN, self.rightClick)
        self.endMenu.SetToolTipString("Loop C-terminus")
        self.loopEnd = -1

        # Minimum loop length (in residues)
        # can't get wx.ComboBox to accept a list of integers as choices
        minmax_length = ['3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19']


        if (platform.system() == "Windows"):
            self.lblMin = wx.StaticText(self, -1, "Minimum length", (10, 190), (140, 20), wx.ALIGN_CENTRE)
            self.lblMin.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        # elif (platform.system() == "Darwin"):
        #     self.lblMin = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/indel/minLength.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(10, 140), size=(140, 20))
        else:
            self.lblMin = wx.StaticText(self, -1, "Minimum length", (20, 190), style=wx.ALIGN_CENTRE)
            self.lblMin.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblEnd, 170, 140)
        self.lblMin.SetForegroundColour("#FFFFFF")
        self.minMenu = wx.ComboBox(self, pos=(10, 210), size=(140, 25), choices=minmax_length, style=wx.CB_READONLY)
        self.minMenu.Bind(wx.EVT_COMBOBOX, self.minMenuSelect)
        self.minMenu.Bind(wx.EVT_RIGHT_DOWN, self.rightClick)
        self.minMenu.SetToolTipString("Minimum length of loop in residues")
        self.minMenu.SetSelection(0)
        self.minLength = int(self.minMenu.GetStringSelection())


        # Max loop length (in residues)

        if (platform.system() == "Windows"):
            self.lblMax = wx.StaticText(self, -1, "Maximum length", (170, 190), (140, 20), wx.ALIGN_CENTRE)
            self.lblMax.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        # elif (platform.system() == "Darwin"):
        #     self.lblMax = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/indel/maxLength.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(10, 140), size=(140, 20))
        else:
            self.lblMax = wx.StaticText(self, -1, "Maximum length", (180, 190), style=wx.ALIGN_CENTRE)
            self.lblMax.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblEnd, 170, 140)
        self.lblMax.SetForegroundColour("#FFFFFF")
        self.maxMenu = wx.ComboBox(self, pos=(170, 210), size=(140, 25), choices=minmax_length, style=wx.CB_READONLY)
        self.maxMenu.Bind(wx.EVT_COMBOBOX, self.maxMenuSelect)
        self.maxMenu.Bind(wx.EVT_RIGHT_DOWN, self.rightClick)
        self.maxMenu.SetToolTipString("Maximum length of loop in residues")
        self.maxMenu.SetSelection(len(minmax_length)-1)
        self.maxLength = int(self.maxMenu.GetStringSelection())


        # Min number of models to evaluate

        if (platform.system() == "Windows"):
            self.lblResultsMin = wx.StaticText(self, -1, "Minumum results", (10, 240), (140, 20), wx.ALIGN_CENTRE)
            self.lblResultsMin.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        # elif (platform.system() == "Darwin"):
        #     self.lblResultsMin = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/indel/minResults.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(10, 240), size=(140, 20))
        else:
            self.lblResultsMin = wx.StaticText(self, -1, "Minimum results", (20, 240), style=wx.ALIGN_CENTRE)
            self.lblResultsMin.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblEnd, 170, 140)
        self.lblResultsMin.SetForegroundColour("#FFFFFF")
        self.ResultsMin = wx.lib.intctrl.IntCtrl(self, pos=(10, 260), size=(140, 25))
        self.ResultsMin.Bind(wx.EVT_COMBOBOX, self.maxMenuSelect)
        self.ResultsMin.Bind(wx.EVT_RIGHT_DOWN, self.rightClick)
        self.ResultsMin.SetToolTipString("Minimum number of models to try.")
        self.ResultsMin.SetValue(10)
        self.minResultsval = self.ResultsMin.GetValue()

        # Max number of models to evaluate

        if (platform.system() == "Windows"):
            self.lblResultsMax = wx.StaticText(self, -1, "Maximum results", (170, 240), (140, 20), wx.ALIGN_CENTRE)
            self.lblResultsMax.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        # elif (platform.system() == "Darwin"):
        #     self.lblResultsMax = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/indel/maxResults.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(10, 240), size=(140, 20))
        else:
            self.lblResultsMax = wx.StaticText(self, -1, "Maximum results", (180, 240), style=wx.ALIGN_CENTRE)
            self.lblResultsMax.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblEnd, 170, 140)
        self.lblResultsMax.SetForegroundColour("#FFFFFF")
        self.ResultsMax = wx.lib.intctrl.IntCtrl(self, pos=(170, 260), size=(140, 25))
        self.ResultsMax.Bind(wx.EVT_COMBOBOX, self.maxMenuSelect)
        self.ResultsMax.Bind(wx.EVT_RIGHT_DOWN, self.rightClick)
        self.ResultsMax.SetToolTipString("Maximum number of models to try.")
        self.ResultsMax.SetValue(25)
        self.maxResultsval = self.ResultsMax.GetValue()



        # if (platform.system() == "Darwin"):
        #     self.btnClear = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/kic/btnClear.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(220, 305), size=(90, 25))
        # else:
        self.btnClear = wx.Button(self, id=-1, label="Clear", pos=(220, 305), size=(90, 25))
        self.btnClear.SetForegroundColour("#000000")
        self.btnClear.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnClear.Bind(wx.EVT_BUTTON, self.clear)
        self.btnClear.SetToolTipString("Clear parameters")

        self.grdLoops = wx.grid.Grid(self)
        self.grdLoops.CreateGrid(0, 2)
        self.grdLoops.SetSize((320, 200))
        self.grdLoops.SetPosition((0, 350))
        self.grdLoops.SetLabelFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.grdLoops.DisableDragColSize()
        self.grdLoops.DisableDragRowSize()
        self.grdLoops.SetColLabelValue(0, "Length")
        self.grdLoops.SetColLabelValue(1, "Score")
        self.grdLoops.SetRowLabelSize(50)
        self.grdLoops.SetColSize(0, 70)
        self.grdLoops.SetColSize(1, 200)
        self.grdLoops.Bind(wx.grid.EVT_GRID_CELL_LEFT_CLICK, self.gridClick)
        self.loops = []
        self.selectedr = -1
        ypos = self.grdLoops.GetPosition()[1] + self.grdLoops.GetSize()[1] + 10
        self.model_selected = ""
        self.previous_model_selected = ""
        self.model_names = []
        self.scores = []
        self.lengths = []


        # if (platform.system() == "Darwin"):
        #     self.btnServerToggle = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/kic/btnServer_Off.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(40, ypos+215), size=(100, 25))
        # else:
        self.btnServerToggle = wx.Button(self, id=-1, label="Server Off", pos=(40, ypos), size=(100, 25))
        self.btnServerToggle.SetForegroundColour("#000000")
        self.btnServerToggle.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.btnServerToggle.Bind(wx.EVT_BUTTON, self.serverToggle)
        self.btnServerToggle.SetToolTipString("Perform KIC simulations locally")
        self.serverOn = False
        self.btnServerToggle.Disable()

        # if (platform.system() == "Darwin"):
        #     self.btnINDEL = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/indel/btnINDEL.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(180, ypos+215), size=(100, 25))
        # else:
        self.btnINDEL = wx.Button(self, id=-1, label="Model!", pos=(180, ypos), size=(100, 25))
        self.btnINDEL.SetForegroundColour("#000000")
        self.btnINDEL.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))

        self.btnINDEL.Bind(wx.EVT_BUTTON, self.INDELClick)
        self.btnINDEL.SetToolTipString("Begin INDEL simulation with selected parameters")
        self.buttonState = "Model!"

        self.scrollh = self.btnINDEL.GetPosition()[1] + self.btnINDEL.GetSize()[1] + 5
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
            browser.open(self.parent.parent.scriptdir + "/help/kic.html")
        else:
            webbrowser.open(self.parent.parent.scriptdir + "/help/kic.html")

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
        # Get the list of all the models in the sequence viewer
        modelList = []
        for r in range(0, self.seqWin.SeqViewer.NumberRows):
            model = self.seqWin.getModelForChain(r)
            if (not(model in modelList)):
                modelList.append(model)
        # Update the combobox list if the list has changed
        if (modelList != self.modelMenu.GetItems()):
            self.modelMenu.Clear()
            self.modelMenu.AppendItems(modelList)
            self.selectedModel = ""
            if (platform.system() == "Windows"):
                self.modelMenu.SetSelection(-1)
            else:
                self.modelMenu.SetSelection(0)
                self.modelMenuSelect(None)
            # Did we lose the model for the data in the loops grid?  If so, clear the loops
            if (len(self.loops) > 0 and not(self.loops[0][2] in modelList)):
                self.loops = []
                self.updateLoops()
        # If the user was deleting things in the sequence window, the specified begin and end positions might
        # not be valid anymore so we should erase them
        poseindx = self.seqWin.getPoseIndexForModel(self.selectedModel)
        if (poseindx >= 0):
            naa = 0
            for ch in self.seqWin.poses[poseindx][0]:
                for residue in ch:
                    if (residue.resname in "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR "):
                        naa = naa + 1
            if (len(self.beginMenu.GetItems()) != naa-1):
                self.selectedModel = ""
                self.modelMenuSelect(None)
        self.Scroll(0, self.winscrollpos)

    def rightClick(self, event):
        # Attempt to fill in loop values from a selection to bypass having to use the ComboBox
        try:
            topLefts = self.seqWin.SeqViewer.GetSelectionBlockTopLeft()
            bottomRights = self.seqWin.SeqViewer.GetSelectionBlockBottomRight()
            row = topLefts[0][0]
            begin = 9999999
            end = 0
            for i in range(0, len(topLefts)):
                for r in range(topLefts[i][0], bottomRights[i][0]+1):
                    if (r != row):
                        continue
                    for c in range(topLefts[i][1], bottomRights[i][1]+1):
                        if (c > end and self.seqWin.sequences[row][c] != "-"):
                            end = c
                        if (c < begin and self.seqWin.sequences[row][c] != "-"):
                            begin = c
            if (begin == end):
                # Have to get at least two residues
                return
            model = self.seqWin.IDs[row]
            chain = model[len(model)-1]
            model = model[:len(model)-2]
            beginres = chain + ":" + self.seqWin.sequences[row][begin] + str(self.seqWin.indxToSeqPos[row][begin][1])
            endres = chain + ":" + self.seqWin.sequences[row][end] + str(self.seqWin.indxToSeqPos[row][end][1])
            mindx = self.modelMenu.GetItems().index(model)
            bindx = self.beginMenu.GetItems().index(beginres)
            eindx = self.endMenu.GetItems().index(endres)
            self.modelMenu.SetSelection(mindx)
            self.beginMenu.SetSelection(bindx)
            self.endMenu.SetSelection(eindx)
            chain = self.beginMenu.GetStringSelection()[0]
            seqpos = self.beginMenu.GetStringSelection()[3:].strip()
            rindx = self.seqWin.getRosettaIndex(self.selectedModel, chain, seqpos)
            self.loopBegin = rindx
            chain = self.endMenu.GetStringSelection()[0]
            seqpos = self.endMenu.GetStringSelection()[3:].strip()
            rindx = self.seqWin.getRosettaIndex(self.selectedModel, chain, seqpos)
            self.loopEnd = rindx
            self.focusView(self.endMenu.GetStringSelection(), self.selectedModel)
            self.populatePivots()
        except:
            pass


    def gridClick(self, event):
        # Set the selected residue's row to blue so it is easy to see what the selection is
        self.selectedr = event.GetRow()
        if (self.selectedr >= self.grdLoops.NumberRows):
            self.selectedr = -1
        if (self.selectedr >= len(self.model_names)):
            event.Skip()
            return
        for r in range(0, self.grdLoops.NumberRows):
            if (r == self.selectedr):
                for c in range(0, self.grdLoops.NumberCols):
                    self.grdLoops.SetCellBackgroundColour(r, c, "light blue")
            else:
                for c in range(0, self.grdLoops.NumberCols):
                    self.grdLoops.SetCellBackgroundColour(r, c, "white")

        self.grdLoops.Refresh()

        # Make sure we're not trying to load an empty row
        if (self.selectedr < len(self.model_names)):
            self.loopEnd = self.begin_seqpos + self.lengths[self.selectedr]
            self.model_selected = self.model_names[self.selectedr]
        else:
            self.selectedr = -1
            event.Skip()

        # Remove the previously selected model from the viewer if it's not the same as the one just selected
        if (self.previous_model_selected != "" or self.previous_model_selected != self.model_selected):
            try:
                self.cmd.remove(self.previous_model_selected)
                self.cmd.delete(self.previous_model_selected)
            except:
                pass

        # Load the model, zoom in on designed loop
        if (self.model_selected != self.previous_model_selected and self.model_selected != ""):
            self.cmd.load(self.model_selected, self.model_selected)
            self.cmd.show()
            self.cmd.show("cartoon")
            self.cmd.hide("lines")
            self.cmd.hide("sticks")
            self.previous_model_selected = self.model_selected
            self.cmd.zoom("resi " + str(self.begin_seqpos) + "-" + str(int(self.begin_seqpos) + 3), 2.0)
            self.cmd.color("blue", self.model_selected)

        event.Skip()

    def modelMenuSelect(self, event):
        # Update the list of positions with the new model
        if (self.selectedModel == self.modelMenu.GetStringSelection()):
            return
        self.selectedModel = self.modelMenu.GetStringSelection()
        logInfo("Selected model " + self.selectedModel)
        # Get the location of the pose
        poseindx = self.seqWin.getPoseIndexForModel(self.selectedModel)
        # Read the positions
        pose = self.seqWin.poses[poseindx]
        positions = []
        for ch in pose[0]:
            for residue in ch:
                if ("ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR ".find(residue.resname) >= 0):
                    chain = ch.id
                    if (len(chain.strip()) == 0):
                        chain = "_"
                    label = chain + ":" + AA3to1(residue.resname) + str(residue.id[1])
                    positions.append(label)
        # Update the beginning and ending positions menus with the available sequence positions
        self.beginMenu.Clear()
        self.beginMenu.AppendItems(positions[0:len(positions)-1])

        if (platform.system() == "Windows"):
            self.beginMenu.SetSelection(-1)
            self.loopBegin = -1
        else:
            self.beginMenu.SetSelection(0)
            self.loopBegin = 1
        self.endMenu.Clear()
        self.endMenu.AppendItems(positions[1:])
        if (platform.system() == "Windows"):
            self.endMenu.SetSelection(-1)
            self.loopEnd = -1
        else:
            self.endMenu.SetSelection(0)
            self.loopEnd = 2
        #self.txtNStruct.Enable()
        #self.populatePivots()

    def changeLoopType(self, event):
        if (self.loopType == "Refine"):
            self.loopType = "Reconstruct"
            # if (platform.system() == "Darwin"):
            #     self.btnLoopType.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/kic/btnLoopType_Reconstruct.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            # else:
            self.btnLoopType.SetLabel(self.loopType)
            self.btnLoopType.SetToolTipString("Reconstruct the current loop using the wildtype sequence")
            self.btnPerturb.Enable()
            self.txtNStruct.Enable()
        elif (self.loopType == "Reconstruct"):
            self.loopType = "De Novo"
            # if (platform.system() == "Darwin"):
            #     self.btnLoopType.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/kic/btnLoopType_DeNovo.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            # else:
            self.btnLoopType.SetLabel(self.loopType)
            self.btnLoopType.SetToolTipString("Construct a new loop with a new sequence")
            self.txtSequence.Enable()
        else:
            self.loopType = "Refine"
            # if (platform.system() == "Darwin"):
            #     self.btnLoopType.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/kic/btnLoopType_Refine.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            # else:
            self.btnLoopType.SetLabel(self.loopType)
            self.btnLoopType.SetToolTipString("Refine a pre-existing loop using the high resolution KIC remodeler only")
            self.txtSequence.Disable()
            self.btnPerturb.Disable()
            self.txtNStruct.Disable()
        logInfo("Changed loop type to " + self.loopType)

    def changePerturbType(self, event):
        if (self.perturbType == "Perturb+Refine"):
            self.perturbType = "Perturb Only, Fullatom"
            # if (platform.system() == "Darwin"):
            #     self.btnPerturb.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/kic/btnPerturb_Fullatom.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            # else:
            self.btnPerturb.SetLabel(self.perturbType)
            self.btnPerturb.SetToolTipString("Perform only KIC coarse perturbations but convert outputted models to repacked fullatom PDBs")
        #elif (self.perturbType == "Perturb Only, Fullatom"):
        #    self.perturbType = "Perturb Only, Centroid"
        #    self.btnPerturb.SetToolTipString("Perform only KIC coarse perturbations and leave outputted PDBs in coarse centroid mode")
        else:
            self.perturbType = "Perturb+Refine"
            # if (platform.system() == "Darwin"):
            #     self.btnPerturb.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/kic/btnPerturb_Refine.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            # else:
            self.btnPerturb.SetLabel(self.perturbType)
            self.btnPerturb.SetToolTipString("Perform KIC coarse perturbation followed by high resolution refinement")
        logInfo("Changed perturbation type to " + self.perturbType)

    def setOutputDir(self, event):
        logInfo("Clicked Output Dir button")
        dlg = wx.DirDialog(
            self, message="Choose a directory",
            defaultPath=self.seqWin.cwd,
            style=wx.DD_DEFAULT_STYLE | wx.DD_DIR_MUST_EXIST)
        if (dlg.ShowModal() == wx.ID_OK):
            path = dlg.GetPath()
            self.outputdir = str(path)
            # Change cwd to the last opened file
            self.seqWin.cwd = self.outputdir
            self.seqWin.saveWindowData(None)
            self.lblDir.SetLabel(self.outputdir)
            self.lblDir.SetForegroundColour("#FFFFFF")
            if (platform.system() == "Linux"):
                resizeTextControlForUNIX(self.lblDir, 130, 190)
            logInfo("Set output directory as " + self.outputdir)
        else:
            logInfo("Cancelled out of Load PDB")

    def populatePivots(self):
        self.menuPivot.Enable()
        # Get the location of the pose
        poseindx = self.seqWin.getPoseIndexForModel(self.selectedModel)
        # Read the positions
        pose = self.seqWin.poses[poseindx]
        positions = []
        ires = 1
        for ch in pose[0]:
            for residue in ch:
                if (ires >= self.loopBegin and ires <= self.loopEnd):
                    if ("ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR ".find(residue.resname) >= 0):
                        chain = ch.id
                        if (len(chain.strip()) == 0):
                            chain = "_"
                        label = chain + ":" + AA3to1(residue.resname) + str(residue.id[1])
                        positions.append(label)
                ires = ires + 1
        self.menuPivot.Clear()
        self.menuPivot.AppendItems(positions)
        self.menuPivot.SetSelection(0)

    def beginMenuSelect(self, event):
        try:
            chain = self.beginMenu.GetStringSelection()[0]
            seqpos = self.beginMenu.GetStringSelection()[3:].strip()
            rindx = self.seqWin.getRosettaIndex(self.selectedModel, chain, seqpos)
            self.loopBegin = rindx
            # If this new loop begin is further down than what is set for loop end, then it needs
            # to be reset and the user should be notified
            if (self.loopEnd >= 0 and self.loopEnd <= rindx):
                if (platform.system() == "Windows"):
                    self.endMenu.SetSelection(-1)
                    self.loopEnd = -1
                else:
                    self.endMenu.SetSelection(self.beginMenu.GetSelection()) # This clears the menu, SetStringSelection/SetValue doesn't seem to work
                    self.endMenuSelect(event)
                #wx.MessageBox("Your selected end loop value is no longer valid.  Please choose an ending position after the one you've selected here.", "Loop End No Longer Valid", wx.OK|wx.ICON_EXCLAMATION)
            self.focusView(self.beginMenu.GetStringSelection(), self.selectedModel)
            logInfo("Selected " + self.beginMenu.GetStringSelection() + " as the beginning of the loop")
        except:
            # Probably the user left the field blank, do nothing
            pass

    def endMenuSelect(self, event):
        try:
            chain = self.endMenu.GetStringSelection()[0]
            seqpos = self.endMenu.GetStringSelection()[3:].strip()
            rindx = self.seqWin.getRosettaIndex(self.selectedModel, chain, seqpos)
            self.loopEnd = rindx
            # If this new loop begin is further up than what is set for loop begin, then it needs
            # to be reset and the user should be notified
            if (self.loopBegin >= 0 and self.loopBegin >= rindx):
                if (platform.system() == "Windows"):
                    self.beginMenu.SetSelection(-1)
                    self.loopBegin = -1
                else:
                    self.beginMenu.SetSelection(self.endMenu.GetSelection()) # This clears the menu, SetStringSelection/SetValue doesn't seem to work
                    self.beginMenuSelect(event)
                wx.MessageBox("Your selected begin loop value is no longer valid.  Please choose a beginning position before the one you've selected here.", "Loop Begin No Longer Valid", wx.OK|wx.ICON_EXCLAMATION)

            self.focusView(self.endMenu.GetStringSelection(), self.selectedModel)
            logInfo("Selected " + self.endMenu.GetStringSelection() + " as the ending of the loop")
        except:
            # Probably the user left the field blank, do nothing
            pass

    def minMenuSelect(self, event):
        self.minLength = int(self.minMenu.GetStringSelection())

    def maxMenuSelect(self,event):
        self.maxLength = int(self.maxMenu.GetStringSelection())

    def updateLoops(self):
        # Redraw the loops grid with current loop information
        scrollpos = self.grdLoops.GetScrollPos(wx.VERTICAL)
        if (self.grdLoops.NumberRows > 0):
            self.grdLoops.DeleteRows(0, self.grdLoops.NumberRows)
        if (len(self.loops) > 0):
            self.grdLoops.AppendRows(len(self.loops))
        row = 0
        for [loopType, sequence, model, begin, pivot, end] in self.loops:
            self.grdLoops.SetRowLabelValue(row, loopType)
            self.grdLoops.SetCellValue(row, 0, sequence)
            chainID, resindx = self.seqWin.getResidueInfo(model, begin)
            if (len(chainID.strip()) == 0):
                chainID = "_"
            self.grdLoops.SetCellValue(row, 1, chainID + "|" + self.seqWin.getResidueTypeFromRosettaIndx(model, begin) + str(resindx))
            chainID, resindx = self.seqWin.getResidueInfo(model, pivot)
            if (len(chainID.strip()) == 0):
                chainID = "_"
            self.grdLoops.SetCellValue(row, 2, chainID + "|" + self.seqWin.getResidueTypeFromRosettaIndx(model, pivot) + str(resindx))
            chainID, resindx = self.seqWin.getResidueInfo(model, end)
            if (len(chainID.strip()) == 0):
                chainID = "_"
            self.grdLoops.SetCellValue(row, 3, chainID + "|" + self.seqWin.getResidueTypeFromRosettaIndx(model, end) + str(resindx))
            readOnly = wx.grid.GridCellAttr()
            readOnly.SetReadOnly(True)
            readOnly.SetAlignment(wx.ALIGN_CENTRE, wx.ALIGN_CENTRE)
            readOnly.SetBackgroundColour("#FFFFFF")
            self.grdLoops.SetRowAttr(row, readOnly)
            row += 1
        self.grdLoops.Scroll(0, scrollpos)

    def add(self, event):
        # Is the loop valid?
        if (self.loopBegin < 0 or self.loopBegin < 0 or self.loopBegin >= self.loopEnd):
            dlg = wx.MessageDialog(self, "You do not have a valid loop specified!", "Loop Not Valid", wx.OK | wx.ICON_ERROR | wx.CENTRE)
            dlg.ShowModal()
            dlg.Destroy()
            return
        # If we're doing a de novo search, is the sequence specified?
        if (self.loopType == "De Novo"):
            sequence = self.txtSequence.GetValue().strip().upper()
            for AA in sequence:
                if (not(AA in "ACDEFGHIKLMNPQRSTVWY")):
                    wx.MessageBox("The sequence you have provided is invalid.  Please only use canonical amino acids.", "Sequence Invalid", wx.OK|wx.ICON_EXCLAMATION)
                    return
            if (len(sequence) == 0):
                wx.MessageBox("You have indicated that you want to design a loop de novo but have not provided the putative sequence of the loop.  Please provide one or switch to use a pre-existing loop.", "No Sequence Indicated", wx.OK|wx.ICON_EXCLAMATION)
                return
        else:
            sequence = ""
        # Did the model change?  If yes, and loops is not empty, then tell the user that this
        # will remove all loops to make room for the new model
        if (len(self.loops) > 0 and self.modelMenu.GetValue() != self.loops[0][2]):
            dlg = wx.MessageDialog(self, "You are attempting to add a loop for a different model.  If you continue, all current loops will be removed.  Is this okay?", "Loop Model Changed", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
            if (dlg.ShowModal() == wx.ID_NO):
                return
            dlg.Destroy()
            self.loops = []
        # Does this loop overlap with a previously-specified loop?  If so, do not add
        i = 1
        for loopType, s, model, begin, pivot, end in self.loops:
            if ((self.loopBegin >= begin and self.loopBegin <= end) or (self.loopEnd >= begin and self.loopEnd <= end)):
                dlg = wx.MessageDialog(self, "The loop you have indicated overlaps with loop " + str(i) + ".  Either change the current loop or remove loop " + str(i) + ".", "Loop Overlap", wx.OK | wx.ICON_ERROR | wx.CENTRE)
                dlg.ShowModal()
                dlg.Destroy()
                return
            i += 1
        # Add this loop to the list of loops currently active
        self.loops.append([self.loopType, sequence, self.modelMenu.GetValue(), self.loopBegin, self.menuPivot.GetSelection() + self.loopBegin, self.loopEnd])
        self.updateLoops()

    def remove(self, event):
        # For this function, remove the indicated loop
        self.activate()
        logInfo("Remove button clicked")
        if (self.selectedr >= 0 and self.selectedr < len(self.loops)):
            self.loops.pop(self.selectedr)
            self.selectedr = -1
        self.updateLoops()

    def clear(self, event):
        logInfo("Clear button clicked")
        # Remove everything
        self.loops = []
        self.updateLoops()

    def viewMenuSelect(self, event):
        try:
            self.focusView(self.viewMenu.GetStringSelection(), self.selectedModel, "kic_view")
            logInfo("Viewing " + self.viewMenu.GetStringSelection())
        except:
            # Probably the user left the field blank, do nothing
            pass

    def focusView(self, posID, origmodel, newmodel=None):
        model = origmodel
        loopEnd = self.loopEnd
        if (posID != "Whole Loop"):
            chain = posID[0]
            seqpos = posID[3:].strip()
            # Loop end needs to be recalculated if this is a view of the de novo loop since the
            # de novo loop may be a different size
            if (newmodel and len(self.txtSequence.GetValue()) > 0):
                loopEnd = self.loopBegin + len(self.txtSequence.GetValue()) + 1 # For the anchor
        else:
            i = 1
            wholeloop_data = []
            for ch in self.KICView[0]:
                for residue in ch:
                    if (i >= self.loopBegin and i <= loopEnd):
                        chain = ch.id
                        seqpos = str(residue.id[1])
                        wholeloop_data.append((chain, seqpos))
                    i = i + 1
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
        # If the loop is validly defined, let's show the whole loop instead of individual residues
        if ((self.loopBegin >= 0 and self.loopEnd >= 0 and not(newmodel)) or posID == "Whole Loop"):
            for i in range(self.loopBegin, loopEnd):
                if (not(newmodel)):
                    (chain, seqpos) = self.seqWin.getResidueInfo(self.selectedModel, i)
                else:
                    (chain, seqpos) = wholeloop_data[i-self.loopBegin]
                if (chain == "_" or len(chain.strip()) == 0):
                    self.cmd.select("viewsele", "viewsele or (resi " + str(seqpos) + " and model " + firstmodel + ")")
                else:
                    self.cmd.select("viewsele", "viewsele or (resi " + str(seqpos) + " and chain " + chain + " and model " + firstmodel + ")")
        self.cmd.select("exviewsele", "model " + firstmodel + " within 12 of viewsele")
        self.cmd.show("cartoon", "exviewsele")
        self.cmd.hide("ribbon", "exviewsele")
        self.cmd.show("sticks", "exviewsele")
        self.cmd.set_bond("stick_radius", 0.1, "exviewsele")
        # Display energy labels for new structures
        if (newmodel):
            relabelEnergies(self.KICView, self.residue_E, newmodel, self.scoretypeMenu.GetStringSelection(), self.cmd, seqpos)
            self.cmd.label("not exviewsele", "")
        self.cmd.zoom("exviewsele")
        #if (chain == " " or chain == "_"):
        #    self.cmd.select("viewsele", "resi " + seqpos + " and model " + firstmodel)
        #else:
        #    self.cmd.select("viewsele", "resi " + seqpos + " and model " + firstmodel + " and chain " + chain)
        self.cmd.show("sticks", "viewsele")
        self.cmd.set_bond("stick_radius", 0.25, "viewsele")
        # Highlight this residue in PyMOL
        self.cmd.select("seqsele", "viewsele")
        if (newmodel):
            # If this is after a protocol, also show the original structure in green for comparison
            self.cmd.select("oldsele", "model " + origmodel + " and symbol c")
            self.cmd.color("green", "oldsele")
            self.cmd.set("cartoon_color", "green", "oldsele")
            #if (chain == " " or chain == "_"):
                #self.cmd.select("viewsele", "resi " + seqpos + " and model " + origmodel)
            #else:
                #self.cmd.select("viewsele", "resi " + seqpos + " and model " + origmodel + " and chain " + chain)
            #self.cmd.select("viewsele", "model " + origmodel + " within 12 of viewsele")
            self.cmd.select("exviewsele", "model " + origmodel + " within 12 of viewsele")
            self.cmd.show("cartoon", "exviewsele")
            self.cmd.hide("ribbon", "exviewsele")
            self.cmd.show("sticks", "exviewsele")
            self.cmd.set_bond("stick_radius", 0.1, "exviewsele")
            self.cmd.zoom("exviewsele")
            self.cmd.delete("oldsele")
            #if (chain == " " or chain == "_"):
                #self.cmd.select("exviewsele", "resi " + seqpos + " and model " + origmodel)
            #else:
                #self.cmd.select("viewsele", "resi " + seqpos + " and model " + origmodel + " and chain " + chain)
            #self.cmd.show("sticks", "viewsele")
            #self.cmd.set_bond("stick_radius", 0.25, "viewsele")
        self.cmd.enable("seqsele")
        self.cmd.delete("viewsele")
        self.cmd.select("exviewsele", "solvent")
        self.cmd.hide("everything", "exviewsele")
        self.cmd.delete("exviewsele")
        self.seqWin.selectUpdate(False)

    def scoretypeMenuSelect(self, event):
        # Make sure there is even a PyMOL_Mover pose loaded
        if (self.selectedModel == ""):
            return
        logInfo("Changed scoretype view to " + self.scoretypeMenu.GetStringSelection())
        recolorEnergies(self.KICView, self.residue_E, "kic_view", self.scoretypeMenu.GetStringSelection(), self.cmd)
        self.viewMenuSelect(event) # To update all the labels

    def serverToggle(self, event):
        if (self.serverOn):
            self.serverOn = False
            # if (platform.system() == "Darwin"):
            #     self.btnServerToggle.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/kic/btnServer_Off.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            # else:
            self.btnServerToggle.SetLabel("Server Off")
            self.btnServerToggle.SetToolTipString("Perform KIC simulations locally")
            logInfo("Turned off KIC server usage")
        else:
            self.serverOn = True
            # if (platform.system() == "Darwin"):
            #     self.btnServerToggle.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/kic/btnServer_On.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            # else:
            self.btnServerToggle.SetLabel("Server On")
            self.btnServerToggle.SetToolTipString("Perform KIC simulations on a remote server")
            logInfo("Turned on KIC server usage")

    def cancelINDEL(self):
        logInfo("Canceled INDEL operation")
        try:
            os.remove("INDELinput")
        except:
            pass
        try:
            os.remove("coarsekicinputtemp")
        except:
            pass
        try:
            os.remove("repacked.pdb")
        except:
            pass
        try:
            os.remove("finekicinput")
        except:
            pass
        self.tmrKIC.Stop()
        self.seqWin.cannotDelete = False
        #self.scoretypeMenu.Disable()
        #self.viewMenu.Disable()
        self.modelMenu.Enable()
        self.beginMenu.Enable()
        self.endMenu.Enable()
        self.minMenu.Enable()
        self.maxMenu.Enable()
        self.ResultsMin.Enable()
        self.ResultsMax.Enable()
        self.btnClear.Enable()

        #self.btnLoopType.Enable()
        #if (self.loopType == "De Novo"):
        #    self.txtSequence.Enable()
        # if (platform.system() == "Darwin"):
        #     self.btnINDEL.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/kic/btnKIC.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
        # else:
        self.btnINDEL.SetLabel("Model!")
        self.buttonState = "Model!"
        self.btnINDEL.SetToolTipString("Perform INDEL simulation with selected parameters")
        deleteInputFiles()
        self.parent.parent.restartDaemon()
        self.parent.GoBtn.Enable()
        # Get rid of the messages
        for i in range(0, len(self.seqWin.msgQueue)):
            if (self.seqWin.msgQueue[i].find("Performing INDEL loop modeling, please be patient...") >= 0):
                self.seqWin.msgQueue.pop(i)
                break
        for i in range(0, len(self.seqWin.msgQueue)):
            if (self.seqWin.msgQueue[i].find("Performing rotamer repacking") >= 0):
                self.seqWin.msgQueue.pop(i)
                break
        for i in range(0, len(self.seqWin.msgQueue)):
            if (self.seqWin.msgQueue[i].find("Performing refined KIC loop modeling") >= 0):
                self.seqWin.msgQueue.pop(i)
                break
        if (len(self.seqWin.msgQueue) > 0):
            self.seqWin.labelMsg.SetLabel(self.seqWin.msgQueue[len(self.seqWin.msgQueue)-1])
        else:
            self.seqWin.labelMsg.SetLabel("")
        self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")

    def INDELClick(self, event):
        # This is also the "Finalize!" button
        if (self.buttonState == "Model!"):
            # Some checking to make sure input parameters make sense
            # TODO make sure that the two end residues aren't selected. AnchoredGraftMover doesn't like grafting at terminal ends
            if (self.minLength > self.maxLength):
                wx.MessageBox("Please choose a maximum length that is greater than or equal to the minimum length.", "Invalid loop lengths" , wx.OK|wx.ICON_EXCLAMATION)
                return
            if (self.maxResultsval < self.minResultsval):
                wx.MessageBox("Please enter a maximum results value that is higher than the minimum results value.", "Invalid parameter", wx.OK|wx.ICON_EXCLAMATION)
            if (self.minResultsval <= 0):
                self.minResultsval = 1
            if (self.maxResultsval <= 0):
                wx.MessageBox("Please enter a maximum results value that is greater than or equal to 1.", "Invalid parameter", wx.OK|wx.ICON_EXCLAMATION)



            self.seqWin.labelMsg.SetLabel("Performing INDEL loop modeling, please be patient...")
            self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
            self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")
            self.seqWin.msgQueue.append("Performing INDEL loop modeling, please be patient...")
            self.seqWin.cannotDelete = True
            self.parent.GoBtn.Disable()
            self.modelMenu.Disable()
            self.beginMenu.Disable()
            self.endMenu.Disable()
            self.minMenu.Disable()
            self.maxMenu.Disable()
            self.ResultsMin.Disable()
            self.ResultsMax.Disable()
            self.btnClear.Disable()



            # if (platform.system() == "Darwin"):
            #     self.btnINDEL.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/kic/btnKIC_Cancel.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            # else:
            self.btnINDEL.SetLabel("Cancel!")
            self.buttonState = "Cancel!"
            self.btnINDEL.SetToolTipString("Cancel the INDEL simulation")
            self.stage = 1
            logInfo("Clicked the INDEL button")
            self.tmrKIC = wx.Timer(self)
            self.Bind(wx.EVT_TIMER, self.threadINDEL, self.tmrKIC)
            self.tmrKIC.Start(1000)

        elif (self.buttonState == "Cancel!"):
            dlg = wx.MessageDialog(self, "Are you sure you want to cancel the INDEL simulation?  All progress will be lost.", "Cancel KIC Simulation", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
            result = dlg.ShowModal()
            if (result == wx.ID_YES):
                self.cancelINDEL()
            dlg.Destroy()
        else:

            # Finalize button, ask whether the changes will be accepted or rejected
            dlg = wx.MessageDialog(self, "Do you want to accept the results of this loop modeling session?", "Accept/Reject Model", wx.YES_NO | wx.CANCEL | wx.ICON_QUESTION | wx.CENTRE)
            result = dlg.ShowModal()
            if (result == wx.ID_YES):
                logInfo("Accepted KIC model")
                accept = True
            elif (result == wx.ID_NO):
                logInfo("Rejected KIC model")
                accept = False
            else:
                logInfo("Cancelled Finalize operation")
                dlg.Destroy()
                return

            # Try to get rid of working loop files in sandbox
            temp_loop_files = glob.glob('loopout_*')
            try:
                for temp_loop in temp_loop_files:
                    os.remove(temp_loop)
            except:
                pass

            # Clear grid of loops
            self.grdLoops.ClearGrid()
            self.grdLoops.DeleteRows(0, self.grdLoops.GetNumberRows())

            # Keep track of the temporary name so we can remove it from the pymol window in a second
            # Rename the chosen model to its final name
            pymol_model_selected = self.model_selected
            os.rename(self.model_selected, self.selectedModel + "_INDEL.pdb")
            self.model_selected = self.selectedModel + "_INDEL.pdb"

            # Try to get rid of the models not chosen
            for i in range(len(self.model_names)):
                try:
                    os.remove(self.model_names[i])
                except:
                    pass

            # Clear internal list of model data
            del self.model_names[:]
            del self.lengths[:]
            del self.scores[:]

            # Re-enable controls
            dlg.Destroy()
            self.modelMenu.Enable()
            self.beginMenu.Enable()
            self.endMenu.Enable()
            self.parent.GoBtn.Enable()
            self.minMenu.Enable()
            self.maxMenu.Enable()
            self.btnClear.Enable()
            self.ResultsMin.Enable()
            self.ResultsMax.Enable()


            #Pop message out of queue
            for i in range(0, len(self.seqWin.msgQueue)):
                if (self.seqWin.msgQueue[i].find("Performing INDEL loop modeling, please be patient...") >= 0):
                    self.seqWin.msgQueue.pop(i)
                    break
            self.seqWin.labelMsg.SetLabel("")



            # if (platform.system() == "Darwin"):
            #     self.btnINDEL.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/kic/btnKIC.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            # else:
            self.btnINDEL.SetLabel("Model!")
            self.buttonState = "Model!"
            self.btnINDEL.SetToolTipString("Perform INDEL simulation with selected parameters")
            self.cmd.label("all", "")
            self.seqWin.cannotDelete = False
            if (not(accept)):
                try:
                    self.cmd.remove(self.model_selected)
                    self.cmd.delete(self.model_selected)
                except:
                    pass
                return
            if (accept and self.selectedr == -1):
                return

            # Get rid of the original pose, save the designed pose, and reload the structure in PyMOL
            poseindx = -1
            for r in range(0, len(self.seqWin.IDs)):
                if (self.seqWin.IDs[r].find(self.selectedModel) >= 0):
                    poseindx = r
                    break

            try:
                self.cmd.load(self.model_selected, self.model_selected)
                # Color final model by ss
                defaultPyMOLView(self.cmd, self.model_selected)
                self.cmd.color('white')
                self.cmd.color('red', 'ss h')
                self.cmd.color('yellow', 'ss s')

                self.cmd.remove(self.selectedModel)
                self.cmd.delete(self.selectedModel)
                self.cmd.remove(pymol_model_selected)
                self.cmd.delete(pymol_model_selected)


                self.seqWin.reloadPose(poseindx, self.model_selected, self.model_selected)
                # IMPORTANT: You have to replace the model in the sandbox with the new designed model
                os.remove(self.selectedModel + ".pdb")


                self.selectedModel = self.model_selected


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
        # TODO re-enable the rest of the controls
        self.seqWin.cannotDelete = False
        self.parent.GoBtn.Enable()
        self.modelMenu.Enable()
        self.btnLoopType.Enable()
        self.beginMenu.Enable()
        self.endMenu.Enable()
        self.txtSequence.Enable()
        self.btnINDEL.Enable()
        # if (platform.system() == "Darwin"):
        #     self.btnINDEL.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/kic/btnKIC.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
        # else:
        self.btnINDEL.SetLabel("KIC!")
        self.buttonState = "KIC!"
        # Get rid of the messages
        for i in range(0, len(self.seqWin.msgQueue)):
            if (self.seqWin.msgQueue[i].find("Performing INDEL loop modeling, please be patient...") >= 0):
                self.seqWin.msgQueue.pop(i)
                break

        if (len(self.seqWin.msgQueue) > 0):
            self.seqWin.labelMsg.SetLabel(self.seqWin.msgQueue[len(self.seqWin.msgQueue)-1])
        else:
            self.seqWin.labelMsg.SetLabel("")
        self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")

    def threadINDEL(self, event):
        # Why am I doing this ridiculous timer thing for this KIC protocol?
        # Because apparently on Linux there's some kind of weird bug that manifests when you
        # attempt to run time.sleep loops looking for files to be generated
        # Pango develops a phobia of periods in strings if you do that????
        # Using this staged timer setup eliminates the error
        # What is the problem?  I don't know.  Why does this fix it?  I don't know
        # The people on StackOverflow said to do it and it fixed it -_-
        # I think it has something to do with Linux not liking things like "time.sleep"
        # and calls to wx in threads
        # Dump a file with the loop modeling parameters for the daemon to pick up
        goToSandbox()
        if (self.stage == 1):
            self.tmrKIC.Stop()
            self.timeoutCount = 0
            #self.nstruct = int(self.txtNStruct.GetValue())
            f = open("INDELinputtemp", "w")
            pdbfile = self.selectedModel + ".pdb"
            # Dump the PDB from PyMOL first in case the coordinates were altered by the user
            self.cmd.save(pdbfile.strip(), "model " + self.selectedModel)
            fixPyMOLSave(pdbfile.strip())
            f.write("PDBFILE\t" + pdbfile.strip() + "\n")
            f2 = open(pdbfile, "r")
            f.write("BEGIN PDB DATA\n")
            for aline in f2:
                f.write(aline.strip() + "\n")
            f.write("END PDB DATA\n")
            f2.close()
            #Grab rosetta indices, write indices and length range
            chain = self.endMenu.GetStringSelection()[0]
            begin_seqpos = self.beginMenu.GetStringSelection()[3:]
            self.begin_seqpos = begin_seqpos
            end_seqpos = self.endMenu.GetStringSelection()[3:]
            begin_index = self.seqWin.getRosettaIndex(self.selectedModel, chain, begin_seqpos)
            end_index = self.seqWin.getRosettaIndex(self.selectedModel, chain, end_seqpos)
            self.maxResultsval = self.ResultsMax.GetValue()
            self.minResultsval = self.ResultsMin.GetValue()
            f.write("LOOP\t" + str(begin_index)  + "\t" + str(end_index) + "\t" + str(self.minLength) + "\t" + str(self.maxLength) + "\t"
                + str(self.minResultsval) + "\t" + str(self.maxResultsval) )
            f.close()


            if (self.serverOn):
                try:
                    self.ID = sendToServer("coarsekicinput")
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
                            if (len(aline.split("\t")) >= 2 and aline.split("\t")[0] == "KIC" and aline.split("\t")[1] == self.ID.strip()):
                                alreadythere = True
                                break
                        f.close()
                    except:
                        pass
                    if (not(alreadythere)):
                        f = open("downloadwatch", "a")
                        f.write("KIC\t" + self.ID.strip() + "\t" + str(datetime.datetime.now().strftime("%A, %B %d - %I:%M:%S %p")) + "\t" + getServerName() + "\t" + desc + "\n")
                        f.close()
                    dlg = wx.MessageDialog(self, "InteractiveROSETTA is now watching the server for job ID " + desc.strip() + ".  You will be notified when the package is available for download.", "Listening for Download", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
                    dlg.ShowModal()
                    dlg.Destroy()
                    # Re-enable everything since we're not waiting for the local daemon to do anything
                    self.scoretypeMenu.Disable()
                    self.viewMenu.Disable()
                    self.modelMenu.Enable()
                    self.beginMenu.Enable()
                    self.endMenu.Enable()
                    self.btnLoopType.Enable()
                    if (self.loopType == "De Novo"):
                        self.txtSequence.Enable()
                    # if (platform.system() == "Darwin"):
                    #     self.btnINDEL.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/kic/btnKIC.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                    # else:
                    self.btnINDEL.SetLabel("KIC!")
                    self.buttonState = "KIC!"
                    self.btnINDEL.SetToolTipString("Perform KIC simulation with selected parameters")
                    self.cmd.label("all", "")
                    self.seqWin.cannotDelete = False
                    self.parent.GoBtn.Enable()
                    # Pop this message out of the queue
                    for i in range(0, len(self.seqWin.msgQueue)):
                        if (self.seqWin.msgQueue[i].find("Performing INDEL loop modeling, please be patient...") >= 0):
                            self.seqWin.msgQueue.pop(i)
                            break
                    if (len(self.seqWin.msgQueue) > 0):
                        self.seqWin.labelMsg.SetLabel(self.seqWin.msgQueue[len(self.seqWin.msgQueue)-1])
                    else:
                        self.seqWin.labelMsg.SetLabel("")
                    logInfo("Coarse KIC input sent to server daemon with ID " + self.ID)
                    return
                except:
                    dlg = wx.MessageDialog(self, "The server could not be reached!  Ensure that you have specified a valid server and that you have an network connection.", "Server Could Not Be Reached", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
                    dlg.ShowModal()
                    dlg.Destroy()
                    return
            else:
                os.rename("INDELinputtemp", "INDELinput")
                self.usingServer = False
                logInfo("INDEL input uploaded locally at INDELinput")
                self.stage = 2
            #if (self.perturbType == "Perturb Only, Centroid"):# or self.loopType == "Refine"):
            #    self.stage = 4
            self.looptimecount = 0
            self.timeout = 18000000
            #self.progress = wx.ProgressDialog("KIC Progress", "Modeling loops in centroid mode...", 100, style=wx.PD_CAN_ABORT | wx.PD_APP_MODAL | wx.PD_ELAPSED_TIME | wx.PD_REMAINING_TIME)
            self.loop_indx = 0
            self.last_progress_indx = 99
            self.tmrKIC.Start(1000)
        elif (self.stage == 2):
            if (os.path.isfile("INDELoutput")):
                self.tmrKIC.Stop()
                # Pop this message out of the queue
                for i in range(0, len(self.seqWin.msgQueue)):
                    if (self.seqWin.msgQueue[i].find("Performing INDEL loop modeling, please be patient...") >= 0):
                        self.seqWin.msgQueue.pop(i)
                        break


                # Parse output file that gives us the filenames, energies, and insertion lengths of all the results
                f = open("INDELoutput")
                self.model_names = []
                self.scores = []
                self.lengths = []
                for line in f:
                    tmp = line.split("\t")
                    self.model_names.append(tmp[0])
                    self.scores.append(tmp[1])
                    self.lengths.append(tmp[2])
                f.close()

                # Clear and populate table
                self.grdLoops.ClearGrid()
                self.grdLoops.AppendRows(len(self.model_names))

                row = 0
                for score, length in zip(self.scores, self.lengths):
                    self.grdLoops.SetCellValue(row, 0, length)
                    self.grdLoops.SetCellValue(row, 1, score)
                    row += 1

                self.btnClear.Disable()
                # We can get rid of the top-level output file now
                try:
                    os.remove("INDELoutput")
                except:
                    pass


                #self.KICView = self.seqWin.pdbreader.get_structure("kic_view", "INDELoutput.pdb")
                self.btnINDEL.Enable()
                #self.enableControls()
                #self.selectedModel = ""
                # if (platform.system() == "Darwin"):
                #     self.btnINDEL.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/kic/btnKIC_Finalize.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                # else:
                self.btnINDEL.SetLabel("Finalize!")
                self.buttonState = "Finalize!"
                self.btnINDEL.SetToolTipString("Accept or reject protocol results")
                #os.remove("INDELoutput.pdb")


            elif (os.path.isfile("errreport")):
                # Something went wrong, tell the user about it (loop sequence probably too short)
                self.tmrKIC.Stop()
                self.parent.parent.restartDaemon() # Has to happen because coarse KIC is threaded
                self.recoverFromError()

            self.looptimecount = self.looptimecount + 1
            if (self.looptimecount > self.timeout):
                # The loop was probably too short and coarse KIC will run forever
                # Kill the daemon and tell the user about it
                self.tmrKIC.Stop()
                # First delete that input file so the new daemon doesn't pick it up right away
                try:
                    os.remove("INDELinput")
                except:
                    pass
                self.parent.parent.restartDaemon() # Has to happen because coarse KIC is threaded
                self.recoverFromError("ERROR: The loop sequence is too short and cannot bridge the endpoint residues!")

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

class PointMutationsPanel(wx.lib.scrolledpanel.ScrolledPanel):
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
            self.lblProt = wx.StaticText(self, -1, "Point Mutations", (25, 15), (270, 25), wx.ALIGN_CENTRE)
            self.lblProt.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblProt = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/pointmutations/lblPointMutations.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(25, 15), size=(270, 25))
        else:
            self.lblProt = wx.StaticText(self, -1, "Point Mutations", (70, 15), style=wx.ALIGN_CENTRE)
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
            self.lblInst = wx.StaticText(self, -1, "Select a model, a position, and a new residue\ntype to manually select optimal rotamers.", (0, 45), (320, 25), wx.ALIGN_CENTRE)
            self.lblInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
        elif (platform.system() == "Darwin"):
            self.lblInst = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/pointmutations/lblInstPointMutations.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 45), size=(270, 25))
        else:
            self.lblInst = wx.StaticText(self, -1, "Select a model, a position, and a new residue\ntype to manually select optimal rotamers.", (5, 45), style=wx.ALIGN_CENTRE)
            self.lblInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
            resizeTextControlForUNIX(self.lblInst, 0, self.GetSize()[0])
        self.lblInst.SetForegroundColour("#FFFFFF")
        
        if (platform.system() == "Windows"):
            self.lblModel = wx.StaticText(self, -1, "Model", (20, 90), (140, 20), wx.ALIGN_CENTRE)
            self.lblModel.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblModel = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/pointmutations/lblModel.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(20, 90), size=(140, 20))
        else:
            self.lblModel = wx.StaticText(self, -1, "Model", (20, 90), style=wx.ALIGN_CENTRE)
            self.lblModel.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblModel, 20, 140)
        self.lblModel.SetForegroundColour("#FFFFFF")
        self.modelMenu = wx.ComboBox(self, pos=(20, 110), size=(140, 25), choices=[], style=wx.CB_READONLY)
        self.modelMenu.Bind(wx.EVT_COMBOBOX, self.modelMenuSelect)
        self.modelMenu.SetToolTipString("Model on which to perform point mutations")
        self.rotamersLoaded = False
        self.selectedModel = ""
        if (platform.system() == "Windows"):
            self.lblPosition = wx.StaticText(self, -1, "Position", (195, 90), (100, 20), wx.ALIGN_CENTRE)
            self.lblPosition.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblPosition = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/pointmutations/lblPosition.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(195, 90), size=(100, 20))
        else:
            self.lblPosition = wx.StaticText(self, -1, "Position", (195, 90), style=wx.ALIGN_CENTRE)
            self.lblPosition.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblPosition, 195, 100)
        self.lblPosition.SetForegroundColour("#FFFFFF")
        self.positionMenu = wx.ComboBox(self, pos=(195, 110), size=(100, 25), choices=[], style=wx.CB_READONLY)
        self.positionMenu.Bind(wx.EVT_COMBOBOX, self.positionMenuSelect)
        self.positionMenu.SetToolTipString("Model residue to mutate")
        self.energiesCalculated = False
        
        if (platform.system() == "Windows"):
            self.lblResidue = wx.StaticText(self, -1, "Residue Type", (20, 140), (140, 20), wx.ALIGN_CENTRE)
            self.lblResidue.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblResidue = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/pointmutations/lblResidueType.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(20, 140), size=(140, 20))
        else:
            self.lblResidue = wx.StaticText(self, -1, "Residue Type", (20, 140), style=wx.ALIGN_CENTRE)
            self.lblResidue.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblResidue, 20, 140)
        self.lblResidue.SetForegroundColour("#FFFFFF")
        self.AAlist = ["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"]
        self.AAlist.sort()
        self.residueMenu = wx.ComboBox(self, pos=(20, 160), size=(140, 25), choices=self.AAlist, style=wx.CB_READONLY)
        self.residueMenu.SetToolTipString("Amino acid type of the rotamers to be searched")
        if (platform.system() == "Darwin"):
            self.btnSearch = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pointmutations/btnSearch.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(195, 160), size=(100, 25))
        else:
            self.btnSearch = wx.Button(self, id=-1, label="Search", pos=(195, 160), size=(100, 25))
            self.btnSearch.SetForegroundColour("#000000")
            self.btnSearch.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnSearch.Bind(wx.EVT_BUTTON, self.searchRotamers)
        self.btnSearch.Disable()
        self.btnSearch.SetToolTipString("Evaluate all rotamers for the chosen amino acid type at the chosen position")

        if (platform.system() == "Darwin"):
            self.scoretypeMenu = wx.ComboBox(self, pos=(7, 190), size=(305, 25), choices=[], style=wx.CB_READONLY)
        else:
            self.scoretypeMenu = wx.ComboBox(self, pos=(7, 190), size=(305, 25), choices=[], style=wx.CB_READONLY | wx.CB_SORT)
        self.scoretypeMenu.Bind(wx.EVT_COMBOBOX, self.scoretypeMenuSelect)
        self.scoretypeMenu.Disable() # Is only enabled after a design and before accepting it
        self.scoretypeMenu.SetToolTipString("Set scoretype by which PyMOL residues will be colored")
        
        self.grdEnergies = wx.grid.Grid(self)
        self.grdEnergies.CreateGrid(0, 1)
        if (winh-270 > 200):
            self.grdEnergies.SetSize((320, winh-270))
        else:
            self.grdEnergies.SetSize((320, 200))
        self.grdEnergies.SetPosition((0, 220))
        self.grdEnergies.SetLabelFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.grdEnergies.DisableDragColSize()
        self.grdEnergies.DisableDragRowSize()
        self.grdEnergies.SetColLabelValue(0, "Energy (REU)")
        self.grdEnergies.SetRowLabelSize(160)
        self.grdEnergies.SetColSize(0, 160)
        self.grdEnergies.Bind(wx.grid.EVT_GRID_CELL_LEFT_CLICK, self.gridClick)
        self.resfile = []
        # Switch for telling the resfile to look for the restype from the sequence window pre-design
        # or from the outputted designed before accepting it
        self.useDesignedSeq = False 
        
        ypos = self.grdEnergies.GetPosition()[1] + self.grdEnergies.GetSize()[1] + 10
        if (platform.system() == "Windows"):
            self.lblTotalELabel = wx.StaticText(self, -1, "Total E:", (0, ypos+2), (160, 20), wx.ALIGN_CENTRE)
            self.lblTotalELabel.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblTotalELabel = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/pointmutations/lblTotalELabel.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, ypos+2), size=(160, 20))
        else:
            self.lblTotalELabel = wx.StaticText(self, -1, "Total E:", (0, ypos+2), style=wx.ALIGN_CENTRE)
            self.lblTotalELabel.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
            resizeTextControlForUNIX(self.lblTotalELabel, 0, 160)
        self.lblTotalELabel.SetForegroundColour("#FFFFFF")
        if (platform.system() != "Linux"):
            self.lblTotalE = wx.StaticText(self, -1, "", (160, ypos+2), (160, 20), wx.ALIGN_CENTRE)
            self.lblTotalE.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        else:
            self.lblTotalE = wx.StaticText(self, -1, "", (160, ypos+2), style=wx.ALIGN_CENTRE)
            self.lblTotalE.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
            resizeTextControlForUNIX(self.lblTotalE, 160, 160)
        self.lblTotalE.SetForegroundColour("#FFFFFF")
        #if (platform.system() == "Darwin"):
        #    self.btnLeft = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pointmutations/btnLeft.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, ypos), size=(25, 25))
        #else:
        #    self.btnLeft = wx.Button(self, id=-1, label="<", pos=(0, ypos), size=(25, 25))
        #    self.btnLeft.SetForegroundColour("#000000")
        #    self.btnLeft.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        #self.btnLeft.Bind(wx.EVT_BUTTON, self.rotamerLeft)
        #self.btnLeft.SetToolTipString("View the previous rotamer")
        #if (platform.system() == "Darwin"):
        #    self.btnRight = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pointmutations/btnRight.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(295, ypos), size=(25, 25))
        #else:
        #    self.btnRight = wx.Button(self, id=-1, label=">", pos=(295, ypos), size=(25, 25))
        #    self.btnRight.SetForegroundColour("#000000")
        #    self.btnRight.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        #self.btnRight.Bind(wx.EVT_BUTTON, self.rotamerRight)
        #self.btnRight.SetToolTipString("View the next rotamer")
        #if (platform.system() == "Windows"):
        #    self.lblSelection = wx.StaticText(self, -1, "Nothing Selected", (30, ypos+2), (265, 15), wx.ALIGN_CENTRE)
        #    self.lblSelection.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        #else:
        #    self.lblSelection = wx.StaticText(self, -1, "Nothing Selected", (30, ypos+2), style=wx.ALIGN_CENTRE)
        #    self.lblSelection.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        #    resizeTextControlForUNIX(self.lblSelection, 30, 265)
        #self.lblSelection.SetForegroundColour("#FFFFFF")
        
        if (platform.system() == "Darwin"):
            self.btnMutate = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/pointmutations/btnMutate.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(100, ypos+30), size=(120, 25))
        else:
            self.btnMutate = wx.Button(self, id=-1, label="Mutate!", pos=(100, ypos+30), size=(120, 25))
            self.btnMutate.SetForegroundColour("#000000")
            self.btnMutate.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.btnMutate.Bind(wx.EVT_BUTTON, self.mutateClick)
        self.btnMutate.Disable()
        self.btnMutate.SetToolTipString("Mutate the current position to the selected rotamer")
        #self.WTType = "NATRO"
        
        self.scrollh = self.btnMutate.GetPosition()[1] + self.btnMutate.GetSize()[1] + 5
        self.SetScrollbars(1, 1, 320, self.scrollh)
        
        self.tmrScore = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, self.threadScore, self.tmrScore)
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
            browser.open(self.parent.parent.scriptdir + "/help/mutations.html")
        else:
            webbrowser.open(self.parent.parent.scriptdir + "/help/mutations.html")
    
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
        if (not(self.rotamersLoaded)):
            self.positionMenu.SetSelection(event.GetRow())
            self.positionMenuSelect(event)
        else:
            self.rotamerClick(event)
        event.Skip()

    def modelMenuSelect(self, event):
        # Did the user already calculate rotamers and not make a mutation?
        # If so, tell them that they haven't done this and that these rotamers will be lost
        # if they change the model now
        if (len(self.selectedModel) > 0 and self.rotamersLoaded):
            dlg = wx.MessageDialog(self, "You have not made a mutation on the previous model.  If you continue, the loaded rotamers will be lost.  Continue?", "Mutation Not Accepted", wx.YES_NO | wx.ICON_EXCLAMATION | wx.CENTRE)
            if (dlg.ShowModal() == wx.ID_NO):
                self.modelMenu.SetStringSelection(self.selectedModel)
                dlg.Destroy()
                return
            dlg.Destroy()
        # Update the list of positions with the new model
        self.selectedModel = self.modelMenu.GetStringSelection()
        logInfo("New model selected: " + self.selectedModel)
        # Get the location of the pose
        poseindx = self.seqWin.getPoseIndexForModel(self.selectedModel)
        # Read the positions
        pose = self.seqWin.poses[poseindx]
        positions = []
        ires = 1
        for ch in pose[0]:
            for residue in ch:
                # Skip NCAAs
                if ("ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR ".find(residue.resname) >= 0):
                    chain = ch.id
                    if (chain == " " or chain == ""):
                        chain = "_"
                    label = chain + ":" + AA3to1(residue.resname) + str(residue.id[1])
                    positions.append(label)
                ires = ires + 1
        oldpos = self.positionMenu.GetSelection()
        self.positionMenu.Clear()
        self.positionMenu.AppendItems(positions)
        if (platform.system() == "Windows"):
            self.positionMenu.SetSelection(-1)
        else:
            if (oldpos < len(positions)):
                self.positionMenu.SetSelection(oldpos)
            else:
                self.positionMenu.SetSelection(0)
        self.energiesCalculated = False
        # Now we're going to calculate the scores of all the residues so the user can see
        # which ones have high energies and so we have access to the PyMOL coloring
        self.grdEnergies.ClearGrid()
        if (self.grdEnergies.NumberRows > 0):
            self.grdEnergies.DeleteRows(0, self.grdEnergies.NumberRows)
        #thrScore = Thread(target=self.threadScore, args=())
        #thrScore.start()
        self.stage = 1
        self.tmrScore.Start(1000)
        logInfo("Attempted to calculate scores for this model")
        self.seqWin.labelMsg.SetLabel("Calculating energies, please be patient...")
        self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")
        self.seqWin.msgQueue.append("Calculating energies, please be patient...")
        
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
        self.Scroll(0, self.winscrollpos)
            
    def relabelEnergies(self, seqpos="none"):
        # This function displays the residue energy on the CA of the selected residue or all visible residues
        # First label everything, then hide the labels of everything not in selstring
        # This version of relabelEnergies is different than the one in tools
        selectedScoretype = str(self.scoretypeMenu.GetStringSelection())
        found = False
        for (scoretypestr, name) in scoretypes.items():
            if (name == selectedScoretype):
                scoretype = scoretypestr
                found = True
                break
        if (not(self.rotamersLoaded)):
            if (found):
                sindx = self.residue_E[0].index(scoretype)
            else:
                sindx = 0 # Default to total_score
            residue_E = []
            for i in range(1, len(self.residue_E)):
                residue_E.append(self.residue_E[i][sindx])
            poseindx = self.seqWin.getPoseIndexForModel(self.selectedModel)
            for i in range(0, len(residue_E)):
                ires = 1
                for ch in self.seqWin.poses[poseindx][0]:
                    for residue in ch:
                        # Skip NCAAs/HETATMs
                        if (not(residue.resname in "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR")):
                            continue
                        if (ires == i+1):
                            chain = ch.id
                            res = residue.id[1]
                            break
                        ires = ires + 1
                if (seqpos != "none" and str(res) != seqpos):
                    continue
                Elabel = str(int(residue_E[i] * 100.0) / 100.0)
                if (chain != "" and chain != " " and chain != "_"):
                    try:
                        self.cmd.select("labelsele", "model " + self.selectedModel + " and chain " + str(chain) + " and resi " + str(res) + " and name ca")
                        self.cmd.label("labelsele", Elabel)
                    except:
                        # Do nothing on a fail (NCAA that has no CA?)
                        pass
                else:
                    try:
                        self.cmd.select("labelsele", "model " + self.selectedModel + " and resi " + str(res) + " and name ca")
                        self.cmd.label("labelsele", Elabel)
                    except:
                        # Do nothing on a fail (NCAA that has no CA?)
                        pass
        else:
            if (found):
                sindx = self.rotamer_E[0].index(scoretype)
            else:
                sindx = 0 # Default to total_score
            Elabel = str(int(self.rotamer_E[self.currentRotamer+1][sindx] * 100.0) / 100.0)
            if (self.rotamerChain != "" and self.rotamerChain != " " and self.rotamerChain != "_"):
                try:
                    self.cmd.select("labelsele", "model rotamer_view and chain " + str(self.rotamerChain) + " and resi " + str(self.rotamerSeqPos) + " and name ca")
                    self.cmd.label("labelsele", Elabel)
                except:
                    # Do nothing on a fail (NCAA that has no CA?)
                    pass
            else:
                try:
                    self.cmd.select("labelsele", "model rotamer_view and resi " + str(self.rotamerSeqPos) + " and name ca")
                    self.cmd.label("labelsele", Elabel)
                except:
                    # Do nothing on a fail (NCAA that has no CA?)
                    pass
        self.cmd.delete("labelsele")
        self.cmd.enable("sele")
            
    def recolorEnergies(self):
        # This version of recolorEnergies is different than the one in tools
        # Color the PyMOL_Mover designed_view object with the selected scoretype
        if (not(self.energiesCalculated)):
            return
        # Now figure out what scoretype it is from the formal user-displayed name
        selectedScoretype = str(self.scoretypeMenu.GetStringSelection())
        found = False
        for (scoretypestr, name) in scoretypes.items():
            if (name == selectedScoretype):
                scoretype = scoretypestr
                found = True
                break
        if (found):
            sindx = self.residue_E[0].index(scoretype)
        else:
            sindx = 0 # Default to total_score
        residue_E = []
        for i in range(1, len(self.residue_E)):
            # indx 0 is the row that has all the scoring types so we skip it
            residue_E.append(self.residue_E[i][sindx])
        residue_E = scale_list(residue_E)
        poseindx = self.seqWin.getPoseIndexForModel(self.selectedModel)
        for i in range(0, len(residue_E)):
            # Update the grdEnergies with the positional energies if rotamers haven't been
            # searched yet
            if (not(self.rotamersLoaded)):
                # Find it, this is necessary because the scoring information might have some
                # NCAAs and HETATMs in it, but grdEnergies skipped all of these
                for j in range(0, self.grdEnergies.NumberRows):
                    if (self.grdEnergies.GetRowLabelValue(j) == self.resIDs[i]):
                        E = float(int(self.residue_E[j+1][sindx] * 100)) / 100.0
                        self.grdEnergies.SetCellValue(j, 0, str(E))
                        break
            ires = 1
            for ch in self.seqWin.poses[poseindx][0]:
                for residue in ch:
                    # Skip NCAAs/HETATMs
                    if (not(residue.resname in "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR")):
                        continue
                    if (ires == i+1):
                        chain = ch.id
                        res = residue.id[1]
                        break
                    ires = ires + 1
            r = residue_E[i]
            b = 255 - r
            g = 0
            if (chain != "" and chain != " " and chain != "_"):
                self.cmd.delete("colorsele")
                self.cmd.select("colorsele", "model " + self.selectedModel + " and chain " + chain + " and resi " + str(res))
                self.cmd.color("0x%02x%02x%02x" % ((r, g, b)), "colorsele")
            else:
                self.cmd.delete("colorsele")
                self.cmd.select("colorsele", "model " + self.selectedModel + " and resi " + str(res))
                self.cmd.color("0x%02x%02x%02x" % ((r, g, b)), "colorsele")
            residue["N"].set_bfactor(r / 255.0 * 100.0)
        self.cmd.delete("colorsele")
        self.cmd.enable("sele")
    
    def positionMenuSelect(self, event):
        # Set the selected residue's row to red so it is easy to see what the selection is
        for r in range(0, self.grdEnergies.NumberRows):
            if (r == self.positionMenu.GetSelection()):
                for c in range(0, self.grdEnergies.NumberCols):
                    self.grdEnergies.SetCellBackgroundColour(r, c, "light blue")
            else:
                for c in range(0, self.grdEnergies.NumberCols):
                    self.grdEnergies.SetCellBackgroundColour(r, c, "white")
        self.grdEnergies.Refresh()
        # Find the neighborhood view on this residue
        try:
            chain = self.positionMenu.GetStringSelection()[0]
        except:
            # There's nothing there, probably because I was called from scoretypeMenuSelect
            return
        seqpos = self.positionMenu.GetStringSelection()[3:].strip() # indx 2 is the AA
        model = self.selectedModel
        logInfo("Position " + self.positionMenu.GetStringSelection() + " was selected")
        self.cmd.hide("all")
        if (chain == " " or chain == "_"):
            self.cmd.select("mutsele", "resi " + seqpos + " and model " + model)
        else:
            self.cmd.select("mutsele", "resi " + seqpos + " and model " + model + " and chain " + chain)
        self.cmd.select("mutsele", "model " + model + " within 12 of mutsele")
        self.cmd.show("cartoon", "mutsele")
        self.cmd.hide("ribbon", "mutsele")
        self.cmd.show("sticks", "mutsele")
        self.cmd.set_bond("stick_radius", 0.1, "mutsele")
        # Display energy labels for designed structures
        if (self.energiesCalculated):
            # When the model is selected, the scores are calculated in the background
            # If they are available now, then use them, if not, use default coloring
            self.relabelEnergies(seqpos)
            self.cmd.label("not mutsele", "")
        self.cmd.zoom("mutsele")
        if (chain == " " or chain == "_"):
            self.cmd.select("mutsele", "resi " + seqpos + " and model " + model)
        else:
            self.cmd.select("mutsele", "resi " + seqpos + " and model " + model + " and chain " + chain)
        self.cmd.show("sticks", "mutsele")
        self.cmd.set_bond("stick_radius", 0.25, "mutsele")
        # Highlight this residue in PyMOL
        self.cmd.select("seqsele", "mutsele")
        self.cmd.enable("seqsele")
        self.cmd.delete("mutsele")
        self.seqWin.selectUpdate(False)

    def rotamerLeft(self, event):
        logInfo("Left button clicked")
        if (self.currentRotamer > 0):
            self.currentRotamer = self.currentRotamer - 1
            #self.lblSelection.SetLabel(self.rotanames[self.currentRotamer])
            #logInfo("Rotamer " + self.lblSelection.GetLabel() + " selected")
            #if (platform.system() != "Windows"):
                #resizeTextControlForUNIX(self.lblSelection, 30, 265)
            self.loadRotamer(self.currentRotamer)
            self.grdEnergies.SetCellBackgroundColour(self.currentRotamer, 0, "light blue")
            self.grdEnergies.SetCellBackgroundColour(self.currentRotamer+1, 0, "white")
            self.grdEnergies.Refresh()
    
    def rotamerRight(self, event):
        logInfo("Right button clicked")
        if (self.currentRotamer < len(self.rotanames)-1):
            self.currentRotamer = self.currentRotamer + 1
            #self.lblSelection.SetLabel(self.rotanames[self.currentRotamer])
            #logInfo("Rotamer " + self.lblSelection.GetLabel() + " selected")
            #if (platform.system() != "Windows"):
                #resizeTextControlForUNIX(self.lblSelection, 30, 265)
            self.loadRotamer(self.currentRotamer)
            self.grdEnergies.SetCellBackgroundColour(self.currentRotamer, 0, "light blue")
            self.grdEnergies.SetCellBackgroundColour(self.currentRotamer-1, 0, "white")
            self.grdEnergies.Refresh()
            
    def rotamerClick(self, event):
        logInfo("Rotamer " + str(event.GetRow()) + " clicked")
        if (self.currentRotamer < len(self.rotanames)):
            self.grdEnergies.SetCellBackgroundColour(self.currentRotamer, 0, "white")
            self.currentRotamer = event.GetRow()
            #self.lblSelection.SetLabel(self.rotanames[self.currentRotamer])
            #logInfo("Rotamer " + self.lblSelection.GetLabel() + " selected")
            #if (platform.system() != "Windows"):
                #resizeTextControlForUNIX(self.lblSelection, 30, 265)
            self.loadRotamer(self.currentRotamer)
            self.grdEnergies.SetCellBackgroundColour(self.currentRotamer, 0, "light blue")
            self.grdEnergies.Refresh()

    def loadRotamer(self, rotindx):
        # This function switches the rotamer in PyMOL to a new one by rotating the chi angles
        # Do nothing for ALA and GLY, which have no rotamers
        if (self.grdEnergies.GetRowLabelValue(rotindx)[0] in "AG"):
            return
        for iatom in range(0, len(self.atompos[self.currentRotamer])):
            [atom, x, y, z] = self.atompos[self.currentRotamer][iatom].split()
            atomname = "model rotamer_view and name " + atom
            self.cmd.alter_state(1, atomname, "x=" + x)
            self.cmd.alter_state(1, atomname, "y=" + y)
            self.cmd.alter_state(1, atomname, "z=" + z)
        #for ichi in range(0, len(self.xatoms)):
        #    atomname1 = "model rotamer_view and name " + self.xatoms[ichi][0]
        #    atomname2 = "model rotamer_view and name " + self.xatoms[ichi][1]
        #    atomname3 = "model rotamer_view and name " + self.xatoms[ichi][2]
        #    atomname4 = "model rotamer_view and name " + self.xatoms[ichi][3]
        #    chival = float(self.chivals[self.currentRotamer][ichi])
        #    print atomname1, atomname2, atomname3, atomname4, chival
        #    self.cmd.set_dihedral(atomname1, atomname2, atomname3, atomname4, chival)
        # The set_dihedral generates a lot of junk that we need to get rid of
        #self.cmd.delete("pk1")
        #self.cmd.delete("pk2")
        #self.cmd.delete("pkbond")
        #self.cmd.delete("pkmol")
        #self.cmd.hide("dihedrals", "model rotamer_view")
        #self.cmd.hide("labels", "model rotamer_view")
        # Get the energy label back
        self.relabelEnergies(self.rotamerSeqPos)
        #self.cmd.align("model rotamer_view and backbone", "model " + self.selectedModel + " and backbone", cycles=0)

    def threadSearch(self, event):
        # Why am I using a Timer?  See the explanation in kic.py
        if (self.stage == 1):
            self.tmrSearch.Stop()
            self.timeoutCount = 0
            # Generate a file for the daemon, the daemon will try all Dunbrack rotamers at this position
            # for this residue type and return a sorted list of these rotamers and their scores
            f = open("rotamerinputtemp", "w")
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
            f.write("RESTYPE\t" + self.residueMenu.GetStringSelection() + "\n")
            f.write("SEQPOS\t" + self.positionMenu.GetStringSelection() + "\n")
            f.close()
            self.rotamerSeqPos = self.positionMenu.GetStringSelection()[3:].strip()
            self.rotamerChain = self.positionMenu.GetStringSelection()[0]
            appendScorefxnParamsInfoToFile("rotamerinputtemp", self.selectWin.weightsfile)
            if (useServer):
                try: 
                    self.ID = sendToServer("rotamerinput")
                    self.usingServer = True
                    logInfo("Rotamer input sent to server daemon with ID " + self.ID)
                except:
                    # Something failed, default to the local daemon
                    os.rename("rotamerinputtemp", "rotamerinput")
                    self.usingServer = False
                    logInfo("Server daemon not available, rotamer input uploaded at rotamerinput")
            else:
                os.rename("rotamerinputtemp", "rotamerinput")
                self.usingServer = False
                logInfo("Rotamer input uploaded locally at rotamerinput")
            self.stage = 2
            self.progress = None
            try:
                os.remove("progress")
            except:
                pass
            self.tmrSearch.Start(1000)
        else:
            if (self.usingServer):
                # See if the file has been uploaded yet and bring it here if so
                queryServerForResults("rotameroutput-" + self.ID)
                self.timeoutCount = self.timeoutCount + 1
            if (self.timeoutCount >= serverTimeout):
                self.tmrSearch.Stop()
                # If this is taking too long, maybe there's something wrong with the server
                # Ask the user if they want to continue waiting or use the local daemon instead
                dlg = wx.MessageDialog(self, "The server is taking a long time to respond.  Continue to wait?  Pressing No will run the calculations locally.", "Delayed Server Response", wx.YES_NO | wx.ICON_EXCLAMATION | wx.CENTRE)
                if (dlg.ShowModal() == wx.ID_YES):
                    # Reset the counter
                    self.timeoutCount = 0
                else:
                    self.usingServer = False
                    self.timeoutCount = 0
                    os.rename("rotamerinputtemp", "rotamerinput")
                    logInfo("Server took too long to respond so the local daemon was used")
                dlg.Destroy()
                self.tmrSearch.Start(1000)
            if (os.path.isfile("rotameroutput")):
                if (self.progress is not None):
                    try:
                        self.progress.Destroy()
                    except:
                        pass
                self.tmrSearch.Stop()
                # Read the output dumped by the child process
                self.rotamer_E = []
                self.rotanames = []
                self.chivals = []
                self.xatoms = []
                self.atompos = []
                f = open("rotameroutput", "r")
                logInfo("Rotamer search output found at rotameroutput")
                for aline in f:
                    if (aline[0:6] == "OUTPUT"):
                        pdbfile = aline.split("\t")[1].strip()
                        self.rotamerView = self.seqWin.pdbreader.get_structure("rotamer_view", pdbfile)
                        self.rotamermodel = pdbfile.split("_R.pdb")[0]
                    elif (aline[0:5] == "INDEX"):
                        self.posepos = int(aline.split("\t")[1])
                        ID = aline.split("\t")[2].strip()
                        if (ID[0] == "_"):
                            self.pymolpos = "model rotamer_view and resi " + ID[1:]
                        else:
                            self.pymolpos = "model rotamer_view and chain " + ID[0] + " and resi " + ID[1:]
                    elif (aline[0:6] == "ENERGY"):
                        if (aline.split()[1] == "total_score"):
                            # This is the scoretype line, row 0 in residue_E
                            self.rotamer_E.append(aline.split()[1:])
                        else:
                            self.rotamer_E.append([])
                            indx = len(self.rotamer_E) - 1
                            for E in aline.split()[1:]:
                                self.rotamer_E[indx].append(float(E))
                    elif (aline[0:4] == "NAME"):
                        self.rotanames.append(aline.split("\t")[1].strip())
                    elif (aline[0:7] == "ROTAMER"):
                        self.atompos.append(aline.split("\t")[1:])
                    #elif (aline[0:8] == "CHIATOMS"):
                        #self.xatoms.append(aline.split("\t")[1].split()[0:4])
                    #elif (aline[0:3] == "CHI"):
                        #self.chivals.append(aline.split("\t")[1:])
                    
                f.close()
                try:
                    self.cmd.remove("rotamer_view")
                    self.cmd.delete("rotamer_view")
                except:
                    pass
                # Load the pdbfile that has the new restype at the selected position
                self.cmd.load(pdbfile, "rotamer_view")
                # Now delete everything except for that new residue, so we can explore the rotamers
                # visually in PyMOL
                self.cmd.hide("labels", "all")
                if (self.rotamerChain != "_"):
                    self.cmd.select("rotsele", "model rotamer_view and chain " + self.rotamerChain + " and resi " + self.rotamerSeqPos)
                    self.cmd.select("rotsele", "model rotamer_view and not rotsele")
                else:
                    self.cmd.select("rotsele", "model rotamer_view and not resi " + self.rotamerSeqPos)
                self.cmd.remove("rotsele")
                self.cmd.select("rotsele", "model rotamer_view")
                self.cmd.center("rotsele")
                self.cmd.show("sticks", "rotsele")
                self.cmd.color("gray", "rotsele and symbol c")
                if (self.rotamerChain != "_"):
                    self.cmd.hide("sticks", "model " + self.rotamermodel + " and chain " + self.rotamerChain + " and resi " + self.rotamerSeqPos)
                    self.cmd.hide("spheres", "model " + self.rotamermodel + " and chain " + self.rotamerChain + " resi " + self.rotamerSeqPos)
                else:
                    self.cmd.hide("sticks", "model " + self.rotamermodel + " and resi " + self.rotamerSeqPos)
                    self.cmd.hide("spheres", "model " + self.rotamermodel + " and resi " + self.rotamerSeqPos)
                self.cmd.select("rotsele", "visible")
                self.cmd.zoom("rotsele")
                self.cmd.delete("rotsele")
                # Now replace the information in the grid with the rotamers and their energies
                self.grdEnergies.ClearGrid()
                self.grdEnergies.DeleteRows(0, self.grdEnergies.NumberRows)
                for i in range(1, len(self.rotamer_E)):
                    self.grdEnergies.AppendRows(1)
                    self.grdEnergies.SetRowLabelValue(i-1, self.rotanames[i-1])
                    E = float(int(self.rotamer_E[i][0] * 100.0) / 100.0)
                    self.grdEnergies.SetCellValue(i-1, 0, str(E))
                    # This needs to happen before the readOnly attr is set otherwise it doesn't apply the center on the
                    # last cell for some bizarre reason
                    self.grdEnergies.SetCellAlignment(i-1, 0, wx.ALIGN_CENTRE, wx.ALIGN_CENTRE)
                    # Very important note: you actually need to create a new GridCellAttr for each new row
                    # You cannot just declare it outside of the loop and use the same one for each row otherwise you
                    # get some pretty nasty crashes when you delete rows
                    readOnly = wx.grid.GridCellAttr()
                    readOnly.SetReadOnly(True)
                    self.grdEnergies.SetRowAttr(i-1, readOnly)
                self.rotamersLoaded = True
                self.currentRotamer = 0
                self.grdEnergies.SetCellBackgroundColour(0, 0, "light blue")
                self.grdEnergies.Refresh()
                #self.lblSelection.SetLabel(self.rotanames[self.currentRotamer])
                #if (platform.system() != "Windows"):
                    #resizeTextControlForUNIX(self.lblSelection, 30, 265)
                # Load the rotamer into view
                self.loadRotamer(self.currentRotamer)
                # Pop this message out of the queue
                for i in range(0, len(self.seqWin.msgQueue)):
                    if (self.seqWin.msgQueue[i].find("Performing rotamer search") >= 0):
                        self.seqWin.msgQueue.pop(i)
                        break
                if (len(self.seqWin.msgQueue) > 0):
                    self.seqWin.labelMsg.SetLabel(self.seqWin.msgQueue[len(self.seqWin.msgQueue)-1])
                else:
                    self.seqWin.labelMsg.SetLabel("")
                self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
                self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")
                self.btnMutate.Enable()
                self.btnSearch.Enable()
                self.parent.GoBtn.Enable()
                os.remove("rotameroutput")
            elif (os.path.isfile("errreport")):
                self.tmrSearch.Stop()
                self.recoverFromError("search")
            elif (os.path.isfile("progress")):
                # The local daemon can output its progress to keep the GUI updated about
                # how far along it is, along with a message
                # This is optional
                # See job/__init__.py for more information
                if (self.progress is None):
                    self.progress = wx.ProgressDialog("Rotamer Search Progress", "Performing rotamer search...", 100, style=wx.PD_APP_MODAL | wx.PD_ELAPSED_TIME | wx.PD_REMAINING_TIME)
                fin = open("progress", "r")
                data = fin.readlines()
                fin.close()
                # First line should be a fraction
                try:
                    num = float(data[0].split("/")[0].strip())
                    den = float(data[0].split("/")[1].strip())
                    # Convert to a percentage
                    percent = int(num / den * 100.0)
                    if (percent > 99):
                        # Let's the appearance of the output file kill the progress bar
                        percent = 99
                except:
                    return
                try:
                    # The optional second line is a new message
                    newmsg = data[1].strip()
                    (keepGoing, skip) = self.progress.Update(percent, newmsg)
                except:
                    (keepGoing, skip) = self.progress.Update(percent)
                if (not(keepGoing)):
                    # User clicked "Cancel" on the progress bar
                    self.cancelSubmit()
                    self.progress.Destroy()
        
    def searchRotamers(self, event):
        logInfo("Search rotamer button clicked")
        # Is there a residue selected?
        if (len(self.residueMenu.GetStringSelection()) < 3 or len(self.positionMenu.GetStringSelection().strip()) == 0):
            dlg = wx.MessageDialog(self, "Please select a residue type and a position.", "Incomplete Input", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
            if (dlg.ShowModal() == wx.ID_OK):
                pass
            dlg.Destroy()
            return
        #thrSearch = Thread(target=self.threadSearch, args=())
        #thrSearch.start()
        self.stage = 1
        self.tmrSearch = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, self.threadSearch, self.tmrSearch)
        self.tmrSearch.Start(1000)
        self.seqWin.labelMsg.SetLabel("Performing rotamer search, please be patient...")
        self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")
        self.seqWin.msgQueue.append("Performing rotamer search, please be patient...")
        self.parent.GoBtn.Disable()
        self.btnSearch.Disable()
        self.btnMutate.Disable()
    
    def recolorGrid(self, allresidue_E, grid, selectedScoretype):
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
            r = residue_E[row]
            b = 255 - r
            g = 0
            for c in range(0, grid.NumberCols):
                grid.SetCellTextColour(row, c, (r, g, b))
        grid.Refresh()
    
    def scoretypeMenuSelect(self, event):
        # Make sure the energies have been calculated, otherwise do nothing
        if (not(self.energiesCalculated)):
            return
        self.recolorEnergies()
        logInfo("Scoretype " + self.scoretypeMenu.GetStringSelection() + " selected")
        if (not(self.rotamersLoaded)):
            self.recolorGrid(self.residue_E, self.grdEnergies, self.scoretypeMenu.GetStringSelection())
            self.positionMenuSelect(event) # To update all the labels
        else:
            # Now we have to replace all the rotamer energies with the energies for the selected
            # scoretype from the menu
            selectedScoretype = str(self.scoretypeMenu.GetStringSelection())
            found = False
            for (scoretypestr, name) in self.scoretypes.items():
                if (name == selectedScoretype):
                    scoretype = scoretypestr
                    found = True
                    break
            if (found):
                sindx = self.rotamer_E[0].index(scoretype)
            else:
                sindx = 0 # Default to total_score
            for i in range(1, len(self.rotamer_E)):
                E = float(int(self.rotamer_E[i][sindx] * 100.0) / 100.0)
                self.grdEnergies.SetCellValue(i-1, 0, str(E))
            self.relabelEnergies(self.rotamerSeqPos)
    
    def enableControls(self, enable=True):
        if (enable):
            self.btnWT.Enable()
            self.btnHydrophobic.Enable()
            self.btnHydrophilic.Enable()
            self.btnAromatic.Enable()
            self.btnAllAA.Enable()
            self.btnAdd.Enable()
            self.btnRemove.Enable()
            self.btnRestrict.Enable()
            self.btnAll.Enable()
            self.btnClear.Enable()
            self.btnApply.Enable()
            self.btnAddRes.Enable()
            self.btnRemoveRes.Enable()
            self.btnLoadResfile.Enable()
            self.btnWTType.Enable()
            self.btnDesign.Enable()
        else:
            self.btnWT.Disable()
            self.btnHydrophobic.Disable()
            self.btnHydrophilic.Disable()
            self.btnAromatic.Disable()
            self.btnAllAA.Disable()
            self.btnAdd.Disable()
            self.btnRemove.Disable()
            self.btnRestrict.Disable()
            self.btnAll.Disable()
            self.btnClear.Disable()
            self.btnApply.Disable()
            self.btnAddRes.Disable()
            self.btnRemoveRes.Disable()
            self.btnLoadResfile.Disable()
            self.btnWTType.Disable()
            self.btnDesign.Disable()
    
    def mutateClick(self, event):
        # Change the chi angles in the pose to the ones selected
        logInfo("Mutate button clicked")
        # Get rid of the original pose, save the designed pose, and reload the structure in PyMOL
        poseindx = self.seqWin.getPoseIndexForModel(self.rotamermodel)
        try:
            self.stored.xyz = []
            self.cmd.iterate_state(1, "model rotamer_view", "stored.xyz.append([name, x, y, z])")
            for ch in self.rotamerView[0]:
                for residue in ch:
                    if (residue.id[1] == int(self.rotamerSeqPos)):
                        for [name, x, y, z] in self.stored.xyz:
                            coords = residue[name].get_coord()
                            coords[0] = x
                            coords[1] = y
                            coords[2] = z
                            residue[name].set_coord(coords)
            self.seqWin.pdbwriter.set_structure(self.rotamerView)
            self.seqWin.pdbwriter.save(self.rotamermodel + "_R.pdb")
            self.cmd.remove(self.rotamermodel)
            self.cmd.delete(self.rotamermodel)
            self.cmd.remove("rotamer_view")
            self.cmd.delete("rotamer_view")
            self.cmd.load(self.rotamermodel + "_R.pdb", self.rotamermodel)
            self.seqWin.reloadPose(poseindx, self.rotamermodel, self.rotamermodel + "_R.pdb")
            defaultPyMOLView(self.cmd, self.rotamermodel)
            del self.rotamerView
            # IMPORTANT: You have to replace the model in the sandbox with the new designed model
            os.remove(self.rotamermodel + ".pdb")
            os.rename(self.rotamermodel + "_R.pdb", self.rotamermodel + ".pdb")
            self.btnMutate.Disable()
            self.rotamersLoaded = False
            self.modelMenuSelect(event) # To regenerate the scores and get the coloring back
        except:
            # Some weird error happened, do nothing instead of crashing
            print "Bug at mutate button click"
            pass        
    
    def recoverFromError(self, stage="score"):
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
        if (stage == "score"):
            # Get rid of the messages
            for i in range(0, len(self.seqWin.msgQueue)):
                if (self.seqWin.msgQueue[i].find("Calculating energies") >= 0):
                    self.seqWin.msgQueue.pop(i)
                    break
            if (len(self.seqWin.msgQueue) > 0):
                self.seqWin.labelMsg.SetLabel(self.seqWin.msgQueue[len(self.seqWin.msgQueue)-1])
            else:
                self.seqWin.labelMsg.SetLabel("")
            self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
            self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")
        else:
            self.seqWin.cannotDelete = False
            self.parent.GoBtn.Enable()
            self.btnSearch.Enable()
            self.btnMutate.Disable()
            # Get rid of the messages
            for i in range(0, len(self.seqWin.msgQueue)):
                if (self.seqWin.msgQueue[i].find("Performing rotamer search") >= 0):
                    self.seqWin.msgQueue.pop(i)
                    break
            if (len(self.seqWin.msgQueue) > 0):
                self.seqWin.labelMsg.SetLabel(self.seqWin.msgQueue[len(self.seqWin.msgQueue)-1])
            else:
                self.seqWin.labelMsg.SetLabel("")
            self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
            self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")
    
    def threadScore(self, event):
        # Why am I using a Timer?  See the explanation in kic.py
        goToSandbox()
        if (self.stage == 1):
            self.tmrScore.Stop()
            self.timeoutCount = 0
            # Dump the input file for the daemon
            f = open("scoreinputtemp", "w")
            pdbfile = self.selectedModel + ".pdb"
            # Dump the PDB from PyMOL first in case the coordinates were altered by the user
            self.scoreModel = self.selectedModel
            self.cmd.save(pdbfile.strip(), "model " + self.selectedModel)
            fixPyMOLSave(pdbfile.strip())
            f.write("PDBFILE\t" + pdbfile.strip() + "\n")
            f2 = open(pdbfile, "r")
            f.write("BEGIN PDB DATA\n")
            for aline in f2:
                f.write(aline.strip() + "\n")
            f.write("END PDB DATA\n")
            f2.close()
            f.close()
            appendScorefxnParamsInfoToFile("scoreinputtemp", self.selectWin.weightsfile)
            if (useServer):
                try: 
                    self.ID = sendToServer("scoreinput")
                    self.usingServer = True
                    logInfo("Scoring input sent to server daemon with ID " + self.ID)
                except:
                    # Something failed, default to the local daemon
                    os.rename("scoreinputtemp", "scoreinput")
                    self.usingServer = False
                    logInfo("Server daemon not available, scoring input uploaded at scoreinput")
            else:
                os.rename("scoreinputtemp", "scoreinput")
                self.usingServer = False
                logInfo("Scoring input uploaded locally at scoreinput")
            self.stage = 2
            self.tmrScore.Start(1000)
        else:
            if (self.usingServer):
                # See if the file has been uploaded yet and bring it here if so
                queryServerForResults("scoreoutput-" + self.ID)
                self.timeoutCount = self.timeoutCount + 1
            if (self.timeoutCount >= serverTimeout):
                self.tmrScore.Stop()
                # If this is taking too long, maybe there's something wrong with the server
                # Ask the user if they want to continue waiting or use the local daemon instead
                dlg = wx.MessageDialog(self, "The server is taking a long time to respond.  Continue to wait?  Pressing No will run the calculations locally.", "Delayed Server Response", wx.YES_NO | wx.ICON_EXCLAMATION | wx.CENTRE)
                if (dlg.ShowModal() == wx.ID_YES):
                    # Reset the counter
                    self.timeoutCount = 0
                else:
                    self.usingServer = False
                    self.timeoutCount = 0
                    os.rename("scoreinputtemp", "scoreinput")
                    logInfo("Server took too long to respond so the local daemon was used")
                dlg.Destroy()
                self.tmrDesign.Start(1000)
            if (os.path.isfile("scoreoutput")):
                logInfo("Found scoring output at scoreoutput")
                self.tmrScore.Stop()
                # Read the output dumped by the child process
                self.residue_E = []
                self.resIDs = []
                f = open("scoreoutput", "r")
                for aline in f:
                    if (aline[0:6] == "OUTPUT"):
                        pdbfile = aline.split("\t")[1].strip()
                    elif (aline[0:7] == "TOTAL_E"):
                        totalE = aline.split("\t")[1].strip()
                    elif (aline[0:6] == "ENERGY"):
                        if (aline.split()[1] == "total_score"):
                            # This is the scoretype line, row 0 in residue_E
                            self.residue_E.append(aline.split()[1:])
                        else:
                            self.residue_E.append([])
                            indx = len(self.residue_E) - 1
                            for E in aline.split()[1:]:
                                self.residue_E[indx].append(float(E))
                    elif (aline[0:2] == "ID"):
                        self.resIDs.append(aline.split("\t")[1].strip())
                f.close()
                self.lblTotalE.SetLabel(totalE)
                self.lblTotalE.SetForegroundColour("#FFFFFF")
                if (platform.system() == "Linux"):
                    resizeTextControlForUNIX(self.lblTotalE, 160, 160)
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
                    if (self.seqWin.msgQueue[i].find("Calculating energies") >= 0):
                        self.seqWin.msgQueue.pop(i)
                        break
                if (len(self.seqWin.msgQueue) > 0):
                    self.seqWin.labelMsg.SetLabel(self.seqWin.msgQueue[len(self.seqWin.msgQueue)-1])
                else:
                    self.seqWin.labelMsg.SetLabel("")
                self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
                self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")
                # Set up the grid for displaying energies
                positions = self.positionMenu.GetItems()
                for i in range(0, len(positions)):
                    self.grdEnergies.AppendRows(1)
                    self.grdEnergies.SetRowLabelValue(i, positions[i])
                    self.grdEnergies.SetCellValue(i, 0, "")
                    # This needs to happen before the readOnly attr is set otherwise it doesn't apply the center on the
                    # last cell for some bizarre reason
                    self.grdEnergies.SetCellAlignment(i, 0, wx.ALIGN_CENTRE, wx.ALIGN_CENTRE)
                    # Very important note: you actually need to create a new GridCellAttr for each new row
                    # You cannot just declare it outside of the loop and use the same one for each row otherwise you
                    # get some pretty nasty crashes when you delete rows
                    readOnly = wx.grid.GridCellAttr()
                    readOnly.SetReadOnly(True)
                    self.grdEnergies.SetRowAttr(i, readOnly)
                self.parent.GoBtn.Enable()
                self.btnSearch.Enable()
                self.cmd.remove(self.scoreModel)
                self.cmd.delete(self.scoreModel)
                self.cmd.load(self.scoreModel + "_S.pdb", self.scoreModel)
                poseindx = self.seqWin.getPoseIndexForModel(self.scoreModel)
                self.seqWin.reloadPose(poseindx, self.scoreModel, self.scoreModel + "_S.pdb")
                defaultPyMOLView(self.cmd, self.scoreModel)
                self.energiesCalculated = True
                self.recolorEnergies()
                self.seqWin.recolorResidues()
                self.recolorGrid(self.residue_E, self.grdEnergies, self.scoretypeMenu.GetStringSelection())
                os.remove("scoreoutput")
            elif (os.path.isfile("errreport")):
                self.tmrScore.Stop()
                self.recoverFromError("score")
import wx
import wx.grid
import wx.lib.scrolledpanel
import os
import os.path
import time
import platform
import multiprocessing
import datetime
#import requests
import webbrowser
from threading import Thread
from tools import *
from Bio.Align.Applications import MuscleCommandline

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
            self.lblProt = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/compmodel/lblCompModel.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(25, 15), size=(270, 25))
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
            self.lblInst = wx.StaticText(self, -1, "Model a structure ab initio using an existing\nstructure as a template.", (0, 45), (320, 25), wx.ALIGN_CENTRE)
            self.lblInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
        elif (platform.system() == "Darwin"):
            self.lblInst = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/compmodel/lblInstCompModel.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 45), size=(320, 80))
        else:
            self.lblInst = wx.StaticText(self, -1, "Model a structure ab initio using an existing\nstructure as a template.", (5, 45), style=wx.ALIGN_CENTRE)
            self.lblInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
            resizeTextControlForUNIX(self.lblInst, 0, self.GetSize()[0]-20)
        self.lblInst.SetForegroundColour("#FFFFFF")
        
        if (platform.system() == "Windows"):
            self.lblFASTA = wx.StaticText(self, -1, "FASTA Sequence:", (10, 90), (320, 20), wx.ALIGN_LEFT)
            self.lblFASTA.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblFASTA = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/compmodel/lblFASTA.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(10, 90), size=(320, 20))
        else:
            self.lblFASTA = wx.StaticText(self, -1, "FASTA Sequence:", (10, 90), style=wx.ALIGN_LEFT)
            self.lblFASTA.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.lblFASTA.SetForegroundColour("#FFFFFF")
        self.txtFASTA = wx.TextCtrl(self, -1, pos=(0, 110), size=(320, 75), style=wx.TE_MULTILINE)
        self.txtFASTA.SetValue("")
        self.txtFASTA.SetToolTipString("FASTA sequence of the protein that will be modeled")
        self.txtFASTA.Bind(wx.EVT_KILL_FOCUS, self.fastaChange)
        self.FASTA = ""
        
        if (platform.system() == "Windows"):
            self.lblProt2 = wx.StaticText(self, -1, "Homologous Templates", (0, 190), (320, 25), wx.ALIGN_CENTRE)
            self.lblProt2.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblProt2 = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/compmodel/lblHomologyTemplates.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(25, 190), size=(270, 25))
        else:
            self.lblProt2 = wx.StaticText(self, -1, "Homologous Templates", (70, 190), style=wx.ALIGN_CENTRE)
            self.lblProt2.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
            resizeTextControlForUNIX(self.lblProt2, 0, self.GetSize()[0]-20)
        self.lblProt2.SetForegroundColour("#FFFFFF")
        
        if (platform.system() == "Windows"):
            self.lblInst2 = wx.StaticText(self, -1, "Provide scaffolds that the above sequence will\nbe threaded onto.  Templates higher in the grid\nare given more priority.  Click on templates in\nthe grid to edit their alignments.", (0, 220), (320, 25), wx.ALIGN_CENTRE)
            self.lblInst2.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
        elif (platform.system() == "Darwin"):
            self.lblInst2 = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/compmodel/lblInstHomologyTemplates.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 220), size=(320, 75))
        else:
            self.lblInst2 = wx.StaticText(self, -1, "Provide scaffolds that the above sequence will\nbe threaded onto.  Templates higher in the grid\nare given more priority.  Click on templates in\nthe grid to edit their alignments.", (5, 220), style=wx.ALIGN_CENTRE)
            self.lblInst2.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
            resizeTextControlForUNIX(self.lblInst2, 0, self.GetSize()[0]-20)
        self.lblInst2.SetForegroundColour("#FFFFFF")
        
        if (platform.system() == "Windows"):
            self.lblModel = wx.StaticText(self, -1, "Templates", (0, 290), (155, 20), wx.ALIGN_CENTRE)
            self.lblModel.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblModel = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/compmodel/lblModelCompModel.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 290), size=(155, 20))
        else:
            self.lblModel = wx.StaticText(self, -1, "Templates", (0, 290), style=wx.ALIGN_CENTRE)
            self.lblModel.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblModel, 0, 155)
        self.lblModel.SetForegroundColour("#FFFFFF")
        self.modelMenu = wx.ComboBox(self, pos=(0, 310), size=(155, 25), choices=[], style=wx.CB_READONLY)
        self.modelMenu.SetToolTipString("Select the models to serve as templates for structure prediction")
        self.templates = []
        self.alignments = []
        
        if (platform.system() == "Darwin"):
            self.btnAddTemplate = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/compmodel/btnAddTemplate.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(165, 310), size=(77, 25))
        else:
            self.btnAddTemplate = wx.Button(self, id=-1, label="Add", pos=(165, 310), size=(77, 25))
            self.btnAddTemplate.SetForegroundColour("#000000")
            self.btnAddTemplate.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnAddTemplate.Bind(wx.EVT_BUTTON, self.addModel)
        self.btnAddTemplate.SetToolTipString("Add model to the list of templates")
        if (platform.system() == "Darwin"):
            self.btnRemoveTemplate = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/compmodel/btnRemoveTemplate.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(242, 310), size=(78, 25))
        else:
            self.btnRemoveTemplate = wx.Button(self, id=-1, label="Remove", pos=(242, 310), size=(78, 25))
            self.btnRemoveTemplate.SetForegroundColour("#000000")
            self.btnRemoveTemplate.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnRemoveTemplate.Bind(wx.EVT_BUTTON, self.removeModel)
        self.btnRemoveTemplate.SetToolTipString("Remove model from the list of templates")
        
        self.grdTemplates = wx.grid.Grid(self)
        self.grdTemplates.CreateGrid(0, 1)
        self.grdTemplates.SetSize((320, 150))
        self.grdTemplates.SetPosition((0, 340))
        self.grdTemplates.SetLabelFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.grdTemplates.DisableDragColSize()
        self.grdTemplates.DisableDragRowSize()
        self.grdTemplates.SetColLabelValue(0, "Templates")
        self.grdTemplates.SetRowLabelSize(160)
        self.grdTemplates.SetColSize(0, 160)
        self.grdTemplates.Bind(wx.grid.EVT_GRID_CELL_LEFT_CLICK, self.gridClick)
        self.selectedr = -1
        ypos = self.grdTemplates.GetPosition()[1] + self.grdTemplates.GetSize()[1] + 10
        
        self.btnMoveUp = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/up-arrow.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, ypos), size=(25, 25))
        self.btnMoveUp.Bind(wx.EVT_BUTTON, self.moveUp)
        self.btnMoveUp.SetToolTipString("Add model to the list of templates")
        self.btnMoveDown = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/down-arrow.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(295, ypos), size=(25, 25))
        self.btnMoveDown.Bind(wx.EVT_BUTTON, self.moveDown)
        self.btnMoveDown.SetToolTipString("Remove model from the list of templates")
        if (platform.system() == "Darwin"):
            self.btnLoadAlign = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/compmodel/btnLoadAlign.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(30, ypos), size=(127, 25))
        else:
            self.btnLoadAlign = wx.Button(self, id=-1, label="Load Alignment", pos=(30, ypos), size=(127, 25))
            self.btnLoadAlign.SetForegroundColour("#000000")
            self.btnLoadAlign.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnLoadAlign.Bind(wx.EVT_BUTTON, self.loadAlign)
        self.btnLoadAlign.SetToolTipString("Load a previously-saved alignment for this template-target pair")
        if (platform.system() == "Darwin"):
            self.btnSaveAlign = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/compmodel/btnSaveAlign.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(163, ypos), size=(127, 25))
        else:
            self.btnSaveAlign = wx.Button(self, id=-1, label="Save Alignment", pos=(163, ypos), size=(127, 25))
            self.btnSaveAlign.SetForegroundColour("#000000")
            self.btnSaveAlign.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnSaveAlign.Bind(wx.EVT_BUTTON, self.saveAlign)
        self.btnSaveAlign.SetToolTipString("Load a previously-saved alignment for this template-target pair")
        
        if (platform.system() == "Windows"):
            self.lblAlignment = wx.StaticText(self, -1, "Alignment:", (10, ypos+30), (320, 20), wx.ALIGN_LEFT)
            self.lblAlignment.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblAlignment = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/compmodel/lblAlignment.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(10, ypos+30), size=(320, 20))
        else:
            self.lblAlignment = wx.StaticText(self, -1, "Alignment:", (10, ypos+30), style=wx.ALIGN_LEFT)
            self.lblAlignment.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.lblAlignment.SetForegroundColour("#FFFFFF")
        self.txtAlignment = wx.TextCtrl(self, -1, pos=(0, ypos+50), size=(320, 125), style=wx.TE_MULTILINE)
        self.txtAlignment.SetValue("")
        self.txtAlignment.SetFont(wx.Font(10, wx.FONTFAMILY_TELETYPE, wx.NORMAL, wx.NORMAL))
        self.txtAlignment.SetToolTipString("Alignment of the currently selected model")
        self.txtAlignment.Bind(wx.EVT_KEY_DOWN, self.alignmentEdit)
        self.maxsize = 30
        
        if (platform.system() == "Windows"):
            self.lblNumModels = wx.StaticText(self, -1, "Models to Generate:", (0, ypos+193), (260, 20), wx.ALIGN_CENTRE)
            self.lblNumModels.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblNumModels = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/compmodel/lblNumModels.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, ypos+193), size=(260, 20))
        else:
            self.lblNumModels = wx.StaticText(self, -1, "Models to Generate:", (0, ypos+193), style=wx.ALIGN_CENTRE)
            self.lblNumModels.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblNumModels, 0, 260)
        self.lblNumModels.SetForegroundColour("#FFFFFF")
        self.txtNumModels = wx.TextCtrl(self, -1, pos=(260, ypos+190), size=(60, 25))
        self.txtNumModels.SetValue("10")
        self.txtNumModels.SetToolTipString("Number of comparative models to generate (1-100)")
        
        if (platform.system() == "Darwin"):
            self.btnThread = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/compmodel/btnThread.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(110, ypos+220), size=(100, 25))
        else:
            self.btnThread = wx.Button(self, id=-1, label="Thread!", pos=(110, ypos+220), size=(100, 25))
            self.btnThread.SetForegroundColour("#000000")
            self.btnThread.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.btnThread.Bind(wx.EVT_BUTTON, self.threadClick)
        self.btnThread.SetToolTipString("Begin threading the sequence onto the template structures")
        self.buttonState = "Thread!"
        
        self.scrollh = self.btnThread.GetPosition()[1] + self.btnThread.GetSize()[1] + 5
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

    def scrolled(self, event):
        self.winscrollpos = self.GetScrollPos(wx.VERTICAL)
        event.Skip()
        
    def activate(self):
        if (self.seqWin):                
            # Select all models
            if (self.seqWin.numChains() > 0):
                rowsSelected = range(0, self.seqWin.numChains())
            else:
                return
            self.models = []
            for r in rowsSelected:
                ID = self.seqWin.IDs[r]
                self.models.append(str(ID))
            if (self.models != self.modelMenu.GetItems()):
                # Update the combo box to have these models as selections
                self.modelMenu.Clear()
                self.modelMenu.AppendItems(self.models)
                if (len(self.models) > 0):
                    if (platform.system() == "Windows"):
                        # Doing this on Linux freezes up the ComboBox, so don't do it on Linux
                        self.modelMenu.SetSelection(0)
                # Take out templates that are not loaded
                self.txtAlignment.SetValue("")
                self.selectedr = -1
                for i in range(len(self.templates)-1, -1, -1):
                    if (self.templates[i] not in self.models):
                        self.templates.pop(i)
                        self.alignments.pop(i)
                self.updateGrid()
        self.Scroll(0, self.winscrollpos)
        
    def recolorModel(self, modelchain):
        # Recolor the model such that aligned portions of the model are shown in purple
        # and unaligned structure is in cyan
        # Purple: Aligned regions, will grab from template in comparative modeling
        # Cyan: Deletions, will not be present at all in the final model
        # Yellow: Insertion, will attempt to add a loop in this region
        addToSelection = False
        model = modelchain[0:len(modelchain)-2]
        chain = modelchain[len(modelchain)-1]
        r = self.seqWin.IDs.index(modelchain)
        c = 1
        # Save the locations of where the insertions will occur
        insertions = []
        for i in range(0, len(self.alignments[self.selectedr][0])):
            if (self.alignments[self.selectedr][0][i] != "-" and self.alignments[self.selectedr][1][i] != "-"):
                self.seqWin.selectNthResidue(r, c, addToSelection, updateSelection=False)
                addToSelection = True
            if (self.alignments[self.selectedr][0][i] != "-" and self.alignments[self.selectedr][1][i] == "-" and c not in insertions):
                insertions.append(c-1)
            if (self.alignments[self.selectedr][1][i] != "-"):
                c += 1
        self.seqWin.selectUpdate(True)
        self.cmd.set("ribbon_color", "purple", "seqsele")
        self.cmd.set("cartoon_color", "purple", "seqsele")
        if (chain != "_"):
            self.cmd.set("ribbon_color", "cyan", "model " + model + " and chain " + chain + " and not seqsele")
            self.cmd.set("cartoon_color", "cyan", "model " + model + " and chain " + chain + " and not seqsele")
        else:
            self.cmd.set("ribbon_color", "cyan", "model " + model + " and not seqsele")
            self.cmd.set("cartoon_color", "cyan", "model " + model + " and not seqsele")
        addToSelection = False
        for c in insertions:
            self.seqWin.selectNthResidue(r, c, addToSelection, updateSelection=False)
            addToSelection = True
        self.seqWin.selectUpdate(True)
        self.cmd.set("ribbon_color", "yellow", "seqsele")
        self.cmd.set("cartoon_color", "yellow", "seqsele")
        self.cmd.disable("seqsele")
        self.cmd.delete("seqsele")
        
    def showModel(self, modelchain):
        model = modelchain[0:len(modelchain)-2]
        chain = modelchain[len(modelchain)-1]
        defaultPyMOLView(self.cmd)
        if (chain != "_"):
            self.cmd.select("sele", "model " + model + " and chain " + chain)
        else:
            self.cmd.select("sele", "model " + model)
        self.cmd.zoom("sele")
        self.cmd.hide("everything", "not sele")
        self.recolorModel(modelchain)
        #self.seqWin.selectUpdate(False)
    
    def modelMenuSelect(self, event):
        logInfo("Selected model " + self.modelMenu.GetValue())
        #self.showModel(self.modelMenu.GetStringSelection())
    
    def updateGrid(self):
        scrollpos = self.grdTemplates.GetScrollPos(wx.VERTICAL)
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
        self.grdTemplates.Scroll(0, scrollpos)
    
    def drawAlignment(self, indx):
        # Renders the alignment in txtAlignment in formal output
        # Columns 1-4: Targ/1st 4 letters of scaffold
        # Column 5: Blank
        # Columns 6+: Alignment
        self.alignmentstr = ""
        # First look to see if we have any gaps aligned to other gaps and remove them
        for i in range(len(self.alignments[indx][0])-1, -1, -1):
            if (self.alignments[indx][0][i] == "-" and self.alignments[indx][1][i] == "-"):
                self.alignments[indx][0] = self.alignments[indx][0][0:i] + self.alignments[indx][0][i+1:]
                self.alignments[indx][1] = self.alignments[indx][1][0:i] + self.alignments[indx][1][i+1:]
        for i in range(0, len(self.alignments[indx][0]), self.maxsize):
            self.alignmentstr += "Targ " + self.alignments[indx][0][i:i+self.maxsize] + "\n"
            self.alignmentstr += self.templates[indx][0:4] + " " + self.alignments[indx][1][i:i+self.maxsize] + "\n\n"
        # On OSX, it will take consecutive "-" and combine them into one, which is super annoying
        # So replace "-" with "~"
        self.txtAlignment.SetValue(self.alignmentstr.replace("-", "~"))
    
    def alignmentEdit(self, event):
        (selpos, end) = self.txtAlignment.GetSelection()
        if (selpos != end):
            return
        # Now we have to do some math and make sure the insertion is in the actual alignment
        # and not in header information or the blank line
        if (platform.system() == "Windows"):
            # Apparently in Windows text boxes newlines register as 2 characters instead of 1 -_-
            linewidth = 16 + 2*self.maxsize # Two alignments, each line has 5 header columns, and 3 newlines
            targstart = 5
            tempstart = 12+self.maxsize
            offset = 1
        else:
            linewidth = 13 + 2*self.maxsize # Two alignments, each line has 5 header columns, and 3 newlines
            targstart = 5
            tempstart = 11+self.maxsize
            offset = 0
        linepos = selpos
        block = 0
        while (linepos >= linewidth):
            linepos -= linewidth
            block += 1
        # Only the SPACEBAR, BACKSPACE, and DELETE keys change the sequence
        # Arrow keys are supported to move the cursor around
        if (int(event.GetKeyCode()) == wx.WXK_SPACE):
            if (linepos < 5 or (linepos > 5+self.maxsize and linepos < 11+self.maxsize) or linepos > 11+2*self.maxsize):
                return
            # Put a gap in this location and add a gap to the end of the other sequence
            shift = True
            if (linepos <= targstart+self.maxsize):
                # Selection in target sequence
                p = linepos - targstart + block * self.maxsize
                seq = 0
            else:
                p = linepos - tempstart + block * self.maxsize
                seq = 1
            self.alignments[self.selectedr][seq] = self.alignments[self.selectedr][seq][0:p] + "-" + self.alignments[self.selectedr][seq][p:]
            if (self.alignments[self.selectedr][seq].endswith("-")):
                self.alignments[self.selectedr][seq] = self.alignments[self.selectedr][seq][0:len(self.alignments[self.selectedr][seq])-1]
            else:
                self.alignments[self.selectedr][1-seq] += "-"
            if (self.alignments[self.selectedr][1-seq][p] == "-"):
                # The doubly-aligned - will be eliminated by drawAlignment
                shift = False
            # Redrawing the alignment resets the scrollbars, this gets the original scroll position back
            vpos = self.txtAlignment.GetScrollPos(wx.VERTICAL)
            self.drawAlignment(self.selectedr)
            self.txtAlignment.SetScrollPos(wx.VERTICAL, vpos)
            if (shift):
                # Move the cursor
                if (linepos == targstart+self.maxsize or linepos == targstart+self.maxsize):
                    # We're at the edge of a line and need to shift by a whole block
                    linepos -= self.maxsize
                    block += 1
                linepos += 1
            selpos = linepos + block * linewidth
            self.txtAlignment.SetSelection(selpos, selpos)
            # Alignment changed, update coloring
            self.recolorModel(self.templates[self.selectedr])
        elif (int(event.GetKeyCode()) == wx.WXK_BACK or int(event.GetKeyCode()) == wx.WXK_DELETE or (platform.system() == "Darwin" and int(event.GetKeyCode()) == 8)):
            # Backspace deletes the gap before the current position and delete deletes the current gap
            if (int(event.GetKeyCode()) == wx.WXK_BACK or (platform.system() == "Darwin" and int(event.GetKeyCode()) == wx.WXK_DELETE)):
                if (linepos < targstart or (linepos > targstart+self.maxsize and linepos < tempstart) or linepos > tempstart+self.maxsize):
                    return
                # Put a gap in this location and add a gap to the end of the other sequence
                shift = True
                if (linepos <= targstart+self.maxsize):
                    # Selection in target sequence
                    p = linepos - targstart + block * self.maxsize
                    seq = 0
                else:
                    p = linepos - tempstart + block * self.maxsize
                    seq = 1
                if (p-1 < 0 or self.alignments[self.selectedr][seq][p-1] != "-"):
                    # You can only delete gaps
                    return
                self.alignments[self.selectedr][seq] = self.alignments[self.selectedr][seq][0:p-1] + self.alignments[self.selectedr][seq][p:] + "-"
            else:
                if (linepos < targstart or (linepos > targstart+self.maxsize and linepos < tempstart) or linepos > tempstart+self.maxsize):
                    return
                # Put a gap in this location and add a gap to the end of the other sequence
                shift = False
                if (linepos <= targstart+self.maxsize):
                    # Selection in target sequence
                    p = linepos - targstart + block * self.maxsize
                    seq = 0
                else:
                    p = linepos - tempstart + block * self.maxsize
                    seq = 1
                if (p >= len(self.alignments[self.selectedr][seq]) or self.alignments[self.selectedr][seq][p] != "-"):
                    # You can only delete gaps
                    return
                self.alignments[self.selectedr][seq] = self.alignments[self.selectedr][seq][0:p] + self.alignments[self.selectedr][seq][p+1:] + "-"
            # Redrawing the alignment resets the scrollbars, this gets the original scroll position back
            vpos = self.txtAlignment.GetScrollPos(wx.VERTICAL)
            self.drawAlignment(self.selectedr)
            self.txtAlignment.SetScrollPos(wx.VERTICAL, vpos)
            if (shift):
                # Move the cursor
                if (linepos == targstart or linepos == tempstart):
                    # We're at the edge of a line and need to shift by a whole block
                    linepos += self.maxsize
                    block -= 1
                linepos -= 1
            selpos = linepos + block * linewidth
            self.txtAlignment.SetSelection(selpos, selpos)
            # Alignment changed, update coloring
            self.recolorModel(self.templates[self.selectedr])
        elif (int(event.GetKeyCode()) == wx.WXK_LEFT):
            if (linepos <= targstart and block == 0):
                linepos = targstart
            elif (linepos <= tempstart and linepos > targstart+self.maxsize and block == 0):
                linepos = tempstart
            elif (linepos < targstart):
                linepos = targstart
            elif (linepos < tempstart and linepos > targstart+self.maxsize):
                linepos = tempstart
            elif (linepos > tempstart+self.maxsize):
                linepos = tempstart+self.maxsize
            elif (linepos == targstart or linepos == tempstart):
                linepos += self.maxsize
                block -= 1
            else:
                linepos -= 1
            selpos = linepos + block * linewidth
            self.txtAlignment.SetSelection(selpos, selpos)
        elif (int(event.GetKeyCode()) == wx.WXK_RIGHT):
            if (linepos < targstart):
                linepos = targstart
            elif (linepos < tempstart and linepos > targstart+self.maxsize):
                linepos = tempstart
            elif (linepos > tempstart+self.maxsize):
                linepos = targstart
                block += 1
            elif (linepos == targstart+self.maxsize or linepos == tempstart+self.maxsize):
                # We're at the edge of a line and need to shift by a whole block
                linepos -= self.maxsize
                block += 1
            else:
                linepos += 1
            # Okay, the last line is weird since it doesn't have to be self.maxsize in length
            # Did we go past the end of a sequence?
            if (linepos <= targstart+self.maxsize):
                # Selection in target sequence
                p = linepos - targstart + block * self.maxsize
                seq = 0
            else:
                p = linepos - tempstart + block * self.maxsize
                seq = 1
            if (p >= len(self.alignments[self.selectedr][seq])):
                p = len(self.alignments[self.selectedr][seq])
                if (seq == 0):
                    linepos = p - (block * self.maxsize) + targstart
                else:
                    linepos = p - (block * self.maxsize) + tempstart
            selpos = linepos + block * linewidth
            if (selpos >= len(self.alignmentstr)):
                # Block is pushing it passed the end of the string
                selpos = linepos + (block-1) * linewidth
            self.txtAlignment.SetSelection(selpos, selpos)
        elif (int(event.GetKeyCode()) == wx.WXK_UP):
            # Again, that last line is annoying since it doesn't have to be maxsize in length
            lastblock = int(len(self.alignments[self.selectedr][0]) / self.maxsize)
            if (block == lastblock):
                maxsize = len(self.alignments[self.selectedr][0]) - lastblock * self.maxsize
                if (linepos < targstart+maxsize):
                    # Still on the top row, so maxsize should not be changed since we're going
                    # up into the previous block which has a full alignment row
                    maxsize = self.maxsize
                else:
                    tempstart = tempstart - self.maxsize + maxsize
            else:
                maxsize = self.maxsize
            # Move up to the next higher sequence if possible
            if (linepos <= targstart and block == 0):
                # In the header region on the first line
                linepos = targstart
            elif (linepos <= tempstart and linepos > targstart+maxsize and block == 0):
                # In the header region on the first alignment, lower sequence
                linepos = targstart
            elif (linepos <= targstart):
                # In header region of target sequence
                linepos = tempstart + offset
                block -= 1
            elif (linepos < tempstart and linepos > targstart+maxsize):
                # In header region of the scaffold sequence
                linepos = targstart
                #block -= 1
            elif (linepos > tempstart+maxsize):
                # In the gap between alignment blocks
                linepos = 11+maxsize + offset
            elif (linepos <= targstart+maxsize and block > 0):
                # In target sequence, move up to scaffold sequence in previous block
                linepos += targstart + maxsize + 1 + offset 
                block -= 1
            elif (linepos >= tempstart and linepos <= tempstart+maxsize):
                # In scaffold sequence, move up to the target in the same block
                linepos -= targstart + maxsize + 1 + offset
            selpos = linepos + block * linewidth
            self.txtAlignment.SetSelection(selpos, selpos)
        elif (int(event.GetKeyCode()) == wx.WXK_DOWN):
            # Again, that last line is annoying since it doesn't have to be maxsize in length
            lastblock = int(len(self.alignments[self.selectedr][0]) / self.maxsize)
            lastmaxsize = len(self.alignments[self.selectedr][0]) - lastblock * self.maxsize
            # Move down to the next lower sequence if possible
            if (linepos <= targstart and block == 0):
                # In the header region on the first line
                linepos = tempstart
            elif (linepos <= targstart):
                # In header region of target sequence
                linepos = tempstart
            elif (linepos < tempstart and linepos > targstart+self.maxsize):
                # In header region of the scaffold sequence
                linepos = targstart
                block += 1
            elif (linepos > tempstart+self.maxsize and block != lastblock):
                # In the gap between alignment blocks, not in the last alignment block
                linepos = targstart
                block += 1
            elif (linepos <= targstart+self.maxsize and block < lastblock):
                # In target sequence, move down to scaffold sequence in same block
                linepos += targstart + self.maxsize + 1 + offset
            elif (linepos <= targstart+self.maxsize and block == lastblock):
                # In target sequence in the last block, shift by the length of the last block
                linepos += targstart + lastmaxsize + 1 + offset
            elif (linepos >= tempstart and linepos <= tempstart+self.maxsize and block < lastblock-1):
                # In scaffold sequence, move down to the target in the next block
                linepos -= targstart + self.maxsize + 1 + offset
                block += 1
            elif (linepos >= tempstart and linepos <= tempstart+self.maxsize and block == lastblock-1):
                # In scaffold sequence but moving to last block, check to make sure we don't go past the edge of the alignment
                linepos -= targstart + self.maxsize + 1 + offset
                block += 1
                if (linepos > targstart+lastmaxsize):
                    linepos = targstart+lastmaxsize
            else:
                # In scaffold of last block, move the cursor to the end of the scaffold
                linepos = tempstart+lastmaxsize
            selpos = linepos + block * linewidth
            self.txtAlignment.SetSelection(selpos, selpos)
        # Select this residue (or the one before if in a gap) in the sequence window
        (selpos, end) = self.txtAlignment.GetSelection()
        linepos = selpos
        block = 0
        while (linepos >= linewidth):
            linepos -= linewidth
            block += 1
        if (linepos <= targstart+self.maxsize):
            # Selection in target sequence
            p = linepos - targstart + block * self.maxsize
        else:
            p = linepos - tempstart + block * self.maxsize
        # Find the row
        for i in range(0, len(self.seqWin.IDs)):
            if (self.templates[self.selectedr] == self.seqWin.IDs[i]):
                r = i
                break
        # Find the residue
        c = 0
        for i in range(0, p+1):
            if (c >= len(self.alignments[self.selectedr][1])):
                break
            if (self.alignments[self.selectedr][1][i] != "-"):
                c += 1
        self.seqWin.selectNthResidue(r, c)
    
    def regenerateAlignments(self, indx=-1):
        # Uses MUSCLE to recalculate the alignments when a change occurs
        # If MUSCLE is not available, a default alignment will be given and the user will
        # have to align manually
        if (indx < 0):
            indx = 0
            l = len(self.alignments)
        else:
            l = indx + 1
        goToSandbox()
        if (platform.system() == "Windows"):
            muscle = self.parent.parent.scriptdir + "\\bin\\muscle_win.exe"
        elif (platform.system() == "Darwin"):
            muscle = self.parent.parent.scriptdir + "/bin/muscle_darwin"
        else:
            muscle = self.parent.parent.scriptdir + "/bin/muscle_unix"
        for i in range(indx, l):
            # Generate the MUSCLE input
            fout = open("muscle.in", "w")
            fout.write("> Target\n")
            fout.write(self.FASTA.strip() + "\n\n")
            fout.write("> Scaffold\n")
            # Get the sequence
            modelchain = self.templates[i]
            seq = ""
            for r in range(0, len(self.seqWin.sequences)):
                if (modelchain == self.seqWin.IDs[r]):
                    seq += self.seqWin.sequences[r]
            fout.write(seq + "\n")
            fout.close()
            # Try to align with MUSCLE
            try:
                muscle_cline = MuscleCommandline(muscle, input="muscle.in", out="muscle.out")
                muscle_cline()
                # Read the results
                targetseq = ""
                scaffoldseq = ""
                onScaffold = False
                fin = open("muscle.out", "r")
                for aline in fin:
                    if (aline.startswith("> Scaffold")):
                        onScaffold = True
                    elif (len(aline.strip()) == 0 or aline.strip()[0] == ">"):
                        continue
                    elif (onScaffold):
                        scaffoldseq += aline.strip()
                    else:
                        targetseq += aline.strip()
                fin.close()
                self.alignments.append([targetseq, scaffoldseq])
            except:
                print "WARNING: MUSCLE could not be found at " + muscle
                argetseq = self.FASTA
                scaffoldseq = seq
                if (len(targetseq) < len(scaffoldseq)):
                    for j in range(0, len(scaffoldseq)-len(targetseq)):
                        targetseq += "-"
                else:
                    for j in range(0, len(targetseq)-len(scaffoldseq)):
                        scaffoldseq += "-"
                self.alignments.append([targetseq, scaffoldseq])
        if (self.selectedr >= 0):
            self.drawAlignment(self.selectedr)
    
    def fastaChange(self, event):
        # Scan the current contents of this text box and take out everything that is invalid
        data = self.txtFASTA.GetValue()#.strip()
        fastaOnly = ""
        for aline in data.split("\n"):
            if (aline.strip().startswith(">")):
                continue
            fastaOnly += aline
        fastaOnly = str(fastaOnly.replace(" ", "").replace("\t", "").upper())
        newseq = ""
        for AA in fastaOnly:
            if (AA not in "ACDEFGHIKLMNPQRSTVWY"):
                dlg = wx.MessageDialog(self, "The amino acid " + AA + " is not valid.  It will be ignored.", "Invalid Amino Acid", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
                dlg.ShowModal()
                dlg.Destroy()
            else:
                newseq += AA
        if (newseq != self.FASTA):
            self.regenerateAlignments()
        self.FASTA = newseq
        self.txtFASTA.SetValue(self.FASTA)
        self.txtFASTA.SetSelection(len(self.FASTA), len(self.FASTA))
        event.Skip()
    
    def addModel(self, event):
        # Add this to the list of template structures
        model = str(self.modelMenu.GetValue())
        if (model not in self.modelMenu.GetItems()):
            return
        if (not(model in self.templates)):
            self.templates.append(model)
            self.regenerateAlignments(len(self.templates)-1)
        # Update the grid with this new information
        self.updateGrid()
        logInfo("Added " + model + " to the list of model structures")
        
    def removeModel(self, event):
        # Just pop out the indicated chain
        model = self.modelMenu.GetValue()
        if (model in self.templates):
            indx = self.templates.index(model)
            self.templates.pop(indx)
            self.alignments.pop(indx)
        # Update the grid with this new information
        self.updateGrid()
        logInfo("Removed " + model + " from the list of static chains")
        
    def selectRow(self, row):
        self.selectedr = row
        if (self.selectedr >= self.grdTemplates.NumberRows):
            self.selectedr = -1
        for r in range(0, self.grdTemplates.NumberRows):
            if (r == self.selectedr):
                for c in range(0, self.grdTemplates.NumberCols):
                    self.grdTemplates.SetCellBackgroundColour(r, c, "light blue")
            else:
                for c in range(0, self.grdTemplates.NumberCols):
                    self.grdTemplates.SetCellBackgroundColour(r, c, "white")
        self.grdTemplates.Refresh()
        self.showModel(self.templates[self.selectedr])
        
    def gridClick(self, event):
        self.selectRow(event.GetRow())
        if (len(self.FASTA.strip()) > 0):
            self.drawAlignment(self.selectedr)
        event.Skip()
        
    def moveUp(self, event):
        # Move the selected template up in the ranking
        if (self.selectedr > 0):
            temp = self.templates[self.selectedr]
            self.templates[self.selectedr] = self.templates[self.selectedr-1]
            self.templates[self.selectedr-1] = temp
            temp = self.alignments[self.selectedr]
            self.alignments[self.selectedr] = self.alignments[self.selectedr-1]
            self.alignments[self.selectedr-1] = temp
            self.updateGrid()
            self.selectedr -= 1
            self.selectRow(self.selectedr)
    
    def moveDown(self, event):
        # Move the selected template down in the ranking
        if (self.selectedr >= 0 and self.selectedr < len(self.templates)-1):
            temp = self.templates[self.selectedr]
            self.templates[self.selectedr] = self.templates[self.selectedr+1]
            self.templates[self.selectedr+1] = temp
            temp = self.alignments[self.selectedr]
            self.alignments[self.selectedr] = self.alignments[self.selectedr+1]
            self.alignments[self.selectedr+1] = temp
            self.updateGrid()
            self.selectedr += 1
            self.selectRow(self.selectedr)
    
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
    
    def loadAlign(self, event):
        logInfo("Load Alignment button clicked")
        # Is a row selected?
        if (len(self.txtAlignment.GetValue().strip()) == 0):
            dlg = wx.MessageDialog(self, "Please select a template from the grid whose alignment is the one attempting to be loaded.", "No Template Selected", wx.OK | wx.ICON_ERROR | wx.CENTRE)
            dlg.ShowModal()
            dlg.Destroy()
            return
        dlg = wx.FileDialog(
            self, message="Choose a File",
            defaultDir=self.seqWin.cwd,
            defaultFile="",
            wildcard="Alignments (*.fasta)|*.fasta",
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
            logInfo("Loaded data from a FASTA alignment", filename)
        # Read the sequences out of the alignment file
        targseq = ""
        readingtarg = False
        readingtemp = False
        tempseq = ""
        fin = open(filename, "r")
        for aline in fin:
            if (aline.startswith("> Target")):
                readingtarg = True
            elif (readingtarg and aline.startswith(">")):
                # Does the ID in here match the one selected?
                ID = aline[2:].strip()
                if (ID != self.templates[self.selectedr]):
                    dlg = wx.MessageDialog(self, "The PDB IDs in the alignment file do not match.\n\nSelected template: " + self.templates[self.selectedr] + "\nFound: " + ID, "PDB IDs Do Not Match", wx.OK | wx.ICON_ERROR | wx.CENTRE)
                    dlg.ShowModal()
                    dlg.Destroy()
                    return
                readingtemp = True
                readingtarg = False
            elif (readingtemp and aline.startswith(">")):
                break
            elif (readingtarg):
                targseq += aline.strip()
            elif (readingtemp):
                tempseq += aline.strip()
        fin.close()
        # Do the sequences match the sequences of the current FASTA and template?
        currtarg = ""
        currtemp = ""
        loadtarg = ""
        loadtemp = ""
        for AA in self.alignments[self.selectedr][0]:
            if (AA != "-"):
                currtarg += AA
        for AA in self.alignments[self.selectedr][1]:
            if (AA != "-"):
                currtemp += AA
        for AA in targseq:
            if (AA != "-"):
                loadtarg += AA
        for AA in tempseq:
            if (AA != "-"):
                loadtemp += AA
        if (currtarg != loadtarg):
            dlg = wx.MessageDialog(self, "The target sequence in the alignment file does not match the current sequence given above!", "Target Sequences Do Not Match", wx.OK | wx.ICON_ERROR | wx.CENTRE)
            dlg.ShowModal()
            dlg.Destroy()
            print currtarg
            print loadtarg
            return
        if (currtemp != loadtemp):
            dlg = wx.MessageDialog(self, "The template sequence in the alignment file does not match the sequence of the selected template!", "Template Sequences Do Not Match", wx.OK | wx.ICON_ERROR | wx.CENTRE)
            dlg.ShowModal()
            dlg.Destroy()
            print currtemp
            print loadtemp
            return
        # If we made it this far, then the sequences should be good
        self.alignments[self.selectedr][0] = targseq
        self.alignments[self.selectedr][1] = tempseq
        self.drawAlignment(self.selectedr)
        self.recolorModel(self.templates[self.selectedr])
    
    def saveAlign(self, event):
        # Is an alignment even there?
        if (len(self.txtAlignment.GetValue().strip()) == 0):
            dlg = wx.MessageDialog(self, "There is no alignment to save!", "No Alignment", wx.OK | wx.ICON_ERROR | wx.CENTRE)
            dlg.ShowModal()
            dlg.Destroy()
            return
        # Okay, now just save the alignment in FASTA format for later use
        dlg = wx.FileDialog(
            self, message="Save an alignment",
            defaultDir=self.seqWin.cwd,
            defaultFile="",
            wildcard="Alignments (*.fasta)|*.fasta",
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
            filename = str(paths[0]).split(".fasta")[0] + ".fasta"
            # Does it exist already?  If so, ask if the user really wants to overwrite it
            if (os.path.isfile(filename)):
                dlg2 = wx.MessageDialog(self, "The file " + filename + " already exists.  Overwrite it?", "Filename Already Exists", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
                if (dlg2.ShowModal() == wx.ID_NO):
                    dlg2.Destroy()
                    logInfo("Canceled save operation due to filename already existing")
            fout = open(filename, "w")
            fout.write("> Target\n\n")
            for i in range(0, len(self.alignments[self.selectedr][0]), 50):
                fout.write(self.alignments[self.selectedr][0][i:i+50] + "\n")
            fout.write("\n> " + self.templates[self.selectedr].strip() + "\n\n")
            for i in range(0, len(self.alignments[self.selectedr][1]), 50):
                fout.write(self.alignments[self.selectedr][1][i:i+50] + "\n")
            fout.close()
    
    def threadClick(self, event):
        logInfo("Clicked the Thread button")
        # This is also the "Finalize!" button
        if (self.buttonState == "Thread!"):
            # Do we have a FASTA sequence?
            if (len(self.FASTA.strip()) == 0):
                wx.MessageBox("You need to input a FASTA sequence whose structure will be predicted!", "FASTA Sequence Required", wx.OK|wx.ICON_ERROR)
                return
            # Do we have at least one template?
            if (len(self.templates) == 0):
                wx.MessageBox("You need to provide at least one homologous structure that will be used for modeling!", "Homologous Structure Required", wx.OK|wx.ICON_ERROR)
                return
            # Do we have a valid number of models?
            try:
                nmodels = int(self.txtNumModels.GetValue())
                if (nmodels < 1 or nmodels > 100):
                    raise Exception()
            except:
                wx.MessageBox("Please enter an integer from 1-100 for the number of models.", "Invalid Number of Models", wx.OK|wx.ICON_ERROR)
                return
            goToSandbox()
            # Write the input file
            fout = open("threadinputtemp", "w")
            fout.write("FASTA\t" + self.FASTA.strip() + "\n")
            fout.write("NMODELS\t" + str(nmodels) + "\n")
            indx = 1
            for template in self.templates:
                chain = template[len(template)-1]
                fout.write("PDBFILE\tPRT" + str(indx) + chain + ".pdb\n")
                indx += 1
                fout.write("BEGIN PDB DATA\n")
                fin = open(template[0:len(template)-2] + ".pdb", "r")
                for aline in fin:
                    if ((aline.startswith("ATOM") or aline.startswith("HETATM")) and aline[21] == template[len(template)-1]):
                        fout.write(aline)
                fin.close()
                fout.write("END PDB DATA\n")
            for i in range(0, len(self.alignments)):
                fout.write("TARGALIGN\t" + self.alignments[i][0].strip() + "\n")
                fout.write("TEMPALIGN\t" + self.alignments[i][1].strip() + "\n")
            fout.close()
            appendScorefxnParamsInfoToFile("threadinputtemp", self.selectWin.weightsfile)
            # Try to send it to the server
            try:
                self.ID = sendToServer("threadinput")
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
                        if (len(aline.split("\t")) >= 2 and aline.split("\t")[0] == "COMPMODEL" and aline.split("\t")[1] == self.ID.strip()):
                            alreadythere = True
                            break
                    f.close()
                except:
                    pass
                if (not(alreadythere)):
                    f = open("downloadwatch", "a")
                    f.write("COMPMODEL\t" + self.ID.strip() + "\t" + str(datetime.datetime.now().strftime("%A, %B %d - %I:%M:%S %p")) + "\t" + getServerName() + "\t" + desc + "\n")
                    f.close()
                dlg = wx.MessageDialog(self, "InteractiveROSETTA is now watching the server for job ID " + desc.strip() + ".  You will be notified when the package is available for download.", "Listening for Download", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
                dlg.ShowModal()
                dlg.Destroy()
            except:
                dlg = wx.MessageDialog(self, "The comparative modeling job could not be sent to the server!  Verify that you have an Internet connection and that the server is running.", "Server Could Not Be Reached", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
                dlg.ShowModal()
                dlg.Destroy()
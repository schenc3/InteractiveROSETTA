import wx
import wx.grid
import wx.lib.scrolledpanel
import os
import os.path
import platform
from tools import *
import webbrowser
from threading import Thread

class SuperimpositionPanel(wx.lib.scrolledpanel.ScrolledPanel):
    def __init__(self, parent, W, H):
        #if (platform.system() == "Windows"):
        wx.lib.scrolledpanel.ScrolledPanel.__init__(self, parent, id=-1, pos=(10, 60), size=(340, H-330), name="ProtSuperimposition")
        winh = H-330
        #else:
        #    wx.lib.scrolledpanel.ScrolledPanel.__init__(self, parent, id=-1, pos=(10, 60), size=(340, H-290), name="ProtSuperimposition")
        #    winh = H-290
        self.SetBackgroundColour("#333333")
        self.parent = parent

        if (platform.system() == "Windows"):
            self.lblProt = wx.StaticText(self, -1, "Superimposition", (25, 15), (270, 25), style=wx.ALIGN_CENTRE)
            self.lblProt.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblProt = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/superimposition/lblSuperimposition.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(25, 15), size=(270, 25))
        else:
            self.lblProt = wx.StaticText(self, -1, "Superimposition", pos=(90, 15), style=wx.ALIGN_CENTRE)
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
            self.lblInst = wx.StaticText(self, -1, "Pre-select models/residues to superimpose", (0, 45), (320, 25), wx.ALIGN_CENTRE)
            self.lblInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
        elif (platform.system() == "Darwin"):
            self.lblInst = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/superimposition/lblInstSuperimposition.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 45), size=(320, 25))
        else:
            self.lblInst = wx.StaticText(self, -1, "Pre-select models/residues to superimpose", pos=(20, 45), style=wx.ALIGN_CENTRE)
            self.lblInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
            resizeTextControlForUNIX(self.lblInst, 0, self.GetSize()[0])
        self.lblInst.SetForegroundColour("#FFFFFF")

        if (platform.system() == "Windows"):
            self.lblBase = wx.StaticText(self, -1, "Superimposition\nBase Model:", (0, 78), (175, 30), wx.ALIGN_CENTRE)
            self.lblBase.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblBase = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/superimposition/lblBase.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 78), size=(175, 40))
        else:
            self.lblBase = wx.StaticText(self, -1, "Superimposition\nBase Model:", pos=(30, 78), style=wx.ALIGN_CENTRE)
            self.lblBase.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblBase, 0, 175)
        self.lblBase.SetForegroundColour("#FFFFFF")
        self.modelMenu = wx.ComboBox(self, pos=(200, 80), size=(110, 25), choices=[], style=wx.CB_READONLY)
        self.modelMenu.SetToolTipString("Model to be used as the immovable reference structure")

        if (platform.system() == "Darwin"):
            self.btnByModel = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/superimposition/btnByModel_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(20, 120), size=(120, 25))
        else:
            self.btnByModel = wx.Button(self, id=-1, label="By Model", pos=(20, 120), size=(120, 25))
            self.btnByModel.SetForegroundColour("#FF0000")
            self.btnByModel.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnByModel.Bind(wx.EVT_BUTTON, self.modelClick)
        self.btnByModel.SetToolTipString("Superimpose whole structures")
        if (platform.system() == "Darwin"):
            self.btnByResidues = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/superimposition/btnByResidues.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(170, 120), size=(120, 25))
        else:
            self.btnByResidues = wx.Button(self, id=-1, label="By Residues", pos=(170, 120), size=(120, 25))
            self.btnByResidues.SetForegroundColour("#000000")
            self.btnByResidues.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnByResidues.Bind(wx.EVT_BUTTON, self.residuesClick)
        self.btnByResidues.SetToolTipString("Superimpose by only the selected residues")
        self.alignBy = "Model"

        if (platform.system() == "Darwin"):
            self.btnCARMSD = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/superimposition/btnCARMSD_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(20, 155), size=(120, 25))
        else:
            self.btnCARMSD = wx.Button(self, id=-1, label="CA RMSD", pos=(20, 155), size=(120, 25))
            self.btnCARMSD.SetForegroundColour("#FF0000")
            self.btnCARMSD.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnCARMSD.Bind(wx.EVT_BUTTON, self.CAclick)
        self.btnCARMSD.SetToolTipString("Superimpose by CA atoms only")
        if (platform.system() == "Darwin"):
            self.btnBBRMSD = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/superimposition/btnBBRMSD.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(170, 155), size=(120, 25))
        else:
            self.btnBBRMSD = wx.Button(self, id=-1, label="BB RMSD", pos=(170, 155), size=(120, 25))
            self.btnBBRMSD.SetForegroundColour("#000000")
            self.btnBBRMSD.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnBBRMSD.Bind(wx.EVT_BUTTON, self.BBclick)
        self.btnBBRMSD.SetToolTipString("Superimpose by all backbone atoms")
        self.RMSDtype = "CA"

        self.grdRMSDResults = wx.grid.Grid(self)
        self.grdRMSDResults.CreateGrid(0, 1)
        if (winh-265 > 200):
            self.grdRMSDResults.SetSize((320, winh-265))
        else:
            self.grdRMSDResults.SetSize((320, 200))
        self.grdRMSDResults.SetPosition((0, 190))
        self.grdRMSDResults.SetLabelFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.grdRMSDResults.DisableDragColSize()
        self.grdRMSDResults.DisableDragRowSize()
        self.grdRMSDResults.SetColSize(0, int(self.grdRMSDResults.GetSize()[0] / 2))
        self.grdRMSDResults.SetRowLabelSize(int(self.grdRMSDResults.GetSize()[0] / 2))
        self.grdRMSDResults.SetColLabelValue(0, "CA RMSD")

        ypos = self.grdRMSDResults.GetPosition()[1] + self.grdRMSDResults.GetSize()[1] + 10
        if (platform.system() == "Darwin"):
            self.btnSuperimpose = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/superimposition/btnSuperimpose.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(80, ypos), size=(150, 25))
        else:
            self.btnSuperimpose = wx.Button(self, id=-1, label="Superimpose!", pos=(80, ypos), size=(150, 25))
            self.btnSuperimpose.SetForegroundColour("#000000")
            self.btnSuperimpose.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.btnSuperimpose.Bind(wx.EVT_BUTTON, self.superimposeClick)
        self.btnSuperimpose.SetToolTipString("Begin superimposition")
        self.models = []

        self.scrollh = self.btnSuperimpose.GetPosition()[1] + self.btnSuperimpose.GetSize()[1] + 5
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
            browser.open(self.parent.parent.scriptdir + "/help/superimposition.html")
        else:
            webbrowser.open(self.parent.parent.scriptdir + "/help/superimposition.html")

    def setSeqWin(self, seqWin):
        self.seqWin = seqWin

    def setPyMOL(self, pymol):
        self.pymol = pymol
        self.cmd = pymol.cmd
        self.stored = pymol.stored

    def scrolled(self, event):
        self.winscrollpos = self.GetScrollPos(wx.VERTICAL)
        event.Skip()

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
                    self.baseModel = ""
                    self.models = []
                    return
            self.models = []
            for r in rowsSelected:
                ID = self.seqWin.IDs[r]
                model = self.seqWin.getModelForChain(r)
                if (not(model in self.models)):
                    self.models.append(str(model))
            # Find out where the poses are for each model
            self.poseindx = []
            for model in self.models:
                self.poseindx.append(self.seqWin.getPoseIndexForModel(model))
            # Update the combo box to have these models as selections
            self.modelMenu.Clear()
            self.modelMenu.AppendItems(self.models)
            if (len(self.models) > 0):
                if (platform.system() == "Windows"):
                    # Doing this on Linux freezes up the ComboBox, so don't do it on Linux
                    self.modelMenu.SetSelection(0)
                self.baseModel = self.models[0]
            else:
                self.baseModel = ""
        self.Scroll(0, self.winscrollpos)

    def modelClick(self, event):
        self.alignBy = "Model"
        if (platform.system() == "Darwin"):
            self.btnByModel.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/superimposition/btnByModel_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            self.btnByResidues.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/superimposition/btnByResidues.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
        else:
            self.btnByModel.SetForegroundColour("#FF0000")
            self.btnByResidues.SetForegroundColour("#000000")
        logInfo("Superimposition type set to By Model")

    def residuesClick(self, event):
        self.alignBy = "Residues"
        if (platform.system() == "Darwin"):
            self.btnByModel.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/superimposition/btnByModel.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            self.btnByResidues.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/superimposition/btnByResidues_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
        else:
            self.btnByModel.SetForegroundColour("#000000")
            self.btnByResidues.SetForegroundColour("#FF0000")
        logInfo("Superimposition type set to By Residues")

    def CAclick(self, event):
        self.RMSDtype = "CA"
        if (platform.system() == "Darwin"):
            self.btnCARMSD.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/superimposition/btnCARMSD_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            self.btnBBRMSD.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/superimposition/btnBBRMSD.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
        else:
            self.btnCARMSD.SetForegroundColour("#FF0000")
            self.btnBBRMSD.SetForegroundColour("#000000")
        logInfo("Superimposition type set to CA RMSD")

    def BBclick(self, event):
        self.RMSDtype = "BB"
        if (platform.system() == "Darwin"):
            self.btnCARMSD.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/superimposition/btnCARMSD.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            self.btnBBRMSD.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/superimposition/btnBBRMSD_Hi.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
        else:
            self.btnCARMSD.SetForegroundColour("#000000")
            self.btnBBRMSD.SetForegroundColour("#FF0000")
        logInfo("Superimposition type set to BB RMSD")

    def superimposeClick(self, event):
        # GetStringSelection does not appear to work right on Linux, so use GetValue instead
        logInfo("Superimposition button clicked")
        self.baseModel = str(self.modelMenu.GetValue())
        if (self.baseModel != "" and len(self.models) > 1):
            self.seqWin.labelMsg.SetLabel("Performing superimposition, please be patient...")
            self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
            self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")
            self.seqWin.msgQueue.append("Performing superimposition, please be patient...")
            self.btnSuperimpose.Disable()
            self.parent.GoBtn.Disable()
            self.seqWin.cannotDelete = True
            self.seqWin.setAlignType("Align")
            #self.seqWin.alignType = "Align"
            if (platform.system() == "Darwin"):
                self.seqWin.AlignBtn.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/sequence/AlignBtn_Align.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.seqWin.AlignBtn.SetLabel("Align")
            self.seqWin.regenerateLookupTable() # To get a MUSCLE alignment in the sequence window
            logInfo("Superimposition started using " + self.baseModel + " as the base model")
            thrSuperimpose = Thread(target=self.doSuperimposition, args=())
            thrSuperimpose.start()
            #self.doSuperimposition()
        elif (self.baseModel == ""):
            dlg = wx.MessageDialog(self, "Please select a base model for the superimposition.", "Superimposition Error", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
            if (dlg.ShowModal() == wx.ID_OK):
                pass
            dlg.Destroy()
        elif (len(self.models) <= 1):
            dlg = wx.MessageDialog(self, "Please select at least two models to superimpose in the sequence window.", "Superimposition Error", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
            if (dlg.ShowModal() == wx.ID_OK):
                pass
            dlg.Destroy()

    def doSuperimposition(self):
        # GetStringSelection does not appear to work right on Linux, so use GetValue instead
        self.baseModel = str(self.modelMenu.GetValue())
        homedir = os.path.expanduser("~")
        if (platform.system() == "Windows"):
            os.chdir(homedir + "\\InteractiveROSETTA")
        else:
            os.chdir(homedir + "/.InteractiveROSETTA")
        if (self.grdRMSDResults.NumberRows > 0):
            self.grdRMSDResults.DeleteRows(0, self.grdRMSDResults.NumberRows)
        if (self.RMSDtype == "CA"):
            self.grdRMSDResults.SetColLabelValue(0, "CA RMSD")
        else:
            self.grdRMSDResults.SetColLabelValue(0, "BB RMSD")
        # For some reason, this is necessary to get the label to update
        w = self.grdRMSDResults.GetColSize(0)
        self.grdRMSDResults.SetColSize(0, w-1)
        self.grdRMSDResults.SetColSize(0, w)
        for m in range(0, len(self.models)):
            model = self.models[m]
            if (model == self.baseModel):
                continue
            # Use PyMOL to do the alignments (very easy to do, and just grab the RMS for reporting afterwards)
            if (self.alignBy == "Model"):
                movingstr = model
                targetstr = self.baseModel
            else:
                movingstr = "seqsele and " + model
                targetstr = "seqsele and " + self.baseModel
            if (self.RMSDtype == "CA"):
                movingstr = movingstr + " and name ca"
                targetstr = targetstr + " and name ca"
            else:
                movingstr = movingstr + " and backbone"
                targetstr = targetstr + " and backbone"
            try:
                RMSD = self.cmd.align(movingstr, targetstr, cycles=0)[0]
            except:
                dlg = wx.MessageDialog(self, "The superimposition for " + model + " could not be performed.  You probably did not select enough residues for this model.", "Superimposition Failed", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
                dlg.ShowModal()
                dlg.Destroy()
                continue
            RMSD = int(RMSD * 1000) / 1000.0
            self.grdRMSDResults.AppendRows(1)
            self.grdRMSDResults.SetRowLabelValue(self.grdRMSDResults.NumberRows-1, self.baseModel + "-" + model)
            self.grdRMSDResults.SetCellValue(self.grdRMSDResults.NumberRows-1, 0, str(RMSD))
            readOnly = wx.grid.GridCellAttr()
            readOnly.SetReadOnly(True)
            self.grdRMSDResults.SetRowAttr(self.grdRMSDResults.NumberRows-1, readOnly)
            self.grdRMSDResults.SetCellAlignment(self.grdRMSDResults.NumberRows-1, 0, wx.ALIGN_CENTRE, wx.ALIGN_CENTRE)
            self.cmd.save(model + ".pdb", model, -1, "pdb")
            fixPyMOLSave(model + ".pdb")
            # Unload old pose and reload new one
            try:
                self.seqWin.poses[self.poseindx[m]] = self.seqWin.pdbreader.get_structure(model, model + ".pdb")
                self.cmd.remove(model)
                self.cmd.delete(model)
                self.cmd.load(model + ".pdb", model)
                defaultPyMOLView(self.cmd, model)
            except:
                # Weird error, pass instead of crash
                pass
            self.seqWin.cannotDelelte = False
            self.cmd.center("all")
        # Pop this message out of the queue
        for i in range(0, len(self.seqWin.msgQueue)):
            if (self.seqWin.msgQueue[i].find("Performing superimposition") >= 0):
                self.seqWin.msgQueue.pop(i)
                break
        if (len(self.seqWin.msgQueue) > 0):
            self.seqWin.labelMsg.SetLabel(self.seqWin.msgQueue[len(self.seqWin.msgQueue)-1])
        else:
            self.seqWin.labelMsg.SetLabel("")
        self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")
        self.btnSuperimpose.Enable()
        self.parent.GoBtn.Enable()
        self.grdRMSDResults.Refresh()
        self.seqWin.cannotDelete = False
        # This will recover the selection
        self.seqWin.selectUpdate(True)
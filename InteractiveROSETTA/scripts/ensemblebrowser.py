import wx
import wx.grid
import wx.lib.scrolledpanel
import os
import os.path
import platform
import webbrowser
from tools import resizeTextControlForUNIX
from tools import logInfo

class EnsembleBrowserPanel(wx.lib.scrolledpanel.ScrolledPanel):
    def __init__(self, parent, W, H):
        if (platform.system() == "Windows"):
            wx.lib.scrolledpanel.ScrolledPanel.__init__(self, parent, id=-1, pos=(10, 60), size=(340, H-330), name="ProtEnsembleBrowse")
            winh = H-330
        else:
            wx.lib.scrolledpanel.ScrolledPanel.__init__(self, parent, id=-1, pos=(10, 60), size=(340, H-290), name="ProtEnsembleBrowse")
            winh = H-290
        self.SetBackgroundColour("#333333")
        self.parent = parent
        
        if (platform.system() == "Windows"):
            self.lblProt = wx.StaticText(self, -1, "Ensemble Browser", (25, 15), (270, 25), style=wx.ALIGN_CENTRE)
            self.lblProt.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblProt = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/ensemblebrowser/lblBrowser.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(25, 15), size=(270, 25))
        else:
            self.lblProt = wx.StaticText(self, -1, "Ensemble Browser", pos=(90, 15), style=wx.ALIGN_CENTRE)
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
            self.lblInst = wx.StaticText(self, -1, "1. Load the desired structures into the\nSequence Window.\n\n2. Superimpose structures if necessary\n\n3. Select regions of the ensembles that\nwill be viewed.", (0, 45), (320, 25), wx.ALIGN_CENTRE)
            self.lblInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
        elif (platform.system() == "Darwin"):
            self.lblInst = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/ensemblebrowser/lblInstBrowser.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 45), size=(320, 120))
        else:
            self.lblInst = wx.StaticText(self, -1, "1. Load the desired structures into the\nSequence Window.\n\n2. Superimpose structures if necessary\n\n3. Select regions of the ensembles that\nwill be viewed.", (40, 45), (320, 25), style=wx.ALIGN_CENTRE)
            self.lblInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
            resizeTextControlForUNIX(self.lblInst, 0, self.GetSize()[0])
        self.lblInst.SetForegroundColour("#FFFFFF")
        
        if (platform.system() == "Windows"):
            self.lblAlways = wx.StaticText(self, -1, "Always Show:", (5, 175), (150, 30), wx.ALIGN_CENTRE)
            self.lblAlways.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblAlways = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/ensemblebrowser/lblAlways.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(25, 175), size=(150, 30))
        else:
            self.lblAlways = wx.StaticText(self, -1, "Always Show:", (25, 175), (150, 30), wx.ALIGN_CENTRE)
            self.lblAlways.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblAlways, 5, 150)
        self.lblAlways.SetForegroundColour("#FFFFFF")
        if (platform.system() == "Darwin"):
            self.modelMenu = wx.ComboBox(self, pos=(160, 175), size=(150, 25), choices=[], style=wx.CB_READONLY)
        else:
            self.modelMenu = wx.ComboBox(self, pos=(160, 170), size=(150, 25), choices=[], style=wx.CB_READONLY)
        self.modelMenu.SetToolTipString("Model to always display with semi-transparency")
        self.selectedModels = []
        
        if (platform.system() == "Darwin"):
            self.btnLeft = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/ensemblebrowser/btnLeft.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(5, 220), size=(25, 25))
        else:
            self.btnLeft = wx.Button(self, id=-1, label="<", pos=(5, 220), size=(25, 25))
            self.btnLeft.SetForegroundColour("#000000")
            self.btnLeft.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnLeft.Bind(wx.EVT_BUTTON, self.leftClick)
        self.btnLeft.SetToolTipString("Cycle to previous model")
        if (platform.system() == "Darwin"):
            self.btnRight = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/ensemblebrowser/btnRight.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(290, 220), size=(25, 25))
        else:
            self.btnRight = wx.Button(self, id=-1, label=">", pos=(290, 220), size=(25, 25))
            self.btnRight.SetForegroundColour("#000000")
            self.btnRight.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnRight.Bind(wx.EVT_BUTTON, self.rightClick)
        self.btnRight.SetToolTipString("Cycle to next model")
        if (platform.system() != "Linux"):
            self.lblSelected = wx.StaticText(self, -1, "Not Selected", (35, 225), (250, 30), wx.ALIGN_CENTRE)
        else:
            self.lblSelected = wx.StaticText(self, -1, "Not Selected", (35, 225), (250, 30), wx.ALIGN_CENTRE)
            resizeTextControlForUNIX(self.lblSelected, 35, 250)
        self.lblSelected.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.lblSelected.SetForegroundColour("#FFFFFF")
        self.currentModel = -1
        
        scrollh = self.btnLeft.GetPosition()[1] + self.btnLeft.GetSize()[1] + 5
        self.SetScrollbars(1, 1, 320, scrollh)
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
            browser.open(self.parent.parent.scriptdir + "/help/browser.html")
        else:
            webbrowser.open(self.parent.parent.scriptdir + "/help/browser.html")

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
        # Get the selected models from the sequence window
        newmodels = self.seqWin.getSelectedModels()
        if (self.selectedModels != newmodels):
            self.selectedModels = newmodels
            # Add the selections to the "Always Show" menu
            self.modelMenu.Clear()
            self.modelMenu.Append("None")
            for model in self.selectedModels:
                self.modelMenu.Append(model)
            if (len(self.selectedModels) == 0):
                self.lblSelected.SetLabel("Not Selected")
            else:
                if (self.currentModel >= len(self.selectedModels) or self.currentModel < 0):
                    self.currentModel = 0
                self.lblSelected.SetLabel(self.selectedModels[self.currentModel])
            if (platform.system() == "Linux"):
                resizeTextControlForUNIX(self.lblSelected, 35, 250)
        self.Scroll(0, self.winscrollpos)

    def updateView(self):
        # Change the PyMOL view to only show the current model and optionally the base model that is
        # always shown, if one is indicated by the modelMenu
        # This is accomplished by making everything outside of the view transparent
        # Why did I use transparency? -> because "hiding" views causes PyMOL to lose the display 
        # information (i.e. it will forget if you had lines/sticks/spheres shown)
        # You don't lose this information by fiddling with transparencies
        if (len(self.selectedModels) == 0):
            return
        model = self.selectedModels[self.currentModel]
        self.lblSelected.SetLabel(model)
        if (platform.system() == "Linux"):
            resizeTextControlForUNIX(self.lblSelected, 35, 250)
        # The selection "seqsele" refers to what the user has indicated should be in the ensemble browsing
        # If it's not in "seqsele", then it will never be shown by the browser
        # First make everything transparent and turn off cartoons (the user will have to deal with the
        # cartoon settings changing)
        self.cmd.hide("cartoon", "all")
        self.cmd.set_bond("stick_transparency", 1, "all")
        self.cmd.set("sphere_transparency", 1, "all") # Use "set" for sphere_transparency
        # Now define a selection for the intersection of the current model atoms and "sele"
        self.cmd.select("ensemblesele", "seqsele and model " + model)
        self.cmd.show("cartoon", "ensemblesele")
        self.cmd.set_bond("stick_transparency", 0, "ensemblesele")
        self.cmd.set("sphere_transparency", 0, "ensemblesele") # Use "set" for sphere_transparency
        # Now show the "always show" model if one if selected and it's not the current model
        alwaysshow = self.modelMenu.GetStringSelection()
        logInfo("Browser view changed to show " + model + " with Always Show set to " + alwaysshow)
        if (alwaysshow != "None" and alwaysshow != "" and alwaysshow != model):
            self.cmd.select("ensemblesele", "seqsele and model " + alwaysshow)
            self.cmd.show("cartoon", "ensemblesele")
            self.cmd.set_bond("stick_transparency", 0.7, "ensemblesele") # So you know which one is the "always shown" one
            self.cmd.set("sphere_transparency", 0, "ensemblesele") # Use "set" for sphere_transparency
        self.cmd.delete("ensemblesele") # Clean up
        self.cmd.enable("seqsele") # This is needed otherwise the sequence window removes the selection since it is not enabled

    def leftClick(self, event):
        self.currentModel = self.currentModel - 1
        if (self.currentModel < 0):
            self.currentModel = 0
        self.updateView()
        logInfo("Left button clicked")

    def rightClick(self, event):
        self.currentModel = self.currentModel + 1
        if (self.currentModel >= len(self.selectedModels)):
            self.currentModel = len(self.selectedModels) - 1
        self.updateView()
        logInfo("Right button clicked")
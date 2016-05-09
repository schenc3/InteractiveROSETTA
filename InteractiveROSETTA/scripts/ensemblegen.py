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
import gzip
from threading import Thread
from tools import *

class EnsembleGenPanel(wx.lib.scrolledpanel.ScrolledPanel):
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
            self.lblProt = wx.StaticText(self, -1, "Ensemble Generation", (25, 15), (270, 25), wx.ALIGN_CENTRE)
            self.lblProt.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblProt = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/ensemblegen/lblEnsembleGeneration.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 15), size=(320, 25))
        else:
            self.lblProt = wx.StaticText(self, -1, "Ensemble Generation", (70, 15), style=wx.ALIGN_CENTRE)
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
            self.lblInst = wx.StaticText(self, -1, "Create structural ensembles for a model.\n\nRelax - Generate structural diversity\nusing minor perturbations.\n\nBackrub - Generate extensive structural diversity\nup to a given RMSD", (0, 45), (320, 25), wx.ALIGN_CENTRE)
            self.lblInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
        elif (platform.system() == "Darwin"):
            self.lblInst = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/ensemblegen/lblInstEnsembleGeneration.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 45), size=(320, 90))
        else:
            self.lblInst = wx.StaticText(self, -1, "Create structural ensembles for a model.\n\nRelax - Generate structural diversity\nusing minor perturbations.\n\nBackrub - Generate extensive structural diversity\nup to a given RMSD", (5, 45), style=wx.ALIGN_CENTRE)
            self.lblInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
            resizeTextControlForUNIX(self.lblInst, 0, self.GetSize()[0])
        self.lblInst.SetForegroundColour("#FFFFFF")
        
        if (platform.system() == "Windows"):
            self.lblLine = wx.StaticText(self, -1, "==========================", (0, 160), (320, 20), wx.ALIGN_CENTRE)
            self.lblLine.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.NORMAL))
        elif (platform.system() == "Darwin"):
            self.lblLine = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/ensemblegen/lblLine.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 160), size=(320, 20))
        else:
            self.lblLine = wx.StaticText(self, -1, "==========================", (0, 160), style=wx.ALIGN_CENTRE)
            self.lblLine.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.NORMAL))
            resizeTextControlForUNIX(self.lblLine, 20, 120)
        self.lblLine.SetForegroundColour("#FFFFFF")
        
        if (platform.system() == "Windows"):
            self.lblRelax = wx.StaticText(self, -1, "Relax Ensemble", (0, 180), (320, 25), wx.ALIGN_CENTRE)
            self.lblRelax.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblRelax = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/ensemblegen/lblRelax.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 180), size=(320, 25))
        else:
            self.lblRelax = wx.StaticText(self, -1, "Relax Ensemble", (0, 180), style=wx.ALIGN_CENTRE)
            self.lblRelax.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
            resizeTextControlForUNIX(self.lblRelax, 0, self.GetSize()[0])
        self.lblRelax.SetForegroundColour("#FFFFFF")
        
        if (platform.system() == "Windows"):
            self.lblModelR = wx.StaticText(self, -1, "Base Model", (0, 210), (230, 20), wx.ALIGN_CENTRE)
            self.lblModelR.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblModelR = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/ensemblegen/lblModelR.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 210), size=(230, 20))
        else:
            self.lblModelR = wx.StaticText(self, -1, "Base Model", (0, 210), style=wx.ALIGN_CENTRE)
            self.lblModelR.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblModelR, 0, 230)
        self.lblModelR.SetForegroundColour("#FFFFFF")
        self.menuModelR = wx.ComboBox(self, pos=(0, 230), size=(230, 25), choices=[], style=wx.CB_READONLY)
        self.menuModelR.SetToolTipString("Select the model that will be diversified")
        
        if (platform.system() == "Windows"):
            self.lblNModelsR = wx.StaticText(self, -1, "# Models", (240, 210), (80, 20), wx.ALIGN_CENTRE)
            self.lblNModelsR.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblNModelsR = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/ensemblegen/lblNModelsR.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(240, 210), size=(80, 20))
        else:
            self.lblNModelsR = wx.StaticText(self, -1, "# Models", (240, 210), style=wx.ALIGN_CENTRE)
            self.lblNModelsR.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblNModelsR, 240, 80)
        self.lblNModelsR.SetForegroundColour("#FFFFFF")
        self.txtNModelsR = wx.TextCtrl(self, -1, pos=(240, 230), size=(80, 25))
        self.txtNModelsR.SetValue("10")
        self.txtNModelsR.SetToolTipString("Number of models to generate")
        
        if (platform.system() == "Darwin"):
            self.btnRelax = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/ensemblegen/btnRelax.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(110, 270), size=(100, 25))
        else:
            self.btnRelax = wx.Button(self, id=-1, label="Relax!", pos=(110, 270), size=(100, 25))
            self.btnRelax.SetForegroundColour("#000000")
            self.btnRelax.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.btnRelax.Bind(wx.EVT_BUTTON, self.relaxClick)
        self.btnRelax.SetToolTipString("Begin relax protocol to generate a structural ensemble")
        self.buttonState = "Ensemble!"
        
        if (platform.system() == "Windows"):
            self.lblLine2 = wx.StaticText(self, -1, "==========================", (0, 300), (320, 20), wx.ALIGN_CENTRE)
            self.lblLine2.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.NORMAL))
        elif (platform.system() == "Darwin"):
            self.lblLine2 = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/ensemblegen/lblLine.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 300), size=(320, 20))
        else:
            self.lblLine2 = wx.StaticText(self, -1, "==========================", (0, 300), style=wx.ALIGN_CENTRE)
            self.lblLine2.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.NORMAL))
            resizeTextControlForUNIX(self.lblLine2, 20, 120)
        self.lblLine2.SetForegroundColour("#FFFFFF")
        
        if (platform.system() == "Windows"):
            self.lblBackrub = wx.StaticText(self, -1, "Backrub Ensemble", (0, 320), (320, 25), wx.ALIGN_CENTRE)
            self.lblBackrub.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblBackrub = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/ensemblegen/lblBackrub.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 320), size=(320, 25))
        else:
            self.lblBackrub = wx.StaticText(self, -1, "Backrub Ensemble", (0, 320), style=wx.ALIGN_CENTRE)
            self.lblBackrub.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
            resizeTextControlForUNIX(self.lblBackrub, 0, self.GetSize()[0])
        self.lblBackrub.SetForegroundColour("#FFFFFF")
        
        if (platform.system() == "Windows"):
            self.lblModelB = wx.StaticText(self, -1, "Base Model", (0, 350), (160, 20), wx.ALIGN_CENTRE)
            self.lblModelB.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblModelB = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/ensemblegen/lblModelR.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 350), size=(230, 20))
        else:
            self.lblModelB = wx.StaticText(self, -1, "Base Model", (0, 350), style=wx.ALIGN_CENTRE)
            self.lblModelB.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblModelB, 0, 160)
        self.lblModelB.SetForegroundColour("#FFFFFF")
        self.menuModelB = wx.ComboBox(self, pos=(0, 370), size=(160, 25), choices=[], style=wx.CB_READONLY)
        self.menuModelB.SetToolTipString("Select the model that will be diversified")
        
        if (platform.system() == "Windows"):
            self.lblNModelsB = wx.StaticText(self, -1, "# Models", (170, 350), (70, 20), wx.ALIGN_CENTRE)
            self.lblNModelsB.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblNModelsB = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/ensemblegen/lblNModelsR.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(170, 350), size=(70, 20))
        else:
            self.lblNModelsB = wx.StaticText(self, -1, "# Models", (170, 350), style=wx.ALIGN_CENTRE)
            self.lblNModelsB.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblNModelsB, 170, 70)
        self.lblNModelsB.SetForegroundColour("#FFFFFF")
        self.txtNModelsB = wx.TextCtrl(self, -1, pos=(170, 370), size=(70, 25))
        self.txtNModelsB.SetValue("10")
        self.txtNModelsB.SetToolTipString("Number of models to generate")
        
        if (platform.system() == "Windows"):
            self.lblMaxRMSD = wx.StaticText(self, -1, "maxRMSD", (250, 350), (70, 20), wx.ALIGN_CENTRE)
            self.lblMaxRMSD.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblMaxRMSD = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/ensemblegen/lblMaxRMSD.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(250, 350), size=(70, 20))
        else:
            self.lblMaxRMSD = wx.StaticText(self, -1, "maxRMSD", (250, 350), style=wx.ALIGN_CENTRE)
            self.lblMaxRMSD.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblMaxRMSD, 250, 70)
        self.lblMaxRMSD.SetForegroundColour("#FFFFFF")
        self.txtMaxRMSD = wx.TextCtrl(self, -1, pos=(250, 370), size=(70, 25))
        self.txtMaxRMSD.SetValue("1")
        self.txtMaxRMSD.SetToolTipString("Number of models to generate")
        
        if (platform.system() == "Darwin"):
            self.btnServerToggle = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/ensemblegen/btnServer_Off.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(40, 410), size=(100, 25))
        else:
            self.btnServerToggle = wx.Button(self, id=-1, label="Server Off", pos=(40, 410), size=(100, 25))
            self.btnServerToggle.SetForegroundColour("#000000")
            self.btnServerToggle.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.btnServerToggle.Bind(wx.EVT_BUTTON, self.serverToggle)
        self.btnServerToggle.SetToolTipString("Perform backrub simulations locally")
        self.serverOn = False
        
        if (platform.system() == "Darwin"):
            self.btnBackrub = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/ensemblegen/btnBackrub.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(180, 410), size=(100, 25))
        else:
            self.btnBackrub = wx.Button(self, id=-1, label="Backrub!", pos=(180, 410), size=(100, 25))
            self.btnBackrub.SetForegroundColour("#000000")
            self.btnBackrub.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.btnBackrub.Bind(wx.EVT_BUTTON, self.backrubClick)
        self.btnBackrub.SetToolTipString("Begin backrub protocol to generate a structural ensemble")
        
        if (platform.system() == "Windows"):
            self.lblLine3 = wx.StaticText(self, -1, "==========================", (0, 440), (320, 20), wx.ALIGN_CENTRE)
            self.lblLine3.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.NORMAL))
        elif (platform.system() == "Darwin"):
            self.lblLine3 = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/ensemblegen/lblLine.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 440), size=(320, 20))
        else:
            self.lblLine3 = wx.StaticText(self, -1, "==========================", (0, 440), style=wx.ALIGN_CENTRE)
            self.lblLine3.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.NORMAL))
            resizeTextControlForUNIX(self.lblLine3, 20, 120)
        self.lblLine3.SetForegroundColour("#FFFFFF")
        
        if (platform.system() == "Windows"):
            self.lblUnpack = wx.StaticText(self, -1, "Unpack Ensemble", (0, 460), (320, 25), wx.ALIGN_CENTRE)
            self.lblUnpack.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblUnpack = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/ensemblegen/lblUnpack.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 460), size=(460, 25))
        else:
            self.lblUnpack = wx.StaticText(self, -1, "Unpack Ensemble", (0, 460), style=wx.ALIGN_CENTRE)
            self.lblUnpack.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
            resizeTextControlForUNIX(self.lblUnpack, 0, self.GetSize()[0])
        self.lblUnpack.SetForegroundColour("#FFFFFF")
        
        if (platform.system() == "Windows"):
            self.lblPrefix = wx.StaticText(self, -1, "Filename Prefix", (0, 490), (230, 20), wx.ALIGN_CENTRE)
            self.lblPrefix.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblPrefix = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/ensemblegen/lblPrefix.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 490), size=(230, 20))
        else:
            self.lblPrefix = wx.StaticText(self, -1, "Filename Prefix", (0, 490), style=wx.ALIGN_CENTRE)
            self.lblPrefix.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblPrefix, 0, 230)
        self.lblPrefix.SetForegroundColour("#FFFFFF")
        self.txtPrefix = wx.TextCtrl(self, -1, pos=(00, 510), size=(210, 25))
        self.txtPrefix.SetValue("")
        self.txtPrefix.SetToolTipString("Filename prefix that will be used to generate the individual PDB files")
        
        if (platform.system() == "Darwin"):
            self.btnUnpack = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/ensemblegen/btnUnpack.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(240, 510), size=(100, 25))
        else:
            self.btnUnpack = wx.Button(self, id=-1, label="Choose File", pos=(220, 510), size=(100, 25))
            self.btnUnpack.SetForegroundColour("#000000")
            self.btnUnpack.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.btnUnpack.Bind(wx.EVT_BUTTON, self.unpackClick)
        self.btnUnpack.SetToolTipString("Choose an ensemble archive to unpack")
        
        self.scrollh = self.btnUnpack.GetPosition()[1] + self.btnUnpack.GetSize()[1] + 5
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
        
    def scrolled(self, event):
        self.winscrollpos = self.GetScrollPos(wx.VERTICAL)
        event.Skip()
        
    def activate(self):
        # Update the menus with all the loaded models
        models = []
        for r in range(0, len(self.seqWin.IDs)):
            if (self.seqWin.poses[r]):
                models.append(self.seqWin.IDs[r][0:len(self.seqWin.IDs[r])-2])
        self.menuModelR.Clear()
        self.menuModelR.AppendItems(models)
        self.menuModelB.Clear()
        self.menuModelB.AppendItems(models)
        self.Scroll(0, self.winscrollpos)
    
    def cancelRelax(self):
        logInfo("Canceled structure relaxation operation")
        try:
            os.remove("relaxinput")
        except:
            pass
        try:
            os.remove("relaxinputtemp")
        except:
            pass
        self.tmrRelax.Stop()
        self.seqWin.cannotDelete = False
        self.btnBackrub.Enable()
        if (platform.system() == "Darwin"):
            self.btnRelax.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/ensemblegen/btnRelax.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
        else:
            self.btnRelax.SetLabel("Relax!")
        self.btnRelax.SetToolTipString("Begin relax protocol to generate a structural ensemble")
        self.parent.parent.restartDaemon()
        self.parent.GoBtn.Enable()
        # Get rid of the messages
        for i in range(0, len(self.seqWin.msgQueue)):
            if (self.seqWin.msgQueue[i].find("Performing structure relaxation") >= 0):
                self.seqWin.msgQueue.pop(i)
                break
        if (len(self.seqWin.msgQueue) > 0):
            self.seqWin.labelMsg.SetLabel(self.seqWin.msgQueue[len(self.seqWin.msgQueue)-1])
        else:
            self.seqWin.labelMsg.SetLabel("")
        self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")
        self.buttonState = "Ensemble!"
    
    def relaxClick(self, event):
        # This is also the "Finalize!" button
        logInfo("Relax button clicked")
        if (self.buttonState == "Ensemble!"):
            # Check to make sure the inputs are valid
            modelOkay = False
            if (len(self.menuModelR.GetValue().strip()) == 0):
                modelOkay = False
            else:
                for ID in self.seqWin.IDs:
                    if (self.menuModelR.GetValue().strip() in ID):
                        modelOkay = True
                        break
            if (not(modelOkay)):
                dlg = wx.MessageDialog(self, "Please select a valid protein model.", "Invalid Model", wx.OK | wx.ICON_ERROR | wx.CENTRE)
                dlg.ShowModal()
                dlg.Destroy()
                return
            try:
                if (int(self.txtNModelsR.GetValue()) < 2):
                    dlg = wx.MessageDialog(self, "Please generate at least 2 models for the ensemble.", "Invalid Ensemble Size", wx.OK | wx.ICON_ERROR | wx.CENTRE)
                    dlg.ShowModal()
                    dlg.Destroy()
                    return
            except:
                dlg = wx.MessageDialog(self, "Please provide an integer greater than 1 for the number of models.", "Invalid Ensemble Size", wx.OK | wx.ICON_ERROR | wx.CENTRE)
                dlg.ShowModal()
                dlg.Destroy()
                return
            self.seqWin.labelMsg.SetLabel("Performing structure relaxation, please be patient...")
            self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
            self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")
            self.seqWin.msgQueue.append("Performing structure relaxation, please be patient...")
            self.parent.GoBtn.Disable()
            self.btnBackrub.Disable()
            self.seqWin.cannotDelete = True
            self.stage = 1
            if (platform.system() == "Darwin"):
                self.btnRelax.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/ensemblegen/btnRelax_Cancel.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnRelax.SetLabel("Cancel!")
            self.buttonState = "Cancel!"
            self.btnRelax.SetToolTipString("Cancel the structure relaxation")
            self.tmrRelax = wx.Timer(self)
            self.Bind(wx.EVT_TIMER, self.threadRelax, self.tmrRelax)
            self.tmrRelax.Start(1000)
        elif (self.buttonState == "Cancel!"):
            dlg = wx.MessageDialog(self, "Are you sure you want to cancel the structure relaxation?  All progress will be lost.", "Cancel Point Mutant Scan", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
            result = dlg.ShowModal()
            if (result == wx.ID_YES):
                self.cancelRelax()
            dlg.Destroy()
    
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
        self.btnRelax.Enable()
        self.btnBackrub.Enable()
        if (platform.system() == "Darwin"):
            self.btnRelax.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/ensemblegen/btnRelax.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
        else:
            self.btnRelax.SetLabel("Scan!")
        if (platform.system() == "Darwin"):
            self.btnBackrub.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/ensemblegen/btnBackrub.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
        else:
            self.btnBackrub.SetLabel("Scan!")
        self.buttonState = "Ensemble!"
        # Get rid of the messages
        for i in range(0, len(self.seqWin.msgQueue)):
            if (self.seqWin.msgQueue[i].find("Performing structure relaxation") >= 0):
                self.seqWin.msgQueue.pop(i)
                break
        for i in range(0, len(self.seqWin.msgQueue)):
            if (self.seqWin.msgQueue[i].find("Performing structure backrub") >= 0):
                self.seqWin.msgQueue.pop(i)
                break
        if (len(self.seqWin.msgQueue) > 0):
            self.seqWin.labelMsg.SetLabel(self.seqWin.msgQueue[len(self.seqWin.msgQueue)-1])
        else:
            self.seqWin.labelMsg.SetLabel("")
        self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")
    
    def threadRelax(self, event):
        # Why am I using a Timer?  See the explanation in kic.py
        goToSandbox()
        if (self.stage == 1):
            self.tmrRelax.Stop()
            self.timeoutCount = 0
            f = open("relaxinputtemp", "w")
            pdbfile = self.menuModelR.GetValue() + ".pdb"
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
            f.write("NMODELS\t" + str(self.txtNModelsR.GetValue()) + "\n")
            f.close()
            appendScorefxnParamsInfoToFile("relaxinputtemp", self.selectWin.weightsfile)
            os.rename("relaxinputtemp", "relaxinput")
            logInfo("Relaxation input uploaded locally at relaxinput")
            self.progress = wx.ProgressDialog("Structure Relaxation", "Performing a structure relaxation...", int(self.txtNModelsR.GetValue()), style=wx.PD_CAN_ABORT | wx.PD_APP_MODAL | wx.PD_ELAPSED_TIME | wx.PD_REMAINING_TIME)
            self.stage = 2
            self.tmrRelax.Start(1000)
        else:
            # Read the output dumped by the child process
            if (os.path.isfile("results.ensb")):
                try:
                    self.progress.Destroy()
                except:
                    pass
                self.tmrRelax.Stop()
                # Save and unpack the resulting ensemble archive
                while (True):
                    dlg = wx.FileDialog(
                            self, message="Save the relaxed ensemble",
                            defaultDir=self.seqWin.cwd,
                            defaultFile="",
                            wildcard="Ensemble Archives (*.ensb)|*.ensb",
                            style=wx.SAVE | wx.CHANGE_DIR)
                    if (dlg.ShowModal() == wx.ID_OK):
                        if (platform.system() == "Darwin"):
                            paths = [dlg.GetPath()]
                        else:
                            paths = dlg.GetPaths()
                        # Change cwd to the last opened file
                        if (platform.system() == "Windows"):
                            lastDirIndx = paths[len(paths)-1].rfind("\\")
                        else:
                            lastDirIndx = paths[len(paths)-1].rfind("/")
                        self.seqWin.cwd = str(paths[len(paths)-1][0:lastDirIndx])
                        self.seqWin.saveWindowData(None)
                        filename = str(paths[0]).split(".ensb")[0] + ".ensb"
                        if (os.path.isfile(filename)):
                            dlg2 = wx.MessageDialog(self, "The file " + filename + " already exists.  Overwrite it?", "Filename Already Exists", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
                            if (dlg2.ShowModal() == wx.ID_NO):
                                dlg2.Destroy()
                                logInfo("Cancelled save operation due to filename already existing")
                                continue
                            else:
                                os.remove(str(filename))
                            dlg2.Destroy()
                        prefix = filename.split(".ensb")[0]
                        goToSandbox()
                        os.rename("results.ensb", filename)
                        self.unpackEnsemble(filename, prefix)
                        break
                    else:
                        dlg2 = wx.MessageDialog(self, "Are you sure you want to cancel?  Your results will be lost!", "Confirm Cancellation", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
                        if (dlg2.ShowModal() == wx.ID_YES):
                            dlg2.Destroy()
                            logInfo("Refused to save relaxation output")
                            break
                        dlg2.Destroy()
                logInfo("Found structure relaxation output at relaxoutput")
                # Pop this message out of the queue
                for i in range(0, len(self.seqWin.msgQueue)):
                    if (self.seqWin.msgQueue[i].find("Performing structure relaxation") >= 0):
                        self.seqWin.msgQueue.pop(i)
                        break
                if (len(self.seqWin.msgQueue) > 0):
                    self.seqWin.labelMsg.SetLabel(self.seqWin.msgQueue[len(self.seqWin.msgQueue)-1])
                else:
                    self.seqWin.labelMsg.SetLabel("")
                self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
                self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")
                self.parent.GoBtn.Enable()
                self.btnBackrub.Enable()
                if (platform.system() == "Darwin"):
                    self.btnRelax.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/ensemblegen/btnRelax.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                else:
                    self.btnRelax.SetLabel("Relax!")
                self.buttonState = "Ensemble!"
                self.btnRelax.SetToolTipString("Begin relax protocol to generate a structural ensemble")
                self.seqWin.cannotDelete = False
            elif (os.path.isfile("errreport")):
                self.tmrRelax.Stop()
                self.recoverFromError()
            elif (os.path.isfile("relaxprogress")):
                f = open("relaxprogress", "r")
                data = f.readlines()
                f.close()
                if (len(data) == 0):
                    indx = 0
                try:
                    indx = int(data[0])
                except:
                    indx = 0
                (keepGoing, skip) = self.progress.Update(indx)
                if (not(keepGoing)):
                    # User clicked "Cancel" on the progress bar
                    self.cancelRelax()
                    self.progress.Destroy()
    
    def cancelBackrub(self):
        logInfo("Canceled structure backrub operation")
        try:
            os.remove("backrubinput")
        except:
            pass
        try:
            os.remove("backrubinputtemp")
        except:
            pass
        self.tmrBackrub.Stop()
        self.seqWin.cannotDelete = False
        self.btnRelax.Enable()
        if (platform.system() == "Darwin"):
            self.btnBackrub.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/ensemblegen/btnBackrub.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
        else:
            self.btnBackrub.SetLabel("Backrub!")
        self.btnRelax.SetToolTipString("Begin backrub protocol to generate a structural ensemble")
        self.parent.parent.restartDaemon()
        self.parent.GoBtn.Enable()
        # Get rid of the messages
        for i in range(0, len(self.seqWin.msgQueue)):
            if (self.seqWin.msgQueue[i].find("Performing structure backrub") >= 0):
                self.seqWin.msgQueue.pop(i)
                break
        if (len(self.seqWin.msgQueue) > 0):
            self.seqWin.labelMsg.SetLabel(self.seqWin.msgQueue[len(self.seqWin.msgQueue)-1])
        else:
            self.seqWin.labelMsg.SetLabel("")
        self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")
        self.buttonState = "Ensemble!"

    def serverToggle(self, event):
        if (self.serverOn):
            self.serverOn = False
            if (platform.system() == "Darwin"):
                self.btnServerToggle.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/ensemblegen/btnServer_Off.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnServerToggle.SetLabel("Server Off")
            self.btnServerToggle.SetToolTipString("Perform backrub simulations locally")
            logInfo("Turned off backrub server usage")
        else:
            self.serverOn = True
            if (platform.system() == "Darwin"):
                self.btnServerToggle.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/ensemblegen/btnServer_On.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnServerToggle.SetLabel("Server On")
            self.btnServerToggle.SetToolTipString("Perform backrub simulations on a remote server")
            logInfo("Turned on backrub server usage")

    def backrubClick(self, event):
        # This is also the "Finalize!" button
        logInfo("Backrub button clicked")
        if (self.buttonState == "Ensemble!"):
            # Check to make sure the inputs are valid
            modelOkay = False
            if (len(self.menuModelB.GetValue().strip()) == 0):
                modelOkay = False
            else:
                for ID in self.seqWin.IDs:
                    if (self.menuModelB.GetValue().strip() in ID):
                        modelOkay = True
                        break
            if (not(modelOkay)):
                dlg = wx.MessageDialog(self, "Please select a valid protein model.", "Invalid Model", wx.OK | wx.ICON_ERROR | wx.CENTRE)
                dlg.ShowModal()
                dlg.Destroy()
                return
            try:
                if (int(self.txtNModelsB.GetValue()) < 2):
                    dlg = wx.MessageDialog(self, "Please generate at least 2 models for the ensemble.", "Invalid Ensemble Size", wx.OK | wx.ICON_ERROR | wx.CENTRE)
                    dlg.ShowModal()
                    dlg.Destroy()
                    return
            except:
                dlg = wx.MessageDialog(self, "Please provide an integer greater than 1 for the number of models.", "Invalid Ensemble Size", wx.OK | wx.ICON_ERROR | wx.CENTRE)
                dlg.ShowModal()
                dlg.Destroy()
                return
            try:
                if (float(self.txtMaxRMSD.GetValue()) <= 0.0):
                    dlg = wx.MessageDialog(self, "Please provide a positive number for the RMSD limit.", "Invalid Max RMSD", wx.OK | wx.ICON_ERROR | wx.CENTRE)
                    dlg.ShowModal()
                    dlg.Destroy()
                    return
            except:
                dlg = wx.MessageDialog(self, "Please provide a positive number for the RMSD limit.", "Invalid Max RMSD", wx.OK | wx.ICON_ERROR | wx.CENTRE)
                dlg.ShowModal()
                dlg.Destroy()
                return
            self.seqWin.labelMsg.SetLabel("Performing structure backrub, please be patient...")
            self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
            self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")
            self.seqWin.msgQueue.append("Performing structure backrub, please be patient...")
            self.parent.GoBtn.Disable()
            self.btnRelax.Disable()
            self.seqWin.cannotDelete = True
            self.stage = 1
            if (platform.system() == "Darwin"):
                self.btnBackrub.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/ensemblegen/btnBackrub_Cancel.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnBackrub.SetLabel("Cancel!")
            self.buttonState = "Cancel!"
            self.btnBackrub.SetToolTipString("Cancel the structure relaxation")
            self.tmrBackrub = wx.Timer(self)
            self.Bind(wx.EVT_TIMER, self.threadBackrub, self.tmrBackrub)
            self.tmrBackrub.Start(1000)
        elif (self.buttonState == "Cancel!"):
            dlg = wx.MessageDialog(self, "Are you sure you want to cancel the structure backrub?  All progress will be lost.", "Cancel Point Mutant Scan", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
            result = dlg.ShowModal()
            if (result == wx.ID_YES):
                self.cancelBackrub()
            dlg.Destroy()
    
    def threadBackrub(self, event):
        # Why am I using a Timer?  See the explanation in kic.py
        goToSandbox()
        if (self.stage == 1):
            self.tmrBackrub.Stop()
            self.timeoutCount = 0
            f = open("backrubinputtemp", "w")
            pdbfile = self.menuModelB.GetValue() + ".pdb"
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
            f.write("NMODELS\t" + str(self.txtNModelsB.GetValue()) + "\n")
            f.write("MAXRMSD\t" + str(self.txtMaxRMSD.GetValue()) + "\n")
            f.close()
            appendScorefxnParamsInfoToFile("backrubinputtemp", self.selectWin.weightsfile)
            if (self.serverOn):
                try: 
                    self.ID = sendToServer("backrubinput")
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
                            if (len(aline.split("\t")) >= 2 and aline.split("\t")[0] == "BACKRUB" and aline.split("\t")[1] == self.ID.strip()):
                                alreadythere = True
                                break
                        f.close()
                    except:
                        pass
                    if (not(alreadythere)):
                        f = open("downloadwatch", "a")
                        f.write("BACKRUB\t" + self.ID.strip() + "\t" + str(datetime.datetime.now().strftime("%A, %B %d - %I:%M:%S %p")) + "\t" + getServerName() + "\t" + desc + "\n")
                        f.close()
                    dlg = wx.MessageDialog(self, "InteractiveROSETTA is now watching the server for job ID " + desc.strip() + ".  You will be notified when the package is available for download.", "Listening for Download", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
                    dlg.ShowModal()
                    dlg.Destroy()
                    # Pop this message out of the queue
                    for i in range(0, len(self.seqWin.msgQueue)):
                        if (self.seqWin.msgQueue[i].find("Performing structure backrub") >= 0):
                            self.seqWin.msgQueue.pop(i)
                            break
                    if (len(self.seqWin.msgQueue) > 0):
                        self.seqWin.labelMsg.SetLabel(self.seqWin.msgQueue[len(self.seqWin.msgQueue)-1])
                    else:
                        self.seqWin.labelMsg.SetLabel("")
                    self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
                    self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")
                    self.parent.GoBtn.Enable()
                    self.btnRelax.Enable()
                    if (platform.system() == "Darwin"):
                        self.btnBackrub.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/ensemblegen/btnBackrub.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                    else:
                        self.btnBackrub.SetLabel("Backrub!")
                    self.buttonState = "Ensemble!"
                    self.btnBackrub.SetToolTipString("Begin backrub protocol to generate a structural ensemble")
                    self.seqWin.cannotDelete = False
                except:
                    dlg = wx.MessageDialog(self, "The server could not be reached!  Ensure that you have specified a valid server and that you have an network connection.", "Server Could Not Be Reached", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
                    dlg.ShowModal()
                    dlg.Destroy()
                    return
            else:
                os.rename("backrubinputtemp", "backrubinput")
                logInfo("Backrub input uploaded locally at backrubinput")
                self.progress = wx.ProgressDialog("Structure Backrub", "Performing a structure backrub...", int(self.txtNModelsB.GetValue()), style=wx.PD_CAN_ABORT | wx.PD_APP_MODAL | wx.PD_ELAPSED_TIME | wx.PD_REMAINING_TIME)
                self.stage = 2
                self.tmrBackrub.Start(1000)
        else:
            # Read the output dumped by the child process
            if (os.path.isfile("results.ensb")):
                try:
                    self.progress.Destroy()
                except:
                    pass
                self.tmrBackrub.Stop()
                # Save and unpack the resulting ensemble archive
                while (True):
                    dlg = wx.FileDialog(
                            self, message="Save the backrubbed ensemble",
                            defaultDir=self.seqWin.cwd,
                            defaultFile="",
                            wildcard="Ensemble Archives (*.ensb)|*.ensb",
                            style=wx.SAVE | wx.CHANGE_DIR)
                    if (dlg.ShowModal() == wx.ID_OK):
                        if (platform.system() == "Darwin"):
                            paths = [dlg.GetPath()]
                        else:
                            paths = dlg.GetPaths()
                        # Change cwd to the last opened file
                        if (platform.system() == "Windows"):
                            lastDirIndx = paths[len(paths)-1].rfind("\\")
                        else:
                            lastDirIndx = paths[len(paths)-1].rfind("/")
                        self.seqWin.cwd = str(paths[len(paths)-1][0:lastDirIndx])
                        self.seqWin.saveWindowData(None)
                        filename = str(paths[0]).split(".ensb")[0] + ".ensb"
                        if (os.path.isfile(filename)):
                            dlg2 = wx.MessageDialog(self, "The file " + filename + " already exists.  Overwrite it?", "Filename Already Exists", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
                            if (dlg2.ShowModal() == wx.ID_NO):
                                dlg2.Destroy()
                                logInfo("Cancelled save operation due to filename already existing")
                                continue
                            else:
                                os.remove(str(filename))
                            dlg2.Destroy()
                        prefix = filename.split(".ensb")[0]
                        goToSandbox()
                        os.rename("results.ensb", filename)
                        self.unpackEnsemble(filename, prefix)
                        break
                    else:
                        dlg2 = wx.MessageDialog(self, "Are you sure you want to cancel?  Your results will be lost!", "Confirm Cancellation", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
                        if (dlg2.ShowModal() == wx.ID_YES):
                            dlg2.Destroy()
                            logInfo("Refused to save backrub output")
                            break
                        dlg2.Destroy()
                logInfo("Found structure backrub output at results.ensb")
                # Pop this message out of the queue
                for i in range(0, len(self.seqWin.msgQueue)):
                    if (self.seqWin.msgQueue[i].find("Performing structure backrub") >= 0):
                        self.seqWin.msgQueue.pop(i)
                        break
                if (len(self.seqWin.msgQueue) > 0):
                    self.seqWin.labelMsg.SetLabel(self.seqWin.msgQueue[len(self.seqWin.msgQueue)-1])
                else:
                    self.seqWin.labelMsg.SetLabel("")
                self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
                self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")
                self.parent.GoBtn.Enable()
                self.btnRelax.Enable()
                if (platform.system() == "Darwin"):
                    self.btnBackrub.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/ensemblegen/btnBackrub.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                else:
                    self.btnBackrub.SetLabel("Backrub!")
                self.buttonState = "Ensemble!"
                self.btnBackrub.SetToolTipString("Begin backrub protocol to generate a structural ensemble")
                self.seqWin.cannotDelete = False
            elif (os.path.isfile("errreport")):
                self.tmrBackrub.Stop()
                self.recoverFromError()
            elif (os.path.isfile("backrubprogress")):
                f = open("backrubprogress", "r")
                data = f.readlines()
                f.close()
                if (len(data) == 0):
                    indx = 0
                try:
                    indx = int(data[0])
                except:
                    indx = 0
                (keepGoing, skip) = self.progress.Update(indx)
                if (not(keepGoing)):
                    # User clicked "Cancel" on the progress bar
                    self.cancelBackrub()
                    self.progress.Destroy()
    
    def unpackEnsemble(self, filename, prefix):
        # Unpack the given archive while using the indicated prefix for saving files
        gzipfile = gzip.open(filename, "rb")
        readingData = False
        for aline in gzipfile:
            if (aline.startswith("BEGIN PDB")):
                if (filename.endswith(".msdar")):
                    if (platform.system() == "Windows"):
                        f = open(self.seqWin.cwd + "\\" + aline.split()[len(aline.split())-1].strip(), "w")
                    else:
                        f = open(self.seqWin.cwd + "/" + aline.split()[len(aline.split())-1].strip(), "w")
                else:
                    indx = aline[aline.rfind("_")+1:].strip()
                    f = open(prefix + "_" + indx, "w")
                readingData = True
            elif (aline.startswith("END PDB")):
                f.close()
                readingData = False
            elif (readingData):
                f.write(aline.strip() + "\n")
        gzipfile.close()
    
    def unpackClick(self, event):
        # Is the prefix valid?
        if (len(self.txtPrefix.GetValue().strip()) == 0 or not(self.txtPrefix.GetValue()[0] in "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789")):
            dlg = wx.MessageDialog(self, "Please provide a valid filename prefix.", "Invalid File Prefix", wx.OK | wx.ICON_ERROR | wx.CENTRE)
            dlg.ShowModal()
            dlg.Destroy()
            return
        dlg = wx.FileDialog(
            self, message="Choose a File",
            defaultDir=self.seqWin.cwd,
            defaultFile="",
            wildcard="Ensemble Archives (*.ensb,*.msdar)|*.ensb;*.msdar",
            style=wx.OPEN | wx.CHANGE_DIR)
        if (dlg.ShowModal() == wx.ID_OK):
            filename = str(dlg.GetPath())
            # Change cwd to the last opened file
            if (platform.system() == "Windows"):
                lastDirIndx = filename.rfind("\\")
            else:
                lastDirIndx = filename.rfind("/")
            self.seqWin.cwd = filename[0:lastDirIndx]
            self.seqWin.saveWindowData(None)
            prefix = self.seqWin.cwd + "/" + self.txtPrefix.GetValue()
            # Unpack the archives contents
            try:
                self.unpackEnsemble(filename, prefix)
            except:
                dlg = wx.MessageDialog(self, "The ensemble archive could not be read properly!", "Ensemble Archive Could Not Be Unpacked", wx.OK | wx.ICON_ERROR | wx.CENTRE)
                dlg.ShowModal()
                dlg.Destroy()
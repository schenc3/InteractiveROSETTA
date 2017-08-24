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

class AntibodyPanel(wx.lib.scrolledpanel.ScrolledPanel):
    def __init__(self, parent, W, H):
        #if (platform.system() == "Windows"):
        wx.lib.scrolledpanel.ScrolledPanel.__init__(self, parent, id=-1, pos=(10, 60), size=(340, H-330), name="ProtAntibody")
        winh = H-330
        #else:
            #wx.lib.scrolledpanel.ScrolledPanel.__init__(self, parent, id=-1, pos=(10, 60), size=(340, H-330), name="ProtMinimization")
            #winh = H-290
        self.SetBackgroundColour("#333333")
        self.parent = parent
        
        if (platform.system() == "Windows"):
            self.lblProt = wx.StaticText(self, -1, "Antibody Modeling", (25, 15), (270, 25), wx.ALIGN_CENTRE)
            self.lblProt.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblProt = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/antibody/lblAntibody.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(25, 15), size=(270, 25))
        else:
            self.lblProt = wx.StaticText(self, -1, "Antibody Modeling", (70, 15), style=wx.ALIGN_CENTRE)
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
            self.lblInst = wx.StaticText(self, -1, "Model antibodies on a remote server", (0, 45), (320, 25), wx.ALIGN_CENTRE)
            self.lblInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
        elif (platform.system() == "Darwin"):
            self.lblInst = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/antibody/lblInstAntibody.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 45), size=(320, 25))
        else:
            self.lblInst = wx.StaticText(self, -1, "Model antibodies on a remote server", (5, 45), style=wx.ALIGN_CENTRE)
            self.lblInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
            resizeTextControlForUNIX(self.lblInst, 0, self.GetSize()[0])
        self.lblInst.SetForegroundColour("#FFFFFF")
        
        if (platform.system() == "Windows"):
            self.lblLightChain = wx.StaticText(self, -1, "Variable Light Chain Sequence:", (10, 90), (320, 20), wx.ALIGN_LEFT)
            self.lblLightChain.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblLightChain = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/antibody/lblLightChain.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(10, 90), size=(320, 20))
        else:
            self.lblLightChain = wx.StaticText(self, -1, "Variable Light Chain Sequence:", (10, 90), style=wx.ALIGN_LEFT)
            self.lblLightChain.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.lblLightChain.SetForegroundColour("#FFFFFF")
        self.txtLightChain = wx.TextCtrl(self, -1, pos=(0, 110), size=(320, 75), style=wx.TE_MULTILINE)
        self.txtLightChain.SetValue("")
        self.txtLightChain.SetToolTipString("Antibody variable light chain protein sequence")
        
        if (platform.system() == "Windows"):
            self.lblHeavyChain = wx.StaticText(self, -1, "Variable Heavy Chain Sequence:", (10, 200), (320, 20), wx.ALIGN_LEFT)
            self.lblHeavyChain.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblHeavyChain = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/antibody/lblHeavyChain.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(10, 200), size=(320, 20))
        else:
            self.lblHeavyChain = wx.StaticText(self, -1, "Variable Heavy Chain Sequence:", (10, 200), style=wx.ALIGN_LEFT)
            self.lblHeavyChain.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.lblHeavyChain.SetForegroundColour("#FFFFFF")
        self.txtHeavyChain = wx.TextCtrl(self, -1, pos=(0, 220), size=(320, 75), style=wx.TE_MULTILINE)
        self.txtHeavyChain.SetValue("")
        self.txtHeavyChain.SetToolTipString("Antibody variable heavy chain protein sequence")
        
        if (platform.system() == "Windows"):
            self.lblNumModels = wx.StaticText(self, -1, "Models to Generate:", (0, 313), (260, 20), wx.ALIGN_CENTRE)
            self.lblNumModels.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblNumModels = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/antibody/lblNumModels.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 313), size=(260, 20))
        else:
            self.lblNumModels = wx.StaticText(self, -1, "Models to Generate:", (0, 313), style=wx.ALIGN_CENTRE)
            self.lblNumModels.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblNumModels, 0, 260)
        self.lblNumModels.SetForegroundColour("#FFFFFF")
        self.txtNumModels = wx.TextCtrl(self, -1, pos=(260, 310), size=(60, 25))
        self.txtNumModels.SetValue("10")
        self.txtNumModels.SetToolTipString("Number of antibody models to generate (1-100)")
        
        if (platform.system() == "Darwin"):
            self.btnImmunize = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/antibody/btnImmunize.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(110, 350), size=(100, 25))
        else:
            self.btnImmunize = wx.Button(self, id=-1, label="Immunize!", pos=(110, 350), size=(100, 25))
            self.btnImmunize.SetForegroundColour("#000000")
            self.btnImmunize.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.btnImmunize.Bind(wx.EVT_BUTTON, self.immunizeClick)
        self.btnImmunize.SetToolTipString("Begin threading the sequence onto the template structures")
        
        self.scrollh = self.btnImmunize.GetPosition()[1] + self.btnImmunize.GetSize()[1] + 5
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
    
    def scrolled(self, event):
        self.winscrollpos = self.GetScrollPos(wx.VERTICAL)
        event.Skip()
    
    def activate(self):
        self.Scroll(0, self.winscrollpos)
    
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
    
    def immunizeClick(self, event):
        # This is also the "Finalize!" button
        logInfo("Immunize button clicked")
        goToSandbox()
        # Are both of the sequences valid?
        LC_seq = self.txtLightChain.GetValue()
        temp = ""
        for line in LC_seq.split("\n"):
            if (len(line.strip()) > 0 and line.strip()[0] != ">"):
                temp += line
        LC_seq = temp
        LC_seq = LC_seq.replace(" ", "").replace("\t", "").replace("\n", "")
        LC_seq = LC_seq.upper()
        for char in LC_seq:
            if (not(char in "ACDEFGHIKLMNPQRSTVWY")):
                dlg = wx.MessageDialog(self, "The sequence for the light chain contains an invalid character: " + char, "Invalid Amino Acid", wx.OK | wx.ICON_ERROR | wx.CENTRE)
                dlg.ShowModal()
                dlg.Destroy()
                return
        if (len(LC_seq) == 0):
            dlg = wx.MessageDialog(self, "You have not provided a sequence for the light chain!", "No Light Chain", wx.OK | wx.ICON_ERROR | wx.CENTRE)
            dlg.ShowModal()
            dlg.Destroy()
            return
        HC_seq = self.txtHeavyChain.GetValue()
        temp = ""
        for line in HC_seq.split("\n"):
            if (len(line.strip()) > 0 and line.strip()[0] != ">"):
                temp += line
        HC_seq = temp
        HC_seq = HC_seq.replace(" ", "").replace("\t", "").replace("\n", "")
        HC_seq = HC_seq.upper()
        for char in HC_seq:
            if (not(char in "ACDEFGHIKLMNPQRSTVWY")):
                dlg = wx.MessageDialog(self, "The sequence for the heavy chain contains an invalid character: " + char, "Invalid Amino Acid", wx.OK | wx.ICON_ERROR | wx.CENTRE)
                dlg.ShowModal()
                dlg.Destroy()
                return
        if (len(HC_seq) == 0):
            dlg = wx.MessageDialog(self, "You have not provided a sequence for the heavy chain!", "No Heavy Chain", wx.OK | wx.ICON_ERROR | wx.CENTRE)
            dlg.ShowModal()
            dlg.Destroy()
            return
        # Is the number of models provided valid?
        try:
            nmodels = int(self.txtNumModels.GetValue())
            if (nmodels < 1 or nmodels > 100):
                raise Exception()
        except:
            dlg = wx.MessageDialog(self, "The number of models is not valid.  Please provide an integer between 1 and 100.", "Number of Models Invalid", wx.OK | wx.ICON_ERROR | wx.CENTRE)
            dlg.ShowModal()
            dlg.Destroy()
            return
        # Write the antibody input file
        fin = open("antibodyinputtemp", "w")
        fin.write("LCSEQ\t" + LC_seq + "\n")
        fin.write("HCSEQ\t" + HC_seq + "\n")
        fin.write("NMODELS\t" + str(nmodels) + "\n")
        fin.close()
        # Send it to the server
        try: 
            self.ID = sendToServer("antibodyinput")
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
                    if (len(aline.split("\t")) >= 2 and aline.split("\t")[0] == "ANTIBODY" and aline.split("\t")[1] == self.ID.strip()):
                        alreadythere = True
                        break
                f.close()
            except:
                pass
            if (not(alreadythere)):
                f = open("downloadwatch", "a")
                f.write("ANTIBODY\t" + self.ID.strip() + "\t" + str(datetime.datetime.now().strftime("%A, %B %d - %I:%M:%S %p")) + "\t" + getServerName() + "\t" + desc + "\n")
                f.close()
            dlg = wx.MessageDialog(self, "InteractiveROSETTA is now watching the server for job ID " + desc.strip() + ".  You will be notified when the package is available for download.", "Listening for Download", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
            dlg.ShowModal()
            dlg.Destroy()
            logInfo("Antibody design input sent to server daemon with ID " + self.ID)
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
### AUTHOR: Evan H. Baugh
### Affiliation: New York University

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
import datetime
from threading import Thread
from tools import *
from io_tools import *

class BioToolsPanel(wx.lib.scrolledpanel.ScrolledPanel):
    ### parent: A handle to the parent of this panel (a panel in the protocols window)
    ###         The protocols frame is the grandparent
    ### W: The width of protocols frame, you should not need to change the sizes of this panel
    ### H: The height of the protocols frame, again do not attempt to change the size osof this panel
    ### W, H are there is case you need the values
    def __init__(self, parent, W, H):
        wx.lib.scrolledpanel.ScrolledPanel.__init__(self, parent, id=-1, pos=(10, 60), size=(340, H-330), name="ProtFixbb")
        winh = H-330
        self.SetBackgroundColour("#333333")
        self.parent = parent

        # This is the message that is displayed underneath the primary sequence
        self.runMsg = "Running an example job..."
        
        # Module Title
        if (platform.system() == "Windows"):
            self.lblProt = wx.StaticText(self, -1, "Biological Tools", (25, 15), (270, 25), wx.ALIGN_CENTRE)
            self.lblProt.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblProt = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/biotools/lblTitle.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(25, 15), size=(270, 25))
        else:
            self.lblProt = wx.StaticText(self, -1, "Template Module", (70, 15), style=wx.ALIGN_CENTRE)
            self.lblProt.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
            resizeTextControlForUNIX(self.lblProt, 0, self.GetSize()[0])
        self.lblProt.SetForegroundColour("#FFFFFF")
        
        # Help button, shows the documentation in the help folder
        if (platform.system() == "Darwin"):
            # You shouldn't need to change the location of this image because it will use the default one
            self.HelpBtn = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/HelpBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(295, 10), size=(25, 25))
        else:
            self.HelpBtn = wx.Button(self, id=-1, label="?", pos=(295, 10), size=(25, 25))
            self.HelpBtn.SetForegroundColour("#0000FF")
            self.HelpBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.HelpBtn.Bind(wx.EVT_BUTTON, self.showHelp)
        self.HelpBtn.SetToolTipString("Display the help file for this window")
        
        # Brief description of the module
        if (platform.system() == "Windows"):
            self.lblInst = wx.StaticText(self, -1, "Here is a short description", (0, 45), (320, 25), wx.ALIGN_CENTRE)
            self.lblInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
        elif (platform.system() == "Darwin"):
            self.lblInst = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/biotools/lblInst.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 45), size=(320, 25))
        else:
            self.lblInst = wx.StaticText(self, -1, "Here is a short description", (5, 45), style=wx.ALIGN_CENTRE)
            self.lblInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
            resizeTextControlForUNIX(self.lblInst, 0, self.GetSize()[0])
        self.lblInst.SetForegroundColour("#FFFFFF")
        
        ### =================================================================================
        ### INSERT YOUR CONTROLS HERE
        ### ...
        ### ...
        ### ...
        ### =================================================================================
        
        ### =================================================================================
        ### HERE ARE SOME EXAMPLES OF CONTROLS FOR YOUR REFERENCE:
        ### DELETE THIS BLOCK WHEN YOU ARE DOING THIS FOR REAL
        ### TEXT LABEL, DOES NOT ALLOW USER INPUT
        if (platform.system() == "Windows"):
            self.lblLabel = wx.StaticText(self, -1, "I am a label", pos=(10, 80), size=(320, 20), style=wx.ALIGN_LEFT)
            self.lblLabel.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblLabel = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/biotools/lblLabel.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(10, 80), size=(320, 20))
        else:
            self.lblLabel = wx.StaticText(self, -1, "I am a label", pos=(10, 80), style=wx.ALIGN_LEFT)
            self.lblLabel.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.lblLabel.SetForegroundColour("#FFFFFF")
        
        ### TEXT ENTRY BOX
        self.txtTextBox = wx.TextCtrl(self, -1, pos=(0, 110), size=(320, 25), style=wx.TE_MULTILINE)
        self.txtTextBox.SetValue("Enter some text here")
        self.txtTextBox.SetToolTipString("I am a Text Box")
        
        ### STANDARD TEXT BUTTON, THIS EXAMPLE SETS UP A FUNCTION TO TOGGLE BETWEEN STATES
        if (platform.system() == "Darwin"):
            self.btnButton = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/biotools/btnButton_1.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(110, 140), size=(100, 25))
        else:
            self.btnButton = wx.Button(self, id=-1, label="State 1", pos=(110, 140), size=(100, 25))
            self.btnButton.SetForegroundColour("#000000")
            self.btnButton.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnButton.Bind(wx.EVT_BUTTON, self.toggleExample)
        self.btnButton.SetToolTipString("I am a button")
        self.toggleState = 1
        
        ### BITMAP BUTTON, DISPLAYS A MESSAGE DIALOG WHEN CLICKED
        self.btnBitmapButton = wx.BitmapButton(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/colorwheel.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(110, 170), size=(100, 25))
        self.btnBitmapButton.Bind(wx.EVT_BUTTON, self.bitmapButtonClick)
        self.btnBitmapButton.SetToolTipString("I am a bitmap button")
        
        ### GRID CONTROL, DISPLAYS VARIOUS SELECTIONS IN THE ACTIVATE FUNCTION
        self.grdGrid = wx.grid.Grid(self)
        self.grdGrid.CreateGrid(0, 3)
        self.grdGrid.SetSize((320, 200))
        self.grdGrid.SetPosition((0, 200))
        self.grdGrid.SetLabelFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.grdGrid.DisableDragColSize()
        self.grdGrid.DisableDragRowSize()
        self.grdGrid.SetColLabelValue(0, "Residues")
        self.grdGrid.SetColLabelValue(1, "Chains")
        self.grdGrid.SetColLabelValue(2, "Models")
        self.grdGrid.SetRowLabelSize(20)
        self.grdGrid.SetColSize(0, 100)
        self.grdGrid.SetColSize(1, 100)
        self.grdGrid.SetColSize(2, 100)
        
        ### STANDARD TEXT BUTTON, THIS EXAMPLE SETS UP A FUNCTION TO TOGGLE BETWEEN STATES
        if (platform.system() == "Darwin"):
            self.btnServer = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/biotools/btnServer_Off.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(110, 430), size=(100, 25))
        else:
            self.btnServer = wx.Button(self, id=-1, label="Server Off", pos=(110, 430), size=(100, 25))
            self.btnServer.SetForegroundColour("#000000")
            self.btnServer.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.btnServer.Bind(wx.EVT_BUTTON, self.toggleServer)
        self.btnServer.SetToolTipString("I am a server toggle button")
        self.serverOn = False
        
        ### END EXAMPLE CONTROLS BLOCK
        ### END DELETION POINT
        ### =================================================================================
        
        # This is the submit button
        # It calls the function "submitClick"
        # You can change its coordinates but make sure it is the farthest thing down the panel
        # so the scrolled panel gets implemented correctly
        if (platform.system() == "Darwin"):
            self.btnSubmit = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/biotools/btnSubmit.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(110, 460), size=(100, 25))
        else:
            self.btnSubmit = wx.Button(self, id=-1, label="Submit!", pos=(110, 460), size=(100, 25))
            self.btnSubmit.SetForegroundColour("#000000")
            self.btnSubmit.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.btnSubmit.Bind(wx.EVT_BUTTON, self.submitClick)
        self.btnSubmit.SetToolTipString("Submit the job")
        self.buttonState = "Submit!"
        
        #lastControl = <insert_your_lowest_control_here>
        # lastControl needs to be specified so the scrollable area gets calculated properly
        lastControl = self.btnSubmit
        
        self.scrollh = lastControl.GetPosition()[1] + lastControl.GetSize()[1] + 5
        self.SetScrollbars(1, 1, 320, self.scrollh)
        self.winscrollpos = 0
        self.Bind(wx.EVT_SCROLLWIN, self.scrolled)
    
    # This function displays the help page for the module
    # It should be located in help/biotools.html
    def showHelp(self, event):
        # Open the help page
        if (platform.system() == "Darwin"):
            try:
                browser = webbrowser.get("Safari")
            except:
                print "Could not load Safari!  The help files are located at " + self.helpdir
                return
            browser.open(self.parent.parent.scriptdir + "/help/index.html")
        else:
            webbrowser.open(self.parent.parent.scriptdir + "/help/index.html")
    
    def scrolled(self, event):
        self.winscrollpos = self.GetScrollPos(wx.VERTICAL)
        event.Skip()
    
    def activate(self):
        # This function gets called everytime the protocol gains focuses
        # Here is some code that stores current selection information
        # You can keep it and use it to your advantage if you so choose
        selectiondata = self.seqWin.getSelectedResidues()
        self.selectedResidues = []
        self.selectedChains = []
        self.selectedModels = []
        topLefts = self.seqWin.SeqViewer.GetSelectionBlockTopLeft()
        bottomRights = self.seqWin.SeqViewer.GetSelectionBlockBottomRight()
        for i in range(0, len(topLefts)):
            for r in range(topLefts[i][0], bottomRights[i][0]+1):
                for c in range(topLefts[i][1], bottomRights[i][1]+1):
                    poseindx = self.seqWin.getPoseIndex(r)
                    resi = self.seqWin.indxToSeqPos[r][c][1] # Integer returned
                    temp = self.seqWin.IDs[r]
                    model = temp[0:len(temp)-2]
                    chain = temp[len(temp)-1]
                    resn = self.seqWin.sequences[r][c]
                    self.selectedResidues.append([resn, resi, model, chain, r, c, poseindx])
                    if (not(temp in self.selectedChains)):
                        self.selectedChains.append(temp)
                    if (not(model in self.selectedModels)):
                        self.selectedModels.append(model)
        ### Each row of self.selectedResidues has the following information:
        ### 0: resn = One-letter amino acid code
        ### 1: resi = PDB index of residue
        ### 2: model = The name of the whole structure the residue is in
        ### 3: chain = The one letter chain ID
        ### 4: r = The row in the Sequence Viewer
        ### 5: c = The column in the Sequence Viewer
        ### 6: poseindx = The location of the BioPython structure object in self.seqWin.poses
        ### =================================================================================
        # THIS UPDATES THE EXAMPLE GRID WITH THE SELECTION INFORMATION
        # YOU MAY REMOVE THIS WHEN YOU REMOVE THE EXAMPLE GRID
        if (self.grdGrid.NumberRows > 0):
            self.grdGrid.DeleteRows(0, self.grdGrid.NumberRows())
        self.grdGrid.AppendRows(len(self.selectedResidues))
        i = 0
        for [resn, resi, model, chain, r, c, poseindx] in self.selectedResidues:
            self.grdGrid.SetCellValue(i, 0, model + "|" + chain + ":" + resn + str(resi))
            i += 1
        i = 0
        for chain in self.selectedChains:
            self.grdGrid.SetCellValue(i, 1, chain)
            i += 1
        i = 0
        for model in self.selectedModels:
            self.grdGrid.SetCellValue(i, 2, model)
            i += 1
        ### END DELETION
        ### =================================================================================
        # You need to keep the following line here
        # Otherwise, whenever you leave the protocol and then come back to it, it will 
        # automatically scroll to include whatever element last had the focus, which gets incredibly
        # annoying for the user.  The following will prevent it from moving
        self.Scroll(0, self.winscrollpos)
    
    # Gives this module a handle on the sequence window
    def setSeqWin(self, seqWin):
        self.seqWin = seqWin
        # So the sequence window knows about what model "designed_view" really is
        self.seqWin.setProtocolPanel(self)
        
    # Gives this module a handle on the PyMOL command line
    def setPyMOL(self, pymol):
        self.pymol = pymol
        self.cmd = pymol.cmd
        self.stored = pymol.stored
        
    # Gives this module a handle on the selection panel
    def setSelectWin(self, selectWin):
        self.selectWin = selectWin
        self.selectWin.setProtPanel(self)
    
    ### =================================================================================
    ### THESE ARE EXAMPLE EVENT FUNCTIONS
    ### YOU MAY DELETE THESE AFTER REMOVING THE EXAMPLE CONTROLS
    def toggleExample(self, event):
        if (self.toggleState == 1):
            self.toggleState = 2
            if (platform.system() == "Darwin"):
                self.btnButton.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/biotools/btnButton_2.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnButton.SetLabel("State 2")
        elif (self.toggleState == 2):
            self.toggleState = 3
            if (platform.system() == "Darwin"):
                self.btnButton.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/biotools/btnButton_3.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnButton.SetLabel("State 3")
        else:
            self.toggleState = 1
            if (platform.system() == "Darwin"):
                self.btnButton.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/biotools/btnButton_1.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnButton.SetLabel("State 1")
    
    def bitmapButtonClick(self, event):
        dlg = wx.MessageDialog(self, "This example is nice isn't it?", "Question", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
        result = dlg.ShowModal()
        if (result == wx.ID_YES):
            dlg2 = wx.MessageDialog(self, "You pressed yes", "Result", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
            dlg2.ShowModal()
            dlg2.Destroy()
        else:
            dlg2 = wx.MessageDialog(self, "You pressed no", "Result", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
            dlg2.ShowModal()
            dlg2.Destroy()
        dlg.Destroy()
        
    def toggleServer(self, event):
        if (self.serverOn):
            self.serverOn = False
            if (platform.system() == "Darwin"):
                self.btnServer.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/biotools/btnServer_Off.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnServer.SetLabel("Server Off")
        else:
            self.serverOn = True
            if (platform.system() == "Darwin"):
                self.btnServer.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/biotools/btnServer_On.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnServer.SetLabel("Server On")
    ### END DELETION
    ### =================================================================================
    
    # This function should be called to officially load a structure into InteractiveROSETTA
    # Temporary PDBs downloaded from the Internet should be downloaded to the "sandbox" first
    # To get to the sandbox, simply execute "goToSandbox()"
    # On Windows, the sandbox is C:\Users\username\InteractiveROSETTA (the folder is hidden)
    # On OSX, the sandbox is /Users/username/.InteractiveROSETTA
    # On Linux, the sandobx is /home/username/.InteractiveROSETTA
    # The sandbox gets cleaned up at the start of every session, so anything that is left there
    # will get deleted the next time InteractiveROSETTA starts up
    def loadPDB(self, pdbfile):
        self.seqWin.PyMOLPDBLoad(1, pdbfile, "Show")
    
    # This function cancels the submission
    def cancelSubmit(self):
        logInfo("Canceled " + self.protocolcode + " operation")
        try:
            os.remove(self.protocolcode + "input")
        except:
            pass
        try:
            os.remove(self.protocolcode + "inputtemp")
        except:
            pass
        self.tmrSubmit.Stop()
        self.seqWin.cannotDelete = False
        self.enableControls()
        if (platform.system() == "Darwin"):
            self.btnSubmit.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/biotools/btnSubmit.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
        else:
            self.btnSubmit.SetLabel("Submit!")
        self.btnSubmit.SetToolTipString("Submit the job")
        deleteInputFiles()
        self.parent.parent.restartDaemon()
        self.parent.GoBtn.Enable()
        self.removeMessage()
        self.buttonState = "Submit!"
    
    # This function adds your message to the label underneath the sequence viewer
    # If you want a new message at any point, set it in "self.runMsg" and call this function
    def updateMessage(self):
        self.seqWin.labelMsg.SetLabel(self.runMsg)
        self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")
        self.seqWin.msgQueue.append(self.runMsg)
        
    # This function removes your message from underneath the sequence viewer
    # Do not change the value of "self.runMsg" before calling this function or it will not
    # be removed properly
    def removeMessage(self):
        # Get rid of the messages
        for i in range(0, len(self.seqWin.msgQueue)):
            if (self.seqWin.msgQueue[i].find(self.runMsg) >= 0):
                self.seqWin.msgQueue.pop(i)
                break
        if (len(self.seqWin.msgQueue) > 0):
            self.seqWin.labelMsg.SetLabel(self.seqWin.msgQueue[len(self.seqWin.msgQueue)-1])
        else:
            self.seqWin.labelMsg.SetLabel("")
        self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")
    
    # This is the function that gets called when a submission occurs
    # It can turn into the Cancel and Finalize buttons if those are appropriate for your module
    def submitClick(self, event):
        # This is also the "Finalize!" button
        logInfo("Submit button clicked")
        if (self.buttonState == "Submit!"):
            ### DO OTHER THINGS AFTER SUBMISSION
            ### e.g. DISABLE SOME CONTROLS, PREPARE INPUTS, CHECK THAT INPUTS ARE VALID
            ### ...
            ### ...
            self.updateMessage()
            # Disable protocol changing and sequence editing while the protocol is active
            self.parent.GoBtn.Disable()
            self.seqWin.cannotDelete = True
            # Uses a Timer to handle interactions with the daemon
            self.stage = 1
            if (platform.system() == "Darwin"):
                self.btnSubmit.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/biotools/btnSubmit_Cancel.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnSubmit.SetLabel("Cancel!")
            self.buttonState = "Cancel!"
            self.btnSubmit.SetToolTipString("Cancel the " + self.protocolcode + " simulation")
            self.tmrSubmit = wx.Timer(self)
            self.Bind(wx.EVT_TIMER, self.threadSubmit, self.tmrSubmit)
            self.tmrSubmit.Start(1000)
        elif (self.buttonState == "Cancel!"):
            dlg = wx.MessageDialog(self, "Are you sure you want to cancel the " + self.protocolcode + " simulation?  All progress will be lost.", "Cancel " + self.protocolcode + " Simulation", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
            result = dlg.ShowModal()
            if (result == wx.ID_YES):
                self.cancelSubmit()
            dlg.Destroy()
        else:
            # Finalize button, ask whether the changes will be accepted or rejected
            dlg = wx.MessageDialog(self, "Do you want to accept the results of this " + self.protocolcode + " simulation?", "Accept/Reject " + self.protocolcode, wx.YES_NO | wx.CANCEL | wx.ICON_QUESTION | wx.CENTRE)
            result = dlg.ShowModal()
            if (result == wx.ID_YES):
                accept = True
                logInfo(self.protocolcode + " accepted")
            elif (result == wx.ID_NO):
                accept = False
                logInfo(self.protocolcode + " rejected")
            else:
                dlg.Destroy()
                logInfo("Finalize operation cancelled")
                return
            dlg.Destroy()
            if (platform.system() == "Darwin"):
                self.btnSubmit.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/biotools/btnSubmit.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnSubmit.SetLabel("Submit!")
            self.buttonState = "Submit!"
            self.btnSubmit.SetToolTipString("Perform " + self.protocolcode)
            self.seqWin.cannotDelete = False
            ### DO THINGS COMMON TO BOTH ACCEPT AND REJECT
            ### ...
            ### ...
            if (not(accept)):
                ### HANDLE A REJECTION
                ### ...
                ### ...
                ### ===========================================================================
                ### EXAMPLE CODE TO DELETE THE OUTPUTTED MODEL FROM PYMOL
                ### DELETE ME LATER
                self.cmd.remove("example_view")
                ### END DELETION
                ### ===========================================================================
                return
            ### HANDLE AN ACCEPT
            ### ...
            ### ...
            ### ===============================================================================
            ### EXAMPLE CODE TO LOAD THE MODEL INTO THE SEQUENCEVIEWER
            ### DELETE ME LATER
            goToSandbox()
            self.cmd.remove("example_view")
            self.seqWin.PyMOLPDBLoad(1, "output.pdb", showDialog="Show")
            ### END DELETION
            ### ===============================================================================
    
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
        self.enableControls(True)
        if (platform.system() == "Darwin"):
            self.btnSubmit.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/biotools/btnSubmit.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
        else:
            self.btnSubmit.SetLabel("Design!")
        self.buttonState = "Submit!"
        self.removeMessage()
    
    def enableControls(self, enabled=True):
        ### Enable or disable controls in this function, using enabled as the boolean to
        ### either enable (True) or disable (False)
        ### Example:
        self.btnSubmit.Enable(enabled)
    
    def threadSubmit(self, event):
        daemon_jobID = self.protocolcode
        inputfile = daemon_jobID + "input"
        outputfile = daemon_jobID + "output"
        # Why am I using a Timer?  See the explanation in kic.py
        goToSandbox()
        if (self.stage == 1):
            self.tmrSubmit.Stop()
            self.timeoutCount = 0
            self.progress = None
            fout = open(inputfile + "temp", "w")
            ### WRITE THE INPUTS INTO A FILE THAT THE DAEMON CAN READ
            ### IF THIS IS GOING TO A REMOTE SERVER YOU NEED TO HAVE ONE OF THE FOLLOWING:
            ### REQUESTED_NODES\t<VALUE>
            ### REQUIRED_NODES\t<VALUE>
            ### IF REQUESTED_NODES IS USED, THEN THE SERVER WILL TRY TO GET AT MOST <VALUE> NODES
            ### IF REQUIRED_NODES IS USED, THEN THE SERVER WILL GET AT LEAST <VALUE> NODES
            ### ...
            ### ...
            ### ================================================================================
            ### EXAMPLE FOR LOCAL DAEMON
            ### DELETE ME LATER
            fout.write("PDBFILE\t1adw_1.pdb\n")
            fout.write("BEGIN PDB DATA\n")
            fin = open(self.datadir + "/1adw_1.pdb")
            for aline in fin:
                fout.write(aline)
            fin.close()
            fout.write("END PDB DATA\n")
            fout.write("RESFILE\t1adw.resfile\n")
            fout.write("BEGIN RESFILE DATA\n")
            fin = open(self.datadir + "/1adw.resfile")
            for aline in fin:
                fout.write(aline)
            fin.close()
            fout.write("END RESFILE DATA\n")
            fout.write("REQUIRED_NODES\t1\n")
            ### ================================================================================
            fout.close()
            ### ADD THE PARAMS AND SCOREFXN DATA TO THE OUTPUTFILE
            ### IT HAS THE FOLLOWING FORMAT:
            ### PARAMS\t<FILENAME>
            ### BEGIN PARAMS DATA
            ### <DATA>
            ### END PARAMS DATA
            ### SCOREFXN\t<FILENAME>
            ### BEGIN SCOREFXN DATA
            ### <DATA>
            ### END SCOREFXN DATA
            appendScorefxnParamsInfoToFile(inputfile + "temp", self.selectWin.weightsfile)
            serverToggleable = False
            try:
                # To find out if the developer defined the serverOn variable
                self.serverOn
                serverToggleable = True
            except:
                pass
            if (self.useServer or (serverToggleable and self.serverOn)):
                try: 
                    ### USE THE CODE BELOW IF YOU WANT TO SEND TO THE USER'S SELECTED SERVER
                    serverName = ""
                    home = os.path.expanduser("~")
                    if (platform.system() == "Windows"):
                        fin = open(home + "/InteractiveROSETTA/seqwindow.cfg", "r")
                    else:
                        fin = open(home + "/.InteractiveROSETTA/seqwindow.cfg", "r")
                    for aline in fin:
                        if (aline.startswith("[SERVER]")):
                            serverName = aline.split("\t")[1].strip()
                    fin.close()
                    if (len(serverName.strip()) == 0):
                        raise Exception("No server specified")
                    self.ID = sendToServer(inputfile)
                    ### ======================================================================
                    ### IF YOU WANT TO FORCE SENDING IT TO A CERTAIN SERVER
                    # serverURL = <URL>
                    # dlg = wx.MessageDialog(self, "Your job will be sent to the following server: " + serverURL + "\n\nIs this okay?  If there is sensitive data in your submission that you don't want on this server, please select No.", "Confirm Server Submission", wx.YES_NO | wx.ICON_EXCLAMATION | wx.CENTRE)
                    # if (dlg.ShowModal() == wx.ID_NO):
                        # self.cancelSubmit()
                        # return
                    # dlg.Destroy()
                    # self.ID = sendToServer(inputfile, serverURL)
                    # serverName = serverURL
                    ### ======================================================================
                    dlg = wx.TextEntryDialog(None, "Enter a description for this submission:", "Job Description", "")
                    if (dlg.ShowModal() == wx.ID_OK):
                        desc = dlg.GetValue()
                        desc = desc.replace("\t", " ").replace("\n", " ").strip()
                    else:
                        desc = self.ID
                    goToSandbox()
                    # First make sure this isn't a duplicate
                    alreadythere = False
                    try:
                        f = open("downloadwatch", "r")
                        for aline in f:
                            if (len(aline.split("\t")) >= 2 and aline.split("\t")[0] == self.protocolcode.upper() and aline.split("\t")[1] == self.ID.strip()):
                                alreadythere = True
                                break
                        f.close()
                    except:
                        pass
                    if (not(alreadythere)):
                        f = open("downloadwatch", "a")
                        f.write(self.protocolcode.upper() + "\t" + self.ID.strip() + "\t" + str(datetime.datetime.now().strftime("%A, %B %d - %I:%M:%S %p")) + "\t" + serverName + "\t" + desc + "\n")
                        f.close()
                    dlg = wx.MessageDialog(self, "InteractiveROSETTA is now watching the server for job ID " + desc.strip() + ".  You will be notified when the package is available for download.", "Listening for Download", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
                    dlg.ShowModal()
                    dlg.Destroy()
                    logInfo(self.protocolcode.upper() + " input sent to server daemon with ID " + self.ID)
                    if (platform.system() == "Darwin"):
                        self.btnSubmit.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/biotools/btnSubmit.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                    else:
                        self.btnSubmit.SetLabel("Submit!")
                    self.buttonState = "Submit!"
                    self.btnSubmit.SetToolTipString("Perform " + self.protocolcode)
                    self.seqWin.cannotDelete = False
                    self.parent.GoBtn.Enable()
                    self.enableControls(True)
                    self.removeMessage()
                    return
                except Exception as e:
                    logInfo("Server daemon not available")
                    if (serverToggleable and serverOn):
                        os.rename(inputfile + "temp", inputfile)
                        logInfo(self.protocolcode + " input uploaded locally at " + inputfile)
                    else:
                        f = open("errreport", "w")
                        f.write(e.message.strip())
                        f.close()
                        self.recoverFromError()
                        return
            else:
                os.rename(inputfile + "temp", inputfile)
                logInfo(self.protocolcode + " input uploaded locally at " + inputfile)
            self.stage = 2
            self.tmrSubmit.Start(1000)
        ### IF YOUR PROTOCOL HAS MORE THAN ONE STEP, YOU CAN IMPLEMENT THAT HERE AS
        ### MULTIPLE STAGES
        else:
            # Read the output dumped by the child process
            if (os.path.isfile(outputfile)):
                if (self.progress is not None):
                    try:
                        self.progress.Destroy()
                    except:
                        pass
                goToSandbox()
                self.tmrSubmit.Stop()
                fin = open(outputfile, "r")
                ### PARSE THE outputfile
                ### ...
                ### ...
                fin.close()
                logInfo("Found " + self.protocolcode + " output at " + outputfile)
                if (platform.system() == "Darwin"):
                    self.btnSubmit.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/biotools/btnSubmit_Finalize.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
                else:
                    self.btnSubmit.SetLabel("Finalize!")
                self.buttonState = "Finalize!"
                self.btnSubmit.SetToolTipString("Accept or reject the results of this protocol")
                os.remove(outputfile)
                self.removeMessage()
                ### DO SOMETHING WITH THE OUTPUT
                ### e.g.) DISPLAY IT IN CONTROLS/PYMOL
                ### ...
                ### ...
                ### ===========================================================================
                ### THIS IS AN EXAMPLE THAT LOADS AN OUTPUT MODEL FOR PYMOL VIEWING
                ### DELETE IT LATER
                self.cmd.load("output.pdb", "example_view")
                defaultPyMOLView(self.cmd, "example_view")
                ### ===========================================================================
            elif (os.path.isfile("errreport")):
                self.tmrSubmit.Stop()
                self.recoverFromError()
                if (self.progress is not None):
                    try:
                        self.progress.Destroy()
                    except:
                        pass
            elif (os.path.isfile("progress")):
                # The local daemon can output its progress to keep the GUI updated about
                # how far along it is, along with a message
                # This is optional
                # See job/__init__.py for more information
                if (self.progress is None):
                    self.progress = wx.ProgressDialog(self.protocolcode.upper() + " Progress", "Performing " + self.protocolcode + "...", 100, style=wx.PD_CAN_ABORT | wx.PD_APP_MODAL | wx.PD_ELAPSED_TIME | wx.PD_REMAINING_TIME)
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
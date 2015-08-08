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

### Replace <ModuleName> below with the name of your module
### This name is what will appear in the drop-down menu
### PROTOCOL NAME: <ModuleName>
### This class must be named "ModulePanel", do not change it
### You may add other classes in this script if you need them
class ModulePanel(wx.lib.scrolledpanel.ScrolledPanel):
    ### parent: A handle to the parent of this panel (a panel in the protocols window)
    ###         The protocols frame is the grandparent
    ### thisfile: The path to this __init__.py file
    ### W: The width of protocols frame, you should not need to change the sizes of this panel
    ### H: The height of the protocols frame, again do not attempt to change the size of this panel
    ### W, H are there is case you need the values
    def __init__(self, parent, thisfile, W, H):
	wx.lib.scrolledpanel.ScrolledPanel.__init__(self, parent, id=-1, pos=(10, 60), size=(340, H-330), name="ProtFixbb")
	winh = H-330
	self.SetBackgroundColour("#333333")
	self.parent = parent
	# Determines the code for this protocol
	# input/output files will be uploaded as <protocolcode>input and <protocolcode>output
	# This is the name the daemon and server use to find the code to execute for this protocol
	self.protocolcode = thisfile.split("__init__.py")[0]
	l = len(self.protocolcode)
	if (platform.system() == "Windows"):
	    indx = self.protocolcode[0:l-1].rfind("\\")
	else:
	    indx = self.protocolcode[0:l-1].rfind("/")
	self.protocolcode = self.protocolcode[indx+1:l-1]
	# The locations of the data and image directories for this module
	self.datadir = thisfile.split("__init__.py")[0] + "/data"
	self.imagedir = thisfile.split("__init__.py")[0] + "/images"
	self.helpdir = thisfile.split("__init__.py")[0] + "/help"
	# This is the message that is displayed underneath the primary sequence
	self.runMsg = "Running special fixbb..."
	### Should this run on the remote server?
	self.useServer = False
	
	# Module Title
	if (platform.system() == "Windows"):
	    self.lblProt = wx.StaticText(self, -1, "Template Module", (25, 15), (270, 25), wx.ALIGN_CENTRE)
	    self.lblProt.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
	elif (platform.system() == "Darwin"):
	    self.lblProt = wx.StaticBitmap(self, -1, wx.Image(self.imagedir + "/lblTitle.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(25, 15), size=(270, 25))
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
	    self.lblInst = wx.StaticText(self, -1, "This module is awesome", (0, 45), (320, 25), wx.ALIGN_CENTRE)
	    self.lblInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
	elif (platform.system() == "Darwin"):
	    self.lblInst = wx.StaticBitmap(self, -1, wx.Image(self.imagedir + "/lblInst.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 45), size=(320, 25))
	else:
	    self.lblInst = wx.StaticText(self, -1, "This module is awesome", (5, 45), style=wx.ALIGN_CENTRE)
	    self.lblInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
	    resizeTextControlForUNIX(self.lblInst, 0, self.GetSize()[0])
	self.lblInst.SetForegroundColour("#FFFFFF")
	
	### =================================================================================
	### INSERT YOUR CONTROLS HERE
	### ...
	### ...
	### ...
	### =================================================================================
	
	# This is the submit button
	# It calls the function "submitClick"
	# You can change its coordinates but make sure it is the farthest thing down the panel
	# so the scrolled panel gets implemented correctly
	if (platform.system() == "Darwin"):
	    self.btnSubmit = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.imagedir + "/btnSubmit.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(110, 60), size=(100, 25))
	else:
	    self.btnSubmit = wx.Button(self, id=-1, label="Submit!", pos=(110, 60), size=(100, 25))
	    self.btnSubmit.SetForegroundColour("#000000")
	    self.btnSubmit.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
	self.btnSubmit.Bind(wx.EVT_BUTTON, self.designClick)
	self.btnSubmit.SetToolTipString("Submit the job")
	self.buttonState = "Submit!"
	
	self.scrollh = self.btnSubmit.GetPosition()[1] + self.btnSubmit.GetSize()[1] + 5
	self.SetScrollbars(1, 1, 320, self.scrollh)
    
    # This function displays the help page for the module
    # It should be located at mymodule/help/index.html
    def showHelp(self, event):
	# Open the help page
	if (platform.system() == "Darwin"):
	    try:
		browser = webbrowser.get("Safari")
	    except:
		print "Could not load Safari!  The help files are located at " + self.helpdir
		return
	    browser.open(self.helpdir + "/index.html")
	else:
	    webbrowser.open(self.helpdir + "/index.html")
    
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
	    self.btnSubmit.SetBitmapLabel(bitmap=wx.Image(self.imagedir + "/btnSubmit.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
	else:
	    self.btnSubmit.SetLabel("Submit!")
	self.btnDesign.SetToolTipString("Submit the job")
	deleteInputFiles()
	self.parent.parent.restartDaemon()
	self.parent.GoBtn.Enable()
	self.removeMessage()
	self.buttonState = "Submit!"
    
    # This function adds your message to the label underneath the sequence viewer
    def updateMessage(self):
	self.seqWin.labelMsg.SetLabel(self.runMsg)
	self.seqWin.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
	self.seqWin.labelMsg.SetForegroundColour("#FFFFFF")
	self.seqWin.msgQueue.append(self.runMsg)
	
    # This function removes your message from underneath the sequence viewer
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
		self.btnSubmit.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/btnDesign_Cancel.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
	    else:
		self.btnSubmit.SetLabel("Cancel!")
	    self.buttonState = "Cancel!"
	    self.btnSubmit.SetToolTipString("Cancel the " + self.protocolcode + " simulation")
	    self.tmrSubmit = wx.Timer(self)
	    self.Bind(wx.EVT_TIMER, self.threadSubmit, self.tmrSubmit)
	    self.tmrSubit.Start(1000)
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
		self.btnSubmit.SetBitmapLabel(bitmap=wx.Image(self.imagedir + "/osx/btnDesign.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
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
		return
	    ### HANDLE AN ACCEPT
	    ### ...
	    ### ...
    
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
	    sessioninfo = os.path.expanduser("~") + "/InteractiveRosetta/sessionlog"
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
	    self.btnSubmit.SetBitmapLabel(bitmap=wx.Image(self.imagedir + "/osx/btnDesign.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
	else:
	    self.btnSubmit.SetLabel("Design!")
	self.buttonState = "Submit!"
	self.removeMessage()
    
    
    def threadSubmit(self, event):
	daemon_jobID = self.protocolcode
	inputfile = daemon_jobID + "input"
	outputfile = daemon_jobID + "output"
	# Why am I using a Timer?  See the explanation in kic.py
	goToSandbox()
	if (self.stage == 1):
	    self.tmrSubmit.Stop()
	    self.timeoutCount = 0
	    fout = open(inputfile + "temp", "w")
	    ### WRITE THE INPUTS INTO A FILE THAT THE DAEMON CAN READ
	    ### IF THIS IS GOING TO A REMOTE SERVER YOU NEED TO HAVE ONE OF THE FOLLOWING:
	    ### REQUESTED_NODES\t<VALUE>
	    ### REQUIRED_NODES\t<VALUE>
	    ### IF REQUESTED_NODES IS USED, THEN THE SERVER WILL TRY TO GET AT MOST <VALUE> NODES
	    ### IF REQUIRED_NODES IS USED, THEN THE SERVER WILL GET AT LEAST <VALUE> NODES
	    ### ...
	    ### ...
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
	    if (self.useServer):
		try: 
		    ### USE THE CODE BELOW IF YOU WANT TO SEND TO THE USER'S SELECTED SERVER
		    serverName = ""
		    home = os.path.expanduser("~")
		    fin = open(home + "/InteractiveROSETTA/seqwindow.cfg", "r")
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
			f.write(self.protocolcode.upper() + "\t" + self.ID.strip() + "\t" + str(datetime.datetime.now().strftime("%A, %B %d - %I:%M:%S %p")) + "\t" + serverName + "\n")
			f.close()
		    dlg = wx.MessageDialog(self, "InteractiveROSETTA is now watching the server for job ID " + self.ID.strip() + ".  You will be notified when the package is available for download.", "Listening for Download", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
		    dlg.ShowModal()
		    dlg.Destroy()
		    logInfo(self.protocolcode.upper() + " design input sent to server daemon with ID " + self.ID)
		except Exception as e:
		    logInfo("Server daemon not available")
		    f = open("errreport", "w")
		    f.write(e.message.strip())
		    f.close()
		    self.recoverFromError()
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
		goToSandbox()
		self.tmrSubmit.Stop()
		fin = open(outputfile, "r")
		### PARSE THE outputfile
		### ...
		### ...
		fin.close()
		logInfo("Found " + self.protocolcode + " output at " + outputfile)
		if (platform.system() == "Darwin"):
		    self.btnSubmit.SetBitmapLabel(bitmap=wx.Image(self.imagedir + "/osx/btnSubmit_Finalize.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
		else:
		    self.btnSubmit.SetLabel("Finalize!")
		self.buttonState = "Finalize!"
		self.btnSubmit.SetToolTipString("Accept or reject the results of this protocol")
		os.remove(outputfile)
		### DO SOMETHING WITH THE OUTPUT
		### e.g.) DISPLAY IT IN CONTROLS/PYMOL
		### ...
		### ...
	    elif (os.path.isfile("errreport")):
		self.tmrSubmit.Stop()
		self.recoverFromError()
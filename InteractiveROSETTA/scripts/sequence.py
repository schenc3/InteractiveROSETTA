import wx
import wx.grid
# This is needed to draw the chain colors on the grid
import wx.lib.mixins.gridlabelrenderer as glr
import os
import os.path
import shutil
import glob
import platform
import urllib
import urllib2
import gzip
import numpy as np
try:
    if (platform.system() == "Windows"):
	import hmmstr4py_win as hmmstr4py
    elif (platform.system() == "Darwin"):
	import hmmstr4py_darwin as hmmstr4py
    else:
	import hmmstr4py
except:
    print "HMMSTR not installed"
import math
import ctypes
import webbrowser
from threading import Thread
import Bio.PDB
import Bio.PDB.DSSP
from Bio.Align.Applications import MuscleCommandline
from tools import *
if (platform.system() == "Windows"):
    # pyperclip and clipboard don't seem to work, use Tk
    from Tkinter import Tk
elif (platform.system() == "Darwin"):
    # Do not use Tkinter on OSX, it crashes the whole program, use pyperclip instead
    try:
	import pyperclip
    except:
	print "pyperclip not installed, sequence copying is disabled"
	print "You can install it by executing: sudo easy_install pyperclip"
else:
    # pyperclip crashes Linux and Tk does not work, use clipboard
    try:
	import clipboard
    except:
	print "clipboard not installed, sequence copying is disabled"
	print "You can install it by executing: sudo easy_install clipboard"

# There is apparently a libc bug on UNIX that caches DNS names that urllib2 tries to open
# If the Internet connection is not valid the first time it tries to establish a connection, it caches this bad
# DNS and never is able to connect in this session, even if the Internet connection comes back
# You have to get a handle on "res_init" in libc and call it before attempting to download more things
if (platform.system() == "Linux"):
    libc = ctypes.cdll.LoadLibrary("libc.so.6")
    res_init = libc.__res_init
elif (platform.system() == "Darwin"):
    libc = ctypes.cdll.LoadLibrary("/usr/lib/libc.dylib")
    res_init = libc.res_init

# ===========================================================================================================
# STRUCTURE LOADER CLASS
# This is the dialog that is displayed when a structure is loaded, asking you what chains you want to load
# and if waters/HETATMs should be loaded

class ProteinDialog(wx.Dialog):
    def __init__(self, parent, PDB, chains, scriptdir):
	if (platform.system() != "Linux"):
	    wx.Dialog.__init__(self, parent, -1, "Protein Loader", size=(305, 330), style=wx.DEFAULT_DIALOG_STYLE | wx.STAY_ON_TOP)
	else:
	    wx.Dialog.__init__(self, parent, -1, "Protein Loader", size=(300, 330), style=wx.DEFAULT_DIALOG_STYLE | wx.STAY_ON_TOP)
	self.scriptdir = scriptdir

	# Structure label
	self.lblPDB = wx.StaticText(self, -1, PDB, (0, 10), (300, 30))
	self.lblPDB.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
	resizeTextControlForUNIX(self.lblPDB, 0, 300)
	midpoint = int(math.ceil(len(chains) / 2.0))

	# Display checkboxes for each chain
	# Space them evenly into two columns
	ygap = int(160 / midpoint)
	self.chkChains = []
	for i in range(0, len(chains)):
	    if (i < midpoint):
		self.chkChains.append(wx.CheckBox(self, -1, "Chain " + chains[i], (20, 40+(ygap*i))))
	    else:
		self.chkChains.append(wx.CheckBox(self, -1, "Chain " + chains[i], (170, 40+(ygap*(i-midpoint)))))
	    self.chkChains[i].SetValue(True)

	# Toggle buttons for loading water and HETATMs
	if (platform.system() == "Darwin"):
	    self.btnWater = wx.BitmapButton(self, -1, wx.Image(self.scriptdir + "/images/osx/sequence/btnWater_Load.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), (5, 230), (140, 30))
	else:
	    self.btnWater = wx.Button(self, id=-1, label="Load Waters", pos=(5, 230), size=(140, 30))
	    self.btnWater.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.btnWater.Bind(wx.EVT_BUTTON, self.waterToggle)
	self.btnWater.SetToolTipString("Load the protein with waters present")
	self.waterState = "Load Waters"
	if (platform.system() == "Darwin"):
	    self.btnHETATM = wx.BitmapButton(self, -1, wx.Image(self.scriptdir + "/images/osx/sequence/btnHETATM_Load.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), (155, 230), (140, 30))
	else:
	    self.btnHETATM = wx.Button(self, id=-1, label="Load HETATMs", pos=(155, 230), size=(140, 30))
	    self.btnHETATM.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.btnHETATM.Bind(wx.EVT_BUTTON, self.HETATMToggle)
	self.btnHETATM.SetToolTipString("Load the protein with heteroatoms present")
	self.hetatmState = "Load HETATMs"

	# OK and Cancel buttons
	if (platform.system() == "Darwin"):
	    self.btnOK = wx.BitmapButton(self, -1, wx.Image(self.scriptdir + "/images/osx/sequence/btnOK_Protein.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), (5, 265), (140, 30))
	else:
	    self.btnOK = wx.Button(self, id=-1, label="OK", pos=(5, 265), size=(140, 30))
	    self.btnOK.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.btnOK.Bind(wx.EVT_BUTTON, self.okDialog)
	self.btnOK.SetToolTipString("Confirm load operation")
	if (platform.system() == "Darwin"):
	    self.btnCancel = wx.BitmapButton(self, -1, wx.Image(self.scriptdir + "/images/osx/sequence/btnCancel_Protein.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), (155, 265), (140, 30))
	else:
	    self.btnCancel = wx.Button(self, id=-1, label="Cancel", pos=(155, 265), size=(140, 30))
	    self.btnCancel.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.btnCancel.Bind(wx.EVT_BUTTON, self.cancelDialog)
	self.btnCancel.SetToolTipString("Cancel this load operation")

	# Center the dialog in the middle of the screen
	self.SetPosition((wx.GetDisplaySize()[0]/2-150, wx.GetDisplaySize()[1]/2-150))

    def waterToggle(self, event):
	# Toggle loading waters on and off
	if (platform.system() == "Darwin"):
	    if (self.waterState == "Load Waters"):
		self.btnWater.SetBitmapLabel(bitmap=wx.Image(self.scriptdir + "/images/osx/sequence/btnWater_No.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
		self.btnWater.SetToolTipString("Load the protein while ignoring waters")
	    else:
		self.btnWater.SetBitmapLabel(bitmap=wx.Image(self.scriptdir + "/images/osx/sequence/btnWater_Load.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
		self.btnWater.SetToolTipString("Load the protein with waters present")
	else:
	    if (self.btnWater.GetLabel() == "Load Waters"):
		self.btnWater.SetLabel("No Waters")
		self.btnWater.SetToolTipString("Load the protein while ignoring waters")
	    else:
		self.btnWater.SetLabel("Load Waters")
		self.btnWater.SetToolTipString("Load the protein with waters present")

    def HETATMToggle(self, event):
	# Toggle loading HETATMs on and off
	if (platform.system() == "Darwin"):
	    if (self.hetatmState == "Load Waters"):
		self.btnHETATM.SetBitmapLabel(bitmap=wx.Image(self.scriptdir + "/images/osx/sequence/btnHETATM_No.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
		self.btnHETATM.SetToolTipString("Load the protein while ignoring heteroatoms")
	    else:
		self.btnHETATM.SetBitmapLabel(bitmap=wx.Image(self.scriptdir + "/images/osx/sequence/btnHETATM_Load.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
		self.btnHETATM.SetToolTipString("Load the protein with heteroatoms present")
	else:
	    if (self.btnHETATM.GetLabel() == "Load HETATMs"):
		self.btnHETATM.SetLabel("No HETATMs")
		self.btnHETATM.SetToolTipString("Load the protein while ignoring heteroatoms")
	    else:
		self.btnHETATM.SetLabel("Load HETATMs")
		self.btnHETATM.SetToolTipString("Load the protein with heteroatoms present")

    # Return codes, the main script knows how to interpret these
    def okDialog(self, event):
	self.SetReturnCode(wx.OK)
	self.EndModal(wx.OK)

    def cancelDialog(self, event):
	self.SetReturnCode(wx.CANCEL)
	self.EndModal(wx.CANCEL)

# ===========================================================================================================
# MODEL NAME AND CHAIN ID DIALOG CLASS
# This is the dialog that is displayed when the user right-clicks on a sequence label and allows them to
# change the model name and/or the chain ID

class ModelDialog(wx.Dialog):
    def __init__(self, parent, model, chain, scriptdir):
	if (platform.system() != "Linux"):
	    wx.Dialog.__init__(self, parent, -1, "Modify Model Data", size=(305, 150))
	else:
	    wx.Dialog.__init__(self, parent, -1, "Modify Model Data", size=(300, 150))
	self.scriptdir = scriptdir

	# We need to have a handle on the parent (the SeqViewer) so we can see the parent's data
	self.parent = parent
	# We also need to know what the original model name and chain ID are
	self.inModel = model
	self.inChain = chain
	if (len(chain.strip()) == 0):
	    self.inputID = model + "|_"
	else:
	    self.inputID = model + "|" + chain

	# Label and text box for the model name
	if (platform.system() == "Darwin"):
	    self.lblModel = wx.StaticBitmap(self, -1, wx.Image(self.scriptdir + "/images/osx/sequence/lblModel.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(5, 13), size=(60, 30))
	else:
	    self.lblModel = wx.StaticText(self, -1, "Model:", (5, 13), (60, 30))
	    self.lblModel.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.txtModel = wx.TextCtrl(self, -1, pos=(70, 10), size=(225, 25))
	self.txtModel.SetValue(model)
	self.txtModel.SetToolTipString("The name of the protein")

	# Label and text box for the chain ID
	if (platform.system() == "Darwin"):
	    self.lblChain = wx.StaticBitmap(self, -1, wx.Image(self.scriptdir + "/images/osx/sequence/lblChain.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(5, 53), size=(140, 30))
	else:
	    self.lblChain = wx.StaticText(self, -1, "Chain Identifier:", (5, 53), (140, 30))
	    self.lblChain.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.txtChain = wx.TextCtrl(self, -1, pos=(150, 50), size=(145, 25))
	self.txtChain.SetValue(chain)
	self.txtChain.SetToolTipString("The identifier of the chain")

	# OK and Cancel buttons
	if (platform.system() == "Darwin"):
	    self.btnOK = wx.BitmapButton(self, -1, wx.Image(self.scriptdir + "/images/osx/sequence/btnOK_Protein.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), (5, 90), (140, 30))
	else:
	    self.btnOK = wx.Button(self, id=-1, label="OK", pos=(5, 90), size=(140, 30))
	    self.btnOK.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.btnOK.Bind(wx.EVT_BUTTON, self.okDialog)
	self.btnOK.SetToolTipString("Confirm model data alteration")
	if (platform.system() == "Darwin"):
	    self.btnCancel = wx.BitmapButton(self, -1, wx.Image(self.scriptdir + "/images/osx/sequence/btnCancel_Protein.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), (155, 90), (140, 30))
	else:
	    self.btnCancel = wx.Button(self, id=-1, label="Cancel", pos=(155, 90), size=(140, 30))
	    self.btnCancel.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.btnCancel.Bind(wx.EVT_BUTTON, self.cancelDialog)
	self.btnCancel.SetToolTipString("Cancel this operation")

	# Center this dialog
	self.SetPosition((wx.GetDisplaySize()[0]/2-150, wx.GetDisplaySize()[1]/2-150))

    def okDialog(self, event):
	# Is the chainID valid?
	chainID = self.txtChain.GetValue().strip().upper()
	# No chain ID, not allowed unless it was empty going into this dialog
	if (len(chainID) == 0 and len(self.inChain.strip()) != 0):
	    dlg = wx.MessageDialog(self, "You have not specified a chain ID!", "No Chain ID", wx.OK | wx.ICON_ERROR | wx.CENTRE)
	    dlg.ShowModal()
	    dlg.Destroy()
	    return
	# Chain ID not valid, only letters are valid
	if (not(chainID in "ABCDEFGHIJKLMNOPQRSTUVWXYZ")):
	    dlg = wx.MessageDialog(self, "The chainID you have supplied is invalid.  Please choose a letter.", "Invalid Chain ID", wx.OK | wx.ICON_ERROR | wx.CENTRE)
	    dlg.ShowModal()
	    dlg.Destroy()
	    return
	# Is the model name okay?
	modelname = self.txtModel.GetValue().strip()
	for char in modelname:
	    if (not(char.upper() in "ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890_-")):
		dlg = wx.MessageDialog(self, "Your model name contains an invalid character: " + char, "Invalid Model Name", wx.OK | wx.ICON_ERROR | wx.CENTRE)
		dlg.ShowModal()
		dlg.Destroy()
		return
	# Is this modelname+chainID combination already taken?
	if (len(chainID.strip()) == 0):
	    chainID = "_"
	myID = modelname + "|" + chainID
	for ID in self.parent.IDs:
	    if (myID == ID and myID != self.inputID):
		dlg = wx.MessageDialog(self, "There is already a chain loaded with this model and chain ID combination.", "Chain Name Taken", wx.OK | wx.ICON_ERROR | wx.CENTRE)
		dlg.ShowModal()
		dlg.Destroy()
		return
	# Everything's good at this point, so return
	self.return_chain = chainID
	self.return_model = modelname
	self.SetReturnCode(wx.OK)
	self.EndModal(wx.OK)

    def cancelDialog(self, event):
	self.SetReturnCode(wx.CANCEL)
	self.EndModal(wx.CANCEL)

# ===========================================================================================================
# DOWNLOAD MANAGER CLASS
# This is the dialog that is displayed when a structure is loaded, asking you what chains you want to load
# and if waters/HETATMs should be loaded

class DownloadManagerDialog(wx.Dialog):
    def __init__(self, parent, serverURL, scriptdir):
	if (platform.system() != "Linux"):
	    wx.Dialog.__init__(self, parent, -1, "Download Manager", size=(405, 330))
	else:
	    wx.Dialog.__init__(self, parent, -1, "Download Manager", size=(400, 330))
	self.scriptdir = scriptdir

	# Save the parent, to access the parent's data
	self.parent = parent
	# The parent will let us know what the saved URL was for the default setting
	self.inURL = serverURL

	# Label, text box, and test button for specifying a remote server URL
	if (platform.system() == "Darwin"):
	    self.lblServerName = wx.StaticBitmap(self, -1, wx.Image(self.scriptdir + "/images/osx/sequence/lblServerName.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(5, 7), size=(395, 23))
	else:
	    self.lblServerName = wx.StaticText(self, -1, "Server Name", (5, 7), (395, 23))
	    self.lblServerName.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
	self.txtServerName = wx.TextCtrl(self, -1, pos=(5, 30), size=(310, 25))
	self.txtServerName.SetValue(serverURL)
	self.txtServerName.SetToolTipString("Enter the HTTP address of the remote server")
	if (platform.system() == "Darwin"):
	    self.btnTest = wx.BitmapButton(self, -1, wx.Image(self.scriptdir + "/images/osx/sequence/btnTest.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), (320, 30), (70, 25))
	else:
	    self.btnTest = wx.Button(self, id=-1, label="Test", pos=(320, 30), size=(70, 25))
	    self.btnTest.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.btnTest.Bind(wx.EVT_BUTTON, self.testConnection)
	self.btnTest.SetToolTipString("Test that InteractiveROSETTA can connect to this server")

	# Set up a scrollable area for displaying all of the jobs that have been sent to the server
	if (platform.system() == "Darwin"):
	    self.lblActive = wx.StaticBitmap(self, -1, wx.Image(self.scriptdir + "/images/osx/sequence/lblActive.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(5, 60), size=(395, 30))
	else:
	    self.lblActive = wx.StaticText(self, -1, "Active Submissions:", (5, 60), (395, 30))
	    self.lblActive.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))

	self.scroll = wx.ScrolledWindow(self, -1, pos=(0, 80))
	self.scroll.SetSize((self.GetSize()[0], 180))

	# Read the downloadwatch file to see the active jobs
	try:
	    jobs = 0
	    self.lblDownloads = []
	    self.btnRemoveJob = []
	    self.realJobIDs = []
	    goToSandbox()
	    fin = open("downloadwatch", "r")
	    for aline in fin:
		if (len(aline.strip()) == 0):
		    continue
		# This is the format for a job submission:
		# jobType <tab> jobID <tab> dateUploaded <tab> remoteServer
		jobType = aline.split("\t")[0].strip()
		if (len(aline.split("\t")) == 5):
		    jobID = aline.split("\t")[4].strip()
		else:
		    jobID = aline.split("\t")[1].strip()
		self.realJobIDs.append(aline.split("\t")[1].strip())
		uploadDate = aline.split("\t")[2].strip()
		server = aline.split("\t")[3].strip()
		status = "Pending"
		# Try to read what the status of the job is, default to "Pending"
		# The server should be posting status updates regularly and all the client does is read them
		URL = server + "/results/" + jobID + "/status"
		try:
		    # Is the server this computer?  If so, we're simply opening a local file
		    # Otherwise, we need urllib2
		    if (URL.lower().startswith("localhost")):
			statusfile = open(URL[URL.find(":")+1:])
		    else:
			statusfile = urllib2.urlopen(URL, timeout=1)
		    status = statusfile.readlines()[0]
		    statusfile.close()
		except urllib2.URLError as e:
		    # Most likely the Internet is down or the server is not available
		    status = "Connection Could Not Be Established"
		except urllib2.HTTPError as e:
		    # Established connection but "status" does not exist on the server
		    status = "Pending"
		except:
		    # Other or local file could not be read for localhost
		    status = "Pending"
		# Generate some labels for the job information
		# First label: jobtype and jobID
		# Second label: The server the job was sent to
		# Third label: The date of the submission
		# Fourth label: The status of the job as reported by the server
		self.lblDownloads.append(wx.StaticText(self.scroll, -1, jobType + " - " + jobID, (5, jobs*100+5), (345, 30)))
		self.lblDownloads[jobs*4].SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
		self.lblDownloads.append(wx.StaticText(self.scroll, -1, "Server: " + server, (15, jobs*100+30), (300, 30)))
		self.lblDownloads[jobs*4+1].SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.NORMAL))
		self.lblDownloads.append(wx.StaticText(self.scroll, -1, "Submitted: " + uploadDate, (15, jobs*100+55), (300, 30)))
		self.lblDownloads[jobs*4+2].SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.NORMAL))
		self.lblDownloads.append(wx.StaticText(self.scroll, -1, "Status: " + status, (15, jobs*100+80), (300, 30)))
		self.lblDownloads[jobs*4+3].SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.NORMAL))
		# Button to send a cancel request for this job
		# The name field is important because it contains the index of this specific cancel button
		# so when one is clicked, we can figure out which one was clicked by looking at the event's name
		if (platform.system() == "Darwin"):
		    self.btnRemoveJob.append(wx.BitmapButton(self, -1, wx.Image(self.scriptdir + "/images/osx/sequence/btnRemoveJob.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), (355, jobs*100), (25, 25)))
		else:
		    self.btnRemoveJob.append(wx.Button(self.scroll, id=-1, label="X", pos=(355, jobs*100), size=(25, 25), name=str(jobs)))
		    self.btnRemoveJob[jobs].SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
		self.btnRemoveJob[jobs].Bind(wx.EVT_BUTTON, self.cancelJob)
		self.btnRemoveJob[jobs].SetToolTipString("Remove this job from the download manager")
		jobs += 1
	    fin.close()
	    self.scroll.SetScrollbars(1, 1, self.scroll.GetSize()[0], max(self.scroll.GetSize()[1], self.lblDownloads[len(self.lblDownloads)-1].GetPosition()[1] + 30))
	except:
	    # The file downloadwatch didn't exist or was corrupt, so we are not aware of any outstanding submissions
	    pass
	# OK and Cancel buttons
	if (platform.system() == "Darwin"):
	    self.btnOK = wx.BitmapButton(self, -1, wx.Image(self.scriptdir + "/images/osx/sequence/btnOK_Protein.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), (40, 265), (140, 30))
	else:
	    self.btnOK = wx.Button(self, id=-1, label="OK", pos=(40, 265), size=(140, 30))
	    self.btnOK.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.btnOK.Bind(wx.EVT_BUTTON, self.okDialog)
	self.btnOK.SetToolTipString("Confirm downloader settings")
	if (platform.system() == "Darwin"):
	    self.btnCancel = wx.BitmapButton(self, -1, wx.Image(self.scriptdir + "/images/osx/sequence/btnCancel_Protein.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), (220, 265), (140, 30))
	else:
	    self.btnCancel = wx.Button(self, id=-1, label="Cancel", pos=(220, 265), size=(140, 30))
	    self.btnCancel.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.btnCancel.Bind(wx.EVT_BUTTON, self.cancelDialog)
	self.btnCancel.SetToolTipString("Cancel this operation")

	# Center this window
	self.SetPosition((wx.GetDisplaySize()[0]/2-150, wx.GetDisplaySize()[1]/2-150))

    def cancelJob(self, event):
	# This function attempts to tell the server to kill this job because the user is no longer interested
	# The server should receive this request and terminate the job, to free up resources for other jobs
	# First figure out which job is being terminated by looking at event's name
	indx = int(event.GetEventObject().GetName())
	if (self.lblDownloads[indx*4+3].GetLabel() == "Status: REMOVED"):
	    # If we already sent the kill request and haven't close this dialog, don't try to send
	    # another kill request
	    return
	# Get the jobID and the server name
	jobID = self.lblDownloads[indx*4].GetLabel().split()[2].strip()
	server = self.lblDownloads[indx*4+1].GetLabel().split()[1].strip()
	# Are you sure, user?  This is irreversible
	dlg = wx.MessageDialog(self, "Are you sure you want to remove job " + jobID + "?  This operation cannot be undone.", "Confirm Job Removal", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
	if (dlg.ShowModal() == wx.ID_YES):
	    # Try to send a kill request to the server
	    try:
		# Usually the character - delimits fields, but the server name might actually have - in it
		# so use | instead
		# sendToServer knows how to deal specifically with kill requests
		sendToServer("kill|" + self.realJobIDs[indx] + "|" + server)
	    except:
		# Don't remove it from our list if the kill request failed
		dlg2 = wx.MessageDialog(self, "InteractiveROSETTA was unable to communicate with the server to remove this job!", "Job Not Terminated", wx.OK | wx.ICON_ERROR | wx.CENTRE)
		dlg2.ShowModal()
		dlg2.Destroy()
		return
	    # Open the downloadwatch file and remove this ID from the watch list
	    goToSandbox()
	    fin = open("downloadwatch", "r")
	    data = []
	    for aline in fin:
		if (len(aline.strip()) == 0):
		    continue
		if (aline.split("\t")[1].strip() == self.realJobIDs[indx]):
		    continue
		data.append(aline)
	    fin.close()
	    fout = open("downloadwatch", "w")
	    for aline in data:
		fout.write(aline)
	    fout.close()
	    self.lblDownloads[indx*4+3].SetLabel("Status: REMOVED")
	dlg.Destroy()

    def testConnection(self, event):
	# Simply send a test request to the server
	# The server receives it and tells us that it was successful, or we catch an error
	setServerName(self.txtServerName.GetValue().strip())
	# Send test input
	try:
	    sendToServer("testinput", remoteServer=self.txtServerName.GetValue().strip())
	    dlg2 = wx.MessageDialog(self, "Server test succeeded!", "Server Successful", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
	    dlg2.ShowModal()
	    dlg2.Destroy()
	except:
	    dlg2 = wx.MessageDialog(self, "Server test failed!  Either the server is not set up to run Rosetta, the server is behind a firewall, or you do not have a network connection.", "Server Failed", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
	    dlg2.ShowModal()
	    dlg2.Destroy()

    def okDialog(self, event):
	# All OK does over Cancel is save the indicated server name
	self.SetReturnCode(wx.OK)
	self.EndModal(wx.OK)

    def cancelDialog(self, event):
	self.SetReturnCode(wx.CANCEL)
	self.EndModal(wx.CANCEL)

# ===========================================================================================================
# CHAIN COLORING CLASSES
# These classes are needed to display the colors of the chains in the SequenceViewer

class SequenceViewer(wx.grid.Grid, glr.GridWithLabelRenderersMixin):
    # This is a special grid that initializes with the custom LabelRenderers to draw the color strips for
    # different chains
    def __init__(self, *args, **kw):
	wx.grid.Grid.__init__(self, *args, **kw)
	glr.GridWithLabelRenderersMixin.__init__(self)

# This is a custom GridLabelRenderer to draw the color strips for chains in the labels
class ChainColorRenderer(glr.GridLabelRenderer):
    def __init__(self, bgcolor):
	# BG color is taken from a list of colors for different chain indices
	self._bgcolor = bgcolor

    def Draw(self, grid, dc, rect, row):
	# This code draws a thin color strip of color _bgcolor in the row labels so the user can see which
	# colors get assigned to which chains
	dc.SetBrush(wx.Brush(self._bgcolor))
	dc.SetPen(wx.TRANSPARENT_PEN)
	saveit = rect[2]
	rect[2] = 6
	dc.DrawRectangleRect(rect)
	rect[2] = saveit
	hAlign, vAlign = grid.GetRowLabelAlignment()
	text = grid.GetRowLabelValue(row)
	self.DrawBorder(grid, dc, rect)
	self.DrawText(grid, dc, rect, text, hAlign, vAlign)

# ===========================================================================================================
# SEQUENCE WINDOW CLASS
# This is the main window frame

class SequenceWin(wx.Frame):
    def __init__(self, W, H, cwd, frozen, poses, sequences, IDs, scriptdir):
	# These are the standard window sizes and positions
	self.stdwinx = 370; self.stdwiny = H-270
	self.stdwinw = W-370; self.stdwinh = 270
	self.screenH = H; self.screenW = W
	winx = self.stdwinx; winy = self.stdwiny
	winw = self.stdwinw; winh = self.stdwinh

	homedir = os.path.expanduser("~")
	# The location of InteractiveROSETTA.py and the last directory InteractiveROSETTA was working in
	self.scriptdir = scriptdir
	self.cwd = cwd
	# Try to get the save values from the cfg file
	try:
	    if (platform.system() == "Windows"):
		f = open(homedir + "\\InteractiveROSETTA\\seqwindow.cfg", "r")
	    else:
		f = open(homedir + "/.InteractiveROSETTA/seqwindow.cfg", "r")
	    for aline in f:
		# We use offsets for the window sizes so that we can attempt to scale the window when the
		# screen resolution changes
		# We also save InteractiveROSETTA's last working directory and the name of the remote server
		if (aline.find("[OFFSET X]") >= 0):
		    winx = winx + int(aline.split()[len(aline.split())-1])
		elif (aline.find("[OFFSET Y]") >= 0):
		    winy = winy + int(aline.split()[len(aline.split())-1])
		elif (aline.find("[OFFSET WIDTH]") >= 0):
		    winw = winw + int(aline.split()[len(aline.split())-1])
		elif (aline.find("[OFFSET HEIGHT]") >= 0):
		    winh = winh + int(aline.split()[len(aline.split())-1])
		elif (aline.find("[CWD]") >= 0):
		    self.cwd = aline.split("\t")[1].strip()
		elif (aline.find("[SERVER]") >= 0):
		    setServerName(aline.split("\t")[1].strip())
	    f.close()
	except:
	    pass
	# Do not allow this window to go off the screen
	if (winx > self.screenW - 100):
	    winx = self.stdwinx
	if (winy > self.screenH - 100):
	    winy = self.stdwiny
	# Catch bad cached sizes
	if (winw < 200):
	    winw = self.stdwinw
	if (winh < 200):
	    winh = self.stdwinh

	wx.Frame.__init__(self, None, 0, "InteractiveROSETTA - Sequence Viewer", size=(winw, winh))
	self.frozen = frozen # This is a legacy Boolean for freezing this window
	self.poses = poses # Contains BioPython structures for loaded PDBs
	self.sequences = sequences # Contains strings of sequences per chain
	self.indxToSeqPos = [] # Contains a mapping of absolute residue index to PDB number index
	self.IDs = IDs # Contains the chainID (our ID: modelname|chainID)
	self.SetPosition((winx, winy))
	self.SetBackgroundColour("#333333")
	self.SetIcon(icon.GetIcon())
	# This window needs a scrollable area in case all the buttons don't fit given the smallest size of
	# the window for certain screen resolutions, so they spill over into a scrolled region
	self.scroll = wx.ScrolledWindow(self, -1)
	self.scroll.SetBackgroundColour("#333333")
	self.scroll.SetSize((winw, winh))

	# Load PDBs button
	if (platform.system() == "Darwin"):
	    self.LoadPDBsBtn = wx.BitmapButton(self.scroll, id=-1, bitmap=wx.Image(self.scriptdir + "/images/osx/sequence/LoadPDBsBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(10, 10), size=(100, 25))
	else:
	    self.LoadPDBsBtn = wx.Button(self.scroll, id=-1, label="Load PDBs", pos=(10, 10), size=(100, 25))
	    self.LoadPDBsBtn.SetForegroundColour("#000000")
	    self.LoadPDBsBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.LoadPDBsBtn.Bind(wx.EVT_BUTTON, self.loadPDBsClick)
	self.LoadPDBsBtn.SetToolTipString("Load new PDBs into PyMOL")

	# Label and text box for searching RCSB for structures
	if (platform.system() == "Darwin"):
	    self.labelRCSB = wx.StaticBitmap(self.scroll, -1, wx.Image(self.scriptdir + "/images/osx/sequence/labelRCSB.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(115, 13), size=(43, 25))
	else:
	    self.labelRCSB = wx.StaticText(self.scroll, -1, "RCSB:", (115, 13), (43, 25))
	    self.labelRCSB.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
	    self.labelRCSB.SetForegroundColour("#FFFFFF")
	self.RCSBTxt = wx.TextCtrl(self.scroll, -1, pos=(158, 10), size=(50, 25))
	self.RCSBTxt.SetValue("")
	self.RCSBTxt.SetToolTipString("Four letter PDB code to search for in the RCSB database")

	# Fetch PDB button to get the structures from RCSB
	if (platform.system() == "Darwin"):
	    self.FetchBtn = wx.BitmapButton(self.scroll, id=-1, bitmap=wx.Image(self.scriptdir + "/images/osx/sequence/FetchPDBBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(210, 10), size=(100, 25))
	else:
	    self.FetchBtn = wx.Button(self.scroll, id=-1, label="Fetch PDB", pos=(210, 10), size=(100, 25))
	    self.FetchBtn.SetForegroundColour("#000000")
	    self.FetchBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.FetchBtn.Bind(wx.EVT_BUTTON, self.fetchClick)
	self.FetchBtn.SetToolTipString("Fetch PDB code from RCSB")

	# Close button to close selected chains, or all chains
	if (platform.system() == "Darwin"):
	    self.CloseBtn = wx.BitmapButton(self.scroll, id=-1, bitmap=wx.Image(self.scriptdir + "/images/osx/sequence/CloseBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(310, 10), size=(100, 25))
	else:
	    self.CloseBtn = wx.Button(self.scroll, id=-1, label="Close", pos=(310, 10), size=(100, 25))
	    self.CloseBtn.SetForegroundColour("#000000")
	    self.CloseBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.CloseBtn.Bind(wx.EVT_BUTTON, self.closeClick)
	self.CloseBtn.SetToolTipString("Close selected/all PDBs")

	# Save button for saving the selected models, or all models
	if (platform.system() == "Darwin"):
	    self.SaveBtn = wx.BitmapButton(self.scroll, id=-1, bitmap=wx.Image(self.scriptdir + "/images/osx/sequence/SaveBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(410, 10), size=(100, 25))
	else:
	    self.SaveBtn = wx.Button(self.scroll, id=-1, label="Save PDB", pos=(410, 10), size=(100, 25))
	    self.SaveBtn.SetForegroundColour("#000000")
	    self.SaveBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.SaveBtn.Bind(wx.EVT_BUTTON, self.saveClick)
	self.SaveBtn.SetToolTipString("Save selected/all loaded models")

	# Save Image button, for taking a picture of the current PyMOL view
	if (platform.system() == "Darwin"):
	    self.SaveImageBtn = wx.BitmapButton(self.scroll, id=-1, bitmap=wx.Image(self.scriptdir + "/images/osx/sequence/SaveImageBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(510, 10), size=(100, 25))
	else:
	    self.SaveImageBtn = wx.Button(self.scroll, id=-1, label="Save Image", pos=(510, 10), size=(100, 25))
	    self.SaveImageBtn.SetForegroundColour("#000000")
	    self.SaveImageBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.SaveImageBtn.Bind(wx.EVT_BUTTON, self.saveImage)
	self.SaveImageBtn.SetToolTipString("Save the current PyMOL view to a PNG image")

	# Join button, for joining two selected chains into one
	if (platform.system() == "Darwin"):
	    self.JoinBtn = wx.BitmapButton(self.scroll, id=-1, bitmap=wx.Image(self.scriptdir + "/images/osx/sequence/JoinChainsBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(610, 10), size=(100, 25))
	else:
	    self.JoinBtn = wx.Button(self.scroll, id=-1, label="Join Chains", pos=(610, 10), size=(100, 25))
	    self.JoinBtn.SetForegroundColour("#000000")
	    self.JoinBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.JoinBtn.Bind(wx.EVT_BUTTON, self.joinChains)
	self.JoinBtn.SetToolTipString("Join the selected chains into a single chain")

	# Renumber button, for renumbering a chain sequence from one starting with the selected residue
	if (platform.system() == "Darwin"):
	    self.RenumberBtn = wx.BitmapButton(self.scroll, id=-1, bitmap=wx.Image(self.scriptdir + "/images/osx/sequence/RenumberBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(710, 10), size=(100, 25))
	else:
	    self.RenumberBtn = wx.Button(self.scroll, id=-1, label="Renumber", pos=(710, 10), size=(100, 25))
	    self.RenumberBtn.SetForegroundColour("#000000")
	    self.RenumberBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.RenumberBtn.Bind(wx.EVT_BUTTON, self.renumber)
	self.RenumberBtn.SetToolTipString("Renumber the selected chain from 1, using the selected residue as the new start")

	# Server button, for displaying the download manager dialog, above
	if (platform.system() == "Darwin"):
	    self.ServerBtn = wx.BitmapButton(self.scroll, id=-1, bitmap=wx.Image(self.scriptdir + "/images/osx/sequence/ServerBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(max(810, W-445), 10), size=(25, 25))
	else:
	    self.ServerBtn = wx.Button(self.scroll, id=-1, label="S", pos=(max(810, W-445), 10), size=(25, 25))
	    self.ServerBtn.SetForegroundColour("#000000")
	    self.ServerBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.ServerBtn.Bind(wx.EVT_BUTTON, self.configureServer)
	self.ServerBtn.SetToolTipString("Connect to a remote server running PyRosetta and Rosetta (this is required for some protocols)")
	# Help button, for opening the Sequence Window's help page
	if (platform.system() == "Darwin"):
	    self.HelpBtn = wx.BitmapButton(self.scroll, id=-1, bitmap=wx.Image(self.scriptdir + "/images/osx/HelpBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(self.ServerBtn.GetPosition()[0]+25, 10), size=(25, 25))
	else:
	    self.HelpBtn = wx.Button(self.scroll, id=-1, label="?", pos=(self.ServerBtn.GetPosition()[0]+25, 10), size=(25, 25))
	    self.HelpBtn.SetForegroundColour("#0000FF")
	    self.HelpBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.HelpBtn.Bind(wx.EVT_BUTTON, self.showHelp)
	self.HelpBtn.SetToolTipString("Display the help file for this window")

	# This is the grid that contains all of the primary sequences
	# Recall that "SequenceViewer" is a special class that inherits from wx.grid.Grid
	# This should be read-only
	self.SeqViewer = SequenceViewer(self.scroll)
	self.SeqViewer.CreateGrid(0, 0)
	self.SeqViewer.SetSize((W-480, 150))
	self.SeqViewer.SetPosition((10, 50))
	self.SeqViewer.SetLabelFont(wx.Font(8, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.SeqViewer.DisableDragColSize()
	self.SeqViewer.DisableDragRowSize()
	# These values are needed to keep track of which cell we are hovering over when displaying PDB numbering
	self.rowpos = -1
	self.colpos = -1
	# This is a dictionary that saves model+chainID+PDB numbers and their associated row/column in the grid
	self.selectLookup = {}
	# This boolean is set to True when a Rosetta protocol is active and we should not be able to tamper with the sequences
	self.cannotDelete = False
	# This boolean is set to True when a protocol view in PyMOL is active, so the selection
	# panel knows it should be changing representations on the protocol_view only
	self.protocol_view_active = False

	# These are the buttons to the right of the SeqViewer
	xpos = self.SeqViewer.GetPosition()[0] + self.SeqViewer.GetSize()[0] + 5
	# BGRecolor button, for changing the color of PyMOL's background
	self.BGRecolorBtn = wx.BitmapButton(self.scroll, -1, wx.Image(self.scriptdir + "/images/colorwheel.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), (xpos, 50), (70, 35))
	self.BGRecolorBtn.Bind(wx.EVT_BUTTON, self.bgRecolor)
	self.BGRecolorBtn.SetToolTipString("Change PyMOL background color")
	# Stereo button, for switching between stereo modes in PyMOL
	if (platform.system() == "Darwin"):
	    self.StereoBtn = wx.BitmapButton(self.scroll, -1, wx.Image(self.scriptdir + "/images/osx/sequence/StereoBtn_Mono.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), (xpos, 90), (70, 35))
	else:
	    self.StereoBtn = wx.Button(self.scroll, id=-1, label="Mono", pos=(xpos, 90), size=(70, 35))
	    self.StereoBtn.SetForegroundColour("#000000")
	    self.StereoBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.StereoBtn.Bind(wx.EVT_BUTTON, self.stereoToggle)
	self.StereoBtn.SetToolTipString("Toggle stereo view on/off")
	# Coloring button, for changing the colors of cells in the SeqViewer
	if (platform.system() == "Darwin"):
	    self.ColoringBtn = wx.BitmapButton(self.scroll, -1, wx.Image(self.scriptdir + "/images/osx/sequence/ColoringBtn_NoColor.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), (xpos, 130), (70, 35))
	else:
	    self.ColoringBtn = wx.Button(self.scroll, id=-1, label="No Color", pos=(xpos, 130), size=(70, 35))
	    self.ColoringBtn.SetForegroundColour("#000000")
	    self.ColoringBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.ColoringBtn.Bind(wx.EVT_BUTTON, self.recolorClick)
	self.ColoringBtn.SetToolTipString("Color primary sequence by secondary structure/B-factor or turn off coloring")
	# Align button, for changing the numbering of the SeqViewer columns
	if (platform.system() == "Darwin"):
	    self.AlignBtn = wx.BitmapButton(self.scroll, -1, wx.Image(self.scriptdir + "/images/osx/sequence/AlignBtn_From1.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), (xpos, 170), (70, 35))
	else:
	    self.AlignBtn = wx.Button(self.scroll, id=-1, label="From 1", pos=(xpos, 170), size=(70, 35))
	    self.AlignBtn.SetForegroundColour("#000000")
	    self.AlignBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
	self.AlignBtn.Bind(wx.EVT_BUTTON, self.realignToggle)
	self.AlignBtn.SetToolTipString("PDB sequences renumbered from 1")
	self.noMUSCLEWarned = False # If MUSCLE is not installed, you get a warning message once per session
	self.viewMode = "Mono"
	self.colorMode = "No Color"
	self.alignType = "From 1"

	# Timers are used for some updates
	# saveTimer is used to periodically save the window dimensions/positions
	# PyMOLUpdateTimer gets fired when a selection happens so the PyMOL window gets updated
	# DownloadTimer is used for periodically checking for completed remote server jobs
	self.saveTimer = wx.Timer(self)
	self.PyMOLUpdateTimer = wx.Timer(self)
	self.DownloadTimer = wx.Timer(self)
	self.activeJobs = []
	self.Bind(wx.EVT_TIMER, self.saveWindowData, self.saveTimer)
	self.Bind(wx.EVT_TIMER, self.updatePyMOLSelection, self.PyMOLUpdateTimer)
	self.Bind(wx.EVT_TIMER, self.downloader, self.DownloadTimer)
	self.DownloadTimer.Start(60000)

	# We need to be aware of when this window gains or loses focus
	self.Bind(wx.EVT_ACTIVATE, self.focusEvent)
	# And when it changes size/position
	self.Bind(wx.EVT_SIZE, self.windowGeometryChange)
	self.Bind(wx.EVT_MOTION, self.windowGeometryChange)

	# Bind a bunch of key and mouse events to the SeqViewer
	self.Bind(wx.grid.EVT_GRID_LABEL_LEFT_CLICK, self.labelClick)
	self.SeqViewer.GetGridWindow().Bind(wx.EVT_LEFT_UP, self.leftRelease)
	self.SeqViewer.GetGridWindow().Bind(wx.EVT_MOTION, self.onMouseOver)
	self.Bind(wx.grid.EVT_GRID_CELL_LEFT_CLICK, self.leftClick)
	self.Bind(wx.grid.EVT_GRID_CELL_RIGHT_CLICK, self.rightClick)
	self.Bind(wx.grid.EVT_GRID_LABEL_LEFT_DCLICK, self.labelDClick)
	self.Bind(wx.grid.EVT_GRID_LABEL_RIGHT_CLICK, self.labelRClick)
	self.Bind(wx.EVT_CHAR_HOOK, self.keyPress)
	#self.SeqViewer.GetGridWindow().Bind(wx.EVT_KEY_DOWN, self.keyPress)
	self.selectedResidues = []

	# This is the label along the bottom of the window, that displays what Rosetta is currently doing
	self.labelMsg = wx.StaticText(self.scroll, -1, "", (10, 205), (winw-50, 25))
	self.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
	self.labelMsg.SetForegroundColour("#FFFFFF")
	# We need a queue in case there are multiple things going on at once
	self.msgQueue = []

	self.scroll.SetScrollbars(1, 1, max(self.StereoBtn.GetPosition()[0] + self.StereoBtn.GetSize()[0] + 5, self.HelpBtn.GetPosition()[0] + self.HelpBtn.GetSize()[0] + 5), winh-40)
	self.Show()

	# These handlers are needed for loading and writing BioPython structures
	self.pdbreader = Bio.PDB.PDBParser()
	self.pdbwriter = Bio.PDB.PDBIO()

	# This is a dictionary that maps DSSP secondary structure codes to cell bg colors and cell text colors
	self.sscolormap = {}
	self.sscolormap["H"] = ("red", "white")
	self.sscolormap["S"] = ("yellow", "black")
	self.sscolormap["L"] = ("white", "black")
	self.sscolormap[""] = ("white", "black")
	self.sscolormap["T"] = ("blue", "white")
	self.sscolormap["B"] = ("blue", "white")
	self.sscolormap["G"] = ("orange", "black")
	self.sscolormap["I"] = ("orange", "black")

	# Try to use DSSP if it is available, to give more accurate secondary structure coloring
	try:
	    if (platform.system() == "Windows"):
		self.dsspexe = glob.glob(self.scriptdir + "\\bin\\dssp_win*")[0]
	    elif (platform.system() == "Darwin"):
		self.dsspexe = glob.glob(self.scriptdir + "/bin/dssp_darwin*")[0]
	    else:
		self.dsspexe = glob.glob(self.scriptdir + "/bin/dssp_unix*")[0]
	except:
	    self.dsspexe = "N/A"
	    print "Note: DSSP could be used to improve secondary structure predictions."
	    if (platform.system() == "Windows"):
		print "      Install it to " + self.scriptdir + "\\bin\\dssp_win to make it available."
	    else:
		print "      Install it to " + self.scriptdir + "/bin/dssp_unix to make it available."

    # ======================================================================================================
    # Functions for setting external handlers

    def setProtWin(self, protWin):
	self.protWin = protWin

    def setProtocolPanel(self, protPanel):
	self.protPanel = protPanel

    def setPyMOL(self, pymol):
	self.pymol = pymol
	self.cmd = pymol.cmd
	self.stored = pymol.stored

    # ======================================================================================================
    # Functions for window frame events

    # Window gains/or loses focus
    def focusEvent(self, event):
	self.inFocus = event.GetActive()
	if (self.inFocus):
	    try:
		self.cmd.enable("sele")
		self.selectUpdate(updatePyMOL=False)
	    except:
		# Immediately after the window is created this function gets called, yet self.cmd is still not defined
		pass
	else:
	    # When they are in PyMOL, everything has to belong to the "sele" selection
	    return
	    try:
		self.cmd.select("sele", "seqsele or sele")
		self.cmd.delete("seqsele")
		self.cmd.enable("sele")
	    except:
		try:
		    self.cmd.select("sele", "seqsele")
		    self.cmd.delete("seqsele")
		    self.cmd.enable("sele")
		except:
		    self.cmd.enable("sele")
	    return
	    # We're leaving the sequence window so update PyMOL with the current selections in the SeqViewer
	    self.inFocus = True
	    self.selectUpdate(self.inFocus)
	    self.inFocus = False
	    #self.selectTimer.Stop()

    # Window geometry changes
    def windowGeometryChange(self, event):
	# This function starts a timer that will write out the size and position of this window to a cfg file
	# so the orientation is saved and can be loaded the next time InteractiveROSETTA is started
	if (not(self.saveTimer.IsRunning())):
	    self.saveTimer.Start(5000)
	event.Skip()

    # Helper function for saving the window data
    def saveWindowData(self, event):
	self.saveTimer.Stop()
	homedir = os.path.expanduser("~")
	data = []
	# Does the config file exist?  It should, but create it if not
	try:
	    if (platform.system() == "Windows"):
		f = open(homedir + "\\InteractiveROSETTA\\seqwindow.cfg", "r")
	    else:
		f = open(homedir + "/.InteractiveROSETTA/seqwindow.cfg", "r")
	    for aline in f:
		data.append(aline)
	    f.close()
	except:
	    pass
	if (platform.system() == "Windows"):
	    f = open(homedir + "\\InteractiveROSETTA\\seqwindow.cfg", "w")
	else:
	    f = open(homedir + "/.InteractiveROSETTA/seqwindow.cfg", "w")
	itemsFound = [False, False, False, False, False, False, False, False, False] # [offX, offY, offW, offH, offpw, offph cwd, serverName, primaryRender]
	# Get the window position, size, and the PyMOL size (can't seem to get the PyMOL position)
	(x, y) = self.GetPosition()
	(w, h) = self.GetSize()
	(pw, ph) = self.cmd.get_viewport()
	# Overwrite the data that was already there
	for aline in data:
	    if (aline.find("[OFFSET X]") >= 0):
		itemsFound[0] = True
		f.write("[OFFSET X] " + str(x-self.stdwinx) + "\n")
	    elif (aline.find("[OFFSET Y]") >= 0):
		itemsFound[1] = True
		f.write("[OFFSET Y] " + str(y-self.stdwiny) + "\n")
	    elif (aline.find("[OFFSET WIDTH]") >= 0):
		itemsFound[2] = True
		f.write("[OFFSET WIDTH] " + str(w-self.stdwinw) + "\n")
	    elif (aline.find("[OFFSET HEIGHT]") >= 0):
		itemsFound[3] = True
		f.write("[OFFSET HEIGHT] " + str(h-self.stdwinh) + "\n")
	    elif (aline.find("[OFFSET PWIDTH]") >= 0):
		itemsFound[4] = True
		f.write("[OFFSET PWIDTH] " + str(pw-self.stdwinw) + "\n")
	    elif (aline.find("[OFFSET PHEIGHT]") >= 0):
		itemsFound[5] = True
		f.write("[OFFSET PHEIGHT] " + str(ph-(self.screenH-340)) + "\n")
	    elif (aline.find("[CWD]") >= 0):
		itemsFound[6] = True
		f.write("[CWD]\t" + str(self.cwd) + "\n")
	    elif (aline.find("[SERVER]") >= 0):
		itemsFound[7] = True
		f.write("[SERVER]\t" + getServerName() + "\n")
	    elif (aline.find("[RENDER]") >= 0):
		itemsFound[8] = True
		f.write("[RENDER]\t" + getPrimaryRender() + "\n")
	    else:
		f.write(aline)
	# Add new data if it was not overwriting a previous entry
	for i in range(0, len(itemsFound)):
	    if (not(itemsFound[i])):
		if (i == 0):
		    f.write("[OFFSET X] " + str(x-self.stdwinx) + "\n")
		elif (i == 1):
		    f.write("[OFFSET Y] " + str(y-self.stdwiny) + "\n")
		elif (i == 2):
		    f.write("[OFFSET WIDTH] " + str(w-self.stdwinw) + "\n")
		elif (i == 3):
		    f.write("[OFFSET HEIGHT] " + str(h-self.stdwinh) + "\n")
		elif (i == 4):
		    f.write("[OFFSET PWIDTH] " + str(pw-self.stdwinw) + "\n")
		elif (i == 5):
		    f.write("[OFFSET PHEIGHT] " + str(ph-self.stdwinh) + "\n")
		elif (i == 6):
		    f.write("[CWD]\t" + str(self.cwd) + "\n")
		elif (i == 7):
		    f.write("[SERVER]\t" + getServerName() + "\n")
		else:
		    f.write("[RENDER]\t" + getPrimaryRender() + "\n")
	f.close()

    # ======================================================================================================
    # Functions for SeqViewer events

    # Mouse moves around in the SeqViewer
    def onMouseOver(self, event):
	# This function is used to dynamically change what SeqViewer's tooltip is so you can see what
	# the residue type and PDB number are easily
	# This was taken from a StackOverflow post

        # Use CalcUnscrolledPosition() to get the mouse position
        # within the entire grid including what's offscreen
        x, y = self.SeqViewer.CalcUnscrolledPosition(event.GetX(),event.GetY())
        coords = self.SeqViewer.XYToCell(x, y)
        # you only need these if you need the value in the cell
        row = coords[0]
        col = coords[1]
        if (row == self.rowpos and col == self.colpos):
	    event.Skip()
	    return
	else:
	    self.rowpos = row
	    self.colpos = col
	if (row < 0 or col < 0):
	    # Nothing if off the grid
	    event.GetEventObject().SetToolTipString("")
	    event.Skip()
	    return
	try:
	    # Display the residue name and PDB number
	    self.SeqViewer.GetCellValue(row, col)
	    resID = self.indxToSeqPos[row][col]
	    chain = self.IDs[row][len(self.IDs[row])-1]
	    if (chain == "_"):
		chain = " "
	    AA3 = self.poses[self.getPoseIndex(row)][0][chain][resID].resname.strip()
	    event.GetEventObject().SetToolTipString(AA3 + str(self.indxToSeqPos[row][col][1]))
	except:
	    event.GetEventObject().SetToolTipString("")
	event.Skip()

    # SeqViewer receives a keyboard press
    def keyPress(self, event):
	# Delete (not backspace) key pressed, or on Macs you need to do fn+DELETE (Mac DELETE == PC BACKSPACE, not PC DELETE)
	if (int(event.GetKeyCode()) == wx.WXK_DELETE or (platform.system() == "Darwin" and int(event.GetKeyCode()) == 8)):
	    if (self.cannotDelete):
		# Active protocol prohibits this action
		wx.MessageBox("You cannot perform deletions while a protocol is active!", "Cannot Delete During Protocol", wx.OK|wx.ICON_EXCLAMATION)
		return
	    logInfo("Pressed the DELETE key")
	    # Are you sure?
	    msg = "Are you sure you want to delete the selected residues?"
	    dlg = wx.MessageDialog(self, msg, "Delete Residues", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
	    if (dlg.ShowModal() == wx.ID_YES):
		# Get the selected residue ranges
		topLefts = self.SeqViewer.GetSelectionBlockTopLeft()
		bottomRights = self.SeqViewer.GetSelectionBlockBottomRight()
		goToSandbox()
		# Get the (r, c) cells to be deleted
		# This list needs to be sorted such that the columns are always decreasing, otherwise the deletions
		# will screw up what the c value refers to since it will alter the later c values
		# if you delete entries in front of them (i.e. if you delete position 6 before 7,
		# then position 7 becomes 6 and old position 8 will be deleted by accident)
		cells = []
		for i in range(0, len(topLefts)):
		    for r in range(topLefts[i][0], bottomRights[i][0]+1):
			for c in range(topLefts[i][1], bottomRights[i][1]+1):
			    cells.append((r, c))
		# Sort the cells
		for i in range(0, len(cells)-1):
		    highest_c = i
		    for j in range(i+1, len(cells)):
			if (cells[j][1] > cells[highest_c][1]):
			    highest_c = j
			elif (cells[j][1] == cells[highest_c][1]):
			    if (cells[j][0] > cells[highest_c][0]):
				highest_c = j
		    temp = cells[i]
		    cells[i] = cells[highest_c]
		    cells[highest_c] = temp
		# First look to see if whole chains have been selected, because then we should use "deleteChain" on these to make it faster
		chain_deletes = []
		for row in range(0, self.SeqViewer.NumberRows):
		    count = 0
		    for (r, c) in cells:
			if (r == row):
			    count = count + 1
		    # Did we get all the residues in this row?
		    if (count == len(self.sequences[row])):
			chain_deletes.append(row)
			# Take these out of cells so we don't delete them one-by-one
			for i in range(len(cells)-1, -1, -1):
			    if (cells[i][0] == row):
				cells.pop(i)
		# Now delete individual cells
		for (r, c) in cells:
		    if (c >= len(self.indxToSeqPos[r])):
			# This occurs if there are blank cells selected at the end of a shorter sequence
			continue
		    logInfo("Deleting " + self.IDs[r] + " position " + str(c+1))
		    # If this is the last residue in a row, then delete the whole row
		    if (len(self.sequences[r]) == 1):
			self.deleteChain(r)
		    else:
			# Remove this residue from the sequence and pop the indxToSeqPos element
			if (self.alignType == "From 1"): # Absolute numbering used
			    # Shift the sequence after the deletion down by 1
			    self.sequences[r] = self.sequences[r][0:c] + self.sequences[r][c+1:]
			    resID = self.indxToSeqPos[r].pop(c)
			else:
			    # If PDB numbering is used, we don't shift the sequence, but change the residue to -
			    self.sequences[r] = self.sequences[r][0:c] + "-" + self.sequences[r][c+1:]
			    resID = self.indxToSeqPos[r][c]
			    self.indxToSeqPos[r][c] = "-"
			pyMOLPos = resID[1] # PDB number as int
			# Take this residue out of the lookup dictionary
			self.selectLookup.pop((self.IDs[r], int(pyMOLPos)))
			# Find the pose that this residue belongs to
			poseloc = r
			while (not(self.poses[poseloc])):
			    poseloc = poseloc - 1
			# Legacy Rosetta pose code, using BioPython now
			# Find the real index of the residue in this pose by summing sequence lengths
			#indx = 0
			#for i in range(poseloc, r):
			#    indx = indx + len(self.sequences[i])
			#indx = indx + c + 1
			chain = self.IDs[r][len(self.IDs[r])-1]
			if (chain == "_"):
			    chain = " "
			# resID is a tuple that came from the BioPython structure
			# Remove the residue from the structure
			self.poses[poseloc][0][chain].detach_child(resID) # You have to give it the tuple otherwise it will screw up on HETATMs
			# IMPORTANT: You have to replace the model in the sandbox with the new model
			currID = self.getModelForChain(poseloc)
			self.pdbwriter.set_structure(self.poses[poseloc])
			self.pdbwriter.save(currID + ".pdb")
			# Update the SeqViewer to reflect this deletion
			for i in range(c, len(self.sequences[r])-1):
			    self.SeqViewer.SetCellValue(r, i, self.SeqViewer.GetCellValue(r, i+1))
			self.SeqViewer.SetCellValue(r, len(self.sequences[r])-1, "")
			# Remove this residue from PyMOL
			# currID is the model name
			fields = self.IDs[r].split("|")
			currID = ""
			for i in range(0, len(fields)-1):
			    currID = currID + fields[i] + "|"
			currID = currID[0:len(currID)-1]
			chainID = fields[len(fields)-1]
			#if (chainID != "_"):
			#    self.cmd.remove("model " + currID + " and chain " + chainID + " and resi " + str(pyMOLPos))
			#else:
			#    self.cmd.remove("model " + currID + " and resi " + str(pyMOLPos))
		# Now delete whole chains that we may have identified
		self.cmd.remove("byres seqsele")
		for i in range(len(chain_deletes)-1, -1, -1):
		    self.deleteChain(chain_deletes[i])
		#self.recolorResidues()
		self.regenerateLookupTable()
	    else:
		logInfo("Canceled DELETE key press")
	    dlg.Destroy()
	    # If the deletion was in the longest sequence, we can delete the ending columns
	    # We should never have a column in the SeqViewer that has only blank spaces in it
	    maxseqlen = 0
	    for i in range(0, len(self.sequences)):
		if (len(self.sequences[i]) > maxseqlen):
		    maxseqlen = len(self.sequences[i])
	    for c in range(self.SeqViewer.NumberCols, maxseqlen, -1):
		self.SeqViewer.DeleteCols(c-1)
	    # Send an activate signal to the protocols panel
	    # Most protocols have a focus event that updates them with sequence information, and the sequences
	    # just changed so we need an update
	    self.protPanel.activate()
	elif (int(event.GetKeyCode()) == wx.WXK_ESCAPE):
	    # ESCAPE key for unselecting everything (right-clicking is easier though)
	    self.SeqViewer.ClearSelection()
	    self.PyMOLUpdateTimer.Stop()
	    self.PyMOLUpdateTimer.Start(100)
	elif ((platform.system() == "Linux" and event.ControlDown() and int(event.GetKeyCode()) == 3) or (event.ControlDown() and int(event.GetKeyCode()) == ord("C"))):
	    # Copy the current selection to the clipboard as FASTA data
	    copystr = ""
	    amISelected = []
	    for r in range(0, len(self.sequences)):
		amISelected.append([])
		for c in range(0, len(self.sequences[r])):
		    amISelected[r].append(0)
	    # Find the selections
	    topLefts = self.SeqViewer.GetSelectionBlockTopLeft()
	    bottomRights = self.SeqViewer.GetSelectionBlockBottomRight()
	    for i in range(0, len(topLefts)):
		for r in range(topLefts[i][0], bottomRights[i][0]+1):
		    for c in range(topLefts[i][1], bottomRights[i][1]+1):
			if (c < len(amISelected[r])):
			    amISelected[r][c] = 1
	    # Do the copying
	    for r in range(0, len(self.sequences)):
		if (sum(amISelected[r]) > 0):
		    copystr += "\n\n> " + self.IDs[r] + "\n"
		    blocklen = 0
		    for c in range(0, len(self.sequences[r])):
			if (amISelected[r][c]):
			    copystr += self.sequences[r][c]
			    blocklen += 1
			    if (blocklen >= 50):
				blocklen = 0
				copystr += "\n"
	    copystr = copystr.strip() + "\n"
	    try:
		if (platform.system() == "Windows"):
		    r = Tk()
		    r.withdraw()
		    r.clipboard_clear()
		    r.clipboard_append(copystr)
		    r.destroy()
		elif (platform.system() == "Darwin"):
		    pyperclip.copy(copystr)
		else:
		    clipboard.copy(copystr)
	    except:
		pass
	event.Skip()

    # Left mouse click on SeqViewer cell
    def leftClick(self, event):
	# This only is used on Windows
	# For some reason CTRL selections on Windows highlight the boxes but don't register as selections so
	# I have to do a selection from the code
	# If CTRL is not held down then the event should be skipped so weird behavior doesn't ensue for normal
	# selections and SHIFT selections
	c = event.GetCol()
	r = event.GetRow()
	# This information is saved for when the mouse is released
	# If the release happens in the same cell, then only this cell is selected
	# Without doing this, you have to click and then drag slightly within the same cell to select only
	# one, which is annoying
	self.clickc = c
	self.clickr = r
	if (platform.system() != "Windows"):
	    # Mac and Windows don't have this CTRL problem, stop here for them
	    event.Skip()
	    return
	if (event.ControlDown()):
	    # Find out if this cell is selected already or not
	    selected = False
	    topLefts = self.SeqViewer.GetSelectionBlockTopLeft()
	    bottomRights = self.SeqViewer.GetSelectionBlockBottomRight()
	    for i in range(0, len(topLefts)):
		if (r >= topLefts[i][0] and r <= bottomRights[i][0] and c >= topLefts[i][1] and c <= bottomRights[i][1]):
		    selected = True
		    break
	    # Toggle selection
	    if (not(selected)):
		self.SeqViewer.SelectBlock(r, c, r, c, True)
	    else:
		self.SeqViewer.DeselectCell(r, c)
	else:
	    # Default behavior
	    event.Skip()

    # Left mouse release on SeqViewer cell
    def leftRelease(self, event):
	# User released the mouse on the grid, probably because they were selecting things
	# Update PyMOL
	# A timer is needed because the selection data is not available at the time this function executes
	# If this was a click and release on the same element, set that element to be the only thing selected
	# This is needed because otherwise you have to click and hold and drag a little to only select one box
	try:
	    if (self.clickc == self.colpos and self.clickr == self.rowpos):
		if (not(event.ControlDown())):
		    self.SeqViewer.ClearSelection()
		self.SeqViewer.SelectBlock(self.clickr, self.clickc, self.clickr, self.clickc, True)
	    self.PyMOLUpdateTimer.Stop()
	    self.PyMOLUpdateTimer.Start(100)
	except:
	    # This first click was on the grid, but not on a cell
	    pass
	event.Skip()

    # Right mouse click on SeqViewer cell
    def rightClick(self, event):
	# Use this to unselect everything because now clicking a single element selects that one element
	# so we need a way to unselect everything without having to do it in PyMOL
	self.SeqViewer.ClearSelection()
	self.PyMOLUpdateTimer.Stop()
	self.PyMOLUpdateTimer.Start(100)
	event.Skip()

    # Single left mouse click on a SeqViewer label
    def labelClick(self, event):
	# This is a dummy function so single clicks on the labels do nothing because I am not skipping the event
	pass

    # Double left mouse click on a SeqViewer label
    def labelDClick(self, event):
	# For some odd reason, normal label clicks, although they look like they select the row/column, do not
	# actually select the row/column in the code so the Timer never picks up the selection
	# Here is an attempt to do it manually
	r = event.GetRow()
	c = event.GetCol()
	if (r >= 0):
	    for c in range(0, self.SeqViewer.NumberCols):
		self.SeqViewer.SelectBlock(r, c, r, c, True)
	elif (c >= 0):
	    for r in range(0, self.SeqViewer.NumberRows):
		self.SeqViewer.SelectBlock(r, c, r, c, True)
	self.selectUpdate(updatePyMOL=True)
	# Do not skip the event, because I don't want anything strange happening

    # Single right mouse click on a SeqViewer label
    def labelRClick(self, event):
	r = event.GetRow()
	ID = self.IDs[r]
	model = ID[0:len(ID)-2]
	chain = ID[len(ID)-1]
	if (chain == "_"):
	    chain = " "
	dlg = ModelDialog(self, model, chain, self.scriptdir)
	if (dlg.ShowModal() == wx.OK):
	    # Now we might need to shuffle some things around
	    # First, model name is new
	    if (model != dlg.return_model):
		# Is this model name already taken?  If so, we need to insert this chain into a PDB file
		modelfound = False
		for i in range(0, len(self.IDs)):
		    if (self.IDs[i][0:len(self.IDs[i])-2] == dlg.return_model):
			# Yes, it is taken already
			# Were we given a chainID?  If not, abort because there needs to be a chain ID in a
			# model with more than one chain
			if (len(dlg.return_chain.strip()) == 0 or dlg.return_chain == "_"):
			    dlg2 = wx.MessageDialog(self, "You are attempting to move a chain to a pre-existing model but have not specified a chain ID.  Please choose a unique chain ID for this chain.", "Unique Chain ID Required", wx.OK | wx.ICON_ERROR | wx.CENTRE)
			    dlg2.ShowModal()
			    dlg2.Destroy()
			    return
			# Find out where the last chain of this model is at
			for j in range(i+1, len(self.IDs)+1):
			    if (j >= len(self.IDs)):
				# You need to do this to get j to go past the end of the loop
				# This happens automatically in other languages but not Python (it will be
				# len-1 instead of len like I want)
				break
			    if (self.IDs[j][0:len(self.IDs[j])-2] != dlg.return_model):
				break
			# j is where we want to insert
			# First, let's see if the renamed chain contains the BioPython structure and if there
			# are many chains on it, because then the pose needs to be moved
			# poseindx = location in self.poses that contains the BioPython structure (either r
			# or some number < r)
			if (not(self.poses[r] is None) and r+1 < len(self.IDs) and self.IDs[r+1][0:len(self.IDs[r+1])-2] == model):
			    self.poses[r+1] = self.poses[r]
			    poseindx = self.getPoseIndex(r+1)
			    c = self.poses[poseindx][0][chain]
			    self.poses[poseindx][0].detach_child(chain)
			else:
			    poseindx = self.getPoseIndex(r)
			    c = self.poses[poseindx][0][chain]
			    self.poses[poseindx][0].detach_child(chain)
			# Pop off the current positions
			if (j > r):
			    j -= 1
			# Take everything out of the current position
			self.poses.pop(r)
			oldsequence = self.sequences.pop(r)
			oldID = self.IDs.pop(r)
			oldindxToSeqPos = self.indxToSeqPos.pop(r)
			# Push them back in the new desired position
			self.poses.insert(j, None)
			self.sequences.insert(j, oldsequence)
			self.IDs.insert(j, oldID)
			self.indxToSeqPos.insert(j, oldindxToSeqPos)
			# Fix the ID
			self.IDs[j] = dlg.return_model + "|" + dlg.return_chain
			if (chain != dlg.return_chain):
			    # If the chain will also be changed, set the chain ID temporarily to $, b/c
			    # when we first add the chain to the new model there might be a conflict with
			    # a pre-existing ID before we actually change the ID
			    c.id = "$"
			    if (len(chain.strip()) == 0):
				self.cmd.alter("model " + model, "chain=\"$\"")
			    else:
				self.cmd.alter("model " + model + " and chain " + chain, "chain=\"$\"")
			    chain = "$"
			# Now change in PyMOL
			if (len(chain.strip()) == 0):
			    self.cmd.create(dlg.return_model, "(model " + dlg.return_model + ") or (model " + model + ")")
			    self.cmd.remove("model " + model)
			else:
			    self.cmd.create(dlg.return_model, "(model " + dlg.return_model + ") or (model " + model + " and chain " + chain + ")")
			    self.cmd.remove("model " + model + " and chain " + chain)
			# Now take the original chain structure from oldpose and add it to the new parent
			poseindx = self.getPoseIndex(j-1)
			self.poses[poseindx][0].add(c)
			modelfound = True
			# r needs to be updated in case we are also renaming the chain, since it is assumed
			# that r refers to the ID being renamed, which is actually j now
			r = j
			break
		# Okay, it's a new model name
		# Is it a chain in the middle of the model it is currently in?
		if (not(modelfound) and r > 0 and r < len(self.IDs) - 1 and self.IDs[r+1][0:len(self.IDs[r+1])-2] == model and self.IDs[r-1][0:len(self.IDs[r-1])-2] == model):
		    # It is, so we have to move things around
		    # Put it before the current model
		    poseindx = self.getPoseIndex(r)
		    j = poseindx
		    self.poses[poseindx][0].detach_child(chain)
		    self.poses.pop(r)
		    oldsequence = self.sequences.pop(r)
		    oldID = self.IDs.pop(r)
		    oldindxToSeqPos = self.indxToSeqPos.pop(r)
		    self.poses.insert(j, None)
		    self.sequences.insert(j, oldsequence)
		    self.IDs.insert(j, oldID)
		    self.indxToSeqPos.insert(j, oldindxToSeqPos)
		    # Fix the ID
		    self.IDs[j] = dlg.return_model + "|" + dlg.return_chain
		    # Now change in PyMOL
		    if (len(chain.strip()) == 0):
			self.cmd.create(dlg.return_model, "model " + model)
			self.cmd.remove("model " + model)
		    else:
			self.cmd.create(dlg.return_model, "model " + model + " and chain " + chain)
			self.cmd.remove("model " + model + " and chain " + chain)
		    # Get the BioPython structure back
		    self.cmd.save("temp.pdb", "model " + dlg.return_model)
		    self.poses[j] = self.pdbreader.get_structure(dlg.return_model, "temp.pdb")
		    # r needs to be updated in case we are also renaming the chain
		    r = j
		elif (not(modelfound)):
		    # Easy, just rename it and fix PyMOL's IDs
		    self.IDs[r] = dlg.return_model + "|" + dlg.return_chain
		    # Detach that chain from the original pose
		    poseindx = self.getPoseIndex(r)
		    if (poseindx != r):
			self.poses[poseindx][0].detach_child(chain)
		    # Now change in PyMOL
		    if (len(chain.strip()) == 0):
			self.cmd.create(dlg.return_model, "model " + model)
			self.cmd.remove("model " + model)
		    else:
			self.cmd.create(dlg.return_model, "model " + model + " and chain " + chain)
			self.cmd.remove("model " + model + " and chain " + chain)
		    # Get the BioPython structure back
		    self.cmd.save("temp.pdb", "model " + dlg.return_model)
		    self.poses[r] = self.pdbreader.get_structure(dlg.return_model, "temp.pdb")
		model = dlg.return_model
	    if (chain != dlg.return_chain):
		# This is easy, just rename the chain in both the ID list and PyMOL
		self.IDs[r] = self.IDs[r][0:len(self.IDs[r])-1] + dlg.return_chain
		self.SeqViewer.SetRowLabelValue(r, self.IDs[r])
		if (len(chain.strip()) == 0):
		    self.cmd.alter("model " + model, "chain=\"" + dlg.return_chain + "\"")
		else:
		    self.cmd.alter("model " + model + " and chain " + chain, "chain=\"" + dlg.return_chain + "\"")
		# Change the IDs in the BioPython structures
		poseindx = self.getPoseIndex(r)
		c = self.poses[poseindx][0][chain]
		self.poses[poseindx][0].detach_child(chain)
		c.detach_parent()
		c.id = dlg.return_chain
		self.poses[poseindx][0].add(c)
	    for k in range(0, len(self.IDs)):
		self.updateSeqViewer(rowToUpdate=k)
	    self.regenerateLookupTable()

    # Helper function for selecting the indicated residue
    # Used by the comparative modeler
    def selectNthResidue(self, r, indx, addToSelection=False, updateSelection=True):
	count = 0
	c = -1
	for i in range(0, len(self.sequences[r])):
	    if (self.sequences[r][i] != "-"):
		count += 1
	    if (count >= indx):
		c = i
		break
	if (c < 0):
	    return
	self.SeqViewer.SelectBlock(r, c, r, c, addToSelection)
	if (updateSelection):
	    self.selectUpdate(updatePyMOL=True)
	# Possibly move the scroll position to bring this residue into view
	scrollpos = self.SeqViewer.GetScrollPos(wx.HORIZONTAL)
	x = self.SeqViewer.CellToRect(r, c)[0]
	w = self.SeqViewer.GetSize()[0] - self.SeqViewer.GetRowLabelSize() - 20
	scrollunitconv = self.SeqViewer.GetScrollPixelsPerUnit()[0]
	scrollpospx = scrollunitconv * scrollpos
	if (x - scrollpospx >= w):
	    shift = int((float(x) - float(scrollpospx)) / float(w))
	    totalscroll = int(scrollpospx + shift * w)
	    scrollunits = int(totalscroll / scrollunitconv)
	    self.SeqViewer.Scroll(scrollunits, 0)
	elif (x - scrollpospx <= 0):
	    shift = int((float(scrollpospx) - float(x)) / float(w)) + 1
	    totalscroll = int(scrollpospx - shift*w)
	    scrollunits = int(totalscroll / scrollunitconv)
	    self.SeqViewer.Scroll(max(0, scrollunits), 0)

    # Helper function for selecting things in PyMOL using PyMOL's cmd API
    def selectInPyMOL(self, r, c, model, chain, res, resend, first):
	self.selectedResidues.append([r, c])
	# If this the first selection, we have to create "sele"
	if (first):
	    # We have to select by model, chain, and residue index using a bunch of "and" operators
	    if ("ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789".find(chain) < 0):
		self.cmd.select("temp", "model " + str(model) + " and resi " + str(res) + "-" + str(resend))
	    else:
		self.cmd.select("temp", "model " + str(model) + " and chain " + str(chain) + " and resi " + str(res) + "-" + str(resend))
	else:
	    # To get PyMOL to add to the existing selection "sele", we used the "or" operator to union "sele"
	    # with the new atoms
	    if ("ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789".find(chain) < 0):
		self.cmd.select("temp", "temp or (model " + str(model) + " and resi " + str(res) + "-" + str(resend) + ")")
	    else:
		self.cmd.select("temp", "temp or (model " + str(model) + " and chain " + str(chain) + " and resi " + str(res) + "-" + str(resend) + ")")

    # Function for updating selections between PyMOL and SeqViewer
    # If updatePyMOL==True, then PyMOL is being updated with SeqViewer information
    # Otherwise, SeqViewer is being updated with PyMOL information
    def selectUpdate(self, updatePyMOL=False):
	# Find out which residues are selected in PyMOL and highlight the corresponding positions in the viewer
	try:
	    self.stored.selected = []
	except:
	    return
	try:
	    # There's this really annoying feature where if you click in the background in PyMOL, it looks like
	    # it deselected everything, but in reality it only disabled the selection, so when you come back here
	    # it sees there's still a selection and reselects everything again
	    # So first make sure the "seqsele" selection is enabled, if it's not then assume nothing is selected
	    # IMPORTANT NOTE: The "enabled" feature does not seem to be recognized properly if the selection is
	    # "sele", which is the default selection name.  You have to name it something besides "sele" to get
	    # this trick to work, hence I used "seqsele"
	    # It's just annoying because the default selection is "sele" so you have to make sure everything ends up
	    # in "seqsele" and not "sele"
	    # If there are issues, try getting rid of the "-i" option in InteractiveROSETTA.py for pymol initialization
	    # because then that internal side window showing the selections will be visible enabling you to see what
	    # is happening with the selections more easily
	    #
	    # We have to try to get everything into the "seqsele" name, but sometimes if the user clicks on things in PyMOL
	    # they end up in sele, such that there is sele and seqsele
	    # If you don't combine the two, then everything in sele ends up getting dropped
	    bothdefined = True
	    try:
		# This will throw an error if there is nothing in seqsele, so maybe only sele is defined
		self.cmd.select("seqsele", "seqsele")
		try:
		    # This will throw an error if sele is not defined, so this determines whether it is [seqsele or seqsele/sele]
		    self.cmd.select("seqsele", "sele or seqsele")
		    defined = True
		    bothdefined = True
		except:
		    # This means only seqsele was defined
		    pass
	    except:
		try:
		    # Seqsele was not defined, so maybe only sele is
		    self.cmd.select("seqsele", "sele")
		except:
		    # Neither was defined, there's nothing selected at all
		    pass
	    defined = True
	    # At this point, everything should be in seqsele unless nothing is selected and seqsele should be enabled
	    # If it is not enabled, then "seqsele" is still defined but nothing is highlighted, in which case nothing is desired
	    # to be selected so we should remove the seqsele name so iterate_state doesn't see these selections
	    # Delete sele first otherwise it will keep coming back
	    try:
		if (bothdefined):
		    # We have to put sele's contents into seqsele
		    self.cmd.select("seqsele", "sele")
		    self.cmd.enable("seqsele")
		self.cmd.delete("sele")
	    except:
		pass
	    if (not("seqsele" in self.cmd.get_names("selections", enabled_only=1))):
		try:
		    self.cmd.delete("seqsele")
		except:
		    pass
		defined = False
	    if (defined):
		try:
		    self.cmd.iterate_state(1, "seqsele", "stored.selected.append([model, chain, resi, resn])")
		    # If nothing is in stored but yet "seqsele" is still defined, delete it because some of the
		    # selecting options assume that if seqsele is empty it doesn't exist
		    if (len(self.stored.selected) == 0):
			self.cmd.delete("seqsele")
		except:
		    pass
	except:
	    pass
	first = True
	amISelected = []
	for r in range(0, self.SeqViewer.NumberRows):
	    amISelected.append([])
	    for c in range(0, len(self.sequences[r])):
		amISelected[r].append(False)
	# If the SequenceWindow is in focus, then the user is fiddling with selecting residues in this window
	# so save the current selected residues BEFORE we add the ones from PyMOL (because otherwise PyMOL will
	# overwrite the selections in this window)
	if (updatePyMOL):
	    self.cmd.delete("seqsele")
	    self.cmd.delete("temp")
	    topLefts = self.SeqViewer.GetSelectionBlockTopLeft()
	    bottomRights = self.SeqViewer.GetSelectionBlockBottomRight()
	    for i in range(0, len(topLefts)):
		for x in range(topLefts[i][0], bottomRights[i][0]+1):
		    for y in range(topLefts[i][1], bottomRights[i][1]+1):
			# If this is a blank space, skip it
			if (y >= len(self.indxToSeqPos[x])):
			    continue
			model = ""
			# The chainID is appended to the ID via a "|"
			# This code recovers the model name if the PDB happened to have "|" as part of the filename
			fields = self.IDs[x].split("|")
			for field in fields[0:len(fields)-1]:
			    model = model + field + "|"
			model = model[0:len(model)-1] # Trim off the last "|"
			chain = fields[len(fields)-1]
			# Save this as a selected residue
			amISelected[x][y] = True
			#self.selectCell(x, y, model, chain, res, first)
	    for r in range(0, len(amISelected)):
		inblock = False
		for c in range(0, len(amISelected[r])):
		    if (self.indxToSeqPos[r][c] == "-"):
			# Skip these, everything on both sides of the dashes should be continuous
			continue
		    elif (amISelected[r][c] and not(inblock)):
			inblock = True
			start = self.indxToSeqPos[r][c][1]
			lastc = c
		    elif (inblock and int(self.indxToSeqPos[r][c][1]) - int(self.indxToSeqPos[r][lastc][1]) < 0):
			# The residues aren't always arranged sequentially in the PDB ordering, so we have to
			# watch out for this
			model = self.IDs[r][0:len(self.IDs[r])-2]
			chain = self.IDs[r][len(self.IDs[r])-1]
			self.selectInPyMOL(r, c, model, chain, start, str(self.indxToSeqPos[r][c-1][1]), first)
			start = self.indxToSeqPos[r][c][1]
			first = False
		    elif (not(amISelected[r][c]) and inblock):
			inblock = False
			model = self.IDs[r][0:len(self.IDs[r])-2]
			chain = self.IDs[r][len(self.IDs[r])-1]
			self.selectInPyMOL(r, c, model, chain, start, str(self.indxToSeqPos[r][lastc][1]), first)
			first = False
		    elif (inblock and amISelected[r][c]):
			lastc = c
		if (inblock):
		    if (self.sequences[r][len(self.sequences[r])-1] != "-"):
			lastc = len(self.sequences[r])-1
		    model = self.IDs[r][0:len(self.IDs[r])-2]
		    chain = self.IDs[r][len(self.IDs[r])-1]
		    self.selectInPyMOL(r, c, model, chain, start, str(self.indxToSeqPos[r][lastc][1]), first)
		    first = False
	    if (not(first)):
		self.cmd.select("seqsele", "temp")
		self.cmd.delete("temp")
		self.cmd.delete("sele")
		# The newest version of PyMOL (1.7.6) does not update the selection until the
		# PyMOL camera is moved, so the following fixes the problem
		self.cmd.turn("x", "1")
		self.cmd.turn("x", "-1")
	else:
	    # Take the selection data from PyMOL and select those residues in the SequenceWindow
	    for atom in self.stored.selected:
		res = atom[2]
		chain = atom[1]
		if (chain == ""):
		    chain = "_"
		model = atom[0]
		if (model == "minimized_view"):
		    # We're selecting a PyMOL_Mover minimize preview so figure out what the real model is to select
		    # the real residues in the Sequence Window
		    model = self.protPanel.selectedModel
		elif (model == "designed_view"):
		    # We're selecting a PyMOL_Mover fixbb preview so figure out what the real model is to select
		    # the real residues in the Sequence Window
		    model = self.protPanel.selectedModel
		try:
		    (r, c) = self.selectLookup[(model + "|" + chain, int(res))]
		except:
		    # Some unrecognized view from a protocol
		    return
		#c = int(res)
		# The last Boolean value is set to "True" if you are adding to the current selection, "False" starts
		# a new selection, so we only want it to be "False" the first time through
		# There doesn't seem to be a SelectCell function that works, so we used SelectBlock instead
		amISelected[r][c] = True
		#self.SeqViewer.SelectBlock(r, c, r, c, not(first))
		#first = False
	    for r in range(0, len(amISelected)):
		inblock = False
		for c in range(0, len(amISelected[r])):
		    if (amISelected[r][c] and not(inblock)):
			inblock = True
			start = c
		    elif (not(amISelected[r][c]) and inblock):
			inblock = False
			self.SeqViewer.SelectBlock(r, start, r, c-1, not(first))
			first = False
		if (inblock):
		    self.SeqViewer.SelectBlock(r, start, r, len(amISelected[r])-1, not(first))
		    first = False
	    if (first):
		# Nothing selected, deselect everything
		for r in range(0, self.SeqViewer.NumberRows):
		    self.SeqViewer.DeselectRow(r)
	self.logSelection()
	self.cmd.enable("seqsele")

    # This function redraws SeqViewer to have the current sequences and/or IDs
    # If rowToUpdate is >= 0, then we are only updating that row (otherwise we are adding new rows)
    # If relabel==True, then the column labels need to be updated (probably because we are labeling by
    # PDB numbering now)
    def updateSeqViewer(self, rowToUpdate=-1, relabel=False):
	# How many extra rows and columns are we getting after this addition?
	newnrows = len(self.sequences)
	newncols = 0
	for i in range(0, newnrows):
	    if (len(self.sequences[i]) > newncols):
		newncols = len(self.sequences[i])
	# Need to add more columns and relabel the columns
	if (newncols > self.SeqViewer.NumberCols or relabel):
	    if (newncols > self.SeqViewer.NumberCols):
		self.SeqViewer.AppendCols(newncols - self.SeqViewer.NumberCols)
	    # Set the column label to A temporarily to get the autosize function to have a good size
	    for i in range(0, self.SeqViewer.NumberCols):
		self.SeqViewer.SetColLabelValue(i, "A")
	    self.SeqViewer.AutoSizeColumns(True)
	    self.SeqViewer.SetColLabelValue(0, "")
	    for i in range(0, self.SeqViewer.NumberCols):
		if (self.alignType != "PDB #"):
		    # Regular numbering, every 5th column gets labeled
		    # Each column can only hold two digits, so extra digits go in the preceding column label
		    if ((i+1) % 5 == 0):
			if (i+1 < 100):
			    self.SeqViewer.SetColLabelValue(i, str(i+1))
			else:
			    # val is the tens+ones digits
			    # val2 is the thousands+hundreds digits
			    val = (i+1) % 100
			    if (val == 0):
				val = "00"
			    elif (val < 10):
				val = "0" + str(val)
			    else:
				val = str(val)
			    val2 = str((i+1) / 100)
			    self.SeqViewer.SetColLabelValue(i, val)
			    self.SeqViewer.SetColLabelValue(i-1, val2)
		    else:
			self.SeqViewer.SetColLabelValue(i, "")
		else:
		    # We have to read the labels from the list generated by realign()
		    indx = self.columnlabels[i]
		    # Is there a gap between this current column and the one before it?
		    # If so, we need to label this column also so the user knows about the break
		    if (i > 0 and indx - self.columnlabels[i-1] != 1):
			foundbreak = True
		    elif (i == 0 and indx != 1):
			foundbreak = True
		    else:
			foundbreak = False
		    if (indx % 5 == 0 or foundbreak):
			# If the position right after a break would overwrite the entry before it, skip this one
			if (i > 0 and self.SeqViewer.GetColLabelValue(i-1).strip() != "" and indx >= 100):
			    self.SeqViewer.SetColLabelValue(i, "")
			    continue
			if (indx < 100):
			    self.SeqViewer.SetColLabelValue(i, str(indx))
			else:
			    val = indx % 100
			    if (val == 0):
				val = "00"
			    elif (val < 10):
				val = "0" + str(val)
			    else:
				val = str(val)
			    val2 = str(indx / 100)
			    self.SeqViewer.SetColLabelValue(i, val)
			    self.SeqViewer.SetColLabelValue(i-1, val2)
		    else:
			self.SeqViewer.SetColLabelValue(i, "")
	elif (newncols < self.SeqViewer.NumberCols):
	    # Delete any extra columns
	    self.SeqViewer.DeleteCols(newncols, self.SeqViewer.NumberCols - newncols)
	if (newnrows > self.SeqViewer.NumberRows and rowToUpdate < 0):
	    self.SeqViewer.AppendRows(newnrows - self.SeqViewer.NumberRows)
	    self.SeqViewer.SetRowLabelValue(newnrows-1, self.IDs[newnrows-1])
	    # To get the color strip, see above
	    self.SeqViewer.SetRowLabelRenderer(newnrows-1, ChainColorRenderer(getChainColor(newnrows-1).replace("0x", "#")))
	    for i in range(0, len(self.sequences[newnrows-1])):
		self.SeqViewer.SetCellValue(newnrows-1, i, self.sequences[newnrows-1][i])
	elif (rowToUpdate >= 0):
	    # Just updating this one row
	    for i in range(0, self.SeqViewer.NumberCols):
		self.SeqViewer.SetCellValue(rowToUpdate, i, "")
		if (i < len(self.sequences[rowToUpdate])):
		    self.SeqViewer.SetCellValue(rowToUpdate, i, self.sequences[rowToUpdate][i])
	    self.SeqViewer.SetRowLabelValue(rowToUpdate, self.IDs[rowToUpdate])
	readOnly = wx.grid.GridCellAttr()
	readOnly.SetReadOnly(True)
	for i in range(newnrows-1, newnrows):
	    self.SeqViewer.SetRowAttr(i, readOnly)
	# Resize the label width to fit long PDB filenames
	font = self.SeqViewer.GetFont()
	dc = wx.WindowDC(self.SeqViewer)
	dc.SetFont(font)
	maxwidth = 40
	for i in range(0, newnrows):
	    (w, h) = dc.GetTextExtent(self.SeqViewer.GetRowLabelValue(i))
	    if (w > maxwidth):
		maxwidth = w
	self.SeqViewer.DisableDragColSize()
	self.SeqViewer.DisableDragRowSize()
	self.SeqViewer.SetRowLabelSize(maxwidth+10)

    def regenerateLookupTable(self):
	# Whenever residues/chains are deleted, the selectLookup table gets all messed up because the rows and columns
	# don't refer to the same things anymore, so we have to recalculate it after a deletion
	# First do a realignment if necessary to get everything in the right places in the grid
	self.realign()
	self.selectLookup = {}
	for r in range(0, self.SeqViewer.NumberRows):
	    for c in range(0, len(self.indxToSeqPos[r])):
		# Don't add dashes
		if (self.indxToSeqPos[r][c] != "-"):
		    self.selectLookup[(self.IDs[r], int(self.indxToSeqPos[r][c][1]))] = (r, c)
	# Recolor the residues
	self.recolorResidues()

    # ======================================================================================================
    # Functions for button events

    def loadPDBsClick(self, event):
	logInfo("Clicked Load PDB button")
	self.labelMsg.SetLabel("Loading PDB into Rosetta, please be patient...")
	self.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
	self.labelMsg.SetForegroundColour("#FFFFFF")
	self.msgQueue.append("Loading PDB into Rosetta, please be patient...")
	# Select one or more files from the standard load dialog
	# Only PDBs are shown to prevent strange data from being loaded into BioPython
	dlg = wx.FileDialog(
	    self, message="Choose a File",
	    defaultDir=self.cwd,
	    defaultFile="",
	    wildcard="PDB Files (*.pdb)|*.pdb",
	    style=wx.OPEN | wx.CHANGE_DIR | wx.MULTIPLE)
	if (dlg.ShowModal() == wx.ID_OK):
	    if (platform.system() == "Darwin"):
		# Macs apparently cannot load multiple files at once, do this to make it compatible with
		# the existing code
		paths = [dlg.GetPath()]
	    else:
		paths = dlg.GetPaths()
	    # Change cwd to the last opened file
	    if (platform.system() == "Windows"):
		lastDirIndx = paths[len(paths)-1].rfind("\\")
	    else:
		lastDirIndx = paths[len(paths)-1].rfind("/")
	    self.cwd = str(paths[len(paths)-1][0:lastDirIndx])
	    # Save this new cwd
	    self.saveWindowData(None)
	    # Load the PDBs into PyMOL
	    for i in range(0, len(paths)):
		filename = str(paths[i])
		filevalid = True
		try:
		    f = open(filename, "r")
		    f.close()
		except:
		    msg = "The file " + filename.strip() + " is invalid."
		    wx.MessageBox(msg, "File Cannot Be Read", wx.OK|wx.ICON_EXCLAMATION)
		    filevalid = False
		if (filevalid):
		    # Freeze the other windows until this is done
		    self.Disable()
		    self.protWin.Disable()
		    # We have to check to see if there are multiple models in this PDB
		    # If there are, then we're going to break them up into separate PDB files
		    # and read them all as distinct models in PyMOL
		    nummodels = 0
		    f = open(filename, "r")
		    for aline in f:
			if (aline[0:5] == "MODEL"):
			    nummodels = nummodels + 1
		    f.close()
		    if (nummodels <= 1):
			filenames = [filename]
		    else:
			# First tell the user this file has multiple models and ask if they really want to load
			# all of the structures, or only one
			dlg2 = wx.MessageDialog(self, "This PDB contains multiple models.  Do you want to load all models? (Yes loads all models, No loads only the first model)", "Structural Ensemble Detected", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
			if (dlg2.ShowModal() == wx.ID_YES):
			    loadAllModels = True
			else:
			    loadAllModels = False
			dlg2.Destroy()
			goToSandbox()
			filenames = []
			f = open(filename, "r")
			if (platform.system() == "Windows"):
			    indx = filename.rfind("\\")
			else:
			    indx = filename.rfind("/")
			if (indx >= 0):
			    modelbasename = filename[indx+1:].split(".pdb")[0]
			else:
			    modelbasename = filename.split(".pdb")[0]
			# All of the models will be the basename from the PDB, plus an extra integer for
			# its index
			writing = False
			aline = "Begin"
			i = 0
			while (aline):
			    aline = f.readline()
			    if (aline[0:5] == "MODEL"):
				# Start of a new model, start writing a PDB file in the sandbox
				writing = True
				i = i + 1
				f2 = open(modelbasename + "_S" + str(i) + ".pdb", "w")
			    elif (aline[0:6] == "ENDMDL"):
				# Done with this model, close it and save it's filename
				writing = False
				f2.close()
				filenames.append(modelbasename + "_S" + str(i) + ".pdb")
				# Break if only loading one model
				if (not(loadAllModels)):
				    break
			    elif (writing):
				# Just writing data to the current file
				f2.write(aline.strip() + "\n")
			f.close()
		    first = True
		    for filename in filenames:
			if (first):
			    first = False
			    # Show means the protein load option dialog will be shown
			    # The 1 is a dummy value from when this function was called by a thread
			    self.PyMOLPDBLoad(1, filename, "Show")
			else:
			    # We use "Reuse" here so that dialog doesn't get displayed again, instead the
			    # options from the first dialog are applied to all subsequent models
			    self.PyMOLPDBLoad(1, filename, "Reuse")
		    logInfo("Loaded PDB file", filename)
		    # Pop this message out of the queue
		    for i in range(0, len(self.msgQueue)):
			if (self.msgQueue[i].find("Loading PDB") >= 0):
			    self.msgQueue.pop(i)
			    break
		    for i in range(0, len(self.msgQueue)):
			if (self.msgQueue[i].find("Fetching PDB") >= 0):
			    self.msgQueue.pop(i)
			    break
		    if (len(self.msgQueue) > 0):
			self.labelMsg.SetLabel(self.msgQueue[len(self.msgQueue)-1])
		    else:
			self.labelMsg.SetLabel("")
		    self.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
		    self.labelMsg.SetForegroundColour("#FFFFFF")
		    self.Enable()
		    self.protWin.Enable()
		    self.recolorResidues()
		    # Some of the protocols recognize changes in the sequence window, so let the protocol panel know something changed
		    self.protPanel.activate()
	else:
	    # Clean up, no PDB loaded
	    # Pop this message out of the queue
	    for i in range(0, len(self.msgQueue)):
		if (self.msgQueue[i].find("Loading PDB") >= 0):
		    self.msgQueue.pop(i)
		    break
	    if (len(self.msgQueue) > 0):
		self.labelMsg.SetLabel(self.msgQueue[len(self.msgQueue)-1])
	    else:
		self.labelMsg.SetLabel("")
	    self.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
	    self.labelMsg.SetForegroundColour("#FFFFFF")
	    logInfo("Cancelled out of Load PDB")

    def fetchClick(self, event):
	# You have to do this on Linux due to the libc bug
	if (platform.system() != "Windows"):
	    res_init()
	# Attempt to get a PDB directly from RCSB so the user doesn't have to download it manually
	pdbCode = str(self.RCSBTxt.GetValue().strip()[0:4])
	pdbUrl = 'http://www.rcsb.org/pdb/downloadFile.do?fileFormat=pdb&compression=NO&structureId='+pdbCode
	try:
	    logInfo("Attempting to fetch PDB " + pdbCode + " from RCSB")
	    try:
		pdbFile = urllib2.urlopen(pdbUrl, timeout=2)
	    except urllib2.URLError as e:
		# Probably no Internet connection
		raise Exception("Could not connect to RSCB.  Do you have an Internet connection?")
	    goToSandbox()
	    validfile = False
	    # Save the file here
	    f = open(pdbCode + ".pdb", "w")
	    for aline in pdbFile:
		if (aline[0:4] == "ATOM"):
		    # Because a failed result returns HTML code, which will crash Rosetta if you try
		    # to load HTML into a structure
		    validfile = True
		f.write(aline.strip() + "\n")
	    f.close()
	    if (validfile):
		self.labelMsg.SetLabel("Fetching PDB from RCSB, please be patient...")
		self.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
		self.labelMsg.SetForegroundColour("#FFFFFF")
		self.msgQueue.append("Fetching PDB from RCSB, please be patient...")
		self.Disable()
		self.protWin.Disable()
		# Ask if the user would like to save a copy before loading into Rosetta
		dlg = wx.MessageDialog(self, "Save a copy of this PDB before loading into Rosetta?", "Save RCSB PDB", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
		if (dlg.ShowModal() == wx.ID_YES):
		    while (True):
			# Get the filename
			# This is in a loop in case the user tries to overwrite a file and then decides
			# not to overwrite it after all, we need to get back to the original save dialog
			dlg2 = wx.FileDialog(
			    self, message="Save a PDB File",
			    defaultDir=self.cwd,
			    defaultFile=pdbCode,
			    wildcard="PDB Files (*.pdb)|*.pdb",
			    style=wx.SAVE | wx.CHANGE_DIR)
			if (dlg2.ShowModal() == wx.ID_OK):
			    path = dlg2.GetPath()
			    # Change cwd to the last opened file
			    if (platform.system() == "Windows"):
				lastDirIndx = path.rfind("\\")
			    else:
				lastDirIndx = path.rfind("/")
			    self.cwd = str(path[0:lastDirIndx])
			    self.saveWindowData(None)
			    # Load the PDBs into PyMOL
			    filename = str(path).split(".pdb")[0] + ".pdb"
			    # Does it exist already?  If so, ask if the user really wants to overwrite it
			    if (os.path.isfile(filename)):
				dlg3 = wx.MessageDialog(self, "The file " + filename + " already exists.  Overwrite it?", "Filename Already Exists", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
				if (dlg3.ShowModal() == wx.ID_NO):
				    dlg3.Destroy()
				    dlg2.Destroy()
				    logInfo("Cancelled save operation due to filename already existing")
				    continue
				dlg3.Destroy()
			    goToSandbox()
			    f = open(pdbCode + ".pdb", "r")
			    f2 = open(str(filename), "w")
			    for aline in f:
				f2.write(aline.strip() + "\n")
			    f2.close()
			    f.close()
			    break
			else:
			    # User cancelled probably, get out of the loop
			    break
		    dlg2.Destroy()
		dlg.Destroy()
		# We have to check to see if there are multiple models in this PDB
		# If there are, then we're going to break them up into separate PDB files
		# and read them all as distinct models in PyMOL
		nummodels = 0
		filename = pdbCode + ".pdb"
		f = open(filename, "r")
		for aline in f:
		    if (aline[0:5] == "MODEL"):
			nummodels = nummodels + 1
		f.close()
		if (nummodels <= 1):
		    filenames = [filename]
		else:
		    # First tell the user this file has multiple models and ask if they really want to load
		    # all of the structures, or only one
		    dlg = wx.MessageDialog(self, "This PDB contains multiple models.  Do you want to load all models? (Yes loads all models, No loads only the first model)", "Structural Ensemble Detected", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
		    if (dlg.ShowModal() == wx.ID_YES):
			loadAllModels = True
		    else:
			loadAllModels = False
		    dlg.Destroy()
		    goToSandbox()
		    filenames = []
		    f = open(filename, "r")
		    if (platform.system() == "Windows"):
			indx = filename.rfind("\\")
		    else:
			indx = filename.rfind("/")
		    if (indx >= 0):
			modelbasename = filename[indx+1:].split(".pdb")[0]
		    else:
			modelbasename = filename.split(".pdb")[0]
		    writing = False
		    aline = "Begin"
		    i = 0
		    while (aline):
			aline = f.readline()
			if (aline[0:5] == "MODEL"):
			    writing = True
			    i = i + 1
			    f2 = open(modelbasename + "_S" + str(i) + ".pdb", "w")
			elif (aline[0:6] == "ENDMDL"):
			    writing = False
			    f2.close()
			    filenames.append(modelbasename + "_S" + str(i) + ".pdb")
			    # Break if only loading one model
			    if (not(loadAllModels)):
				break
			elif (writing):
			    f2.write(aline.strip() + "\n")
		    f.close()
		first = True
		for filename in filenames:
		    if (first):
			first = False
			# Show the protein load dialog the first time
			self.PyMOLPDBLoad(1, filename, "Show")
		    else:
			# Reuse the options from the first load dialog
			self.PyMOLPDBLoad(1, filename, "Reuse")
		#self.PyMOLPDBLoad(1, pdbCode + ".pdb")
		logInfo("Fetched PDB file", pdbCode + ".pdb")
		# Pop this message out of the queue
		for i in range(0, len(self.msgQueue)):
		    if (self.msgQueue[i].find("Loading PDB") >= 0):
			self.msgQueue.pop(i)
			break
		for i in range(0, len(self.msgQueue)):
		    if (self.msgQueue[i].find("Fetching PDB") >= 0):
			self.msgQueue.pop(i)
			break
		if (len(self.msgQueue) > 0):
		    self.labelMsg.SetLabel(self.msgQueue[len(self.msgQueue)-1])
		else:
		    self.labelMsg.SetLabel("")
		self.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
		self.labelMsg.SetForegroundColour("#FFFFFF")
		self.Enable()
		self.protWin.Enable()
		self.recolorResidues()
		# Some of the protocols recognize changes in the sequence window, so let the protocol panel know something changed
		self.protPanel.activate()
	    else:
		raise Exception("The PDB code " + pdbCode + " could not be found in RCSB!")
	except Exception as e:
	    wx.MessageBox(e.message, "File Cannot Be Read", wx.OK|wx.ICON_EXCLAMATION)
	    try:
		os.remove(pdbCode + ".pdb")
	    except:
		self.Enable()
		self.protWin.Enable()
		# Pop this message out of the queue
		for i in range(0, len(self.msgQueue)):
		    if (self.msgQueue[i].find("Loading PDB") >= 0):
			self.msgQueue.pop(i)
			break
		for i in range(0, len(self.msgQueue)):
		    if (self.msgQueue[i].find("Fetching PDB") >= 0):
			self.msgQueue.pop(i)
			break
		if (len(self.msgQueue) > 0):
		    self.labelMsg.SetLabel(self.msgQueue[len(self.msgQueue)-1])
		else:
		    self.labelMsg.SetLabel("")
		self.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
		self.labelMsg.SetForegroundColour("#FFFFFF")
		pass

    # Helper function for loading a PDBfile into PyMOL
    # dummy is an unused variable for when this function ran in a thread
    # pdbfile is the filename to load
    # showDialog is "Show" if the protein load dialog should be shown
    # showDialog should be "Reuse" for extra models in an ensemble pdbfile
    # FlexiblePeptide is specific for when FlexPepDock attempts to create the peptide
    #    Usually the name "flexpeptide" is reserved but this tells the loader to use the
    #    flexpeptide keyword for this peptide
    def PyMOLPDBLoad(self, dummy, pdbfile, showDialog="NoShow", flexiblePeptide=False):
	pdbfile = str(pdbfile)
	# Check the PDB for duplicate atoms
	# We want to rename dupliates otherwise BioPython will drop them
	#cleanPDB(pdbfile)
	# Set the ID of this PDB as the filename without the .pdb extension
	if (platform.system() == "Windows"):
	    startindx = pdbfile.rfind("\\") + 1
	else:
	    startindx = pdbfile.rfind("/") + 1
	endindx = pdbfile.rfind(".pdb")
	data = []
	# This is also really annoying...PyRosetta apparently cannot handle file paths with whitespace, so we have
	# to create a directory for InteractiveROSETTA in the user's home directory and cd there
	# All of the uploaded files will be stored here, so we can use relative filenames and avoid whitespace,
	# unless there's whitespace in the PDB name itself, in which case the user can be told to rename it
	goToSandbox()
	newID = pdbfile[startindx:endindx]
	# Special characters in the file name will mess up PyMOL, so take them out
	for i in range(0, len(newID)):
	    if (newID[i] == "@"):
		newID = newID[0:i] + "A" + newID[i+1:]
	    elif (newID[i] == "$"):
		newID = newID[0:i] + "S" + newID[i+1:]
	    elif (newID[i] == "&"):
		newID = newID[0:i] + "A" + newID[i+1:]
	    elif (newID[i] == "*"):
		newID = newID[0:i] + "o" + newID[i+1:]
	    elif (newID[i] == "("):
		newID = newID[0:i] + "_" + newID[i+1:]
	    elif (newID[i] == ")"):
		newID = newID[0:i] + "_" + newID[i+1:]
	pdbID = newID
	# The PyMOL ID is normally just the PDB filename without the ".pdb"
	# We have to check to make sure that the same PDB has not been loaded twice because if it is an extra
	# copy, the ID needs to be changed
	taken = False
	for ID in self.IDs:
	    # This is fancy stuff for parsing the actual model name out of the ID, since the ID is the model
	    # + "|" + the chainID
	    fields = ID.split("|")
	    currID = ""
	    for i in range(0, len(fields)-1):
		currID = currID + fields[i] + "|"
	    currID = currID[0:len(currID)-1]
	    if (currID == newID):
		taken = True
		break
	# Sometimes I use special names for selections, and if a modelname is the same as
	# these selections it can screw things up
	if (newID in ["temp", "sele", "seqsele", "params", "designed_view", "minimized_view"]):
	    taken = True
	if (newID == "flexpeptide" and not(flexiblePeptide)):
	    taken = True
	# Replace whitespace with _ to avoid PyMOL issues
	newID = newID.replace(" ", "_")
	newID = newID.replace("\n", "_")
	newID = newID.replace("\t", "_")
	if (taken):
	    for i in range(2, 9999):
		taken = False
		for ID in self.IDs:
		    fields = ID.split("|")
		    '''currID = ""
		    for j in range(0, len(fields)-1):
			currID = currID + fields[j] + "|"
		    currID = currID[0:len(currID)-1'''
		    currID = '|'.join(fields[:len(fields)-1])
		    if (currID == newID + "_" + str(i)):
			taken = True
			break
		if (not(taken)):
		    newID = newID + "_" + str(i)
		    break
	chainIDs = []
	prev_res = "0000"
	f = open(pdbfile, "r")
	for aline in f:
	    data.append(aline)
	f.close
	try:
	    if (platform.system() == "Windows"):
		shutil.copyfile(pdbfile, os.path.expanduser("~") + "\\InteractiveROSETTA\\" + newID + ".pdb")
	    else:
		shutil.copyfile(pdbfile, os.path.expanduser("~") + "/.InteractiveROSETTA/" + newID + ".pdb")
	except:
	    # This happens when you attempt to fetch the PDB and a copy is already downloaded
	    # to the sandbox
	    pass
	#offset = 0
	#for aline in data:
	    #if (aline[0:3] != "TER"):
	    # Take out the unrecognized residues
	#    if ((aline[0:4] == "ATOM" or aline[0:6] == "HETATM") and not(aline[17:20].strip() in getRecognizedTypes())):
	#	offset = offset + 1
	#    elif (aline[0:4] == "ATOM" or aline[0:6] == "HETATM" or aline[0:3] == "TER"):
	#	try:
	#	    atomno = int(aline[7:11])
	#	    atomno = atomno - offset
	#	    aline = aline[0:7] + ("%4i" % atomno) + aline[11:]
	#	except:
	#	    pass
	#	f.write(aline)
	#    else:
	#	f.write(aline)
	#f.close()
	# Check the PDB for duplicate atoms
	# We want to rename dupliates otherwise BioPython will drop them
	# This will only modify the sandbox PDB and not the actual PDB
	# Of course, if the users save their PDBfile it will reflect these changes
	if (platform.system() == "Windows"):
	    pdbdata = cleanPDB(os.path.expanduser("~") + "/InteractiveROSETTA/" + newID + ".pdb")
	else:
	    pdbdata = cleanPDB(os.path.expanduser("~") + "/.InteractiveROSETTA/" + newID + ".pdb")
	# Separate out the chain sequences, but keep the whole pose intact for later Rosetta operations
	try:
	    biopdb = self.pdbreader.get_structure(newID, newID + ".pdb")
	except:
	    # Tell the user the file failed to load and undo the additions to self.indxToSeqPos
	    msg = "There was an error loading " + pdbfile + "\n\n"
	    msg = msg + "Things to check for:\n"
	    msg = msg + "Is the file a valid PDB file?\n"
	    msg = msg + "Are some residues missing a large number of atoms?"
	    dlg = wx.MessageDialog(self, msg, "Close PDBs", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
	    if (dlg.ShowModal() == wx.ID_OK):
		self.labelMsg.SetLabel("")
		self.labelMsg.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
		self.labelMsg.SetForegroundColour("#FFFFFF")
		self.Enable()
		self.protWin.Enable()
		return
	doNotLoad = False
	# Before showing the chain selection dialog, let's first see if there's either multiple chains
	# or there are HETATMs present, otherwise there's no need to show the dialog
	if (len(biopdb[0]) > 1):
	    # Multiple chains, show it
	    pass
	else:
	    # Are there HETATMs?
	    foundHETATMs = False
	    for ch in biopdb[0]:
		for residue in ch:
		    if (not(residue.resname in "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR ADE THY CYT GUA")):
			logInfo('residue not found: %s'%(residue.resname))
			foundHETATMs = True
			break
		if (foundHETATMs):
		    break
	    if (not(foundHETATMs)):
		showDialog = "NoShow"
	if (showDialog != "NoShow"):
	    # Get the chains, then display a dialog that will let the user uncheck chains they do not want to load
	    # and possibly opt to not load waters/hetatms
	    chains = []
	    for ch in biopdb[0]:
		if (ch.id == " "):
		    chains.append("_")
		else:
		    chains.append(ch.id)
	    if (showDialog == "Show"): # Otherwise this is an NMR ensemble and the options from the first model should be applied to all
		dlg = ProteinDialog(self, pdbID + ".pdb", chains, self.scriptdir)
		if (dlg.ShowModal() != wx.OK):
		    # Canceled or closed, do nothing
		    doNotLoad = True
		# Now check the status of those checkboxes
		self.chainsToDelete = []
		for i in range(0, len(chains)):
		    if (not(dlg.chkChains[i].GetValue())):
			self.chainsToDelete.append(chains[i])
		self.keepWaters = dlg.btnWater.GetLabel() == "Load Waters"
		self.keepHETATMs = dlg.btnHETATM.GetLabel() == "Load HETATMs"
		dlg.Destroy()
	    if (len(chains) == len(self.chainsToDelete)):
		# Load nothing
		doNotLoad = True
	    if (not(doNotLoad)):
		# Get rid of the indicated chains
		for chain in self.chainsToDelete:
		    if (chain == "_"):
			biopdb[0].detach_child(" ")
		    else:
			biopdb[0].detach_child(chain)
		# Now take out HETATMs/waters if requested
		toRemove = []
		for ch in biopdb[0]:
		    for residue in ch:
			if (not(self.keepWaters) and residue.resname == "HOH"):
			    toRemove.append((ch.id, residue.id))
			elif (not(self.keepHETATMs) and not(residue.resname in "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR ADE THY CYT GUA")):
			    toRemove.append((ch.id, residue.id))
		for (ch, res) in toRemove:
		    biopdb[0][ch].detach_child(res)
		# Now we have to look for empty chains incase there were HETATM only chains and remove them
		for ch in biopdb[0]:
		    length = 0
		    for residue in ch:
			length = length + 1
		    if (length == 0):
			biopdb[0].detach_child(ch.id)
		# Save it for PyMOL, since atoms and/or chains may have been removed
		self.pdbwriter.set_structure(biopdb)
		self.pdbwriter.save(newID + ".pdb")
		# Take the waters/HETATMs out of the data if requested
		for i in range(len(pdbdata)-1, -1, -1):
		    if (len(pdbdata[i]) >= 20 and pdbdata[i][17:20] == "HOH" and not(self.keepWaters)):
			pdbdata.pop(i)
		    elif (pdbdata[i].startswith("HETATM") and not(self.keepHETATMs)):
			pdbdata.pop(i)
		    elif (pdbdata[i].startswith("HETATM") or pdbdata[i].startswith("ATOM")):
			if (pdbdata[i][21] == " " and "_" in self.chainsToDelete):
			    pdbdata.pop(i)
			elif (pdbdata[i][21] in self.chainsToDelete):
			    pdbdata.pop(i)
	if (not(doNotLoad)):
	    self.poses.append(biopdb)
	    # Use DSSP if available to characterize secondary structure
	    if (self.dsspexe != "N/A"):
		try:
		    dsspdata = Bio.PDB.DSSP(biopdb[0], newID + ".pdb", dssp=self.dsspexe)
		except:
		    self.dsspexe = "N/A"
		    print "WARNING: DSSP could not be run.  Is its binary executable?"
	    first = True
	    i = 0
	    # Things get messed up if loading a PDB when in PDB# or Align mode, so turn these off
	    # Turn off Align for good and require the user to use superimposition again to get it back
	    # For PDB#, turn it back on after loading this whole structure to get the right alignment
	    if (self.alignType == "Align"):
		retoggle = True
	    elif (self.alignType == "PDB #"):
		retoggle = True
		self.alignType = "From 1"
		if (platform.system() == "Darwin"):
		    self.AlignBtn.SetBitmapLabel(bitmap=wx.Image(self.scriptdir + "/images/osx/sequence/AlignBtn_From1.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
		else:
		    self.AlignBtn.SetLabel(self.alignType)
	    else:
		retoggle = False
	    ssdata = {}
	    for c in biopdb[0]:
		if (not(first)):
		    self.poses.append(None) # Just to keep all the lists the same size
		    # Use getPoseIndex to find out where the real structure is given a row
		first = False
		chain = c.id
		if (len(chain.strip()) == 0):
		    chain = "_"
		seq = ""
		for r in c:
		    seq = seq + AA3to1(r.resname)
		self.sequences.append(seq)
		self.IDs.append(newID + "|" + str(chain))
		self.indxToSeqPos.append([])
		ires = 0
		for r in c:
		    self.indxToSeqPos[len(self.indxToSeqPos)-1].append(r.id) # This is actually a 3-tuple, the seqpos is element 1
		    self.selectLookup[(newID + "|" + chain, r.id[1])] = (len(self.sequences)-1, ires)
		    if (self.dsspexe != "N/A"):
			# Alter the states in PyMOL with the DSSP secondary structure predictions
			try:
			    ss = dsspdata[(c.id, r.id[1])][1]
			    # Not all the DSSP codes match up with PyMOL, so make them match
			    if (ss == "E"):
				# Beta sheet
				ss = "S"
			    elif (ss == "S"):
				# Bend
				ss = "B"
			    elif (ss == "-"):
				# Loop
				ss = "L"
			    elif (ss == "T"):
				ss = "B"
			    elif (ss == "I"):
				ss = "G"
			    #if (c.id != " "):
				#self.cmd.alter("model " + newID + " and chain " + c.id + " and resi " + str(r.id[1]), "ss=\"" + ss + "\"")
			    #else:
				#self.cmd.alter("model " + newID + " and resi " + str(r.id[1]), "ss=\"" + ss + "\"")
			except:
			    ss = "L"
			ssdata[c.id + ("%4i" % r.id[1]) + r.id[2]] = ss
		    ires = ires + 1
		i = i + 1
		self.updateSeqViewer()
	    # Rewrite the PDB with B-factors that encode secondary structure coloring
	    fout = open(newID + ".pdb", "w")
	    for aline in pdbdata:
		if (aline.startswith("ATOM") or aline.startswith("HETATM")):
		    ID = aline[21] + aline[22:27]
		    if (aline[12:16] == " CA "):
			if (ssdata[ID] == "H"):
			    fout.write(aline[0:60] + "  30.0" + aline[66:])
			elif (ssdata[ID] == "S"):
			    fout.write(aline[0:60] + "  70.0" + aline[66:])
			else:
			    fout.write(aline[0:60] + "  50.0" + aline[66:])
		    elif (aline[12:16] == " C  "):
			if (ssdata[ID] == "B"):
			    fout.write(aline[0:60] + "  30.0" + aline[66:])
			elif (ssdata[ID] == "G"):
			    fout.write(aline[0:60] + "  70.0" + aline[66:])
			else:
			    fout.write(aline[0:60] + "  50.0" + aline[66:])
		    else:
			fout.write(aline)
		else:
		    fout.write(aline)
	    fout.close()
	    self.cmd.load(newID + ".pdb", newID)
	    self.cmd.alter("byres model " + newID + " and b < 50 and name ca", "ss=\"H\"")
	    self.cmd.alter("byres model " + newID + " and b > 50 and name ca", "ss=\"S\"")
	    self.cmd.alter("byres model " + newID + " and b < 50 and name c", "ss=\"B\"")
	    self.cmd.alter("byres model " + newID + " and b > 50 and name c", "ss=\"G\"")
	    if (retoggle):
		self.realignToggle(None)
	    # Have PyMOL load from this duplicated PDB file so when protocols are run on the PDB we can easily dump
	    # the PDB from Rosetta and then reload it in PyMOL without touching the original PDB unless the user
	    # explicitly tells us to save the new PDB
	    if (self.dsspexe != "N/A"):
		self.cmd.rebuild()
	    defaultPyMOLView(self.cmd, newID)

    def closeClick(self, event):
	# Find out which rows have a selection
	logInfo("Clicked on Close PDB button")
	selectedRows = []
	topLefts = self.SeqViewer.GetSelectionBlockTopLeft()
	bottomRights = self.SeqViewer.GetSelectionBlockBottomRight()
	for i in range(0, len(topLefts)):
	    for j in range(topLefts[i][0], bottomRights[i][0]+1):
		if (not(j in selectedRows)):
		    selectedRows.append(j)
	selectedRows.sort()
	# If no rows were selected, default to closing everything
	if (len(selectedRows) == 0):
	    for i in range(0, self.SeqViewer.NumberRows):
		selectedRows.append(i)
	if (len(selectedRows) == self.SeqViewer.NumberRows):
	    msg = "Are you sure you want to close all chains?"
	else:
	    msg = "Are you sure you want to close the selected chains?"
	dlg = wx.MessageDialog(self, msg, "Close PDBs", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
	if (dlg.ShowModal() == wx.ID_YES):
	    # Go back to renumbering from 1
	    self.alignType = "From 1"
	    if (platform.system() == "Darwin"):
		self.AlignBtn.SetBitmapLabel(bitmap=wx.Image(self.scriptdir + "/images/osx/sequence/AlignBtn_From1.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
	    else:
		self.AlignBtn.SetLabel(self.alignType)
	    self.regenerateLookupTable()
	    # Take off the rows backwards so the renumbering doesn't screw things up
	    while (len(selectedRows) > 0):
		currRow = selectedRows.pop()
		self.deleteChain(currRow)
	    # Shave off any extra columns
	    maxseqlen = 0
	    for i in range(0, len(self.sequences)):
		if (len(self.sequences[i]) > maxseqlen):
		    maxseqlen = len(self.sequences[i])
	    for c in range(self.SeqViewer.NumberCols, maxseqlen, -1):
		self.SeqViewer.DeleteCols(c-1)
	    self.regenerateLookupTable()
	else:
	    logInfo("Cancelled out of Close PDB")
	dlg.Destroy()
	#self.recolorResidues()
	self.SeqViewer.ClearSelection()
	self.protPanel.activate()

    def deleteChain(self, currRow):
	if (self.cannotDelete):
	    wx.MessageBox("You cannot perform deletions while a protocol is active!", "Cannot Delete During Protocol", wx.OK|wx.ICON_EXCLAMATION)
	    return
	logInfo("Deleted row " + self.IDs[currRow] + " in sequence viewer")
	goToSandbox()
	self.SeqViewer.DeleteRows(currRow)
	# Pop the sequences and IDs out that have been deleted
	seq = self.sequences.pop(currRow)
	ID = self.IDs.pop(currRow)
	indexes = self.indxToSeqPos.pop(currRow)
	# Unload the model from PyMOL
	fields = ID.split("|")
	currID = ""
	for i in range(0, len(fields)-1):
	    currID = currID + fields[i] + "|"
	currID = currID[0:len(currID)-1]
	chainID = fields[len(fields)-1]
	# Update the dictionary
	for indx in indexes:
	    if (indx != "-"):
		self.selectLookup.pop((ID, int(indx[1])))
	if (chainID != "_"):
	    self.cmd.remove("model " + currID + " and chain " + chainID)
	else:
	    self.cmd.remove("model " + currID)
	# Update the labels of everything after the deleted row
	for i in range(0, self.SeqViewer.NumberRows):
	    self.SeqViewer.SetRowLabelValue(i, self.IDs[i])
	# Remove the residues from the PyRosetta pose
	# This is a little tricky since the whole pose is only stored in the first chain and there are Nones
	# for the other chains
	# So we have to find this pose, delete all the residues from the deleted chain, and make sure the
	# pose ends up in the right spot in the list
	if (self.poses[currRow] and (currRow == len(self.poses) - 1 or self.poses[currRow+1])):
	    # Easy, only one chain, just delete this one pose and we're done
	    self.poses.pop(currRow) # This is popping a real pose
	    # If this is the last chain, we need to delete the model name from PyMOL, otherwise if you try
	    # to reload this PDB PyMOL will not load it because there will be a namespace clash
	    self.cmd.delete(currID)
	elif (self.poses[currRow]):
	    # We have the structure, but we're deleting the first chain, so we have to shift where the structure
	    # is stored in self.poses
	    s = self.poses.pop(currRow)
	    chain = ID[len(ID)-1]
	    if (chain == "_"):
		chain = " "
	    s[0].detach_child(chain)
	    self.poses[currRow] = s # Place the structure in the spot for the new first chain
	    # IMPORTANT: You have to replace the model in the sandbox with the new model
	    self.pdbwriter.set_structure(s)
	    self.pdbwriter.save(currID + ".pdb")
	else:
	    # We're on a later chain, which means we have to backtrack to find the structure
	    poseloc = currRow - 1
	    while (not(self.poses[poseloc])):
		poseloc = poseloc - 1
	    # Find the true delete window by adding the lengths of the sequences leading up to this chain
	    start = 0
	    for i in range(poseloc, currRow):
		start = start + len(self.sequences[i])
	    end = start + len(seq)
	    chain = ID[len(ID)-1]
	    if (chain == "_"):
		chain = " "
	    self.poses[poseloc][0].detach_child(chain)
	    self.poses.pop(currRow) # This is popping a blank entry, not a real structure
	    # IMPORTANT: You have to replace the model in the sandbox with the new model
	    self.pdbwriter.set_structure(self.poses[poseloc])
	    self.pdbwriter.save(currID + ".pdb")
	self.protWin.Selection.recolorSavedChainColors()

    def saveClick(self, event):
	logInfo("Clicked on Save PDB button")
	models = self.getSelectedModels()
	# Default to save all of them if none selected
	if (len(models) == 0):
	    for r in range(0, self.SeqViewer.NumberRows):
		currID = self.getModelForChain(r)
		if (not(currID in models)):
		    models.append(currID)
	# Find the pose for this model
	for model in models:
	    poseindx = self.getPoseIndexForModel(model)
	    while (True):
		dlg = wx.FileDialog(
		    self, message="Save a PDB File",
		    defaultDir=self.cwd,
		    defaultFile=model,
		    wildcard="PDB Files (*.pdb)|*.pdb",
		    style=wx.SAVE | wx.CHANGE_DIR)
		if (dlg.ShowModal() == wx.ID_OK):
		    path = dlg.GetPath()
		    # Change cwd to the last opened file
		    if (platform.system() == "Windows"):
			lastDirIndx = path.rfind("\\")
		    else:
			lastDirIndx = path.rfind("/")
		    self.cwd = str(path[0:lastDirIndx])
		    self.saveWindowData(None)
		    # Load the PDBs into PyMOL
		    filename = str(path).split(".pdb")[0] + ".pdb"
		    # Does it exist already?  If so, ask if the user really wants to overwrite it
		    if (os.path.isfile(filename)):
			dlg2 = wx.MessageDialog(self, "The file " + filename + " already exists.  Overwrite it?", "Filename Already Exists", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
			if (dlg2.ShowModal() == wx.ID_NO):
			    dlg2.Destroy()
			    logInfo("Cancelled save operation due to filename already existing")
			    continue
			dlg2.Destroy()
		    goToSandbox()
		    # These shenanigans are necessary because Rosetta cannot dump PDBs to filenames with
		    # whitespace in them apparently
		    logInfo("Saved a PDB to " + filename.strip())
		    #self.pdbwriter.set_structure(self.poses[poseindx])
		    #self.pdbwriter.save(filename)
		    # Use PyMOL to save it because the coordinates could have changed using the manipulation panel
		    # and these changes will not be reflected in BioPython
		    self.cmd.save(filename.strip(), "model " + model)
		    # PyMOL can do weird things with adding TER lines around NCAAs, which crashes Rosetta
		    # Take those out
		    fixPyMOLSave(filename.strip())
		    # Update the IDs of this model in both the sequence window and in PyMOL
		    newmodel = str(path[lastDirIndx+1:]).split(".pdb")[0]
		    taken = False
		    for ID in self.IDs:
			if (newmodel == ID.split("|")[0] and newmodel != model):
			    taken = True
			    break
		    # Sometimes I use special names for selections, and if a modelname is the same as
		    # these selections it can screw things up
		    if (newmodel in ["temp", "sele", "seqsele", "params", "designed_view", "minimized_view"]):
			taken = True
		    if (newmodel == "flexpeptide" and not(flexiblePeptide)):
			taken = True
		    # Replace whitespace with _ to avoid PyMOL issues
		    newmodel = newmodel.replace(" ", "_")
		    newmodel = newmodel.replace("\n", "_")
		    newmodel = newmodel.replace("\t", "_")
		    if (taken):
			# Find a new ID that is not taken
			for i in range(2, 1000):
			    for ID in self.IDs:
				if (newmodel + "_" + str(i) == ID.split("|")[0]):
				    continue
			    break
			newmodel = newmodel + "_" + str(i)
		    self.cmd.set_name(model, newmodel)
		    # Now we have to save it in the sandbox otherwise the protocols will be looking for a file with the new
		    # name and will not find it
		    self.cmd.save(newmodel + ".pdb", "model " + newmodel)
		    for i in range(0, len(self.IDs)):
			if (self.IDs[i][0:len(self.IDs[i])-2] == model):
			    self.IDs[i] = newmodel + "|" + self.IDs[i][len(self.IDs[i])-1]
			    self.SeqViewer.SetRowLabelValue(i, self.IDs[i])
		    self.regenerateLookupTable() # Have to do this because the old dictionary is based off the old IDs
		else:
		    logInfo("Cancelled out of Save PDB")
		break
	    dlg.Destroy()

    def saveImage(self, event):
	logInfo("Clicked on Save Image button")
	dlg = wx.FileDialog(
	    self, message="Save an Image",
	    defaultDir=self.cwd,
	    defaultFile="",
	    wildcard="PNG Images (*.png)|*.png",
	    style=wx.SAVE | wx.CHANGE_DIR)
	if (dlg.ShowModal() == wx.ID_OK):
	    if (platform.system() == "Darwin"):
		filename = str(dlg.GetPath())
	    else:
		filename = str(dlg.GetPaths()[0])
	    if (platform.system() == "Windows"):
		# Have to double the number of "\" otherwise the filename can get garbled when
		# PyMOL interprets things like \t as tabs...
		for i in range(len(filename)-1, -1, -1):
		    if (filename[i] == "\\"):
			filename = filename[0:i] + "\\" + filename[i:]
	    self.cmd.disable("seqsele")
	    self.cmd.disable("sele")
	    # Does it exist already?  If so, ask if the user really wants to overwrite it
	    if (os.path.isfile(filename)):
		dlg2 = wx.MessageDialog(self, "The file " + filename + " already exists.  Overwrite it?", "Filename Already Exists", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
		if (dlg2.ShowModal() == wx.ID_NO):
		    dlg2.Destroy()
		    logInfo("Canceled save operation due to filename already existing")
		dlg2.Destroy()
	    self.cmd.png(filename)
	    self.cmd.enable("seqsele")
	    self.cmd.enable("sele")
	    logInfo("Saved an image to " + filename)
	else:
	    logInfo("Cancelled out of Save Image")

    def joinChains(self, event):
	# Do we have more than one chain selected, if not, do nothing
	selectedrows = self.getSelectedChains()
	if (len(selectedrows) <= 1):
	    return
	# Make sure the user really wants to do this
	dlg = wx.MessageDialog(self, "Are you sure you want to join these chains together?", "Confirm Chain Join", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
	if (dlg.ShowModal() == wx.ID_NO):
	    dlg.Destroy()
	    return
	dlg.Destroy()
	# Are multiple models selected?  If so, tell the user and ask if they really want to move a chain to a different model
	if (len(self.getSelectedModels()) > 1):
	    dlg = wx.MessageDialog(self, "You have multiple models selected.  The chain from the lower model will be moved to the higher model.  Is this okay?", "Confirm Chain Join", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
	    if (dlg.ShowModal() == wx.ID_NO):
		dlg.Destroy()
		return
	    dlg.Destroy()
	base_r = selectedrows.pop(0)
	basechain = self.IDs[base_r][len(self.IDs[base_r])-1]
	baseposeindx = self.getPoseIndex(base_r)
	basemodel = self.getModelForChain(base_r)
	# Renumber the base model from 1 to make sure there won't be any numbering conflicts
	ires = 1
	for residue in self.poses[baseposeindx][0][basechain]:
	    thisID = residue.id
	    residue.id = (thisID[0], ires, thisID[2])
	    ires = ires + 1
	nrowsdelete = len(selectedrows)
	for r in selectedrows:
	    chain = self.IDs[r][len(self.IDs[r])-1]
	    poseindx = self.getPoseIndex(r)
	    model = self.getModelForChain(r)
	    for residue in self.poses[poseindx][0][chain]:
		thisID = residue.id
		residue.id = (thisID[0], ires, thisID[2])
		ires = ires + 1
		residue.detach_parent()
		self.poses[baseposeindx][0][basechain].add(residue)
	    self.poses[poseindx][0].detach_child(chain)
	    self.sequences[r] = "RemoveMe"
	    if (r != len(self.sequences) - 1 and not(self.poses[r+1])):
		self.poses[r+1] = self.poses[r]
	    # Take it out of PyMOL
	    if (len(chain.strip()) > 0):
		self.cmd.remove("model " + model + " and chain " + chain)
	    else:
		self.cmd.remove("model " + model)
	# Delete the unoccupied rows
	for i in range(len(selectedrows)-1, -1, -1):
	    r = selectedrows[i]
	    self.SeqViewer.DeleteRows(pos=r)
	    self.sequences.pop(r)
	    self.indxToSeqPos.pop(r)
	    self.poses.pop(r)
	    self.IDs.pop(r)
	# Change back to numbering from one so we can easily get back the proper sequences and indexes
	if (self.alignType != "From 1"):
	    self.realignToggle(None)
	else:
	    # Regenerating the lookup table calls the realigner which will regenerate the sequences properly
	    self.regenerateLookupTable()
	# Now we have to reload the structure in PyMOL
	self.pdbwriter.set_structure(self.poses[baseposeindx])
	self.pdbwriter.save(basemodel + ".pdb")
	self.poses[baseposeindx] = self.pdbreader.get_structure(basemodel, basemodel + ".pdb")
	self.cmd.remove(basemodel)
	self.cmd.delete(basemodel)
	self.cmd.load(basemodel + ".pdb", basemodel)
	# Redo DSSP if available
	if (self.dsspexe != "N/A"):
	    try:
		dsspdata = Bio.PDB.DSSP(self.poses[baseposeindx][0], basemodel + ".pdb", dssp=self.dsspexe)
	    except:
		self.dsspexe = "N/A"
		print "WARNING: DSSP could not be run.  Is its binary executable?"
	    for c in self.poses[baseposeindx][0]:
		for r in c:
		    if (self.dsspexe != "N/A"):
			# Alter the states in PyMOL with the DSSP secondary structure predictions
			try:
			    ss = dsspdata[(c.id, r.id[1])][1]
			    # Not all the DSSP codes match up with PyMOL, so make them match
			    if (ss == "E"):
				# Beta sheet
				ss = "S"
			    elif (ss == "S"):
				# Bend
				ss = "B"
			    elif (ss == "-"):
				# Loop
				ss = "L"
			    if (c.id != " "):
				self.cmd.alter("model " + basemodel + " and chain " + c.id + " and resi " + str(r.id[1]), "ss=\"" + ss + "\"")
			    else:
				self.cmd.alter("model " + basemodel + " and resi " + str(r.id[1]), "ss=\"" + ss + "\"")
			except:
			    ss = "S"
	defaultPyMOLView(self.cmd, basemodel)
	# Clear the selection otherwise deleted rows are still considered selected
	self.SeqViewer.ClearSelection()

    def renumber(self, event):
	# Get the selected elements and make sure there is only one
	topLefts = self.SeqViewer.GetSelectionBlockTopLeft()
	bottomRights = self.SeqViewer.GetSelectionBlockBottomRight()
	cells = []
	for i in range(0, len(topLefts)):
	    for r in range(topLefts[i][0], bottomRights[i][0]+1):
		for c in range(topLefts[i][1], bottomRights[i][1]+1):
		    cells.append((r, c))
	if (len(cells) == 0):
	    return
	elif (len(cells) > 1):
	    dlg = wx.MessageDialog(self, "Please select only one residue to indicate the chain's new starting residue.", "Select Only One Residue", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
	    dlg.ShowModal()
	    dlg.Destroy()
	    return
	(r, c) = cells[0]
	poseindx = self.getPoseIndex(r)
	model = self.getModelForChain(r)
	ires = 1
	chain = self.IDs[r][len(self.IDs[r])-1]
	if (chain == "_"):
	    chain = " "
	leadingres = []
	# Take out the residues leading up to the new start
	for residue in self.poses[poseindx][0][chain]:
	    if (residue.id[1] == self.indxToSeqPos[r][c][1]):
		break
	    leadingres.append(residue)
	for residue in leadingres:
	    self.poses[poseindx][0][chain].detach_child(residue.id)
	# Add back the original leading residues at the end
	for residue in leadingres:
	    #thisID = residue.id
	    #residue.id = (thisID[0], ires, thisID[2])
	    #ires = ires + 1
	    self.poses[poseindx][0][chain].add(residue)
	# Renumber everything still in the chain
	for residue in self.poses[poseindx][0][chain]:
	    thisID = residue.id
	    residue.id = (thisID[0], ires, thisID[2])
	    ires = ires + 1
	# Change back to numbering from one so we can easily get back the proper sequences and indexes
	if (self.alignType != "From 1"):
	    self.realignToggle(None)
	else:
	    # Regenerating the lookup table calls the realigner which will regenerate the sequences properly
	    self.regenerateLookupTable()
	# Now we have to reload the structure in PyMOL
	self.pdbwriter.set_structure(self.poses[poseindx])
	self.pdbwriter.save(model + ".pdb")
	self.poses[poseindx] = self.pdbreader.get_structure(model, model + ".pdb")
	self.cmd.remove(model)
	self.cmd.delete(model)
	self.cmd.load(model + ".pdb", model)
	# Redo DSSP if available
	if (self.dsspexe != "N/A"):
	    try:
		dsspdata = Bio.PDB.DSSP(self.poses[poseindx][0], model + ".pdb", dssp=self.dsspexe)
	    except:
		self.dsspexe = "N/A"
		print "WARNING: DSSP could not be run.  Is its binary executable?"
	    for c in self.poses[poseindx][0]:
		for r in c:
		    if (self.dsspexe != "N/A"):
			# Alter the states in PyMOL with the DSSP secondary structure predictions
			try:
			    ss = dsspdata[(c.id, r.id[1])][1]
			    # Not all the DSSP codes match up with PyMOL, so make them match
			    if (ss == "E"):
				# Beta sheet
				ss = "S"
			    elif (ss == "S"):
				# Bend
				ss = "B"
			    elif (ss == "-"):
				# Loop
				ss = "L"
			    if (c.id != " "):
				self.cmd.alter("model " + model + " and chain " + c.id + " and resi " + str(r.id[1]), "ss=\"" + ss + "\"")
			    else:
				self.cmd.alter("model " + model + " and resi " + str(r.id[1]), "ss=\"" + ss + "\"")
			except:
			    ss = "S"
	defaultPyMOLView(self.cmd, model)
	# Clear the selection otherwise deleted rows are still considered selected
	self.SeqViewer.ClearSelection()

    def configureServer(self, event):
	# This button allows the user to give a name for the remote server
	# After submission, a test will be made and the user will be notified if the test was sucessful
	dlg = DownloadManagerDialog(self, getServerName(), self.scriptdir)
	#dlg = wx.TextEntryDialog(self, "Enter the URL of the remote server:", "Configure Remote Server", "", style=wx.OK | wx.CANCEL)
        #dlg.SetValue(getServerName())
        if (dlg.ShowModal() == wx.OK):
	    setServerName(dlg.txtServerName.GetValue().strip())
	    self.saveWindowData(None)
        dlg.Destroy()

    def showHelp(self, event):
	# Open the help page
	if (platform.system() == "Darwin"):
	    try:
		browser = webbrowser.get("Safari")
	    except:
		print "Could not load Safari!  The help files are located at " + self.scriptdir + "/help"
		return
	    browser.open(self.scriptdir + "/help/sequence.html")
	else:
	    webbrowser.open(self.scriptdir + "/help/sequence.html")

    def bgRecolor(self, event):
	logInfo("Clicked on the background color button")
	dlg = wx.ColourDialog(self)
	dlg.GetColourData().SetChooseFull(True)
	if (dlg.ShowModal() == wx.ID_OK):
	    data = dlg.GetColourData()
	    mycolor = "0x%02x%02x%02x" % data.GetColour().Get()
	    self.cmd.bg_color(mycolor)
	    logInfo("Set background color to " + mycolor)
	else:
	    logInfo("Cancelled out of background color selection")
	dlg.Destroy()

    def stereoToggle(self, event):
	if (self.viewMode == "Mono"):
	    self.viewMode = "Crosseye"
	    if (platform.system() == "Darwin"):
		self.StereoBtn.SetBitmapLabel(bitmap=wx.Image(self.scriptdir + "/images/osx/sequence/StereoBtn_Crosseye.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
	    else:
		self.StereoBtn.SetLabel(self.viewMode)
	    self.cmd.stereo("crosseye")
	    logInfo("Turned stereo crosseye on")
	elif (self.viewMode == "Crosseye"):
	    self.viewMode = "Walleye"
	    if (platform.system() == "Darwin"):
		self.StereoBtn.SetBitmapLabel(bitmap=wx.Image(self.scriptdir + "/images/osx/sequence/StereoBtn_Walleye.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
	    else:
		self.StereoBtn.SetLabel(self.viewMode)
	    self.cmd.stereo("walleye")
	    logInfo("Turned stereo walleye on")
	elif (self.viewMode == "Walleye"):
	    self.viewMode = "Quadbuffer"
	    if (platform.system() == "Darwin"):
		self.StereoBtn.SetBitmapLabel(bitmap=wx.Image(self.scriptdir + "/images/osx/sequence/StereoBtn_Quadbuffer.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
	    else:
		self.StereoBtn.SetLabel(self.viewMode)
	    self.cmd.stereo("quadbuffer")
	    logInfo("Turned stereo quadbuffer on")
	else:
	    self.viewMode = "Mono"
	    if (platform.system() == "Darwin"):
		self.StereoBtn.SetBitmapLabel(bitmap=wx.Image(self.scriptdir + "/images/osx/sequence/StereoBtn_Mono.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
	    else:
		self.StereoBtn.SetLabel(self.viewMode)
	    self.cmd.stereo("off")
	    logInfo("Turned stereo off")

    def recolorClick(self, event):
	if (self.colorMode == "No Color"):
	    self.colorMode = "SS"
	    if (platform.system() == "Darwin"):
		self.ColoringBtn.SetBitmapLabel(bitmap=wx.Image(self.scriptdir + "/images/osx/sequence/ColoringBtn_SS.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
	    else:
		self.ColoringBtn.SetLabel(self.colorMode)
	    logInfo("Colored residues by SS")
	elif (self.colorMode == "SS"):
	    self.colorMode = "BFactor"
	    if (platform.system() == "Darwin"):
		self.ColoringBtn.SetBitmapLabel(bitmap=wx.Image(self.scriptdir + "/images/osx/sequence/ColoringBtn_BFactor.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
	    else:
		self.ColoringBtn.SetLabel(self.colorMode)
	    logInfo("Colored residues by B-factor")
	elif (self.colorMode == "BFactor"):
	    try:
		hmmstr4py
		self.colorMode = "HMMSTR"
		if (platform.system() == "Darwin"):
		    self.ColoringBtn.SetBitmapLabel(bitmap=wx.Image(self.scriptdir + "/images/osx/sequence/ColoringBtn_HMMSTR.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
		else:
		    self.ColoringBtn.SetLabel(self.colorMode)
		logInfo("Colored residues by HMMSTR")
	    except:
		self.colorMode = "No Color"
		if (platform.system() == "Darwin"):
		    self.ColoringBtn.SetBitmapLabel(bitmap=wx.Image(self.scriptdir + "/images/osx/sequence/ColoringBtn_NoColor.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
		else:
		    self.ColoringBtn.SetLabel(self.colorMode)
		logInfo("Turned residue coloring off")
	else:
	    self.colorMode = "No Color"
	    if (platform.system() == "Darwin"):
		self.ColoringBtn.SetBitmapLabel(bitmap=wx.Image(self.scriptdir + "/images/osx/sequence/ColoringBtn_NoColor.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
	    else:
		self.ColoringBtn.SetLabel(self.colorMode)
	    logInfo("Turned residue coloring off")
	self.recolorResidues()

    def recolorResidues(self):
	# This function will recolor the cells in the SeqViewer to either a default white mode,
	# a secondary structure coloring, or a B-factor coloring
	if (self.colorMode == "SS"):
	    # Default everything to white
	    for r in range(0, self.SeqViewer.NumberRows):
		for c in range(0, self.SeqViewer.NumberCols):
		    self.SeqViewer.SetCellBackgroundColour(r, c, "white")
		    self.SeqViewer.SetCellTextColour(r, c, "black")
	    # For each row in the SeqViewer, select the helices and B-sheets, get the selection
	    # in here, and color appropriately
	    for r in range(0, self.SeqViewer.NumberRows):
		model = self.getModelForChain(r)
		chain = self.IDs[r][len(self.IDs[r])-1]
		if (chain != "_" and chain != " " and chain != ""):
		    #self.cmd.select("colorsele", "model " + model + " and chain " + chain)
		    colorsele = "model " + model + " and chain " + chain
		else:
		    #self.cmd.select("colorsele", "model " + model)
		    colorsele = "colorsele", "model " + model
		self.stored.selected = []
		self.cmd.iterate_state(1, colorsele + " and name ca", "stored.selected.append(str(resi)+\"|\"+ss)")
		#self.stored.selected = list(set(self.stored.selected))
		for code in self.stored.selected:
		    resi = int(code.split("|")[0])
		    ss = code.split("|")[1].upper()
		    (color, tcolor) = self.sscolormap[ss]
		    for x in range(0, len(self.indxToSeqPos[r])):
			if (self.indxToSeqPos[r][x] != "-" and self.indxToSeqPos[r][x][1] == int(resi)):
			    c = x
			    break
		    self.SeqViewer.SetCellBackgroundColour(r, c, color)
		    self.SeqViewer.SetCellTextColour(r, c, tcolor)
	    #self.cmd.delete("colorsele")
	elif (self.colorMode == "BFactor"):
	    # A B-Factor of 0 is blue, 100 is red, and everything else is a gradient in between
	    # Default everything to white
	    for r in range(0, self.SeqViewer.NumberRows):
		for c in range(0, self.SeqViewer.NumberCols):
		    self.SeqViewer.SetCellBackgroundColour(r, c, "white")
		    self.SeqViewer.SetCellTextColour(r, c, "black")
	    r = 0
	    while (r < self.SeqViewer.NumberRows):
		if (self.poses[r]):
		    # Skip over the blank entries
		    for ch in self.poses[r][0]:
			c = 0
			for residue in ch:
			    try:
				bfactor = residue["N"].get_bfactor()
			    except:
				# This NCAA does not have an "N" named atom
				continue
			    if (bfactor < 0):
				bfactor = 0.0
			    elif (bfactor > 100):
				bfactor = 100.0
			    red = int(2.55 * bfactor)
			    blue = 255 - red
			    color = "#%02x%02x%02x" % (red, 0, blue)
			    (row, col) = self.selectLookup[(self.IDs[r], residue.id[1])]
			    self.SeqViewer.SetCellBackgroundColour(row, col, color.upper())
			    self.SeqViewer.SetCellTextColour(row, col, "white")
			    c = c + 1
			r = r + 1
		else:
		    r = r + 1
	elif (self.colorMode == "HMMSTR"):
	    ### =================================================================================
	    ### Author: Oluwadamilola Lawal
	    ### PI: Prof. Chris Bystroff
	    ### Rosetta Summer Intern, 2015
	    ### =================================================================================
	    # A HMMSTR of 0 is blue, 1 is red, and everything else is a gradient in between
	    # Default everything to white
	    for r in range(0, self.SeqViewer.NumberRows):
		for c in range(0, self.SeqViewer.NumberCols):
		    self.SeqViewer.SetCellBackgroundColour(r, c, "white")
		    self.SeqViewer.SetCellTextColour(r, c, "black")

	    AA_list2 = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

	    AA_list = ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU',
	'MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']

	    background = [0.0828,0.0194,0.0586,0.0599,0.0401,0.0809,0.0227,0.0555,0.0596,0.0802,0.0207,0.0473,0.0460,0.0373,0.0464,0.0625,0.0589,0.0687,0.0151,0.0376]

	    r = 0
	    while (r < self.SeqViewer.NumberRows):
		if (self.poses[r]):
		    # Skip over the blank entries
		    #Get all non-water residue as indexed in AA_list into a list
		    for ch in self.poses[r][0]:
			seq = []
			for res in ch:
			    if (res.resname != 'HOH') and (res.resname in AA_list):
				seq.append(AA_list.index(res.resname)+1)
			    elif (res.resname != 'HOH'):
				seq.append(21)
			nres = len(seq)
			np.asarray(seq)
			#Get the xyz coordinates of 'N', 'CA', 'C', 'O'
			backbone = np.zeros([3,(4*nres)],dtype=np.float32,order='F')
			backbone[:,:] = 999.0
			n = 0
			for res in ch:
			    if res.resname in AA_list:
				backbone[:,(n)] = res["N"].get_coord()
				backbone[:,(n+1)] = res["CA"].get_coord()
				backbone[:,(n+2)] = res["C"].get_coord()
				backbone[:,(n+3)] = res["O"].get_coord()
				n = n + 4
			    else:
				n = n + 4
			#Send seq and backbone information to HMMSTR, and return gamma sequence profile
			if(platform.system() == "Windows"):
			    gprofile = hmmstr4py.hmmstr.hmmstr_gamma(self.scriptdir + "\\data\\model_R.hmm",seq,backbone)
			else:
			    gprofile = hmmstr4py.hmmstr.hmmstr_gamma(self.scriptdir + "/data/model_R.hmm",seq,backbone)
			#get the corresponding gamma value for each residue in the sequence
			llratio = []
			for pos in xrange(len(seq)):
			    if seq[pos] != 21:
				gp = gprofile[pos,(seq[pos]-1)]
				bg = background[(seq[pos]-1)]
				llrt = np.log(gp/bg)
				if llrt > 3.0:
				    llratio.append(3.0)
				elif llrt < -3.0:
				    llratio.append(-3.0)
				else:
				    llratio.append(llrt)
			    else:
				gamma_seq.append(21)
			#Use HMMSTR gamma value to generate a color
			c = 0
			for residue in ch:
			    if residue.resname != 'HOH':
				try:
				    bfactor = residue["CA"].get_bfactor()
				except:
				    c = c + 1
				    continue

				green = int(128 + (llratio[c] * 42))
				red = 254 - green
				color = "#%02x%02x%02x" % (red, green, 0)
				(row, col) = self.selectLookup[(self.IDs[r], residue.id[1])]
				self.SeqViewer.SetCellBackgroundColour(row, col, color.upper())
				self.SeqViewer.SetCellTextColour(row, col, "white")
				c = c + 1
			r = r + 1
		else:
		    r = r + 1
	else:
	    # Default everything to white
	    for r in range(0, self.SeqViewer.NumberRows):
		for c in range(0, self.SeqViewer.NumberCols):
		    self.SeqViewer.SetCellBackgroundColour(r, c, "white")
		    self.SeqViewer.SetCellTextColour(r, c, "black")
	self.cmd.enable("sele")
	self.SeqViewer.Refresh() # Needed to repaint otherwise the changes don't occur

    def realignToggle(self, event):
	if (self.alignType == "From 1"):
	    self.alignType = "PDB #"
	    if (platform.system() == "Darwin"):
		self.AlignBtn.SetBitmapLabel(bitmap=wx.Image(self.scriptdir + "/images/osx/sequence/AlignBtn_PDB.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
	    else:
		self.AlignBtn.SetLabel(self.alignType)
	    self.regenerateLookupTable()
	#elif (self.alignType == "PDB #"):
	#    self.AlignBtn.SetLabel("Align")
	#    self.regenerateLookupTable()
	else:
	    self.alignType = "From 1"
	    if (platform.system() == "Darwin"):
		self.AlignBtn.SetBitmapLabel(bitmap=wx.Image(self.scriptdir + "/images/osx/sequence/AlignBtn_From1.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
	    else:
		self.AlignBtn.SetLabel(self.alignType)
	    self.regenerateLookupTable()

    def setAlignType(self, val):
	self.alignType = val

    def realign(self):
	if (self.alignType == "Align"):
	    # Use MUSCLE to perform the multiple sequence alignment using BioPython's tools
	    # First, does MUSCLE exist?  If not, tell the user this feature is disabled until it is download
	    try:
		if (platform.system() == "Windows"):
		    muscle = glob.glob(self.scriptdir + "\\bin\\muscle_win*")[0]
		elif (platform.system() == "Darwin"):
		    muscle = glob.glob(self.scriptdir + "/bin/muscle_darwin*")[0]
		else:
		    muscle = glob.glob(self.scriptdir + "/bin/muscle_unix*")[0]
	    except:
		if (not(self.noMUSCLEWarned)):
		    if (platform.system() == "Windows"):
			msg = "MUSCLE installation not found.  If you want to use the alignment feature, then download a copy of MUSCLE to " + self.scriptdir + "\\bin\\muscle_win."
		    else:
			msg = "MUSCLE installation not found.  If you want to use the alignment feature, then download a copy of MUSCLE to " + self.scriptdir + "/bin/muscle_unix."
		    dlg = wx.MessageDialog(self, msg, "MUSCLE Not Found", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
		    dlg.ShowModal()
		    dlg.Destroy()
		    self.noMUSCLEWarned = True
		self.realignToggle(None) # To get back to normal #
		return
	    # Get all the models
	    models = []
	    for r in range(0, len(self.IDs)):
		if (not(self.IDs[r][0:len(self.IDs[r])-2] in models)):
		    models.append(self.IDs[r][0:len(self.IDs[r])-2])
	    # Get all the sequences
	    nseqs = 0
	    if (len(self.IDs) == 0):
		self.realignToggle(None) # To get back to normal #
		return
	    else:
		lastmodel = self.IDs[0][0:len(self.IDs[0])-2]
	    seq = ""
	    seqs = []
	    # Generate a FASTA file with each chain as a different entry
	    # They need to be marked with an index, because MUSCLE can shuffle their order and the order
	    # needs to be reconstructed so the sequences match up in the window
	    # This is not a totally ideal MSA because the chains are all treated separately
	    # I may add a chain joining button at some point to make this easier so I don't have to mess with the MSA
	    f = open("alignin.fasta", "w")
	    indx = 0
	    for r in range(0, len(self.IDs)):
		model =  self.IDs[r][0:len(self.IDs[r])-2]
		if (model != lastmodel):
		    nseqs = nseqs + 1
		    lastmodel = model
		f.write("> " + str(indx) + "\n")
		f.write(self.sequences[r].replace("-", "").strip() + "\n")
		indx = indx + 1
	    f.close()
	    nseqs = nseqs + 1
	    if (nseqs <= 1):
		self.realignToggle(None) # To get back to normal #
		return
	    # Do alignment with MUSCLE
	    muscle_cline = MuscleCommandline(muscle, input="alignin.fasta", out="alignout.fasta")
	    muscle_cline()
	    alignment = range(0, indx)
	    # Read the alignment into a list
	    f = open("alignout.fasta", "r")
	    seq = ""
	    first = True
	    for aline in f:
		if (">" in aline and first):
		    first = False
		    indx = int(aline[aline.index(">")+1:])
		elif (">" in aline):
		    alignment[indx] = seq
		    indx = int(aline[aline.index(">")+1:])
		    seq = ""
		else:
		    seq = seq + aline.strip()
	    alignment[indx] = seq
	    f.close()
	    os.remove("alignin.fasta")
	    os.remove("alignout.fasta")
	    # Now update the data structures
	    thisseq = alignment.pop(0)
	    offset = 0
	    for r in range(0, len(self.sequences)):
		# Find the length of the aligned sequence that corresponds to this chain
		lenchain = 0
		for char in self.sequences[r]:
		    if (char != "-" and char != " "):
			lenchain = lenchain + 1
		thislen = 0
		for i in range(0, len(thisseq)):
		    if (thisseq[i] != "-" and thisseq[i] != " "):
			thislen = thislen + 1
		    if (thislen == lenchain):
			break
		self.sequences[r] = ""
		for j in range(0, offset):
		    self.sequences[r] = self.sequences[r] + "-"
		self.sequences[r] = self.sequences[r] + thisseq[0:i+1]
		charsleft = 0
		for char in thisseq[i+1:].strip():
		    if (char != "-" and len(char.strip()) != 0):
			charsleft = charsleft + 1
		if (charsleft == 0):
		    try:
			thisseq = alignment.pop(0)
			offset = 0
		    except:
			# Last one, the loop will exit after this iteration anyway
			pass
		else:
		    offset = offset + len(self.sequences[r])
		    thisseq = thisseq[i+1:]
	    # Now we have to fix indxToSeqPos
	    # It should already be in PDB order since it was just in PDB # mode, so let's take all the
	    # dashes out, then iterate down each of the sequences and insert dashes when appropriate
	    for r in range(0, len(self.indxToSeqPos)):
		for c in range(len(self.indxToSeqPos[r])-1, -1, -1):
		    if (self.indxToSeqPos[r][c] == "-"):
			self.indxToSeqPos[r].pop(c)
	    for r in range(0, len(self.sequences)):
		for c in range(0, len(self.sequences[r])):
		    if (self.sequences[r][c] == "-"):
			self.indxToSeqPos[r].insert(c, "-")
	    # MUSCLE likes to rename all the HETATMs as "X", so change them back to Os and Zs
	    for r in range(0, len(self.sequences)):
		poseindx = self.getPoseIndex(r)
		for c in range(0, len(self.sequences[r])):
		    if (self.sequences[r][c] == "-"):
			continue
		    resID = self.indxToSeqPos[r][c]
		    chain = self.IDs[r][len(self.IDs[r])-1]
		    if (chain == "_"):
			chain = " "
		    if (self.poses[poseindx][0][chain][resID].resname == "HOH"):
			self.sequences[r] = self.sequences[r][0:c] + "O" + self.sequences[r][c+1:]
		    elif (self.sequences[r][c] == "X"):
			self.sequences[r] = self.sequences[r][0:c] + "Z" + self.sequences[r][c+1:]
	    # Now update the SequenceWindow
	    for r in range(0, len(self.sequences)):
		if (r < len(self.sequences) - 1):
		    self.updateSeqViewer(r, relabel=False)
		else:
		    self.updateSeqViewer(r, relabel=True)
	elif (self.alignType == "PDB #"):
	    # This shouldn't be too hard at first: all you have to do is find the max PDB number for each chain
	    # then create a sequence of the same number of - and then iterate down the residues and replace each
	    # position with the appropriate letter
	    ncols = 0
	    for r in range(0, len(self.sequences)):
		poseindx = self.getPoseIndex(r)
		chain = self.IDs[r][len(self.IDs[r])-1]
		if (chain == "_"):
		    chain = " "
		nres = 0
		for residue in self.poses[poseindx][0][chain]:
		    nres = max(nres, residue.id[1])
		self.sequences[r] = ""
		self.indxToSeqPos[r] = []
		for i in range(0, nres):
		    self.sequences[r] = self.sequences[r] + "-"
		    self.indxToSeqPos[r].append("-")
		for residue in self.poses[poseindx][0][chain]:
		    ires = residue.id[1] - 1
		    self.sequences[r] = self.sequences[r][0:ires] + AA3to1(residue.resname) + self.sequences[r][ires+1:]
		    self.indxToSeqPos[r][ires] = residue.id
		ncols = max(ncols, len(self.sequences[r]))
	    # Here's where things get dicey...sometimes the PDB numbering is such that there are very large ranges
	    # of empty space (like HOHs being at position 5,000+, which generates a ton of columns)
	    # So now go through the sequences and remove columns that are only filled with -
	    # But then we're going to have to keep track of these deletions so the column labels can reflect what
	    # is actually at the position whenever there is a break like this
	    self.columnlabels = []
	    columns_to_delete = []
	    for c in range(0, ncols):
		columnbad = True
		for r in range(0, len(self.sequences)):
		    if (c < len(self.sequences[r])):
			if (self.sequences[r][c] != "-"):
			    columnbad = False
			    break
		    if (not(columnbad)):
			break
		if (columnbad):
		    columns_to_delete.append(c)
		else:
		    self.columnlabels.append(c+1)
	    # Now take out these columns
	    for i in range(len(columns_to_delete)-1, -1, -1):
		col = columns_to_delete[i]
		for r in range(0, len(self.sequences)):
		    if (col < len(self.sequences[r])):
			self.indxToSeqPos[r].pop(col)
			self.sequences[r] = self.sequences[r][0:col] + self.sequences[r][col+1:].strip()
	    # Now update the SequenceWindow
	    for r in range(0, len(self.sequences)):
		if (r < len(self.sequences) - 1):
		    self.updateSeqViewer(r, relabel=False)
		else:
		    self.updateSeqViewer(r, relabel=True)
	else:
	    # Simple, just make everything numbered from one again
	    r = 0
	    for pose in self.poses:
		if (pose):
		    for i in range(0, len(pose[0])):
			chain = self.IDs[r][len(self.IDs[r])-1]
			if (chain == "_"):
			    chain = " "
			ch = pose[0][chain]
			self.indxToSeqPos[r] = []
			self.sequences[r] = ""
			for residue in ch:
			    self.indxToSeqPos[r].append(residue.id)
			    self.sequences[r] = self.sequences[r] + AA3to1(residue.resname)
			r = r + 1
	    # Now update the SequenceWindow
	    for r in range(0, len(self.sequences)):
		if (r < len(self.sequences) - 1):
		    self.updateSeqViewer(r, relabel=False)
		else:
		    self.updateSeqViewer(r, relabel=True)

    # ======================================================================================================
    # Download manager functions

    def recoverFromError(self, protocol):
	# This function tells the user what the error was and tries to revert the protocol
	# back to the pre-daemon state so the main GUI can continue to be used
	f = open("errreport", "r")
	errmsg = "An error was encountered during a submitted " + protocol + " protocol:\n\n"
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

    def downloader(self, event):
	self.DownloadTimer.Stop()
	# You have to do this on Linux due to the libc bug
	if (platform.system() == "Linux"):
	    res_init()
	# Read the jobs that are active
	self.activeJobs = []
	removeJob = []
	changed = False
	goToSandbox()
	try:
	    f = open("downloadwatch", "r")
	except:
	    return
	for aline in f:
	    self.activeJobs.append(aline.strip())
	    removeJob.append(False)
	f.close()
	for i in range(0, len(self.activeJobs)):
	    job = self.activeJobs[i]
	    #elif (job.startswith("MSD") or job.startswith("ANTIBODY") or job.startswith("DOCK") or job.startswith("PMUTSCAN") or job.startswith("BACKRUB") or job.startswith("KIC") or job.startswith("FLEXPEP")):
	    jobtype = job.split("\t")[0]
	    ID = job.split("\t")[1].strip()
	    jobURL = job.split("\t")[3].strip()
	    if (len(job.split("\t")) == 5):
		desc = job.split("\t")[4].strip()
	    else:
		desc = ID
	    if (job.startswith("MSD")):
		packageext = ".msdar"
	    elif (job.startswith("PMUTSCAN")):
		packageext = ".scan"
	    else:
		packageext = ".ensb"
	    if (jobURL.lower().startswith("localhost")):
		serverlocation = getServerName()[jobURL.find(":")+1:]
		if (platform.system() == "Windows"):
		    resultsdir = serverlocation + "\\results\\" + ID + "\\"
		    filepath = resultsdir + "results" + packageext
		else:
		    resultsdir = serverlocation + "/results/" + ID + "/"
		    filepath = resultsdir + "results" + packageext
	    else:
		serverlocation = None
		URL = jobURL + "/results/" + ID + "/results" + packageext
	    try:
		if (serverlocation):
		    if (not(os.path.isfile(filepath)) and not(os.path.isfile(filepath.split(packageext)[0] + ".gz"))):
			raise Exception()
		    if (jobtype == "MSD"):
			dlg = wx.MessageDialog(self, "Your MSD package job ID " + desc + " is ready.", "MSD Ready", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
		    elif (jobtype == "ANTIBODY"):
			dlg = wx.MessageDialog(self, "Your antibody package job ID " + desc + " is ready.", "Antibodies Ready", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
		    elif (jobtype == "COMPMODEL"):
			dlg = wx.MessageDialog(self, "Your comparative modeling package job ID " + desc + " is ready.", "Structures Ready", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
		    elif (jobtype == "DOCK"):
			dlg = wx.MessageDialog(self, "Your docking package job ID " + desc + " is ready.", "Docking Ready", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
		    elif (jobtype == "PMUTSCAN"):
			dlg = wx.MessageDialog(self, "Your point mutant scanning job ID " + desc + " is ready.", "Point Mutants Ready", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
		    elif (jobtype == "BACKRUB"):
			dlg = wx.MessageDialog(self, "Your backrub ensemble job ID " + desc + " is ready.", "Backrub Ensemble Ready", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
		    elif (jobtype == "KIC"):
			dlg = wx.MessageDialog(self, "Your KIC ensemble job ID " + desc + " is ready.", "KIC Ensemble Ready", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
		    elif (jobtype == "FLEXPEP"):
			dlg = wx.MessageDialog(self, "Your flexible peptide docking package job ID " + desc + " is ready.", "Flexible Peptide Docking Ready", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
		    else:
			# CUSTOM MODULE
			dlg = wx.MessageDialog(self, "Your " + jobtype + " package job ID " + desc + " is ready.", jobtype + " Ready", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
		    dlg.ShowModal()
		    dlg.Destroy()
		    try:
			os.rename(filepath, "results" + packageext)
		    except:
			# Not there, maybe it's a custom .gz file
			os.rename(filepath.split(packageext)[0] + ".gz", "results.gz")
			filepath = filepath.split(packageext)[0] + ".gz"
			packageext = ".gz"
		else:
		    try:
			downloadpage = urllib2.urlopen(URL, timeout=1) # To make sure its there before display the dialog
			downloadpage.close()
		    except:
			downloadpage = urllib2.urlopen(URL.split(packageext)[0] + ".gz", timeout=1) # To make sure its there before display the dialog
			downloadpage.close()
			URL = URL.split(packageext)[0] + ".gz"
			packageext = ".gz"
			# Maybe it's a custom module that uses a .gz extension
		    if (jobtype == "MSD"):
			dlg = wx.MessageDialog(self, "Your MSD package job ID " + desc + " is ready.", "MSD Download Ready", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
		    elif (jobtype == "ANTIBODY"):
			dlg = wx.MessageDialog(self, "Your antibody package job ID " + desc + " is ready.", "Antibody Download Ready", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
		    elif (jobtype == "COMPMODEL"):
			dlg = wx.MessageDialog(self, "Your comparative modeling package job ID " + desc + " is ready.", "Structure Download Ready", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
		    elif (jobtype == "DOCK"):
			dlg = wx.MessageDialog(self, "Your docking package job ID " + desc + " is ready.", "Docking Download Ready", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
		    elif (jobtype == "PMUTSCAN"):
			dlg = wx.MessageDialog(self, "Your point mutant scanning job ID " + desc + " is ready.", "Point Mutants Download Ready", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
		    elif (jobtype == "BACKRUB"):
			dlg = wx.MessageDialog(self, "Your backrub ensemble job ID " + desc + " is ready.", "Backrub Ensemble Download Ready", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
		    elif (jobtype == "KIC"):
			dlg = wx.MessageDialog(self, "Your KIC ensemble job ID " + desc + " is ready.", "KIC Ensemble Download Ready", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
		    elif (jobtype == "FLEXPEP"):
			dlg = wx.MessageDialog(self, "Your flexible peptide docking package job ID " + desc + " is ready.", "Flexible Peptide Docking Download Ready", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
		    else:
			# CUSTOM MODULE
			dlg = wx.MessageDialog(self, "Your " + jobtype + " package job ID " + desc + " is ready.", jobtype + " Download Ready", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
		    dlg.ShowModal()
		    dlg.Destroy()
		    if (jobtype == "MSD"):
			busyDlg = wx.BusyInfo("Downloading MSD archive, please wait...")
		    elif (jobtype == "ANTIBODY"):
			busyDlg = wx.BusyInfo("Downloading antibody archive, please wait...")
		    elif (jobtype == "COMPMODEL"):
			busyDlg = wx.BusyInfo("Downloading structure archive, please wait...")
		    elif (jobtype == "DOCK"):
			busyDlg = wx.BusyInfo("Downloading docking archive, please wait...")
		    elif (jobtype == "PMUTSCAN"):
			busyDlg = wx.BusyInfo("Downloading point mutant scan report, please wait...")
		    elif (jobtype == "BACKRUB"):
			busyDlg = wx.BusyInfo("Downloading backrub archive, please wait...")
		    elif (jobtype == "KIC"):
			busyDlg = wx.BusyInfo("Downloading KIC archive, please wait...")
		    elif (jobtype == "FLEXPEP"):
			busyDlg = wx.BusyInfo("Downloading flexpep archive, please wait...")
		    else:
			# CUSTOM MODULE
			busyDlg = wx.BusyInfo("Downloading " + jobtype + " archive, please wait...")
		    (oldfilename, info) = urllib.urlretrieve(URL, "results" + packageext)
		    busyDlg.Destroy()
		    del busyDlg
		while (True):
		    custom_module = False
		    if (jobtype == "MSD"):
			dlg = wx.FileDialog(
			    self, message="Save the MSD package",
			    defaultDir=self.cwd,
			    defaultFile=ID,
			    wildcard="MSD Archives (*.msdar)|*.msdar",
			    style=wx.SAVE | wx.CHANGE_DIR)
		    elif (jobtype == "ANTIBODY"):
			dlg = wx.FileDialog(
			    self, message="Save the antibody package",
			    defaultDir=self.cwd,
			    defaultFile=ID,
			    wildcard="Ensemble Archives (*.ensb)|*.ensb",
			    style=wx.SAVE | wx.CHANGE_DIR)
		    elif (jobtype == "COMPMODEL"):
			dlg = wx.FileDialog(
			    self, message="Save the structure package",
			    defaultDir=self.cwd,
			    defaultFile=ID,
			    wildcard="Ensemble Archives (*.ensb)|*.ensb",
			    style=wx.SAVE | wx.CHANGE_DIR)
		    elif (jobtype == "DOCK"):
			dlg = wx.FileDialog(
			    self, message="Save the docking package",
			    defaultDir=self.cwd,
			    defaultFile=ID,
			    wildcard="Ensemble Archives (*.ensb)|*.ensb",
			    style=wx.SAVE | wx.CHANGE_DIR)
		    elif (jobtype == "PMUTSCAN"):
			dlg = wx.FileDialog(
			    self, message="Save the scan file",
			    defaultDir=self.cwd,
			    defaultFile=ID,
			    wildcard="Point Mutant Scan (*.scan)|*.scan",
			    style=wx.SAVE | wx.CHANGE_DIR)
		    elif (jobtype == "BACKRUB"):
			dlg = wx.FileDialog(
			    self, message="Save the backrub package",
			    defaultDir=self.cwd,
			    defaultFile=ID,
			    wildcard="Ensemble Archives (*.ensb)|*.ensb",
			    style=wx.SAVE | wx.CHANGE_DIR)
		    elif (jobtype == "KIC"):
			dlg = wx.FileDialog(
			    self, message="Save the KIC package",
			    defaultDir=self.cwd,
			    defaultFile=ID,
			    wildcard="Ensemble Archives (*.ensb)|*.ensb",
			    style=wx.SAVE | wx.CHANGE_DIR)
		    elif (jobtype == "FLEXPEP"):
			dlg = wx.FileDialog(
			    self, message="Save the flexpep package",
			    defaultDir=self.cwd,
			    defaultFile=ID,
			    wildcard="Ensemble Archives (*.ensb)|*.ensb",
			    style=wx.SAVE | wx.CHANGE_DIR)
		    else:
			# CUSTOM MODULE
			custom_module = True
			if (packageext == ".ensb"):
			    dlg = wx.FileDialog(
				self, message="Save the results package",
				defaultDir=self.cwd,
				defaultFile=ID,
				wildcard="Ensemble Archives (*.ensb)|*.ensb",
				style=wx.SAVE | wx.CHANGE_DIR)
			else:
			    dlg = wx.FileDialog(
				self, message="Save the results package",
				defaultDir=self.cwd,
				defaultFile=ID,
				wildcard="GZipped Archives (*.gz)|*.gz",
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
			self.cwd = str(paths[len(paths)-1][0:lastDirIndx])
			self.saveWindowData(None)
			filename = str(paths[0]).split(packageext)[0] + packageext
			# Does it exist already?  If so, ask if the user really wants to overwrite it
			if (os.path.isfile(filename)):
			    dlg2 = wx.MessageDialog(self, "The file " + filename + " already exists.  Overwrite it?", "Filename Already Exists", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
			    if (dlg2.ShowModal() == wx.ID_NO):
				dlg2.Destroy()
				logInfo("Cancelled save operation due to filename already existing")
				continue
			    else:
				os.remove(str(filename))
			    dlg2.Destroy()
			break
		goToSandbox()
		os.rename("results" + packageext, str(filename))
		if (jobtype != "PMUTSCAN"):
		    # Unpackage the files
		    gzipfile = gzip.open(str(filename), "rb")
		    prefix = filename.split(packageext)[0]
		    readingData = False
		    for aline in gzipfile:
			if (aline.startswith("BEGIN PDB")):
			    if (jobtype == "MSD"):
				if (platform.system() == "Windows"):
				    f = open(self.cwd + "\\" + aline.split()[len(aline.split())-1].strip(), "w")
				else:
				    f = open(self.cwd + "/" + aline.split()[len(aline.split())-1].strip(), "w")
			    elif (jobtype in ["ANTIBODY", "DOCK", "BACKRUB", "KIC", "FLEXPEP", "COMPMODEL"] or custom_module):
				indx = aline[aline.rfind("_")+1:].strip()
				f = open(prefix + "_" + indx, "w")
			    readingData = True
			elif (aline.startswith("BEGIN FILE")):
			    f = open(self.cwd + "/" + aline.split("BEGIN FILE")[1].strip(), "w")
			    readingData = True
			elif (aline.startswith("END PDB")):
			    f.close()
			    readingData = False
			elif (aline.startswith("END FILE")):
			    f.close()
			    readingData = False
			elif (readingData):
			    f.write(aline.strip() + "\n")
		    gzipfile.close()
		changed = True
		removeJob[i] = True
	    except:
		# Not there yet
		pass
	    URL = jobURL + "/results/" + ID + "/errreport"
	    try:
		# Look for an error file
		if (serverlocation):
		    downloadpage = open(resultsdir + "errreport")
		else:
		    downloadpage = urllib2.urlopen(URL, timeout=1) # To make sure its there before display the dialog
		f = open("errreport", "w")
		for aline in downloadpage:
		    f.write(aline.strip() + "\n")
		downloadpage.close()
		f.close()
		self.recoverFromError(jobtype)
		changed = True
		removeJob[i] = True
	    except:
		pass
	# If a change happened, update the "downloadwatch" file to take the completed jobs out
	if (changed):
	    goToSandbox()
	    try:
		os.remove("downloadwatch")
	    except:
		pass
	    f = open("downloadwatch", "w")
	    for i in range(0, len(self.activeJobs)):
		if (not(removeJob[i])):
		    f.write(self.activeJobs[i].strip() + "\n")
	    f.close()
	self.DownloadTimer.Start(60000)

    # ======================================================================================================
    # Helper functions used by other protocols

    def doesResidueExist(self, model, chain, seqpos):
	# Useful function for determining whether a residue exists on a model
	# You can easily figure this out by trying to access it in the lookup table
	# Return its coordinates so resfile data can be updated easily
	try:
	    return self.selectLookup[(model + "|" + chain, int(seqpos))]
	except:
	    return False

    def doesAtomExist(self, model, chain, seqpos, atomname):
	# Useful function for determining whether a residue exists on a model
	# You can easily figure this out by trying to access it in the lookup table
	# Return its coordinates so constraints data can be updated easily
	try:
	    (r, c) = self.selectLookup[(model + "|" + chain, int(seqpos))]
	    if (chain == "_"):
		chain = " "
	    poseindx = self.getPoseIndex(r)
	    for residue in self.poses[poseindx][0][chain]:
		if (residue.id[1] == int(seqpos)):
		    # Now see if the atom name is there
		    for atom in residue:
			if (atom.id.upper() == atomname.upper()):
			    return (r, c)
	    raise Exception()
	except:
	    return False

    def updatePyMOLSelection(self, event):
	self.PyMOLUpdateTimer.Stop()
	self.selectUpdate(updatePyMOL=True)

    def getSelectedResidues(self):
	# Useful function for getting selected information
	# First we have to update the selection in case the user was fiddling in PyMOL and
	# then went straight to a protocol.  The sequence window doesn't know about the PyMOL
	# changes yet but the protocols read selections from it
	self.selectUpdate(None)
	topLefts = self.SeqViewer.GetSelectionBlockTopLeft()
	bottomRights = self.SeqViewer.GetSelectionBlockBottomRight()
	selectionkey = []
	for i in range(0, len(topLefts)):
	    for r in range(topLefts[i][0], bottomRights[i][0]+1):
		if (r >= len(self.poses)):
		    continue
		offset = 0
		for j in range(r, -1, -1):
		    if (self.poses[j]):
			poseindx = j
			break
		    offset = offset + len(self.sequences[j-1])
		for c in range(topLefts[i][1], bottomRights[i][1]+1):
		    if (c >= len(self.indxToSeqPos[r])):
			continue
		    if (self.indxToSeqPos[r][c] != "-"):
			selectionkey.append(str(r) + ":" + str(poseindx) + ":" + str(c+offset))
	# Sort
	selectionkey.sort()
	# Now change it back to meaningful integers
	lastr = -1
	selection = []
	for key in selectionkey:
	    r = int(key.split(":")[0])
	    poseindx = int(key.split(":")[1])
	    resi = int(key.split(":")[2])
	    if (lastr != r):
		selection.append([r, poseindx, []])
	    selection[len(selection)-1][2].append(resi)
	return selection

    def numChains(self):
	return self.SeqViewer.NumberRows

    def getChainPose(self, row):
	# Gets a pose object for only one chain given either by the row in the sequence window or the ID
	if (not("int" in str(type(row)))):
	    row = self.IDs.index(row)
	poseindx = self.getPoseIndex(row)
	chain = self.IDs[row][len(self.IDs[row])-1]
	if (chain == "_"):
	    chain = " "
	return self.poses[poseindx][0][chain]

    def getModelForChain(self, r):
	# Useful function to figure out what the name of the model is from an ID in the SeqViewer
	ID = self.IDs[r]
	model = ""
	fields = ID.split("|")
	for field in fields[0:len(fields)-1]:
	    model = model + field + "|"
	model = model[0:len(model)-1]
	return model

    def getPoseIndex(self, r):
	# Given a row in the SeqViewer, find out the index of the pose for this sequence
	poseindx = -1
	for i in range(r, -1, -1):
	    if (self.poses[i]):
		poseindx = i
		break
	return poseindx

    def getPoseIndexForModel(self, model):
	# Given a model, find the poseindx
	poseindx = -1
	for i in range(0, self.SeqViewer.NumberRows):
	    if (self.poses[i] and model.strip() == self.IDs[i][0:len(self.IDs[i].strip())-2]):
		poseindx = i
		break
	return poseindx

    def getRosettaIndex(self, model, chain, seqpos):
	# Useful function for finding what the Rosetta index of a model/chain/seqpos is
	poseindx = self.getPoseIndexForModel(model)
	ires = 1
	for ch in self.poses[poseindx][0]:
	    for residue in ch:
		thischain = ch.id
		thisseqpos = residue.id[1]
		if (chain == "_" and (thischain == "_" or len(thischain.strip()) == 0) and int(thisseqpos) == int(seqpos)):
		    return ires
		elif (chain == thischain and int(thisseqpos) == int(seqpos)):
		    return ires
		ires = ires + 1
	return -1

    def getResidueInfo(self, model, rosettaindx):
	# Useful function for getting the chain and seqpos given a model and rosetta index
	ires = 1
	poseindx = self.getPoseIndexForModel(model)
	for ch in self.poses[poseindx][0]:
	    for residue in ch:
		if (ires == rosettaindx):
		    return (ch.id, residue.id[1])
		ires = ires + 1
	return None

    def getResidueTypeFromRosettaIndx(self, model, rosettaindx):
	# Given the absolute Rosetta index of a residue, get its residue type
	poseindx = self.getPoseIndexForModel(model)
	ires = rosettaindx
	for r in range(poseindx, len(self.sequences)):
	    if (ires - len(self.sequences[r]) <= 0):
		return self.sequences[r][ires-1]
	    ires -= len(self.sequences[r])

    def getIsCanonicalAA(self, r, c):
	# Find out if this is a canonical CAA or not
	poseindx = self.getPoseIndex(r)
	offset = 1
	for i in range(poseindx, r):
	    offset = offset + len(self.sequences[i])
	chain = self.IDs[r][len(self.IDs[r])-1]
	if (chain == "_"):
	    chain = " "
	restype = self.poses[poseindx][0][chain][self.indxToSeqPos[r][c]].resname
	if (not(restype in "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR ")):
	    return False
	return True

    def getSelectedModels(self):
	# Gets the names of the models that are currently selected for easy PyMOL access
	topLefts = self.SeqViewer.GetSelectionBlockTopLeft()
	bottomRights = self.SeqViewer.GetSelectionBlockBottomRight()
	models = []
	rowsSelected = []
	for i in range(0, len(topLefts)):
	    for r in range(topLefts[i][0], bottomRights[i][0]+1):
		if (not(r in rowsSelected)):
		    rowsSelected.append(r)
	rowsSelected.sort()
	for r in rowsSelected:
	    newmodel = self.getModelForChain(r)
	    if (not(newmodel in models)):
		models.append(newmodel)
	return models

    def getSelectedChains(self):
	# Gets all the selected chains
	topLefts = self.SeqViewer.GetSelectionBlockTopLeft()
	bottomRights = self.SeqViewer.GetSelectionBlockBottomRight()
	models = []
	rowsSelected = []
	for i in range(0, len(topLefts)):
	    for r in range(topLefts[i][0], bottomRights[i][0]+1):
		if (not(r in rowsSelected)):
		    rowsSelected.append(r)
	rowsSelected.sort()
	return rowsSelected

    def areResiduesSelected(self):
	# Used by the selection panel to see if residues are selected
	# If not, the button usually defaults to applying the operation to everything
	return (len(self.SeqViewer.GetSelectionBlockTopLeft()) > 0)

    def reloadPose(self, poseindx, model, newpdb):
	# This function reloads a pose from a protocol and reloads the sequence, which may
	# have changed
	# Create a copy of the inputted pose
	self.poses[poseindx] = self.pdbreader.get_structure(model, newpdb)
	# Now reload the sequences
	i = 0
	for ch in self.poses[poseindx][0]:
	    r = poseindx + i
	    seq = ""
	    self.indxToSeqPos[r] = []
	    for residue in ch:
		seq = seq + AA3to1(residue.resname)
		self.indxToSeqPos[r].append(residue.id)
	    self.sequences[r] = seq
	    # Some of the protocols change the chain IDs, so we have to update the IDs and labels
	    chain = ch.id
	    if (len(chain.strip()) == 0):
		chain = "_"
	    self.IDs[r] = model + "|" + chain
	    self.updateSeqViewer(r)
	    i = i + 1
	# Rerun DSSP
	self.rerunDSSPForModel(model, newpdb)
	self.regenerateLookupTable()

    def rerunDSSPForModel(self, model, pdbmodel="None"):
	try:
	    if (self.dsspexe == "N/A"):
		return
	    if (pdbmodel == "None"):
		pdbmodel = model + ".pdb"
	    poseindx = self.getPoseIndexForModel(model)
	    dsspdata = Bio.PDB.DSSP(self.poses[poseindx][0], pdbmodel, dssp=self.dsspexe)
	    for c in self.poses[poseindx][0]:
		chain = c.id
		if (len(chain.strip()) == 0):
		    chain = "_"
		for r in c:
		    # Alter the states in PyMOL with the DSSP secondary structure predictions
		    try:
			ss = dsspdata[(c.id, r.id[1])][1]
			# Not all the DSSP codes match up with PyMOL, so make them match
			if (ss == "E"):
			    # Beta sheet
			    ss = "S"
			elif (ss == "S"):
			    # Bend
			    ss = "B"
			elif (ss == "-"):
			    # Loop
			    ss = "L"
			if (c.id != " "):
			    self.cmd.alter("model " + model + " and chain " + c.id + " and resi " + str(r.id[1]), "ss=\"" + ss + "\"")
			else:
			    self.cmd.alter("model " + model + " and resi " + str(r.id[1]), "ss=\"" + ss + "\"")
		    except:
			ss = "L"
	except:
	    # In case something goes wrong, just use PyMOL's SS coloring and don't crash
	    pass

    def reloadFromPyMOL(self, poseindx):
	# This function reloads the BioPython data dumped from PyMOL
	# You need this for docking in case the user was rotating the binding partners in PyMOL before the
	# docking simulation
	model = self.IDs[poseindx][0:len(self.IDs[poseindx])-2]
	self.cmd.save(model + ".pdb", model)
	fixPyMOLSave(model + ".pdb")
	self.poses[poseindx] = self.pdbreader.get_structure(model, model + ".pdb")

    def logSelection(self):
	# Periodically log which residues are selected
	topLefts = self.SeqViewer.GetSelectionBlockTopLeft()
	bottomRights = self.SeqViewer.GetSelectionBlockBottomRight()
	for i in range(0, len(topLefts)):
	    if (topLefts[i][0] >= len(self.IDs)):
		continue
	    topLeftPos = self.IDs[topLefts[i][0]] + " position " + str(topLefts[i][1]+1)
	    bottomRightPos = self.IDs[bottomRights[i][0]] + " position " + str(bottomRights[i][1]+1)
	    logInfo("New Selection Entry: ")
	    logInfo("Selected " + topLeftPos + " to " + bottomRightPos)
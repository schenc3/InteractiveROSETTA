import wx
import os
import os.path
import sys
import platform
import subprocess
import multiprocessing
import psutil
from tools import *
# If this fails, it is probably because PyRosetta is not being imported
try:
    from daemon import daemonLoop
except:
    if (platform.system() == "Windows"):
        addon = "\n\nDefault Install Location: C:\Program Files\PyRosetta"
    else:
        addon = ""
    dlg = wx.MessageDialog(None, "InteractiveROSETTA cannot import PyRosetta!\n\nDid you install PyRosetta already?" + addon, "PyRosetta Cannot Be Imported", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
    while (True):
        if (dlg.ShowModal() == wx.ID_YES):
            # Ask the user to show us where they installed it
            dlg2 = wx.DirDialog(None, "Please navigate to and enter the PyRosetta folder",
                                defaultPath=os.path.expanduser("~"),
                                style=wx.DD_DEFAULT_STYLE
                                | wx.DD_DIR_MUST_EXIST
                                | wx.DD_CHANGE_DIR)
            if (dlg2.ShowModal() == wx.ID_OK):
                # Let's try to open import it from here
                rosettapath = str(dlg2.GetPath().strip())
                # Look for the database
                if (os.path.exists(os.path.join(rosettapath, "database"))):
                    rosettadb = os.path.join(rosettapath, "database")
                else:
                    rosettadb = os.path.join(rosettapath, "rosetta_database")
                # Save these settings so when we attempt to import from the daemon again
                # it will look in these directories
                if (platform.system() == "Windows"):
                    cfgfile = os.path.expanduser("~") + "\\InteractiveROSETTA\\seqwindow.cfg"
                else:
                    cfgfile = os.path.expanduser("~") + "/.InteractiveROSETTA/seqwindow.cfg"
                # Read our current settings
                data = []
                try:
                    f = open(cfgfile, "r")
                    for aline in f:
                        if (not("[ROSETTAPATH]" in aline) and not("[ROSETTADB]" in aline)):
                            data.append(aline.strip())
                    f.close()
                except:
                    pass
                # Write the selected settings back out
                f = open(cfgfile, "w")
                for aline in data:
                    f.write(aline + "\n")
                f.write("[ROSETTAPATH]\t" + rosettapath.strip() + "\n")
                f.write("[ROSETTADB]\t" + rosettadb.strip() + "\n")
                f.close()
                try:
                    # Try again
                    from daemon import daemonLoop
                    # Success, get out of here
                    break
                except:
                    # Still failed, tell the user and go back to the directory selection
                    dlg3 = wx.MessageDialog(None, "PyRosetta still could not be imported.", "PyRosetta Cannot Be Imported", wx.OK | wx.ICON_ERROR | wx.CENTRE)
                    dlg3.ShowModal()
                    dlg3.Destroy()
                    continue
            else:
                # Giving up, the user cannot find PyRosetta
                pass
            dlg2.Destroy()
        # If no, did they at least download it?  We can try to install it on Windows, or untar it on OSX/Linux
        if (platform.system() == "Windows"):
            msg = "Did you download the Windows installer for PyRosetta?"
        elif (platform.system() == "Darwin"):
            msg = "Did you download the Mac package for PyRosetta?"
        else:
            msg = "Did you download the Linux package for PyRosetta?"
        dlg_dl = wx.MessageDialog(None, msg, "PyRosetta Cannot Be Imported", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
        if (dlg_dl.ShowModal() == wx.ID_YES):
            # Show us where it is
            gotPyRosetta = False
            while (True):
                if (platform.system() == "Windows"):
                    msg = "Select the PyRosetta .exe Installer"
                    wc = "Installer (*.exe)|*.exe"
                else:
                    msg = "Select the PyRosetta .tar.bz2 Package"
                    wc = "Package (*.bz2;*.tar.gz;*.tgz)|*.bz2;*.tar.gz;*.tgz"
                dlg2 = wx.FileDialog(
                    None, message=msg,
                    defaultDir=os.path.expanduser("~"),
                    defaultFile="",
                    wildcard=wc,
                    style=wx.OPEN | wx.CHANGE_DIR)
                if (dlg2.ShowModal() == wx.ID_OK):
                    filename = str(dlg2.GetPath())
                    # Unpack it
                    try:
                        olddir = os.getcwd()
                        os.chdir(os.path.expanduser("~"))
                        if platform.system() == 'Darwin':
                                os.chdir('/Applications/InteractiveROSETTA.app')
                        if (platform.system() == "Windows"):
                            # Simply execute the file
                            import subprocess
                            subprocess.call([filename])
                            # Go back up and ask the user to find the folder, since they just installed it
                            continue
                        else:
                            import shutil
                            import tarfile
                            # Move it to the user's home directory
                            print "Installing PyRosetta, please wait..."
                            packagename = filename[filename.rfind("/")+1:].strip()
                            try:
                                shutil.copy(filename, packagename) #shutil.move
                            except:
                                # It's probably already in the home directory
                                pass
                            # Unpack it
                            # What is the file extension so we can figure out the reading format
                            if (packagename.endswith(".bz2")):
                                readformat = "r:bz2"
                            else:
                                readformat = "r:gz"
                            progress = wx.ProgressDialog("PyRosetta Unpack Progress", "", 7, style=wx.PD_ELAPSED_TIME | wx.PD_REMAINING_TIME)
                            tar = tarfile.open(packagename, readformat)
                            folders = ["database", "rosetta/basic", "rosetta/core", "rosetta/utility", "rosetta/numeric", "rosetta/protocols"]
                            directory = None
                            for item in tar:
                                print "Extracting " + str(item)
                                for i in range(0, len(folders)):
                                    folder = folders[i]
                                    if (folder in str(item)):
                                        folders.pop(i)
                                        progress.Update(6-len(folders), "Unpacking " + folder)
                                        break
                                # Figure out the name of the folder we are extracting into
                                if (directory is None):
                                    try:
                                        if ("/" not in str(item)):
                                            raise Exception()
                                        directory = str(item)[str(item).find("'")+1:min(str(item).find("/"), str(item).rfind("'"))]
                                    except:
                                        pass
                                tar.extract(item)
                            tar.close()
                            progress.Destroy()
                            # Move the directory to "PyRosetta"
                            print "Unpacked directory " + directory
                            try:
                                print "Setting up PyRosetta import..."
                                # Now let's try to import it, since we know where it is
                                rosettapath = os.path.expanduser("~") + "/" + directory
                                #rename if OS X
                                if platform.system() == 'Darwin':
                                        import os; os.remove(packagename)
                                        import commands; print commands.getstatusoutput('mv -v /Applications/InteractiveROSETTA.app/%s /Applications/InteractiveROSETTA.app/PyRosetta'%(directory))[1]
                                        directory = 'PyRosetta'
                                        rosettapath = '/Applications/InteractiveROSETTA.app/%s'%(directory)
                                # Look for the database
                                if (os.path.exists(os.path.join(rosettapath, "database"))):
                                    rosettadb = os.path.join(rosettapath, "database")
                                else:
                                    rosettadb = os.path.join(rosettapath, "rosetta_database")
                                # Save these settings so when we attempt to import from the daemon again
                                # it will look in these directories
                                cfgfile = os.path.expanduser("~") + "/.InteractiveROSETTA/seqwindow.cfg"
                                data = []
                                try:
                                    f = open(cfgfile, "r")
                                    for aline in f:
                                        if (not("[ROSETTAPATH]" in aline) and not("[ROSETTADB]" in aline)):
                                            data.append(aline.strip())
                                    f.close()
                                except:
                                    pass
                                # Write the selected settings back out
                                f = open(cfgfile, "w")
                                for aline in data:
                                    f.write(aline + "\n")
                                f.write("[ROSETTAPATH]\t" + rosettapath.strip() + "\n")
                                f.write("[ROSETTADB]\t" + rosettadb.strip() + "\n")
                                f.close()
                                from daemon import daemonLoop
                                os.chdir(olddir)
                                gotPyRosetta = True
                                break
                            except:
                                # Failed for some reason, but it did unpack
                                if platform.system() == 'Darwin':
                                        dlg3 = wx.MessageDialog(None, "InteractiveROSETTA unpacked to /Applications/InteractiveROSETTA.app/PyRosetta", "Please restart InteractiveROSETTA", wx.OK | wx.ICON_ERROR | wx.CENTRE)
                                        dlg3.ShowModal()
                                        dlg3.Destroy()
                                        exit()
                                else:
                                        dlg3 = wx.MessageDialog(None, "InteractiveROSETTA could not import PyRosetta, but it was unpacked to " + os.path.expanduser("~") + "/%s.  Please try starting InteractiveROSETTA again."%(directory), "PyRosetta Cannot Be Imported", wx.OK | wx.ICON_ERROR | wx.CENTRE)
                                        dlg3.ShowModal()
                                        dlg3.Destroy()
                                        exit()
                    except:
                        # Failed for some reason
                        dlg3 = wx.MessageDialog(None, "InteractiveROSETTA could not install PyRosetta.  Please install it manually.", "PyRosetta Cannot Be Installed", wx.OK | wx.ICON_ERROR | wx.CENTRE)
                        dlg3.ShowModal()
                        dlg3.Destroy()
                        exit()
                else:
                    break
            if (gotPyRosetta):
                break
        dlg_dl.Destroy()
        # Tell the user to install PyRosetta before continuing
        if (platform.system() == "Windows"):
            dlg2 = wx.MessageDialog(None, "You need to download the Windows installer and double click on the downloaded exe file to install PyRosetta.\n\nDownload: http://www.pyrosetta.org/dow", "Please Install PyRosetta", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
        elif (platform.system() == "Darwin"):
            dlg2 = wx.MessageDialog(None, "You need to download the Mac package and do not untar it. \n\nDownload: http://www.pyrosetta.org/dow", "Please Install PyRosetta", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
        else:
            dlg2 = wx.MessageDialog(None, "You need to download the Linux package and untar the package.  To do that, open a terminal and execute \"tar xjvf [package]\" where [package] is the name of the downloaded file.\n\nDownload: http://www.pyrosetta.org/dow", "Please Install PyRosetta", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
        dlg2.ShowModal()
        dlg2.Destroy()
        # We're done until PyRosetta is correctly installed and imported
        exit()
    dlg.Destroy()
# Is that molfile2params stuff unpacked?
# Check for it
# I used to unpack it into the scripts directory during installation, but it could end up failing
# Now I am going to unpack it to the sandbox, since it does not require root access to write to that
# directory and can be performed automatically here without having the user execute untar commands
molfile_unpacked = True
if (platform.system() == "Windows"):
    if (not(os.path.isfile(os.path.expanduser("~") + "\\InteractiveROSETTA\\molfile2params\\molfile_to_params.py"))):
        molfile_unpacked = False
    elif (not(os.path.isfile(os.path.expanduser("~") + "\\InteractiveROSETTA\\molfile2params\\rosetta_py\\__init__.py"))):
        molfile_unpacked = False
    elif (not(os.path.isfile(os.path.expanduser("~") + "\\InteractiveROSETTA\\molfile2params\\rosetta_py\\io\\__init__.py"))):
        molfile_unpacked = False
    elif (not(os.path.isfile(os.path.expanduser("~") + "\\InteractiveROSETTA\\molfile2params\\rosetta_py\\io\\mdl_molfile.py"))):
        molfile_unpacked = False
    elif (not(os.path.isfile(os.path.expanduser("~") + "\\InteractiveROSETTA\\molfile2params\\rosetta_py\\io\\pdb.py"))):
        molfile_unpacked = False
    elif (not(os.path.isfile(os.path.expanduser("~") + "\\InteractiveROSETTA\\molfile2params\\rosetta_py\\utility\\__init__.py"))):
        molfile_unpacked = False
    elif (not(os.path.isfile(os.path.expanduser("~") + "\\InteractiveROSETTA\\molfile2params\\rosetta_py\\utility\\r3.py"))):
        molfile_unpacked = False
    elif (not(os.path.isfile(os.path.expanduser("~") + "\\InteractiveROSETTA\\molfile2params\\rosetta_py\\utility\\rankorder.py"))):
        molfile_unpacked = False
else:
    if (not(os.path.isfile(os.path.expanduser("~") + "/.InteractiveROSETTA/molfile2params/molfile_to_params.py"))):
        molfile_unpacked = False
    elif (not(os.path.isfile(os.path.expanduser("~") + "/.InteractiveROSETTA/molfile2params/rosetta_py/__init__.py"))):
        molfile_unpacked = False
    elif (not(os.path.isfile(os.path.expanduser("~") + "/.InteractiveROSETTA/molfile2params/rosetta_py/io/__init__.py"))):
        molfile_unpacked = False
    elif (not(os.path.isfile(os.path.expanduser("~") + "/.InteractiveROSETTA/molfile2params/rosetta_py/io/mdl_molfile.py"))):
        molfile_unpacked = False
    elif (not(os.path.isfile(os.path.expanduser("~") + "/.InteractiveROSETTA/molfile2params/rosetta_py/io/pdb.py"))):
        molfile_unpacked = False
    elif (not(os.path.isfile(os.path.expanduser("~") + "/.InteractiveROSETTA/molfile2params/rosetta_py/utility/__init__.py"))):
        molfile_unpacked = False
    elif (not(os.path.isfile(os.path.expanduser("~") + "/.InteractiveROSETTA/molfile2params/rosetta_py/utility/r3.py"))):
        molfile_unpacked = False
    elif (not(os.path.isfile(os.path.expanduser("~") + "/.InteractiveROSETTA/molfile2params/rosetta_py/utility/rankorder.py"))):
        molfile_unpacked = False
if (not(molfile_unpacked)):
    import tarfile
    import shutil
    if (platform.system() == "Windows"):
        cfgfile = os.path.expanduser("~") + "\\InteractiveROSETTA\\seqwindow.cfg"
    else:
        cfgfile = os.path.expanduser("~") + "/.InteractiveROSETTA/seqwindow.cfg"
    # Find PyRosetta
    if (platform.system() == "Windows"):
        filetofind = "\\toolbox\\molfile2params.tar.gz"
    else:
        filetofind = "/toolbox/molfile2params.tar.gz"
    fin = open(cfgfile, "r")
    rosettapath = ""
    for aline in fin:
        if (aline.startswith("[ROSETTAPATH]")):
            rosettapath = aline.split("\t")[1].strip()
        break
    fin.close()
    if (not(os.path.isfile(rosettapath + filetofind))):
        # Was PyRosetta imported from the environment variable correctly?
        if (platform.system() == "Windows"):
            pythonpaths = os.environ["PYTHONPATH"].split(";")
        else:
            pythonpaths = os.environ["PYTHONPATH"].split(":")
        for pythonpath in pythonpaths:
            if (os.path.isfile(pythonpath + filetofind)):
                rosettapath = pythonpath
                break
    # Let's copy this package to the sandbox, because then we can unpack it without root access
    if (platform.system() == "Windows"):
        molfiletgz = rosettapath + "\\toolbox\\molfile2params.tar.gz"
        localmolfiletgz = os.path.expanduser("~") + "\\InteractiveROSETTA\\molfile.tgz"
    else:
        molfiletgz = rosettapath + "/toolbox/molfile2params.tar.gz"
        localmolfiletgz = os.path.expanduser("~") + "/.InteractiveROSETTA/molfile.tgz"
    shutil.copy(molfiletgz, localmolfiletgz)
    # Unpack it
    olddir = os.getcwd()
    goToSandbox()
    tar = tarfile.open(localmolfiletgz, "r:gz")
    for item in tar:
        tar.extract(item)
    tar.close()
    os.chdir(olddir)
from selection import SelectPanel
from superimposition import SuperimpositionPanel
from minimization import MinimizationPanel
from ensemblebrowser import EnsembleBrowserPanel
from compmodel import CompModelPanel
from fixbb import FixbbPanel
from dagview import DagViewPanel
#try:
    # This uses openbabel to convert PDBs to MOL2 files, which makes it way easier
    # to add new parameters files since we can just read the PDB that the user is trying
    # to load directly
from advresiduecreator import ResidueCreatorPanel
#except:
    # Maybe the user couldn't compile openbabel, so we'll default to the original version
    # that requires the user to generate their own MOL2 files
#    print "A Python OpenBabel installation was not detected on your system."
#    print "Although OpenBabel is not required, it can greatly simplify parameterizing new residues"
#    print "On Windows, you need to install the main OpenBabel as well: http://openbabel.org/wiki/Category:Installation"
#    print "On Debian, Python OpenBabel can be installed using apt-get: sudo apt-get install python-openbabel"
#    print "On Mac and RedHat, you need to compile from source: http://open-babel.readthedocs.org/en/latest/Installation/install.html#compiling-open-babel"
#    from residuecreator import ResidueCreatorPanel
from pointmutations import PointMutationsPanel
from kic import KICPanel
from docking import DockingPanel
from msd import MSDPanel
from pmutscan import PointMutantScanPanel
from surfaces import SurfacesPanel
from antibody import AntibodyPanel
from ensemblegen import EnsembleGenPanel
from flexpepdock import FlexPepDockPanel
from modulemanager import ModuleManagerPanel
from biotools import BioToolsPanel
from indel import INDELmodelPanel
import glob

class ProtocolsPanel(wx.Panel):
    def __init__(self, parent, W, H):
        wx.Panel.__init__(self, parent, id=-1, pos=(0, 0), size=(350, H-265), name="ProtocolsPanel")
        self.W = W
        self.H = H
        self.parent = parent
        self.SetBackgroundColour("#333333")
        self.Show()

        if (platform.system() == "Windows"):
            self.label = wx.StaticText(self, -1, "Protocols", (5, 5), (340, 25), wx.ALIGN_CENTRE)
            self.label.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.label = wx.StaticBitmap(self, -1, wx.Image(self.parent.scriptdir + "/images/osx/protocols/lblProtocols.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(5, 5), size=(340, 25))
        else:
            self.label = wx.StaticText(self, -1, "Protocols", style=wx.ALIGN_CENTRE)
            self.label.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
            resizeTextControlForUNIX(self.label, 0, self.GetSize()[0])
        self.label.SetForegroundColour("#FFFFFF")

        if (platform.system() == "Darwin"):
            self.protMenu = wx.ComboBox(self, pos=(5, 30), size=(230, 25), choices=[], style=wx.CB_READONLY)
        else:
            self.protMenu = wx.ComboBox(self, pos=(5, 30), size=(230, 25), choices=[], style=wx.CB_READONLY | wx.CB_SORT)
        self.protMenu.SetToolTipString("List of currently available protocols")

        self.default_protocols = ["Antibody Modeling",
                   "Biological Tools",
                   "Docking",
                   "Energy Minimization",
                   "Ensemble Browser",
                   "Ensemble Generation",
                   "Flexible Peptide Docking",
                   "INDEL Loop Modeler",
                   "Loop Modeling (KIC)",
                   "Module Manager",
                   "Molecular Surfaces",
                   "Pathway Visualization (GeoFold)",
                   "Point Mutant Scan",
                   "Point Mutations",
                   "Protein Design (Fixbb)",
                   "Protein Design (MSD)",
                   "Residue/Ligand Creator",
                   "Structure Prediction (Comparative Modeling)",
                   "Superimposition"]
        self.readModules()
        self.protMenu.SetSelection(self.protocols.index("Superimposition"))

        if (platform.system() == "Darwin"):
            self.GoBtn = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.scriptdir + "/images/osx/protocols/GoBtn.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(240, 30), size=(100, 25))
        else:
            self.GoBtn = wx.Button(self, id=-1, label="Go!", pos=(240, 30), size=(100, 25))
            self.GoBtn.SetForegroundColour("#000000")
            self.GoBtn.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.GoBtn.Bind(wx.EVT_BUTTON, self.changeProtocol)
        self.GoBtn.SetToolTipString("Change to the selected protocol")
        self.currentProtocol = "Superimposition"

        self.protPanel = SuperimpositionPanel(self, W, H)

    def readModules(self):
        self.protocols = self.default_protocols[:]
        # Let's see if there are any custom modules to import
        # Everything needs to be in a try/except block to prevent corrupted modules from
        # crashing the whole program
        # Ignore anything that does not import properly and write a message out to the terminal
        olddir = os.getcwd()
        goToSandbox("modules")
        sys.path.append(os.getcwd())
        moduledirs = glob.glob("*")
        self.modules = {}
        for moduledir in moduledirs:
            # Ignore the template
            if (moduledir == "template"):
                continue
            try:
                module = __import__(moduledir)
            except:
                print "Failed to import " + moduledir + ", does " + moduledir + "/__init__.py exist?"
                continue
            try:
                # Does the panel exist?
                module.ModulePanel
            except:
                print "Failed to import " + moduledir + ", could not find the class \"ModulePanel\""
                continue
            # Look for the name of the protocol in the comments of the module's script
            # If the name does not exist, just use the directory name
            modulename = moduledir
            fin = open(moduledir + "/__init__.py", "r")
            for aline in fin:
                if (aline.startswith("### PROTOCOL NAME:")):
                    modulename = aline.split("### PROTOCOL NAME:")[1].strip()
                    break
            fin.close()
            # Does this name exist already
            nameTaken = False
            for name in self.protocols:
                if (name.strip() == modulename.strip()):
                    nameTaken = True
                    break
            if (nameTaken):
                i = 2
                while (True):
                    nameTaken = False
                    for name in self.protocols:
                        if (name.strip() == modulename.strip() + " (" + str(i) + ")"):
                            nameTaken = True
                            break
                    if (not(nameTaken)):
                        modulename = modulename.strip() + " (" + str(i) + ")"
                        break
                    i += 1
            # Add it into the list of protocols
            addedIn = False
            for i in range(0, len(self.protocols)):
                if (modulename < self.protocols[i]):
                    addedIn = True
                    self.protocols.insert(i, modulename)
                    break
            if (not(addedIn)):
                self.protocols.append(modulename)
            self.modules[modulename] = module
        os.chdir(olddir)
        self.protMenu.Clear()
        self.protMenu.AppendItems(self.protocols)

    def changeProtocol(self, event):
        logInfo("Go button clicked")
        selectedProtocol = self.protMenu.GetStringSelection()
        if (selectedProtocol != self.currentProtocol and selectedProtocol != ""):
            # Destroy the old panel so we can free memory up and create the new protocol panel
            if (self.protPanel is not None):
                if (self.currentProtocol == "Superimposition"):
                    self.protPanel.Destroy()
                    del self.protPanel
                elif (self.currentProtocol == "Energy Minimization"):
                    # Check to see if the user has accepted the minimization
                    # If not, ask if they really want to proceed
                    if (self.protPanel.buttonState == "Finalize!"):
                        dlg = wx.MessageDialog(self, "You have not accepted your minimization and the results will be lost if you proceed.  Proceed anyway?", "Minimization Not Accepted", wx.YES_NO | wx.ICON_EXCLAMATION | wx.CENTRE)
                        if (dlg.ShowModal() == wx.ID_NO):
                            logInfo("Go cancelled due to unaccepted minimization")
                            dlg.Destroy()
                            return
                        dlg.Destroy()
                        logInfo("Minimization job rejected")
                    # Try to delete the minimized model if the user doesn't want to accept it
                    try:
                        self.cmd.remove("minimized_view")
                        self.cmd.delete("minimized_view")
                    except:
                        pass
                    self.cmd.label("all", "")
                    self.protPanel.Destroy()
                    del self.protPanel
                elif (self.currentProtocol == "Ensemble Browser"):
                    self.protPanel.Destroy()
                    del self.protPanel
                elif (self.currentProtocol == "Biological Tools"):
                    self.protPanel.Destroy()
                    del self.protPanel
                elif (self.currentProtocol == "Module Manager"):
                    self.protPanel.Destroy()
                    del self.protPanel
                elif (self.currentProtocol == "Ensemble Generation"):
                    # Check to see if the user has accepted the design
                    # If not, ask if they really want to proceed
                    if (self.protPanel.buttonState == "Cancel!"):
                        dlg = wx.MessageDialog(self, "Your ensemble is not finished and the results will be lost if you proceed.  Proceed anyway?", "Design Not Accepted", wx.YES_NO | wx.ICON_EXCLAMATION | wx.CENTRE)
                        if (dlg.ShowModal() == wx.ID_NO):
                            logInfo("Go cancelled due to unfinished ensemblegen job")
                            dlg.Destroy()
                            return
                        logInfo("Ensemblegen job rejected")
                        dlg.Destroy()
                    self.protPanel.Destroy()
                    del self.protPanel
                elif (self.currentProtocol == "Protein Design (Fixbb)"):
                    # Check to see if the user has accepted the design
                    # If not, ask if they really want to proceed
                    if (self.protPanel.buttonState == "Finalize!"):
                        dlg = wx.MessageDialog(self, "You have not accepted your design and the results will be lost if you proceed.  Proceed anyway?", "Design Not Accepted", wx.YES_NO | wx.ICON_EXCLAMATION | wx.CENTRE)
                        if (dlg.ShowModal() == wx.ID_NO):
                            logInfo("Go cancelled due to unaccepted fixbb job")
                            dlg.Destroy()
                            return
                        logInfo("Fixbb job rejected")
                        dlg.Destroy()
                    # Try to delete the design model if the user doesn't want to accept it
                    try:
                        self.cmd.remove("designed_view")
                        self.cmd.delete("designed_view")
                    except:
                        pass
                    self.cmd.label("all", "")
                    self.protPanel.Destroy()
                    del self.protPanel
                elif (self.currentProtocol == "Residue/Ligand Creator"):
                    try:
                        self.cmd.remove("params")
                        self.cmd.delete("params")
                    except:
                        pass
                    self.protPanel.Destroy()
                    del self.protPanel
                elif (self.currentProtocol == "Point Mutations"):
                    try:
                        self.cmd.remove("rotamer_view")
                        self.cmd.delete("rotamer_view")
                    except:
                        pass
                    self.protPanel.Destroy()
                    del self.protPanel
                elif (self.currentProtocol == "Point Mutant Scan"):
                    self.protPanel.Destroy()
                    del self.protPanel
                elif (self.currentProtocol == "Loop Modeling (KIC)"):
                    # Check to see if the user has accepted the loop model
                    # If not, ask if they really want to proceed
                    if (self.protPanel.buttonState == "Finalize!"):
                        dlg = wx.MessageDialog(self, "You have not accepted your loop model and the results will be lost if you proceed.  Proceed anyway?", "KIC Not Accepted", wx.YES_NO | wx.ICON_EXCLAMATION | wx.CENTRE)
                        if (dlg.ShowModal() == wx.ID_NO):
                            logInfo("Go cancelled due to unaccepted KIC job")
                            dlg.Destroy()
                            return
                        logInfo("KIC job rejected")
                        dlg.Destroy()
                    # Try to delete the design model if the user doesn't want to accept it
                    try:
                        self.cmd.remove("kic_view")
                        self.cmd.delete("kic_view")
                    except:
                        pass
                    self.cmd.label("all", "")
                    self.protPanel.Destroy()
                elif (self.currentProtocol == "Docking"):
                    # Check to see if the user has accepted the docking model
                    # If not, ask if they really want to proceed
                    if (self.protPanel.buttonState != "Dock!"):
                        dlg = wx.MessageDialog(self, "You have not accepted your docking model and the results will be lost if you proceed.  Proceed anyway?", "Docking Not Accepted", wx.YES_NO | wx.ICON_EXCLAMATION | wx.CENTRE)
                        if (dlg.ShowModal() == wx.ID_NO):
                            logInfo("Go cancelled due to unaccepted docking job")
                            dlg.Destroy()
                            return
                        logInfo("Docking job rejected")
                        dlg.Destroy()
                    # Try to delete the design model if the user doesn't want to accept it
                    try:
                        self.cmd.remove("dock_view")
                        self.cmd.delete("dock_view")
                    except:
                        pass
                    self.cmd.label("all", "")
                    self.protPanel.Destroy()
                elif (self.currentProtocol == "Structure Prediction (Comparative Modeling)"):
                    # Check to see if the user has accepted the docking model
                    # If not, ask if they really want to proceed
                    if (self.protPanel.buttonState == "Finalize!"):
                        dlg = wx.MessageDialog(self, "You have not accepted your comparative modeling structure and the results will be lost if you proceed.  Proceed anyway?", "Docking Not Accepted", wx.YES_NO | wx.ICON_EXCLAMATION | wx.CENTRE)
                        if (dlg.ShowModal() == wx.ID_NO):
                            logInfo("Go cancelled due to unaccepted comparative modeling job")
                            dlg.Destroy()
                            return
                        logInfo("Comparative modeling job rejected")
                        dlg.Destroy()
                    # Try to delete the design model if the user doesn't want to accept it
                    try:
                        self.cmd.remove("thread_view")
                        self.cmd.delete("thread_view")
                    except:
                        pass
                    self.cmd.label("all", "")
                    self.protPanel.Destroy()
                elif (self.currentProtocol == "Protein Design (MSD)"):
                    self.protPanel.Destroy()
                    del self.protPanel
                elif (self.currentProtocol == "Antibody Modeling"):
                    self.protPanel.Destroy()
                    del self.protPanel
                elif (self.currentProtocol == "Molecular Surfaces"):
                    try:
                        self.cmd.delete("curr_surf_recp")
                    except:
                        pass
                    try:
                        self.cmd.delete("curr_surf_lig")
                    except:
                        pass
                    self.cmd.hide("surface", "all")
                    if (self.parent.Selection.showSurf):
                        self.currentProtocol = "n/a"
                        self.parent.Selection.displaySurfaces()
                    self.protPanel.Destroy()
                    del self.protPanel
                elif (self.currentProtocol == "Flexible Peptide Docking"):
                    self.protPanel.Destroy()
                    del self.protPanel
                else:
                    # Must be a custom module
                    self.protPanel.Destroy()
                    del self.protPanel
            self.currentProtocol = selectedProtocol
            self.seqWin.cannotDelete = False
            # Restart the Rosetta daemon to clear its memory up
            self.parent.restartDaemon()
            logInfo("Changed protocol to " + selectedProtocol)
            try:
                if (selectedProtocol == "Superimposition"):
                    self.protPanel = SuperimpositionPanel(self, self.W, self.H)
                elif (selectedProtocol == "Energy Minimization"):
                    self.protPanel = MinimizationPanel(self, self.W, self.H)
                    self.protPanel.setSelectWin(self.selectWin)
                elif (selectedProtocol == "Ensemble Browser"):
                    self.protPanel = EnsembleBrowserPanel(self, self.W, self.H)
                elif (selectedProtocol == "Biological Tools"):
                    self.protPanel = BioToolsPanel(self, self.W, self.H)
                elif (selectedProtocol == "Module Manager"):
                    self.protPanel = ModuleManagerPanel(self, self.W, self.H)
                elif (selectedProtocol == "Ensemble Generation"):
                    self.protPanel = EnsembleGenPanel(self, self.W, self.H)
                    self.protPanel.setSelectWin(self.selectWin)
                elif (selectedProtocol == "Protein Design (Fixbb)"):
                    self.protPanel = FixbbPanel(self, self.W, self.H)
                    self.protPanel.setSelectWin(self.selectWin)
                elif (selectedProtocol == "Residue/Ligand Creator"):
                    self.protPanel = ResidueCreatorPanel(self, self.W, self.H)
                elif (selectedProtocol == "Point Mutations"):
                    self.protPanel = PointMutationsPanel(self, self.W, self.H)
                    self.protPanel.setSelectWin(self.selectWin)
                elif (selectedProtocol == "Point Mutant Scan"):
                    self.protPanel = PointMutantScanPanel(self, self.W, self.H)
                    self.protPanel.setSelectWin(self.selectWin)
                elif (selectedProtocol == "Loop Modeling (KIC)"):
                    self.protPanel = KICPanel(self, self.W, self.H)
                    self.protPanel.setSelectWin(self.selectWin)
                elif (selectedProtocol == "Docking"):
                    self.protPanel = DockingPanel(self, self.W, self.H)
                    self.protPanel.setSelectWin(self.selectWin)
                elif (selectedProtocol == "Structure Prediction (Comparative Modeling)"):
                    self.protPanel = CompModelPanel(self, self.W, self.H)
                    self.protPanel.setSelectWin(self.selectWin)
                elif (selectedProtocol == "Protein Design (MSD)"):
                    self.protPanel = MSDPanel(self, self.W, self.H)
                    self.protPanel.setSelectWin(self.selectWin)
                elif (selectedProtocol == "Antibody Modeling"):
                    self.protPanel = AntibodyPanel(self, self.W, self.H)
                elif (selectedProtocol == "Molecular Surfaces"):
                    self.protPanel = SurfacesPanel(self, self.W, self.H)
                    self.cmd.hide("surface", "all")
                elif (selectedProtocol == "Flexible Peptide Docking"):
                    self.protPanel = FlexPepDockPanel(self, self.W, self.H)
                    self.protPanel.setSelectWin(self.selectWin)
                elif selectedProtocol == "Pathway Visualization (GeoFold)":
                    self.protPanel = DagViewPanel(self,self.W,self.H)
                elif selectedProtocol == "INDEL Loop Modeler":
                    self.protPanel = INDELmodelPanel(self, self.W, self.H)
                    self.protPanel.setSelectWin(self.selectWin)
                else:
                    # Custom module
                    # Custom modules are aware of PyMOL, the Sequence Window, and the SelectionPanel
                    # They also know about the Protocols Window through the parent variable
                    # They also will receive activation events
                    self.protPanel = self.modules[selectedProtocol].ModulePanel(self, self.modules[selectedProtocol].__file__, self.W, self.H)
                    self.protPanel.setSelectWin(self.selectWin)
                self.protPanel.setSeqWin(self.seqWin)
                self.protPanel.setPyMOL(self.pymol)
                self.protPanel.activate()
                self.seqWin.setProtocolPanel(self.protPanel)
            except:
                self.protPanel = None

    def setSeqWin(self, seqWin):
        self.seqWin = seqWin
        self.protPanel.setSeqWin(seqWin)
        self.seqWin.setProtocolPanel(self.protPanel)

    def setPyMOL(self, pymol):
        self.pymol = pymol
        self.cmd = pymol.cmd
        self.stored = pymol.stored
        self.protPanel.setPyMOL(pymol)

    def setSelectWin(self, selectWin):
        self.selectWin = selectWin

    def activate(self):
        self.cmd.enable("seqsele")
        try:
            self.protPanel.activate()
        except:
            # Probably a custom module that does not have the activate function implemented
            print "WARNING: Could not activate the ProtocolPanel for " + self.currentProtocol

class ProtocolsWin(wx.Frame):
    def __init__(self, W, H, scriptdir):
        if (platform.system() == "Darwin"):
            self.stdwinx = 0; self.stdwiny = 24
        else:
            self.stdwinx = 0; self.stdwiny = 0
        self.stdwinw = 370; self.stdwinh = H-40
        self.screenH = H; self.screenW = W
        winx = self.stdwinx; winy = self.stdwiny
        winw = self.stdwinw; winh = self.stdwinh
        self.scriptdir = scriptdir
        homedir = os.path.expanduser("~")
        # Try to get the save values from the cfg file
        try:
            if (platform.system() == "Windows"):
                f = open(homedir + "\\InteractiveROSETTA\\protwindow.cfg", "r")
            else:
                f = open(homedir + "/.InteractiveROSETTA/protwindow.cfg", "r")
            for aline in f:
                if (aline.find("[OFFSET X]") >= 0):
                    winx = winx + int(aline.split()[len(aline.split())-1])
                elif (aline.find("[OFFSET Y]") >= 0):
                    winy = winy + int(aline.split()[len(aline.split())-1])
                elif (aline.find("[OFFSET WIDTH]") >= 0):
                    winw = winw + int(aline.split()[len(aline.split())-1])
                elif (aline.find("[OFFSET HEIGHT]") >= 0):
                    winh = winh + int(aline.split()[len(aline.split())-1])
            f.close()
        except:
            pass
        if (winx > self.screenW - 100):
            winx = self.stdwinx
        if (winy > self.screenH - 100):
            winy = self.stdwiny
        # Catch bad cached sizes
        if (winw < 200):
            winw = self.stdwinw
        if (winh < 200):
            winh = self.stdwinh
        # Maybe the screen resolution has changed and the saved dimensions put the windows in
        # weird places, so default them to better positions and the user can change them later
        #if (winw < 350):
        #    winw = 370
        #elif (winw > W):
        #    winw = 370
        #if (winx < 0):
        #    winx = 0
        #elif (winx > W-winw):
        #    winx = W-winw
        #if (winh > H - 40):
        #    winh = H - 40
        #if (winy < 0):
        #    winy = 0
        #elif (winy > H-winh):
        #    winh = H-40
        wx.Frame.__init__(self, None, -1, "InteractiveROSETTA - Protocols", size=(winw, winh))
        self.SetPosition((winx, winy))
        self.SetBackgroundColour("#333333")
        self.SetSizeHints(330, 560, 370, H)
        self.SetIcon(icon.GetIcon())
        self.Show()

        self.Protocols = ProtocolsPanel(self, winw, winh)
        self.Selection = SelectPanel(self, winw, winh)
        self.Protocols.setSelectWin(self.Selection)
        self.Selection.setProtPanel(self.Protocols)

        self.saveTimer = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, self.saveWindowData, self.saveTimer)
        self.Bind(wx.EVT_SIZE, self.windowGeometryChange)
        self.Bind(wx.EVT_MOTION, self.windowGeometryChange)
        self.Bind(wx.EVT_ACTIVATE, self.focusEvent)

        # Start the Rosetta daemon that will run in the background looking for job input files generated
        # by the main GUI
        # It could be the case that the user was in the middle of a protocol, then quits suddenly so the
        # daemon doesn't terminate itself because it's in the middle of a protocol
        # The the user restarts the GUI before the daemon has a chance to finish the protocol
        # This checks to see if the daemon is already active and doesn't spawn a new one if one is already
        # running
        stillrunning = False
        count = 0
        for proc in psutil.process_iter():
            try:
                if (platform.system() == "Windows"):
                    if (len(proc.cmdline()) >= 3 and proc.cmdline()[0].find("python") >= 0 and proc.cmdline()[2].find("from multiprocessing.forking") >= 0):
                        stillrunning = True
                        break
                else:
                    # On Unix systems you just have to make sure two instances of python are running
                    # because there isn't any forking information in the daemon's instance of python
                    if (len(proc.cmdline()) >= 2 and proc.cmdline()[0].find("python") >= 0 and proc.cmdline()[1].find("InteractiveROSETTA.py") >= 0):
                        count = count + 1
            except:
                # In Windows it will crash if you try to read process information for the Administrator
                # Doesn't matter though since InteractiveROSETTA is run by a non-Administrator
                # But we need to catch these errors since we don't know which processes are admin ones
                pass
        if (platform.system() != "Windows" and count == 2):
            stillrunning = True
        if (not(stillrunning)):
            print "Starting Rosetta protocol daemon..."
            #if (platform.system() == "Linux"):
                #self.daemon_process = subprocess.Popen(args=["cd " + self.scriptdir + "/scripts" + "; python", "daemon.py"], shell=False)
            #else:
            self.daemon_process = multiprocessing.Process(target=daemonLoop)
            self.daemon_process.start()

    def restartDaemon(self):
        # This function kills the current daemon process and starts a new one
        # This function is called whenever the protocol panel changes to a different protocol
        # This function is necessary because on Windows PyRosetta can start to use well over
        # 4GB of memory if the daemon has loaded memory from multiple protocols
        # Killing the daemon and restarting it clears up that memory so the user's computer
        # doesn't slow down as it moves a lot of stuff into swap space
        savedir = os.getcwd()
        os.chdir(self.scriptdir)
        self.daemon_process.terminate()
        print "Restarting Rosetta protocol daemon..."
        self.daemon_process = multiprocessing.Process(target=daemonLoop)
        self.daemon_process.start()
        os.chdir(savedir)

    def focusEvent(self, event):
        if (event.GetActive()):
            # Protocols read selection information from the sequence window, so update the sequence window
            # If we're going from PyMOL->Protocols, since updates usually happen on the sequence window focus event
            self.seqWin.selectUpdate(False) # Update PyMOL changes in sequence
            self.Protocols.activate()
        event.Skip()

    def setSeqWin(self, seqWin):
        self.seqWin = seqWin
        self.Protocols.setSeqWin(seqWin)
        self.Selection.setSeqWin(seqWin)

    def windowGeometryChange(self, event):
        # This function starts a timer that will write out the size and position of this window to a cfg file
        # so the orientation is saved and can be loaded the next time InteractiveROSETTA is started
        if (not(self.saveTimer.IsRunning())):
            self.saveTimer.Start(5000)
        # We have to do some finagling if this window gets resized because the minimum height
        # of the specific protocol panel needs to be at least 300
        # If the window is shrunk, the protocol panel will shrink down to 300 but then the select
        # panel needs to shrink
        (w, h) = self.GetSize()
        return
        if (h > 560):
            try:
                self.Protocols.protPanel.SetSize((w, h-330))
                self.Protocols.protPanel.SetScrollbars(1, 1, 320, 800)
            except:
                pass
            (x, y) = self.Selection.GetPosition()
            self.Selection.SetPosition((x, h-270))
            self.Selection.SetSize((w-20, self.Selection.GetSize()[1]))
        event.Skip()

    def saveWindowData(self, event):
        homedir = os.path.expanduser("~")
        data = []
        try:
            if (platform.system() == "Windows"):
                f = open(homedir + "\\InteractiveROSETTA\\protwindow.cfg", "r")
            else:
                f = open(homedir + "/.InteractiveROSETTA/protwindow.cfg", "r")
            for aline in f:
                data.append(aline)
            f.close()
        except:
            pass
        if (platform.system() == "Windows"):
            f = open(homedir + "\\InteractiveROSETTA\\protwindow.cfg", "w")
        else:
            f = open(homedir + "/.InteractiveROSETTA/protwindow.cfg", "w")
        itemsFound = [False, False, False, False] # [offX, offY, offW, offH]
        (x, y) = self.GetPosition()
        (w, h) = self.GetSize()
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
            else:
                f.write(aline)
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
        f.close()

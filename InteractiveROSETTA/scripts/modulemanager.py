import wx
import wx.grid
import wx.lib.scrolledpanel
import wx.lib.dialogs
import os
import os.path
import platform
import webbrowser
import shutil
from contextlib import closing
from zipfile import ZipFile, ZIP_DEFLATED
from tools import resizeTextControlForUNIX
from tools import logInfo

class ModuleManagerPanel(wx.lib.scrolledpanel.ScrolledPanel):
    def __init__(self, parent, W, H):
        if (platform.system() == "Windows"):
            wx.lib.scrolledpanel.ScrolledPanel.__init__(self, parent, id=-1, pos=(10, 60), size=(340, H-330), name="ProtManageModules")
            winh = H-330
        else:
            wx.lib.scrolledpanel.ScrolledPanel.__init__(self, parent, id=-1, pos=(10, 60), size=(340, H-290), name="ProtManageModules")
            winh = H-290
        self.SetBackgroundColour("#333333")
        self.parent = parent
        
        if (platform.system() == "Windows"):
            self.lblProt = wx.StaticText(self, -1, "Module Manager", (25, 15), (270, 25), style=wx.ALIGN_CENTRE)
            self.lblProt.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblProt = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/modulemanager/lblModule.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 15), size=(320, 25))
        else:
            self.lblProt = wx.StaticText(self, -1, "Module Manager", pos=(90, 15), style=wx.ALIGN_CENTRE)
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
            self.lblInst = wx.StaticText(self, -1, "Create module packages for distribution\nand install pre-made packages.\nSee the documentation for how to set up a package.", (0, 45), (320, 25), wx.ALIGN_CENTRE)
            self.lblInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
        elif (platform.system() == "Darwin"):
            self.lblInst = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/modulemanager/lblInstModule.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 45), size=(320, 60))
        else:
            self.lblInst = wx.StaticText(self, -1, "Create module packages for distribution\nand install pre-made packages.\nSee the documentation for how to set up a package.", (40, 45), (320, 25), style=wx.ALIGN_CENTRE)
            self.lblInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
            resizeTextControlForUNIX(self.lblInst, 0, self.GetSize()[0])
        self.lblInst.SetForegroundColour("#FFFFFF")
        
        if (platform.system() == "Windows"):
            self.lblInstall = wx.StaticText(self, -1, "Install Module:", (10, 103), (195, 20), wx.ALIGN_LEFT)
            self.lblInstall.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblInstall = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/modulemanager/lblInstallModule.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(10, 103), size=(195, 20))
        else:
            self.lblInstall = wx.StaticText(self, -1, "Install Module:", (10, 103), style=wx.ALIGN_LEFT)
            self.lblInstall.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.lblInstall.SetForegroundColour("#FFFFFF")
        if (platform.system() == "Darwin"):
            self.btnInstall = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/modulemanager/btnInstall.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(210, 100), size=(100, 25))
        else:
            self.btnInstall = wx.Button(self, id=-1, label="Install", pos=(210, 100), size=(100, 25))
            self.btnInstall.SetForegroundColour("#000000")
            self.btnInstall.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnInstall.Bind(wx.EVT_BUTTON, self.installClick)
        self.btnInstall.SetToolTipString("Install an IRM module package")
        
        if (platform.system() == "Windows"):
            self.lblUninstall = wx.StaticText(self, -1, "Uninstall Module:", (10, 133), (195, 20), wx.ALIGN_LEFT)
            self.lblUninstall.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblUninstall = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/modulemanager/lblUninstallModule.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(10, 133), size=(195, 20))
        else:
            self.lblUninstall = wx.StaticText(self, -1, "Uninstall Module:", (10, 133), style=wx.ALIGN_LEFT)
            self.lblUninstall.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.lblUninstall.SetForegroundColour("#FFFFFF")
        if (platform.system() == "Darwin"):
            self.btnUninstall = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/modulemanager/btnUninstall.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(210, 130), size=(100, 25))
        else:
            self.btnUninstall = wx.Button(self, id=-1, label="Uninstall", pos=(210, 130), size=(100, 25))
            self.btnUninstall.SetForegroundColour("#000000")
            self.btnUninstall.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnUninstall.Bind(wx.EVT_BUTTON, self.uninstallClick)
        self.btnUninstall.SetToolTipString("Uninstall a module")
        
        if (platform.system() == "Windows"):
            self.lblCreate = wx.StaticText(self, -1, "Create Module:", (10, 163), (195, 20), wx.ALIGN_LEFT)
            self.lblCreate.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblCreate = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/modulemanager/lblCreateModule.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(10, 163), size=(195, 20))
        else:
            self.lblCreate = wx.StaticText(self, -1, "Create Module:", (10, 163), style=wx.ALIGN_LEFT)
            self.lblCreate.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.lblCreate.SetForegroundColour("#FFFFFF")
        if (platform.system() == "Darwin"):
            self.btnCreate = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/modulemanager/btnCreate.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(210, 160), size=(100, 25))
        else:
            self.btnCreate = wx.Button(self, id=-1, label="Create", pos=(210, 160), size=(100, 25))
            self.btnCreate.SetForegroundColour("#000000")
            self.btnCreate.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnCreate.Bind(wx.EVT_BUTTON, self.createClick)
        self.btnCreate.SetToolTipString("Create an IRM module for distribution")
        
        scrollh = self.btnCreate.GetPosition()[1] + self.btnCreate.GetSize()[1] + 5
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

    def scrolled(self, event):
        self.winscrollpos = self.GetScrollPos(wx.VERTICAL)
        event.Skip()

    def activate(self):
        self.Scroll(0, self.winscrollpos)
    
    def setSeqWin(self, seqWin):
        self.seqWin = seqWin
    
    def setPyMOL(self, cmd):
        pass
    
    def installClick(self, event):
        # Browse for the package and install it
        # The package is an IRM file, but it is basically just a zip file
        dlg = wx.FileDialog(
            self, message="Choose a File",
            defaultDir=self.seqWin.cwd,
            defaultFile="",
            wildcard="InteractiveROSETTA Modules (*.irm)|*.irm",
            style=wx.OPEN | wx.CHANGE_DIR)
        if (dlg.ShowModal() != wx.ID_OK):
            return
        module = dlg.GetPath()
        # Prevent the user from trying to package up the template
        if (module.split(".irm")[0].endswith("template")):
            dlg = wx.MessageDialog(self, "The template module is reserved, please do not attempt to overwrite it.", "Operation Forbidden", wx.OK | wx.ICON_ERROR | wx.CENTRE)
            dlg.ShowModal()
            dlg.Destroy()
            return
        if (platform.system() == "Windows"):
            packagecode = module[module.rfind("\\")+1:].split(".irm")[0]
        else:
            packagecode = module[module.rfind("/")+1:].split(".irm")[0]
        # Let's see if this module is already installed
        home = os.path.expanduser("~")
        if (platform.system() == "Windows"):
            if (os.path.exists(home + "/InteractiveROSETTA/modules/" + packagecode)):
                dlg = wx.MessageDialog(self, "This module is already installed.  Do you want to overwrite it?", "Module Already Installed", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
                if (dlg.ShowModal() == wx.ID_NO):
                    return
                dlg.Destroy()
                # Delete it
                shutil.rmtree(home + "/InteractiveROSETTA/modules/" + packagecode, ignore_errors=True)
            os.mkdir(home + "/InteractiveROSETTA/modules/" + packagecode)
        else:
            if (os.path.exists(home + "/.InteractiveROSETTA/modules/" + packagecode)):
                dlg = wx.MessageDialog(self, "This module is already installed.  Do you want to overwrite it?", "Module Already Installed", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
                if (dlg.ShowModal() == wx.ID_NO):
                    return
                dlg.Destroy()
                # Delete it
                shutil.rmtree(home + "/.InteractiveROSETTA/modules/" + packagecode, ignore_errors=True)
            os.mkdir(home + "/.InteractiveROSETTA/modules/" + packagecode)
        # Unpack the irm package to the ~/InteractiveROSETTA/modules directory
        fin = open(module, 'rb')
        z = ZipFile(fin)
        for name in z.namelist():
            if (platform.system() == "Windows"):
                outpath = home + "/InteractiveROSETTA/modules/" + packagecode
            else:
                outpath = home + "/.InteractiveROSETTA/modules/" + packagecode
            z.extract(name, outpath)
        fin.close()
        # Remove the "server" directory if it exists since that belongs on the server, not the client
        try:
            if (platform.system() == "Windows"):
                shutil.rmtree(home + "/InteractiveROSETTA/modules/" + packagecode + "/server", ignore_errors=True)
            else:
                shutil.rmtree(home + "/.InteractiveROSETTA/modules/" + packagecode + "/server", ignore_errors=True)
        except:
            pass
        # Is there a license?  If so, display it.
        if (platform.system() == "Windows"):
            if (os.path.isfile(home + "/InteractiveROSETTA/modules/" + packagecode + "/license")):
                fin = open(home + "/InteractiveROSETTA/modules/" + packagecode + "/license")
                licensetext = ""
                for aline in fin:
                    licensetext += aline
                fin.close()
                dlg1 = wx.lib.dialogs.ScrolledMessageDialog(None, licensetext, packagecode + " License")
                dlg1.ShowModal()
                dlg1.Destroy()
                dlg = wx.MessageDialog(None, "Do you accept the license agreement?", packagecode + " License", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
                if (dlg.ShowModal() == wx.ID_NO):
                    # Delete everything
                    shutil.rmtree(home + "/InteractiveROSETTA/modules/" + packagecode)
                    return
        else:
            if (os.path.isfile(home + "/.InteractiveROSETTA/modules/" + packagecode + "/license")):
                fin = open(home + "/.InteractiveROSETTA/modules/" + packagecode + "/license")
                licensetext = ""
                for aline in fin:
                    licensetext += aline
                fin.close()
                dlg1 = wx.lib.dialogs.ScrolledMessageDialog(None, licensetext, packagecode + " License")
                dlg1.ShowModal()
                dlg1.Destroy()
                dlg = wx.MessageDialog(None, "Do you accept the license agreement?", packagecode + " License", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
                if (dlg.ShowModal() == wx.ID_NO):
                    # Delete everything
                    shutil.rmtree(home + "/.InteractiveROSETTA/modules/" + packagecode)
                    return
        # Is there a message from the developer?  If so, display it
        if (platform.system() == "Windows"):
            if (os.path.isfile(home + "/InteractiveROSETTA/modules/" + packagecode + "/message")):
                fin = open(home + "/InteractiveROSETTA/modules/" + packagecode + "/message")
                msgtext = ""
                for aline in fin:
                    msgtext += aline
                fin.close()
                dlg = wx.MessageDialog(self, msgtext, "Message From the Developer", wx.OK | wx.ICON_ERROR | wx.CENTRE)
                dlg.ShowModal()
                dlg.Destroy()
        else:
            if (os.path.isfile(home + "/.InteractiveROSETTA/modules/" + packagecode + "/message")):
                fin = open(home + "/.InteractiveROSETTA/modules/" + packagecode + "/message")
                msgtext = ""
                for aline in fin:
                    msgtext += aline
                fin.close()
                dlg = wx.MessageDialog(self, msgtext, "Message From the Developer", wx.OK | wx.ICON_ERROR | wx.CENTRE)
                dlg.ShowModal()
                dlg.Destroy()
        # Reload the modules to get the updated list of protocols
        self.parent.readModules()
    
    def uninstallClick(self, event):
        # Ask the user which module they want to uninstall
        if (platform.system() == "Windows"):
            dlg = wx.DirDialog(self, "Select the module to uninstall",
                            defaultPath=os.path.expanduser("~") + "/InteractiveROSETTA/modules",
                            style=wx.DD_DEFAULT_STYLE 
                            | wx.DD_DIR_MUST_EXIST)
        else:
            dlg = wx.DirDialog(self, "Select the module to uninstall",
                            defaultPath=os.path.expanduser("~") + "/.InteractiveROSETTA/modules",
                            style=wx.DD_DEFAULT_STYLE 
                            | wx.DD_DIR_MUST_EXIST)
        if (dlg.ShowModal() != wx.ID_OK):
            return
        module = dlg.GetPath().strip()
        # Make sure that the folder selected was an actual module and not some random folder
        if (platform.system() == "Windows"):
            path = module[0:module.rfind("\\")]
            expectedloc = os.path.expanduser("~") + "\\InteractiveROSETTA\\modules"
        else:
            path = module[0:module.rfind("/")]
            expectedloc = os.path.expanduser("~") + "/.InteractiveROSETTA/modules"
        if (path.strip() != expectedloc.strip()):
            dlg = wx.MessageDialog(self, "That is not an InteractiveROSETTA module.  Select a folder in the InteractiveROSETTA/modules directory.", "Invalid Directory Chosen", wx.OK | wx.ICON_ERROR | wx.CENTRE)
            dlg.ShowModal()
            dlg.Destroy()
            return
        # Prevent the user from uninstalling the template set
        if (module.strip().endswith("template")):
            dlg = wx.MessageDialog(self, "The template module is reserved, please do not attempt to remove it.", "Operation Forbidden", wx.OK | wx.ICON_ERROR | wx.CENTRE)
            dlg.ShowModal()
            dlg.Destroy()
            return
        # Make sure the user wants to do this
        dlg = wx.MessageDialog(self, "Are you sure you want to remove this module?", "Confirm Uninstallation", wx.YES_NO | wx.ICON_QUESTION | wx.CENTRE)
        if (dlg.ShowModal() == wx.ID_NO):
            return
        dlg.Destroy()
        # Remove it
        shutil.rmtree(module, ignore_errors=True)
        # Update the menu
        self.parent.readModules()
    
    def createClick(self, event):
        # Ask the user for the directory containing the panel __init__.py script
        if (platform.system() == "Windows"):
            dlg = wx.DirDialog(self, "Select the directory containing __init__.py for your ModulePanel",
                            defaultPath=os.path.expanduser("~") + "/InteractiveROSETTA/modules",
                            style=wx.DD_DEFAULT_STYLE 
                            | wx.DD_DIR_MUST_EXIST
                            | wx.DD_CHANGE_DIR)
        else:
            dlg = wx.DirDialog(self, "Select the directory containing __init__.py for your ModulePanel",
                    defaultPath=os.path.expanduser("~") + "/.InteractiveROSETTA/modules",
                    style=wx.DD_DEFAULT_STYLE 
                    | wx.DD_DIR_MUST_EXIST
                    | wx.DD_CHANGE_DIR)
        if (dlg.ShowModal() != wx.ID_OK):
            return
        packageloc = dlg.GetPath().strip()
        # Prevent the user from trying to package up the template
        if (packageloc.strip().endswith("template")):
            dlg = wx.MessageDialog(self, "The template module is reserved, please do not attempt to package it.", "Operation Forbidden", wx.OK | wx.ICON_ERROR | wx.CENTRE)
            dlg.ShowModal()
            dlg.Destroy()
            return
        if (platform.system() == "Windows"):
            packagecode = packageloc[packageloc.rfind("\\")+1:].strip()
        else:
            packagecode = packageloc[packageloc.rfind("/")+1:].strip()
        # Let's make sure the __init__.py file is there
        if (not(os.path.isfile(packageloc + "/__init__.py"))):
            dlg = wx.MessageDialog(self, "There is no __init__.py file in that directory!", "Invalid Directory Chosen", wx.OK | wx.ICON_ERROR | wx.CENTRE)
            dlg.ShowModal()
            dlg.Destroy()
            return
        # Ask the user for the directory to save the package in
        dlg = wx.DirDialog(self, "Specify the save location",
                           defaultPath=self.seqWin.cwd,
                           style=wx.DD_DEFAULT_STYLE 
                           | wx.DD_CHANGE_DIR)
        if (dlg.ShowModal() != wx.ID_OK):
            return
        # The IRM package is really just a zip file containing the specified folder and all its
        # contents.  It is a zip file so Windows user's can examine its contents easily
        with closing(ZipFile(dlg.GetPath() + "/" + packagecode + ".irm", "w", ZIP_DEFLATED)) as z:
            for root, dirs, files in os.walk(packageloc):
                #NOTE: ignore empty directories
                for fn in files:
                    absfn = os.path.join(root, fn)
                    zfn = absfn[len(packageloc)+len(os.sep):] #XXX: relative path
                    z.write(absfn, zfn)
        dlg = wx.MessageDialog(self, "Your " + packagecode + " module has been created successfully!", "Module Successfully Created", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
        dlg.ShowModal()
        dlg.Destroy()
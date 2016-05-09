import wx
import wx.grid
import wx.lib.scrolledpanel
import os
import os.path
import sys
import platform
import glob
import webbrowser
try:
    # It may have been installed locally in the sandbox
    if (platform.system() == "Windows"):
        sys.path.append(os.path.expanduser("~") + "\\InteractiveROSETTA\\molfile2params")
    else:
        sys.path.append(os.path.expanduser("~") + "/.InteractiveROSETTA/molfile2params")
    import molfile_to_params
except:
    try:
        if (platform.system() == "Windows"):
            shutil.rmtree(os.path.expanduser("~") + "\\InteractiveROSETTA\\molfile2params")
        else:
            shutil.rmtree(os.path.expanduser("~") + "/.InteractiveROSETTA/molfile2params")
    except:
        pass
    dlg = wx.MessageDialog(None, "InteractiveROSETTA was unable to import molfile2params!\nInteractiveROSETTA attempted to remove the previous molfile2params install.  Run InteractiveROSETTA again to reinstall it automatically.  If it still does not import, please contact the developer.", "Molfile2Params Missing", wx.OK | wx.ICON_ERROR | wx.CENTRE)
    dlg.ShowModal()
    dlg.Destroy()
    exit()
from tools import *

class ResidueCreatorPanel(wx.lib.scrolledpanel.ScrolledPanel):
    def __init__(self, parent, W, H):
        #if (platform.system() == "Windows"):
        wx.lib.scrolledpanel.ScrolledPanel.__init__(self, parent, id=-1, pos=(10, 60), size=(340, H-330), name="ProtSuperimposition")
        winh = H-330
        #else:
            #wx.lib.scrolledpanel.ScrolledPanel.__init__(self, parent, id=-1, pos=(10, 60), size=(340, H-290), name="ProtSuperimposition")
            #winh = H-290
        self.SetBackgroundColour("#333333")
        self.parent = parent
        self.sizer = wx.GridBagSizer(5, 5)
        #self.sb = wx.ScrollBar(self)
        #self.sb.SetScrollRange(10)
        #print self.sb.GetScrollRange(0)
        
        if (platform.system() == "Windows"):
            self.lblProt = wx.StaticText(self, -1, "Residue/Ligand Creator", (25, 15), (270, 25), style=wx.ALIGN_CENTRE)
            self.lblProt.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblProt = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/advresiduecreator/lblLigand.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(25, 15), size=(270, 25))
        else:
            self.lblProt = wx.StaticText(self, -1, "Residue/Ligand Creator", pos=(60, 15), style=wx.ALIGN_CENTRE)
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
            self.lblInst = wx.StaticText(self, -1, "Upload a PDB file containing unrecognized atoms.\nThen select unrecognized residue types to parameterize.", (0, 45), (320, 25), wx.ALIGN_CENTRE)
            self.lblInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
        elif (platform.system() == "Darwin"):
            self.lblInst = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/advresiduecreator/lblInstLigand.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 45), size=(320, 50))
        else:
            self.lblInst = wx.StaticText(self, -1, "Upload a PDB file containing unrecognized atoms.\nThen select unrecognized residue types to parameterize.", pos=(0, 45), style=wx.ALIGN_CENTRE)
            self.lblInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
            resizeTextControlForUNIX(self.lblInst, 0, self.GetSize()[0])
        self.lblInst.SetForegroundColour("#FFFFFF")
        
        self.lblMOL2 = wx.StaticText(self, -1, "None Uploaded", pos=(10, 103), size=(180, 25), style=wx.ALIGN_CENTRE)
        self.lblMOL2.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        if (platform.system() == "Linux"):
            resizeTextControlForUNIX(self.lblMOL2, 10, 180)
        self.lblMOL2.SetForegroundColour("#FFFFFF")
        if (platform.system() == "Darwin"):
            self.btnLoad = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/advresiduecreator/btnLoadMOL2.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(200, 100), size=(110, 25))
        else:
            self.btnLoad = wx.Button(self, id=-1, label="Load .MOL2", pos=(200, 100), size=(110, 25))
            self.btnLoad.SetForegroundColour("#000000")
            self.btnLoad.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnLoad.Bind(wx.EVT_BUTTON, self.loadMOL2)
        self.btnLoad.SetToolTipString("Load a .mol2 file containing a ligand/NCAA")
        
        if (platform.system() == "Darwin"):
            self.btnType = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/advresiduecreator/btnType_Ligand.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(40, 140), size=(100, 25))
        else:
            self.btnType = wx.Button(self, id=-1, label="Ligand", pos=(40, 140), size=(100, 25))
            self.btnType.SetForegroundColour("#000000")
            self.btnType.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnType.Bind(wx.EVT_BUTTON, self.typeToggle)
        self.btnType.SetToolTipString("Uploaded .mol2 file represents a ligand")
        if (platform.system() == "Darwin"):
            self.btnCreate = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/advresiduecreator/btnCreate.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(180, 140), size=(100, 25))
        else:
            self.btnCreate = wx.Button(self, id=-1, label="Create!", pos=(180, 140), size=(100, 25))
            self.btnCreate.SetForegroundColour("#000000")
            self.btnCreate.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnCreate.Bind(wx.EVT_BUTTON, self.createParams)
        self.btnCreate.Disable()
        self.btnCreate.SetToolTipString("Parameterize the uploaded .mol2 structure")
        self.paramsType = "Ligand"
        self.paramsAtoms = []
        
        self.grdParamsAtoms = wx.grid.Grid(self)
        self.grdParamsAtoms.CreateGrid(0, 1)
        if (winh-265 > 200):
            self.grdParamsAtoms.SetSize((320, winh-265))
        else:
            self.grdParamsAtoms.SetSize((320, 200))
        self.grdParamsAtoms.SetPosition((0, 175))
        self.grdParamsAtoms.SetLabelFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.grdParamsAtoms.DisableDragColSize()
        self.grdParamsAtoms.DisableDragRowSize()
        self.grdParamsAtoms.SetColLabelValue(0, "Assigned Type")
        self.grdParamsAtoms.Bind(wx.grid.EVT_GRID_CELL_LEFT_CLICK, self.gridClick)
        
        ypos = self.grdParamsAtoms.GetPosition()[1] + self.grdParamsAtoms.GetSize()[1] + 10
        if (platform.system() == "Darwin"):
            self.atomMenu = wx.ComboBox(self, pos=(10, ypos), size=(140, 25), choices=[], style=wx.CB_READONLY)
        else:
            self.atomMenu = wx.ComboBox(self, pos=(10, ypos), size=(140, 25), choices=[], style=wx.CB_READONLY | wx.CB_SORT)
        self.atomMenu.Bind(wx.EVT_COMBOBOX, self.atomMenuSelect)
        self.atomMenu.SetToolTipString("Select .mol2 atoms for editing")
        if (platform.system() == "Darwin"):
            self.typeMenu = wx.ComboBox(self, pos=(170, ypos), size=(140, 25), choices=[], style=wx.CB_READONLY)
        else:
            self.typeMenu = wx.ComboBox(self, pos=(170, ypos), size=(140, 25), choices=[], style=wx.CB_READONLY | wx.CB_SORT)
        self.typeMenu.Bind(wx.EVT_COMBOBOX, self.typeMenuSelect)
        self.typeMenu.SetToolTipString("Change the parameterization type for the selected atom")
        # Useful dictionary of atom elements mapped to available types
        self.atomtypes = {}
        self.atomtypes["H"] = ["Hpol", "Haro", "Hapo"]
        self.atomtypes["C"] = ["CH3", "CH2", "CH1", "aroC", "CNH2", "COO", "CAbb", "CObb"]
        self.atomtypes["N"] = ["Nlys", "NH2O", "Ntrp", "Nhis", "Npro", "Nbb"]
        self.atomtypes["O"] = ["OH", "Oaro", "OOC", "OCbb", "ONH2"]
        self.atomtypes["S"] = ["S"]
        self.atomtypes["P"] = ["Phos"]
        self.atomtypes["F"] = ["F"]
        self.atomtypes["CL"] = ["Cl"]
        self.atomtypes["BR"] = ["Br"]
        self.atomtypes["I"] = ["I"]
        self.atomtypes["NA"] = ["Na1p"]
        self.atomtypes["K"] = ["K1p"]
        self.atomtypes["MG"] = ["Mg2p"]
        self.atomtypes["FE"] = ["Fe3p"]
        self.atomtypes["CA"] = ["Ca2p"]
        self.atomtypes["ZN"] = ["Zn2p"]
        
        if (platform.system() == "Darwin"):
            self.lblNterm = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/advresiduecreator/lblNterm.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(10, ypos+35), size=(40, 25))
        else:
            self.lblNterm = wx.StaticText(self, -1, "Nterm:", pos=(10, ypos+35), size=(40, 25))
            self.lblNterm.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
            self.lblNterm.SetForegroundColour("#FFFFFF")
        self.NtermMenu = wx.ComboBox(self, pos=(60, ypos+30), size=(90, 25), choices=[], style=wx.CB_READONLY)
        self.NtermMenu.Bind(wx.EVT_COMBOBOX, self.atomMenuSelect)
        self.NtermMenu.SetToolTipString("N-terminus atom for NCAAs")
        if (platform.system() == "Darwin"):
            self.lblCterm = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/advresiduecreator/lblCterm.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(170, ypos+35), size=(40, 25))
        else:
            self.lblCterm = wx.StaticText(self, -1, "Cterm:", pos=(170, ypos+35), size=(40, 25))
            self.lblCterm.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
            self.lblCterm.SetForegroundColour("#FFFFFF")
        if (platform.system() == "Darwin"):
            self.CtermMenu = wx.ComboBox(self, pos=(220, ypos+30), size=(90, 25), choices=[], style=wx.CB_READONLY)
        else:
            self.CtermMenu = wx.ComboBox(self, pos=(220, ypos+30), size=(90, 25), choices=[], style=wx.CB_READONLY | wx.CB_SORT)
        self.CtermMenu.Bind(wx.EVT_COMBOBOX, self.typeMenuSelect)
        self.CtermMenu.SetToolTipString("C-terminus atom for NCAAs")
        self.NtermMenu.Disable()
        self.CtermMenu.Disable()
        
        if (platform.system() == "Darwin"):
            self.lblCode = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/advresiduecreator/lblCode.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(10, ypos+65), size=(40, 25))
        else:
            self.lblCode = wx.StaticText(self, -1, "Code:", pos=(10, ypos+65), size=(40, 25))
            self.lblCode.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
            self.lblCode.SetForegroundColour("#FFFFFF")
        self.txtCode = wx.TextCtrl(self, -1, pos=(50, ypos+60), size=(40, 25))
        self.txtCode.SetValue("UNK")
        self.txtCode.SetToolTipString("Three-letter amino acid code for the ligand/NCAA")
        if (platform.system() == "Darwin"):
            self.btnAdd = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/advresiduecreator/btnAddDB.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(110, ypos+60), size=(200, 25))
        else:
            self.btnAdd = wx.Button(self, id=-1, label="Add to Database", pos=(110, ypos+60), size=(200, 25))
            self.btnAdd.SetForegroundColour("#000000")
            self.btnAdd.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.btnAdd.Bind(wx.EVT_BUTTON, self.addToDB)
        self.btnAdd.Disable()
        self.btnAdd.SetToolTipString("Add the ligand/NCAA to the Rosetta database with the selected parameters")
        
        if (platform.system() == "Windows"):
            self.lblLine = wx.StaticText(self, -1, "==========================", (0, ypos+90), (320, 20), wx.ALIGN_CENTRE)
            self.lblLine.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.NORMAL))
        elif (platform.system() == "Darwin"):
            self.lblLine = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/advresiduecreator/lblLine.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, ypos+90), size=(320, 20))
        else:
            self.lblLine = wx.StaticText(self, -1, "==========================", (0, ypos+90), style=wx.ALIGN_CENTRE)
            self.lblLine.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.NORMAL))
            resizeTextControlForUNIX(self.lblLine, 20, 120)
        self.lblLine.SetForegroundColour("#FFFFFF")
        
        self.removeMenu = wx.ComboBox(self, pos=(10, ypos+110), size=(90, 25), choices=[], style=wx.CB_READONLY)
        self.removeMenu.Bind(wx.EVT_COMBOBOX, self.resMenuSelect)
        self.removeMenu.SetToolTipString("Select residues already parameterized for removal")
        if (platform.system() == "Darwin"):
            self.btnRemove = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/advresiduecreator/btnRemoveDB.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(110, ypos+110), size=(200, 25))
        else:
            self.btnRemove = wx.Button(self, id=-1, label="Remove from DB", pos=(110, ypos+110), size=(200, 25))
            self.btnRemove.SetForegroundColour("#000000")
            self.btnRemove.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnRemove.Bind(wx.EVT_BUTTON, self.removeParams)
        self.btnRemove.SetToolTipString("Uploaded .mol2 file represents a ligand")
        
        goToSandbox("params")
        paramsfiles = glob.glob("*.fa.params")
        paramslist = []
        for param in paramsfiles:
            if (param != "HOH.fa.params"):
                paramslist.append(param.split(".fa.params")[0])
        self.removeMenu.AppendItems(paramslist)
        
        #self.SetSizerAndFit(self.sizer)
        #self.SetupScrolling()
        scrollh = self.btnRemove.GetPosition()[1] + self.btnRemove.GetSize()[1] + 5
        self.SetScrollbars(1, 1, 320, scrollh)
        self.grdParamsAtoms.SetColSize(0, int(self.grdParamsAtoms.GetSize()[0] / 2))
        self.grdParamsAtoms.SetRowLabelSize(int(self.grdParamsAtoms.GetSize()[0] / 2))
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
            browser.open(self.parent.parent.scriptdir + "/help/ligand.html")
        else:
            webbrowser.open(self.parent.parent.scriptdir + "/help/ligand.html")

    def setSeqWin(self, seqWin):
        self.seqWin = seqWin
        
    def setPyMOL(self, pymol):
        self.pymol = pymol
        self.cmd = pymol.cmd
        self.stored = pymol.stored

    def gridClick(self, event):
        self.atomMenu.SetSelection(event.GetRow())
        self.atomMenuSelect(event)
        event.Skip()

    def scrolled(self, event):
        self.winscrollpos = self.GetScrollPos(wx.VERTICAL)
        event.Skip()

    def activate(self):
        self.Scroll(0, self.winscrollpos)

    def loadMOL2(self, event):
        # Get the structure from a MOL2 file and load it into PyMOL
        logInfo("Load MOL2 button clicked")
        dlg = wx.FileDialog(
            self, message="Choose a File",
            defaultDir=self.seqWin.cwd,
            defaultFile="",
            wildcard="MOL2 Files (*.mol2)|*.mol2",
            style=wx.OPEN | wx.CHANGE_DIR)
        if (dlg.ShowModal() == wx.ID_OK):
            paths = dlg.GetPaths()
            # Change cwd to the last opened file
            if (platform.system() == "Windows"):
                lastDirIndx = paths[len(paths)-1].rfind("\\")
            else:
                lastDirIndx = paths[len(paths)-1].rfind("/")
            self.seqWin.cwd = str(paths[len(paths)-1][0:lastDirIndx])
            filename = str(paths[0])
            self.loadedfile = filename
            localfilename = filename[lastDirIndx+1:]
            # Delete a file if we're loading a new one
            try:
                self.cmd.remove("params")
                self.cmd.delete("params")
            except:
                pass
            try:
                self.cmd.load(filename, "params")
            except:
                wx.MessageBox("The file " + filename + " could not be read!", "File Cannot Be Read", wx.OK|wx.ICON_EXCLAMATION)
                return
            logInfo("MOL2 file loaded", filename)
            self.cmd.select("paramssele", "model params")
            self.cmd.hide("everything", "paramssele")
            self.cmd.show("sticks", "paramssele")
            self.cmd.color("gray", "paramssele and symbol c")
            self.cmd.show("spheres", "metal and model params")
            self.btnCreate.Enable()
            self.lblMOL2.SetLabel(localfilename)
            self.lblMOL2.SetForegroundColour("#FFFFFF")
            if (platform.system() == "Linux"):
                resizeTextControlForUNIX(self.lblMOL2, 10, 180)
        else:
            logInfo("Load MOL2 operation cancelled")
            
    def typeToggle(self, event):
        if (self.paramsType == "Ligand"):
            self.paramsType = "Polymer"
            self.NtermMenu.Enable()
            self.CtermMenu.Enable()
            if (platform.system() == "Darwin"):
                self.btnType.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/advresiduecreator/btnType_Polymer.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnType.SetLabel(self.paramsType)
            self.btnType.SetToolTipString("Uploaded .mol2 file represents an NCAA embedded as part of a polypeptide sequence")
        else:
            self.paramsType = "Ligand"
            self.NtermMenu.Disable()
            self.CtermMenu.Disable()
            if (platform.system() == "Darwin"):
                self.btnType.SetBitmapLabel(bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/advresiduecreator/btnType_Ligand.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap())
            else:
                self.btnType.SetLabel(self.paramsType)
            self.btnType.SetToolTipString("Uploaded .mol2 file represents a ligand")
        logInfo("Params type changed to " + self.paramsType)
        
    def createParams(self, event):
        # Change to the sandbox location
        goToSandbox()
        # Attempt to generate the params file
        try:
            if (os.path.isfile("LG.params")):
                os.remove("LG.params")
            if (os.path.isfile("LG.fa.params")):
                os.remove("LG.fa.params")
            if (os.path.isfile("LG.cen.params")):
                os.remove("LG.cen.params")
            molfile_to_params.main([self.loadedfile, "--no-pdb", "--keep-names", "-c"])
        except:
            wx.MessageBox("The file " + self.loadedfile + " could not be converted to a params file!", "File Cannot Be Processed", wx.OK|wx.ICON_EXCLAMATION)
            return
        logInfo("Params file created successfully")
        # Now read the LG.params file and grab out the atom names and their assigned types
        # so the user can see them and modify them if desired
        f = open("LG.fa.params", "r")
        if (self.grdParamsAtoms.NumberRows > 0):
            self.grdParamsAtoms.DeleteRows(0, self.grdParamsAtoms.NumberRows)
        self.atomnames = []
        atomtypes = []
        for aline in f:
            if (aline[0:4] == "ATOM"):
                atomname = aline.split()[1]
                atomtype = aline.split()[2]
                self.atomnames.append(atomname)
                atomtypes.append(atomtype)
        f.close()
        # Sort the atomnames to make it easier for the user to find things
        for i in range(0, len(self.atomnames)-1):
            lowest = i
            for j in range(i+1, len(self.atomnames)):
                if (self.atomnames[j] < self.atomnames[lowest]):
                    lowest = j
            temp = self.atomnames[i]
            self.atomnames[i] = self.atomnames[lowest]
            self.atomnames[lowest] = temp
            temp = atomtypes[i]
            atomtypes[i] = atomtypes[lowest]
            atomtypes[lowest] = temp
        # Now add things to the grid
        for i in range(0, len(self.atomnames)):
            self.grdParamsAtoms.AppendRows(1)
            self.grdParamsAtoms.SetRowLabelValue(i, self.atomnames[i])
            self.grdParamsAtoms.SetCellValue(i, 0, atomtypes[i])
            self.grdParamsAtoms.SetCellAlignment(i, 0, wx.ALIGN_CENTRE, wx.ALIGN_CENTRE)
            readOnly = wx.grid.GridCellAttr()
            readOnly.SetReadOnly(True)
            self.grdParamsAtoms.SetRowAttr(i, readOnly)
        # Update some of the atom selection menus with the list of atomnames
        self.atomMenu.Clear()
        self.atomMenu.AppendItems(self.atomnames)
        self.NtermMenu.Clear()
        self.NtermMenu.AppendItems(self.atomnames)
        self.CtermMenu.Clear()
        self.CtermMenu.AppendItems(self.atomnames)
        self.btnAdd.Enable()

    def atomMenuSelect(self, event):
        # PyMOL should make the selected atom stand out, so it will generate a small sphere
        # for the selected atom
        logInfo("Atom " + self.atomMenu.GetStringSelection() + " selected")
        # Set the selected residue's row to red so it is easy to see what the selection is
        for r in range(0, self.grdParamsAtoms.NumberRows):
            if (r == self.atomMenu.GetSelection()):
                for c in range(0, self.grdParamsAtoms.NumberCols):
                    self.grdParamsAtoms.SetCellBackgroundColour(r, c, "light blue")
            else:
                for c in range(0, self.grdParamsAtoms.NumberCols):
                    self.grdParamsAtoms.SetCellBackgroundColour(r, c, "white")
        self.grdParamsAtoms.Refresh()
        self.cmd.select("paramssele", "model params")
        self.cmd.hide("spheres", "paramssele")
        try:
            self.cmd.select("paramssele", "model params and name " + self.atomMenu.GetStringSelection())
            self.cmd.show("spheres", "paramssele")
            self.cmd.set("sphere_scale", 0.3, "paramssele")
            self.cmd.select("sele", "paramssele")
            self.cmd.enable("sele")
            self.cmd.delete("paramssele")
            # Get the element symbol from PyMOL
            self.pymol.stored.element = ""
            self.cmd.iterate_state(1, "sele", "stored.element = elem")
            # Put the available choices into the typeMenu
            self.typeMenu.Clear()
            self.typeMenu.AppendItems(self.atomtypes[self.pymol.stored.element.upper()])
            for i in range(0, self.grdParamsAtoms.NumberRows):
                if (self.atomMenu.GetStringSelection().strip() == self.grdParamsAtoms.GetRowLabelValue(i).strip()):
                    currentType = self.grdParamsAtoms.GetCellValue(i, 0).strip()
                    break
            indx = self.typeMenu.GetItems().index(currentType)
            self.typeMenu.SetSelection(indx)
        except:
            pass

    def typeMenuSelect(self, event):
        # Replace the type in the grid with the selected type
        atomname = self.atomMenu.GetStringSelection()
        atomtype = self.typeMenu.GetStringSelection()
        logInfo("Atom type " + atomtype + " selected")
        for r in range(0, self.grdParamsAtoms.NumberRows):
            if (atomname.strip() == self.grdParamsAtoms.GetRowLabelValue(r).strip()):
                self.grdParamsAtoms.SetCellValue(r, 0, atomtype)
                break
        
    def addToDB(self, event):
        # First check to make sure that a valid code is given
        code = self.txtCode.GetValue().strip().upper()
        if (len(code) > 3):
            wx.MessageBox("You have not entered a valid 3-letter code.  Please enter a valid 3-letter code.", "Bad Code", wx.OK|wx.ICON_EXCLAMATION)
            return
        # If this is a polymer, make sure both an N and C termini are specified
        if (self.paramsType == "Polymer"):
            Nterm = self.NtermMenu.GetStringSelection().strip()
            Cterm = self.CtermMenu.GetStringSelection().strip()
            if (len(Nterm) == 0 or len(Cterm) == 0):
                wx.MessageBox("Please choose an N and C terminus atom for your polymer residue.", "Termini Not Specified", wx.OK|wx.ICON_EXCLAMATION)
                return
        # Now make sure this parameters file isn't already in our database
        if (platform.system() == "Windows"):
            paramsfile = "params\\" + code + ".params"
            paramslist = glob.glob("params\\*.params")
        else:
            paramsfile = "params/" + code + ".params"
            paramslist = glob.glob("params/*.params")
        if (paramsfile in paramslist):
            dlg = wx.MessageDialog(self, "There is already a parameters file for this code in the database.  Do you want to overwrite the previous entry?", "Duplicate Parameters", wx.YES_NO | wx.ICON_EXCLAMATION | wx.CENTRE)
            if (dlg.ShowModal() == wx.ID_NO):
                dlg.Destroy()
                return
            dlg.Destroy()
        logInfo("Saved params file with the code " + self.txtCode.GetValue().strip().upper())
        if (self.paramsType == "Polymer"):
            logInfo("The N-terminus was " + Nterm + " and the C-terminus was " + Cterm)
        # Now we have to read the LG.params file and replace the user uploaded data
        for (origfile, cen_fa) in [("LG.fa.params", "fa"), ("LG.cen.params", "cen")]:
            f = open(origfile, "r")
            f2 = open(paramsfile + "." + cen_fa + ".params", "w")
            wroteTermini = False
            for aline in f:
                if (aline[0:4] == "NAME"):
                    f2.write("NAME " + code + "\n")
                elif (aline[0:9] == "IO_STRING"):
                    f2.write("IO_STRING " + code + " Z\n")
                elif (aline[0:4] == "TYPE"):
                    f2.write("TYPE " + self.paramsType.upper() + "\n")
                elif (aline[0:4] == "ATOM"):
                    # Now we have to find the type from our graph
                    thisatomname = aline[5:9].strip()
                    for r in range(0, self.grdParamsAtoms.NumberRows):
                        if (thisatomname == self.grdParamsAtoms.GetRowLabelValue(r).strip()):
                            thisatomtype = self.grdParamsAtoms.GetCellValue(r, 0).strip()
                            break
                    # Get the right amount of whitespace on all sides
                    if (len(thisatomtype) == 1):
                        thisatomtype = thisatomtype + "   "
                    elif (len(thisatomtype) == 2):
                        thisatomtype = thisatomtype + "  "
                    elif (len(thisatomtype) == 3):
                        thisatomtype = thisatomtype + " "
                    # Replace the old entry
                    aline = aline[0:10] + thisatomtype + aline[14:]
                    f2.write(aline)
                elif (aline[0:4] == "BOND" and not(wroteTermini) and self.paramsType == "Polymer"):
                    f2.write("LOWER_CONNECT " + Nterm + "\n")
                    f2.write("UPPER_CONNECT " + Cterm + "\n")
                    f2.write(aline)
                    wroteTermini = True
                else:
                    f2.write(aline)
            f.close()
            f2.close()
        # Delete the LG.params file
        os.remove("LG.fa.params")
        os.remove("LG.cen.params")
        wx.MessageBox("Your parameters file was created successfully!  InteractiveROSETTA will now recognize " + self.txtCode.GetValue() + " entries.", "Params Creation Successful", wx.OK|wx.ICON_EXCLAMATION)
        self.grdParamsAtoms.ClearGrid()
        self.atomMenu.Clear()
        self.typeMenu.Clear()
        self.btnAdd.Disable()
        self.removeMenu.Append(self.selectedType)
        self.cmd.remove("params")
        
    def removeParams(self, event):
        # Take the indicated parameters file out of the database
        paramsToRemove = self.removeMenu.GetStringSelection().strip()
        if (len(paramsToRemove) == 0):
            return
        dlg = wx.MessageDialog(self, "This operation will remove " + paramsToRemove + " from the database.  Are you sure you want to proceed?", "Parameters Removal", wx.YES_NO | wx.ICON_EXCLAMATION | wx.CENTRE)
        if (dlg.ShowModal() == wx.ID_NO):
            dlg.Destroy()
            return
        dlg.Destroy()
        goToSandbox("params")
        try:
            os.remove(paramsToRemove + ".fa.params")
        except:
            pass
        goToSandbox()
        self.removeMenu.Delete(self.removeMenu.GetItems().index(paramsToRemove))
        dlg = wx.MessageDialog(self, "If any loaded models contained " + paramsToRemove + ", you need to unload and reload them to prevent unexpected behavior.", "Parameters Removal", wx.OK | wx.ICON_EXCLAMATION | wx.CENTRE)
        dlg.ShowModal()
        dlg.Destroy()
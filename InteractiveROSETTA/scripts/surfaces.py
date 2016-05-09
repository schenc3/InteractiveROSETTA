import wx
import wx.grid
import wx.lib.scrolledpanel
import os
import os.path
import platform
from tools import *
import webbrowser
from threading import Thread

class SurfacesPanel(wx.lib.scrolledpanel.ScrolledPanel):
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
            self.lblProt = wx.StaticText(self, -1, "Molecular Surfaces", (25, 15), (270, 25), style=wx.ALIGN_CENTRE)
            self.lblProt.SetFont(wx.Font(12, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblProt = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/surfaces/lblSurfaces.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 15), size=(320, 25))
        else:
            self.lblProt = wx.StaticText(self, -1, "Molecular Surfaces", pos=(90, 15), style=wx.ALIGN_CENTRE)
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
            self.lblInst = wx.StaticText(self, -1, "Draw and configure structural surfaces", (0, 45), (320, 25), wx.ALIGN_CENTRE)
            self.lblInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
        elif (platform.system() == "Darwin"):
            self.lblInst = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/surfaces/lblInstSurfaces.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 45), size=(320, 25))
        else:
            self.lblInst = wx.StaticText(self, -1, "Draw and configure structural surfaces", pos=(20, 45), style=wx.ALIGN_CENTRE)
            self.lblInst.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.NORMAL))
            resizeTextControlForUNIX(self.lblInst, 0, self.GetSize()[0])
        self.lblInst.SetForegroundColour("#FFFFFF")
        
        if (platform.system() == "Windows"):
            self.lblReceptor = wx.StaticText(self, -1, "Receptor", (0, 90), (155, 20), wx.ALIGN_CENTRE)
            self.lblReceptor.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblReceptor = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/surfaces/lblReceptor.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 90), size=(155, 20))
        else:
            self.lblReceptor = wx.StaticText(self, -1, "Receptor", (0, 90), style=wx.ALIGN_CENTRE)
            self.lblReceptor.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblReceptor, 0, 155)
        self.lblReceptor.SetForegroundColour("#FFFFFF")
        if (platform.system() == "Windows"):
            self.lblLigand = wx.StaticText(self, -1, "Ligand", (160, 90), (155, 20), wx.ALIGN_CENTRE)
            self.lblLigand.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblLigand = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/surfaces/lblLigand.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(160, 90), size=(155, 20))
        else:
            self.lblLigand = wx.StaticText(self, -1, "Ligand", (160, 90), style=wx.ALIGN_CENTRE)
            self.lblLigand.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblLigand, 160, 155)
        self.lblLigand.SetForegroundColour("#FFFFFF")
        
        if (platform.system() == "Darwin"):
            self.btnRecpSelect = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/surfaces/btnRecpSelect.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 120), size=(155, 25))
        else:
            self.btnRecpSelect = wx.Button(self, id=-1, label="Set Selection", pos=(0, 120), size=(155, 25))
            self.btnRecpSelect.SetForegroundColour("#000000")
            self.btnRecpSelect.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnRecpSelect.Bind(wx.EVT_BUTTON, self.recpSetSelection)
        self.btnRecpSelect.SetToolTipString("Set the current selection as the receptor.  Use only the receptor if you simply want a surface and not an interface.")
        if (platform.system() == "Darwin"):
            self.btnLigSelect = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/surfaces/btnRecpSelect.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(160, 120), size=(155, 25))
        else:
            self.btnLigSelect = wx.Button(self, id=-1, label="Set Selection", pos=(160, 120), size=(155, 25))
            self.btnLigSelect.SetForegroundColour("#000000")
            self.btnLigSelect.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnLigSelect.Bind(wx.EVT_BUTTON, self.ligSetSelection)
        self.btnLigSelect.SetToolTipString("Set the current selection as the ligand.  If you have ligand and receptor specified on the same model, then the interface between the receptor and ligand will be drawn.")
        
        if (platform.system() == "Darwin"):
            self.btnRecpHighlight = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/surfaces/btnRecpHighlight.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 150), size=(155, 25))
        else:
            self.btnRecpHighlight = wx.Button(self, id=-1, label="Highlight", pos=(0, 150), size=(155, 25))
            self.btnRecpHighlight.SetForegroundColour("#000000")
            self.btnRecpHighlight.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnRecpHighlight.Bind(wx.EVT_BUTTON, self.recpHighlight)
        self.btnRecpHighlight.SetToolTipString("Show the residues currently selected for the receptor.")
        if (platform.system() == "Darwin"):
            self.btnLigHighlight = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/surfaces/btnRecpHighlight.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(160, 150), size=(155, 25))
        else:
            self.btnLigHighlight = wx.Button(self, id=-1, label="Highlight", pos=(160, 150), size=(155, 25))
            self.btnLigHighlight.SetForegroundColour("#000000")
            self.btnLigHighlight.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnLigHighlight.Bind(wx.EVT_BUTTON, self.ligHighlight)
        self.btnLigHighlight.SetToolTipString("Show the residues currently selected for the ligand.")
        
        if (platform.system() == "Windows"):
            self.lblColoring = wx.StaticText(self, -1, "Coloring:", (0, 190), (80, 20), wx.ALIGN_CENTRE)
            self.lblColoring.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblColoring = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/surfaces/lblColoring.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 183), size=(80, 20))
        else:
            self.lblColoring = wx.StaticText(self, -1, "Coloring:", (0, 190), style=wx.ALIGN_CENTRE)
            self.lblColoring.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblColoring, 0, 80)
        self.lblColoring.SetForegroundColour("#FFFFFF")
        if (platform.system() == "Darwin"):
            self.btnCPKColor = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/surfaces/btnCPKColor.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(80, 180), size=(80, 35))
        else:
            self.btnCPKColor = wx.Button(self, id=-1, label="CPK", pos=(80, 180), size=(80, 35))
            self.btnCPKColor.SetForegroundColour("#000000")
            self.btnCPKColor.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnCPKColor.Bind(wx.EVT_BUTTON, self.CPKColor)
        self.btnCPKColor.SetToolTipString("Color the surface according to CPK coloring")
        if (platform.system() == "Darwin"):
            self.btnElecColor = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/surfaces/btnElecColor.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(160, 180), size=(80, 35))
        else:
            self.btnElecColor = wx.Button(self, id=-1, label="Elec", pos=(160, 180), size=(80, 35))
            self.btnElecColor.SetForegroundColour("#000000")
            self.btnElecColor.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnElecColor.Bind(wx.EVT_BUTTON, self.elecColor)
        self.btnElecColor.SetToolTipString("Color the surface according to electrostatic coloring")
        self.btnCustomColor = wx.BitmapButton(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/colorwheel.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), (240, 180), (80, 35))
        self.btnCustomColor.Bind(wx.EVT_BUTTON, self.customColor)
        self.btnCustomColor.SetToolTipString("Color the surface according to a user-defined color")
        
        if (platform.system() == "Windows"):
            self.lblTransparency = wx.StaticText(self, -1, "Transparency:", (0, 223), (100, 20), wx.ALIGN_CENTRE)
            self.lblTransparency.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblTransparency = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/surfaces/lblTransparency.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 223), size=(100, 20))
        else:
            self.lblTransparency = wx.StaticText(self, -1, "Transparency:", (0, 223), style=wx.ALIGN_CENTRE)
            self.lblTransparency.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblTransparency, 0, 100)
        self.lblTransparency.SetForegroundColour("#FFFFFF")
        self.sldTransparency = wx.Slider(self, -1, 0, 0, 100, (110, 220), (210, 25), wx.SL_HORIZONTAL)
        self.sldTransparency.SetTickFreq(5, 1)
        self.sldTransparency.Bind(wx.EVT_SLIDER, self.transSlide)
        
        if (platform.system() == "Windows"):
            self.lblSurfaceName = wx.StaticText(self, -1, "Surface Name:", (0, 253), (140, 20), wx.ALIGN_CENTRE)
            self.lblSurfaceName.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        elif (platform.system() == "Darwin"):
            self.lblSurfaceName = wx.StaticBitmap(self, -1, wx.Image(self.parent.parent.scriptdir + "/images/osx/surfaces/lblSurfaceName.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(0, 253), size=(140, 20))
        else:
            self.lblSurfaceName = wx.StaticText(self, -1, "Surface Name:", (0, 253), style=wx.ALIGN_CENTRE)
            self.lblSurfaceName.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
            resizeTextControlForUNIX(self.lblSurfaceName, 0, 140)
        self.lblSurfaceName.SetForegroundColour("#FFFFFF")
        self.txtSurfaceName = wx.TextCtrl(self, -1, pos=(140, 250), size=(180, 25))
        self.txtSurfaceName.SetValue("")
        self.txtSurfaceName.SetToolTipString("Give the surface a name")
        
        if (platform.system() == "Darwin"):
            self.menuSurfaces = wx.ComboBox(self, pos=(0, 280), size=(220, 25), choices=[], style=wx.CB_READONLY)
        else:
            self.menuSurfaces = wx.ComboBox(self, pos=(0, 280), size=(220, 25), choices=[], style=wx.CB_READONLY | wx.CB_SORT)
        self.menuSurfaces.Bind(wx.EVT_COMBOBOX, self.surfaceSelect)
        self.menuSurfaces.SetToolTipString("List of currently active surfaces")
        if (platform.system() == "Darwin"):
            self.btnDelete = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/surfaces/btnDelete.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(220, 280), size=(100, 25))
        else:
            self.btnDelete = wx.Button(self, id=-1, label="Delete", pos=(220, 280), size=(100, 25))
            self.btnDelete.SetForegroundColour("#000000")
            self.btnDelete.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnDelete.Bind(wx.EVT_BUTTON, self.deleteSurface)
        self.btnDelete.SetToolTipString("Delete the selected surface")
        
        if (platform.system() == "Darwin"):
            self.btnDraw = wx.BitmapButton(self, id=-1, bitmap=wx.Image(self.parent.parent.scriptdir + "/images/osx/surfaces/btnDraw.png", wx.BITMAP_TYPE_PNG).ConvertToBitmap(), pos=(110, 310), size=(100, 25))
        else:
            self.btnDraw = wx.Button(self, id=-1, label="Draw!", pos=(110, 310), size=(100, 25))
            self.btnDraw.SetForegroundColour("#000000")
            self.btnDraw.SetFont(wx.Font(10, wx.DEFAULT, wx.ITALIC, wx.BOLD))
        self.btnDraw.Bind(wx.EVT_BUTTON, self.drawClick)
        self.btnDraw.SetToolTipString("Draw the indicated surface")
        
        self.scrollh = self.btnDraw.GetPosition()[1] + self.btnDraw.GetSize()[1] + 5
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
        try:
            self.cmd.delete("curr_surf_recp")
        except:
            pass
        try:
            self.cmd.delete("curr_surf_lig")
        except:
            pass
        self.cmd.hide("surface", "all")
        items = []
        for name in self.cmd.get_names("selections"):
            if (name.startswith("surf_recp_")):
                items.append(name[10:])
        self.menuSurfaces.AppendItems(items)

    def scrolled(self, event):
        self.winscrollpos = self.GetScrollPos(wx.VERTICAL)
        event.Skip()

    def activate(self):
        self.Scroll(0, self.winscrollpos)

    def recpSetSelection(self, event):
        self.seqWin.selectUpdate()
        try:
            self.cmd.hide("surface", "curr_surf_recp")
        except:
            pass
        try:
            self.cmd.show("surface", "seqsele")
            self.cmd.select("curr_surf_recp", "seqsele")
            self.cmd.disable("curr_surf_recp")
            self.cmd.enable("seqsele")
        except:
            self.cmd.show("surface", "all")
            self.cmd.select("curr_surf_recp", "all")
            self.cmd.disable("curr_surf_recp")
            self.cmd.enable("seqsele")
    
    def ligSetSelection(self, event):
        self.seqWin.selectUpdate()
        try:
            self.cmd.flag("ignore", "curr_surf_lig", "clear")
        except:
            pass
        try:
            self.cmd.flag("ignore", "seqsele", "set")
            self.cmd.select("curr_surf_lig", "seqsele")
            self.cmd.disable("curr_surf_lig")
            self.cmd.enable("seqsele")
            self.cmd.rebuild()
        except:
            pass
    
    def recpHighlight(self, event):
        try:
            self.cmd.select("seqsele", "curr_surf_recp")
            self.seqWin.selectUpdate()
        except:
            pass
    
    def ligHighlight(self, event):
        try:
            self.cmd.select("seqsele", "curr_surf_lig")
            self.seqWin.selectUpdate()
        except:
            pass
    
    def CPKColor(self, event):
        try:
            self.cmd.set("surface_color", "gray", "curr_surf_recp")
            self.cmd.set("surface_color", "blue", "curr_surf_recp and symbol n")
            self.cmd.set("surface_color", "red", "curr_surf_recp and symbol o")
            self.cmd.set("surface_color", "yellow", "curr_surf_recp and symbol s")
            self.cmd.set("surface_color", "purple", "curr_surf_recp and symbol p")
        except:
            pass
    
    def elecColor(self, event):
        try:
            total_charge = 0
            natoms = 0
            self.stored.bfactors = []
            self.cmd.iterate_state(1, "curr_surf_recp", "stored.bfactors.append([model, chain, resi, name, formal_charge])")
            for model, chain, resi, name, c in self.stored.bfactors:
                natoms += 1
                total_charge += c / 4
                if (len(chain.strip()) == 0):
                    chain = " "
                    selstring = "model " + model + " and resi " + str(resi) + " and name " + name
                else:
                    selstring = "model " + model + " and chain " + chain + " and resi " + str(resi) + " and name " + name
                if (c <= 0):
                    red = 255
                    blue = int(255 * float(4+c)/4.0)
                    green = blue
                else:
                    blue = 255
                    red = int(255 * float(4-c)/4.0)
                    green = red
                color = "0x%02x%02x%02x" % (red, green, blue)
                self.cmd.set("surface_color", color, selstring)
            if (total_charge < natoms / -2.0 or total_charge > natoms / 2):
                self.cmd.set("surface_color", "white", "curr_surf_recp")
                self.cmd.set("surface_color", "blue", "curr_surf_recp and (resn lys or resn arg)")
                self.cmd.set("surface_color", "red", "curr_surf_recp and (resn asp or resn glu)")
        except:
            pass
    
    def customColor(self, event):
        dlg = wx.ColourDialog(self)
        dlg.GetColourData().SetChooseFull(True)
        if (dlg.ShowModal() == wx.ID_OK):
            data = dlg.GetColourData()
            mycolor = "0x%02x%02x%02x" % data.GetColour().Get()
            try:
                self.cmd.set("surface_color", mycolor, "curr_surf_recp")
            except:
                pass
            logInfo("Set surface color to " + mycolor)
        dlg.Destroy()
    
    def transSlide(self, event):
        self.cmd.set("transparency", str(float(self.sldTransparency.GetValue()) / 100.0))
    
    def surfaceSelect(self, event):
        # Hide everything except for the current surface
        self.cmd.hide("surface", "all")
        self.cmd.flag("ignore", "all", "clear")
        try:
            self.cmd.select("curr_surf_lig", "surf_lig_" + self.menuSurfaces.GetValue())
            self.cmd.disable("curr_surf_lig")
        except:
            pass
        self.cmd.select("curr_surf_recp", "surf_recp_" + self.menuSurfaces.GetValue())
        self.cmd.disable("curr_surf_recp")
        try:
            self.cmd.flag("ignore", "curr_surf_lig", "set")
        except:
            pass
        self.cmd.show("surface", "curr_surf_recp")
        self.cmd.enable("seqsele")
    
    def deleteSurface(self, event):
        # Delete this surface from the selections
        if (len(self.menuSurfaces.GetValue()) == 0):
            return
        try:
            self.cmd.delete("surf_lig_" + self.menuSurfaces.GetValue())
        except:
            pass
        self.cmd.delete("surf_recp_" + self.menuSurfaces.GetValue())
        items = self.menuSurfaces.GetItems()
        indx = items.index(self.menuSurfaces.GetValue())
        items.pop(indx)
        self.menuSurfaces.Clear()
        self.menuSurfaces.AppendItems(items)
    
    def drawClick(self, event):
        # Is there a name in the text box?
        if (len(self.txtSurfaceName.GetValue().strip()) == 0):
            dlg = wx.MessageDialog(self, "You have not provided a name for this surface!  Please provide a name.", "No Surface Name", wx.OK | wx.ICON_ERROR | wx.CENTRE)
            dlg.ShowModal()
            dlg.Destroy()
            return
        elif ("surf_recp_" + self.txtSurfaceName.GetValue().strip() in self.menuSurfaces.GetItems()):
            dlg = wx.MessageDialog(self, "The name " + self.txtSurfaceName.GetValue().strip() + " is already taken.  Please provide a unique surface name or delete the surface with this name.", "Name Already Taken", wx.OK | wx.ICON_ERROR | wx.CENTRE)
            dlg.ShowModal()
            dlg.Destroy()
            return
        # Save this setup as a selection name in PyMOL and add it to the list
        try:
            self.cmd.select("surf_recp_" + self.txtSurfaceName.GetValue().strip(), "curr_surf_recp")
            self.cmd.disable("surf_recp_" + self.txtSurfaceName.GetValue().strip())
        except:
            dlg = wx.MessageDialog(self, "You have not selected residues to be in the receptor.  A surface cannot be drawn without this information.", "No Surface Receptor", wx.OK | wx.ICON_ERROR | wx.CENTRE)
            dlg.ShowModal()
            dlg.Destroy()
            return
        try:
            self.cmd.select("surf_lig_" + self.txtSurfaceName.GetValue().strip(), "curr_surf_lig")
            self.cmd.disable("surf_lig_" + self.txtSurfaceName.GetValue().strip())
        except:
            # It's okay if this is not defined, it's optional
            pass
        self.cmd.enable("seqsele")
        self.menuSurfaces.Append(self.txtSurfaceName.GetValue().strip())
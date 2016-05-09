#!/usr/bin/python
import wx
import os
import commands

### This is a simple widget I wrote to make it easy to create those stupid PNG images that
### are needed for text buttons on OSX
### It depends upon Linux Imagemagick
### Don't know if OSX also has Imagemagick, but the command you need is the "convert" binary
class MainFrame(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self, None, wx.ID_ANY, "InteractiveROSETTA OSX Text Creator", size=(400, 400))
        
        self.lblWidth = wx.StaticText(self, -1, "Button Width", size=(190, 20), style=wx.ALIGN_LEFT)
        self.lblWidth.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.lblHeight = wx.StaticText(self, -1, "Button Height", size=(190, 20), style=wx.ALIGN_LEFT)
        self.lblHeight.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        
        self.txtWidth = wx.TextCtrl(self, -1, size=(190, 25))
        self.txtWidth.SetValue("")
        self.txtWidth.SetToolTipString("The width of the button in pixels")
        self.txtHeight = wx.TextCtrl(self, -1, size=(190, 25))
        self.txtHeight.SetValue("")
        self.txtHeight.SetToolTipString("The height of the button in pixels")
        
        self.lblText = wx.StaticText(self, -1, "Button Text", size=(190, 20), style=wx.ALIGN_LEFT)
        self.lblText.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.txtText = wx.TextCtrl(self, -1, size=(190, 25))
        self.txtText.SetValue("")
        self.txtText.SetToolTipString("The text to be display on the button")
        
        self.chkBold = wx.CheckBox(self, -1, "Bold", size=(190, 25))
        self.chkItalic = wx.CheckBox(self, -1, "Italic", size=(190, 25))
        
        self.outputDir = os.getcwd()
        self.lblOutputDir = wx.StaticText(self, -1, "Output Directory", size=(390, 20), style=wx.ALIGN_LEFT)
        self.lblOutputDir.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnOutputDir = wx.Button(self, id=-1, label=self.outputDir, size=(390, 30))
        self.btnOutputDir.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnOutputDir.Bind(wx.EVT_BUTTON, self.selectOutputDir)
        self.btnOutputDir.SetToolTipString("The directory to which images will be outputted")
        
        self.lblName = wx.StaticText(self, -1, "Image Name", size=(190, 20), style=wx.ALIGN_LEFT)
        self.lblName.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.txtName = wx.TextCtrl(self, -1, size=(190, 25))
        self.txtName.SetValue("")
        self.txtName.SetToolTipString("The name of the image (without an extension)")
        
        self.chkRed = wx.CheckBox(self, -1, "Red Text", size=(190, 25))
        self.btnCreate = wx.Button(self, id=-1, label="Create!", size=(190, 30))
        self.btnCreate.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        self.btnCreate.Bind(wx.EVT_BUTTON, self.createImage)
        self.btnCreate.SetToolTipString("Create the image")
        
        sizer = wx.GridBagSizer(hgap=5, vgap=5)
        sizer.Add(self.lblWidth, pos=(0, 0), flag=wx.ALL, border=5)
        sizer.Add(self.lblHeight, pos=(0, 1), flag=wx.ALL, border=5)
        sizer.Add(self.txtWidth, pos=(1, 0), flag=wx.ALL, border=5)
        sizer.Add(self.txtHeight, pos=(1, 1), flag=wx.ALL, border=5)
        sizer.Add(self.lblText, pos=(2, 0), flag=wx.ALL, border=5)
        sizer.Add(self.txtText, pos=(2, 1), flag=wx.ALL, border=5)
        sizer.Add(self.chkBold, pos=(3, 0), flag=wx.ALL, border=5)
        sizer.Add(self.chkItalic, pos=(3, 1), flag=wx.ALL, border=5)
        sizer.Add(self.lblOutputDir, pos=(4, 0), span=(1, 2), flag=wx.ALL, border=5)
        sizer.Add(self.btnOutputDir, pos=(5, 0), span=(1, 2), flag=wx.ALL, border=5)
        sizer.Add(self.lblName, pos=(6, 0), flag=wx.ALL, border=5)
        sizer.Add(self.chkRed, pos=(6, 1), flag=wx.ALL, border=5)
        sizer.Add(self.txtName, pos=(7, 0), flag=wx.ALL, border=5)
        sizer.Add(self.btnCreate, pos=(7, 1), flag=wx.ALL, border=5)
        self.SetSizer(sizer)
        self.Layout()
        self.Show()
        
    def selectOutputDir(self, event):
        dlg = wx.DirDialog(self, "Select the output directory",
                            defaultPath=self.outputDir,
                            style=wx.DD_DEFAULT_STYLE 
                            | wx.DD_DIR_MUST_EXIST)
        if (dlg.ShowModal() != wx.ID_OK):
            return
        path = dlg.GetPath().strip()
        self.outputDir = str(path)
        self.btnOutputDir.SetLabel(self.outputDir)
    
    def createImage(self, event):
        # Do we have all the input data?
        try:
            width = int(self.txtWidth.GetValue())
        except:
            dlg = wx.MessageDialog(self, "Please give an integer width in pixels.", "Invalid Width", wx.OK | wx.ICON_ERROR | wx.CENTRE)
            dlg.ShowModal()
            dlg.Destroy()
            return
        try:
            height = int(self.txtHeight.GetValue())
        except:
            dlg = wx.MessageDialog(self, "Please give an integer height in pixels.", "Invalid Width", wx.OK | wx.ICON_ERROR | wx.CENTRE)
            dlg.ShowModal()
            dlg.Destroy()
            return
        textstring = self.txtText.GetValue().strip()
        if (len(textstring) == 0):
            dlg = wx.MessageDialog(self, "Please provide text for the button.", "Invalid Label", wx.OK | wx.ICON_ERROR | wx.CENTRE)
            dlg.ShowModal()
            dlg.Destroy()
            return
        imagename = self.txtName.GetValue().strip()
        if (len(imagename) == 0):
            dlg = wx.MessageDialog(self, "Please provide a name for the image file.", "Invalid Filename", wx.OK | wx.ICON_ERROR | wx.CENTRE)
            dlg.ShowModal()
            dlg.Destroy()
            return
        # Calculate a font size that will fit the size of the button
        fontsize = 48
        fsize = 10
        dc = wx.ScreenDC()
        if (self.chkBold.GetValue()):
            bold = wx.BOLD
        else:
            bold = wx.NORMAL
        if (self.chkItalic.GetValue()):
            italic = wx.ITALIC
        else:
            italic = wx.NORMAL
        for fontsize in range(48, 14, -2):
            dc.SetFont(wx.Font(fontsize, wx.DEFAULT, italic, bold, face="FreeSans"))
            (w, h) = dc.GetTextExtent(textstring)
            if (w <= width and h <= height):
                fsize = fontsize + 4
                tw = w
                th = h
                break
        print tw, th, wx.SystemSettings.GetFont(0).GetFamily(), wx.SystemSettings.GetFont(0).GetFaceName()
        commandline = "convert -size " + str(width) + "x" + str(height) + " xc:transparent -font FreeSans"
        if (self.chkBold.GetValue()):
            commandline += "Bold"
        if (self.chkItalic.GetValue()):
            commandline += "Oblique"
        commandline += " -pointsize " + str(fsize) + " -fill "
        if (self.chkRed.GetValue()):
            commandline += "red"
        else:
            commandline += "black"
        commandline += " -draw \"text " + str((width-tw)/2) + "," + str(height-((height-th)/2)) + " '" + textstring + "'\" \"" + self.outputDir + "/" + imagename + ".png\""
        print commandline
        res, output = commands.getstatusoutput(commandline)

if (__name__ == "__main__"):    
    app = wx.App()
    frame = MainFrame()
    app.MainLoop()

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

class ConstraintPanel(wx.lib.scrolledpanel.ScrolledPanel):
    '''A panel allowing for the selection of constraints.  Easily added into any
    pre-existing module.  Note: Doesn't need to be a ScrolledPanel because its
    parent will be a ScrolledPanel'''

    def __init__(self, parent):
      print 'creating constraint panel'
      wx.lib.scrolledpanel.ScrolledPanel.__init__(self,parent,-1)
      print 'Panel initialized'
      #sizer
      vbox = wx.BoxSizer(wx.VERTICAL)
      hbox = wx.BoxSizer(wx.HORIZONTAL)
      #Add constraint button
      self.ConstraintBtn = wx.Button(self,-1,label="Add Constraint",pos=(10,10),size=(120,25))
      self.ConstraintBtn.SetForegroundColour("#000000")
      self.ConstraintBtn.SetFont(wx.Font(10,wx.DEFAULT,wx.ITALIC,wx.BOLD))
      self.ConstraintBtn.Bind(wx.EVT_BUTTON,self.addConstraint)
      self.ConstraintBtn.SetToolTipString("Add a new constraint to the simulation")
      hbox.Add(self.ConstraintBtn,1)
      print 'constraint button done'
      #Save Constraints button
      self.SaveConstraintsBtn = wx.Button(self,-1,label="Save Constraints",pos=(90,10),size=(120,25))
      self.SaveConstraintsBtn.SetForegroundColour("#000000")
      self.SaveConstraintsBtn.SetFont(wx.Font(10,wx.DEFAULT,wx.ITALIC,wx.BOLD))
      self.SaveConstraintsBtn.Bind(wx.EVT_BUTTON,self.saveConstraints)
      self.SaveConstraintsBtn.SetToolTipString("Save current constraints to a file for use in later simulations")
      hbox.Add(self.SaveConstraintsBtn,1)
      print 'save button'
      #Load Constraints file button
      self.LoadConstraintsBtn = wx.Button(self,-1,label="Load Constraints",pos=(180,10),size=(120,25))
      self.LoadConstraintsBtn.SetForegroundColour("#000000")
      self.LoadConstraintsBtn.SetFont(wx.Font(10,wx.DEFAULT,wx.ITALIC,wx.BOLD))
      self.LoadConstraintsBtn.Bind(wx.EVT_BUTTON,self.loadConstraints)
      self.LoadConstraintsBtn.SetToolTipString("Load constraints from a file")
      hbox.Add(self.LoadConstraintsBtn,1)
      print 'load button'

      vbox.Add(hbox)

      self.SetSizer(vbox)
      self.SetupScrolling()
      print 'scrolling'

    def setSeqWin(self,seqWin):
      self.seqWin = seqWin

    def setPyMOL(self, pymol):
      self.pymol = pymol
      self.cmd = pymol.cmd
      self.stored = pymol.stored

    def setSelectWin(self, selectWin):
      self.selectWin = selectWin
      self.selectWin.setProtPanel(self)

    #Event Listeners

    def loadConstraints(self,event):
      '''Loads a set of constraints from a file to be used for this simulation'''
      logInfo('Load Constraints button clicked!')

    def saveConstraints(self,event):
      '''Saves the current set of constraints to a file that can be used in
      later simulations'''
      logInfo('Save Constraints button clicked!')

    def addConstraint(self,event):
      '''Adds a new constraint to the set of constraints'''
      logInfo('Add Constraint button clicked!')

    def gridClick(self,event):
      logInfo('Constraints Grid Clicked!')

#Figure out what wx Class this should inherit from
#wait until you finish writing the panel class to work on this....
class ConstraintPopUp():
    '''This class will be the pop-up menu where the user actually selects the
    constraints to add'''

    def __init__(self,parent,pos,W,H):
      logInfo('Initializing Contraint Pop-up!')
      #inherited class __init__
      #custom options

    def setSeqWin(self,seqWin):
        self.seqWin = seqWin

    def setPyMOL(self, pymol):
        self.pymol = pymol
        self.cmd = pymol.cmd
        self.stored = pymol.stored

    def setSelectWin(self, selectWin):
        self.selectWin = selectWin
        self.selectWin.setProtPanel(self)
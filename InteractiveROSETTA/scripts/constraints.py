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

class ConstraintPanel(wxlib.scrolledpanel.ScrolledPanel):
    '''A panel allowing for the selection of constraints.  Easily added into any
    pre-existing module.  Note: Doesn't need to be a ScrolledPanel because its
    parent will be a ScrolledPanel'''

    def __init__(self, parent, pos, W, H):
      wx.lib.scrolledpanel.ScrolledPanel.__init__(self,parent,id=-1,pos=pos,size=(W,H),name='Constraints')
      print 'Panel initialized'
      #Add constraint button
      x,y = pos
      self.ConstraintBtn = wx.Button(self,-1,label="Add Constraint",pos=(x+10,y+10),size=(57,25))
      self.ConstraintBtn.SetForegroundColour("#000000")
      self.ConstraintBtn.SetFont(wx.Font(10,wx.DEFAULT,wx.ITALIC,wx.BOLD))
      self.ConstraintBtn.Bind(wx.EVT_BUTTON,self.addConstraint)
      self.ConstraintBtn.SetToolTipString("Add a new constraint to the simulation")
      print 'constraint button done'
#==============================================================================
#       #Constraint Grid Type Atoms Function
#       self.ConstraintsGrid = wx.grid.Grid(self)
#       self.ConstraintsGrid.CreateGrid(0,3)
#       if H-235 > 200:
#         self.ConstraintsGrid.SetSize((320,H-235))
#       else:
#         self.ConstraintsGrid.SetSize((320,200))
#       self.ConstraintsGrid.SetPosition((0,y+15))
#       self.ConstraintsGrid.setLabelFont(wx.Font(10,wx.DEFAULT,wx.NORMAL,wx.BOLD))
#       self.ConstraintsGrid.DisableDragColSize()
#       self.ConstraintsGrid.DisableDragRowSize()
#       self.ConstraintsGrid.SetColLabelValue(0,"Type")
#       self.ConstraintsGrid.SetColLabelValue(1,"Atoms")
#       self.ConstraintsGrid.SetColLabelValue(2,"Function")
#       self.ConstraintsGrid.SetRowLabelSize(80)
#       self.ConstraintsGrid.SetColSize(0,100)
#       self.ConstraintsGrid.SetColSize(1,100)
#       self.ConstraintsGrid.SetColSize(2,100)
#       self.ConstraintsGrid.Bind(wx.grid.EVT_GRID_CELL_LEFT_CLICK,self.gridClick)
#       self.selectedr = -1
#       print 'grid set up'
#==============================================================================
      #Save Constraints button
      self.SaveConstraintsBtn = wx.Button(self,-1,label="Save Constraints",pos=(x+77,y+10),size=(57,25))
      self.SaveConstraintsBtn.SetForegroundColour("#000000")
      self.SaveConstraintsBtn.SetFont(wx.Font(10,wx.DEFAULT,wx.ITALIC,wx.BOLD))
      self.SaveConstraintsBtn.Bind(wx.EVT_BUTTON,self.saveConstraints)
      self.SaveConstraintsBtn.SetToolTipString("Save current constraints to a file for use in later simulations")
      print 'save button'
      #Load Constraints file button
      self.LoadConstraintsBtn = wx.Button(self,-1,label="Load Constraints",pos=(x+77+67,y+10),size=(57,25))
      self.LoadConstraintsBtn.SetForegroundColour("#000000")
      self.LoadConstraintsBtn.SetFont(wx.Font(10,wx.DEFAULT,wx.ITALIC,wx.BOLD))
      self.LoadConstraintsBtn.Bind(wx.EVT_BUTTON,self.loadConstraints)
      self.LoadConstraintsBtn.SetToolTipString("Load constraints from a file")
      print 'load button'

      self.scrollh = self.btnMinimize.GetPosition()[1] + self.btnMinimize.GetSize()[1] + 5
      self.SetScrollbars(1, 1, 320, self.scrollh)
      self.winscrollpos = 0
      self.Bind(wx.EVT_SCROLLWIN, self.scrolled)
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

    def scrolled(self,event):
      self.winscrollpos = self.getScrollPos(wx.VERTICAL)
      event.Skip()

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
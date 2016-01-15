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

    def __init__(self, parent,minPanel):
      print 'creating constraint panel'
      wx.lib.scrolledpanel.ScrolledPanel.__init__(self,parent,-1)
      self.minPanel = minPanel
      print 'Panel initialized'
      #sizer
      sizer = wx.GridBagSizer(10,10)
      self.SetSizer(sizer)

      #Add constraint button
      self.ConstraintBtn = wx.Button(self,-1,label="Add Constraint",size=(120,25))
      self.ConstraintBtn.SetForegroundColour("#000000")
      self.ConstraintBtn.SetFont(wx.Font(10,wx.DEFAULT,wx.ITALIC,wx.BOLD))
      self.ConstraintBtn.Bind(wx.EVT_BUTTON,self.addConstraint)
      self.ConstraintBtn.SetToolTipString("Add a new constraint to the simulation")
      self.Sizer.Add(self.ConstraintBtn,(0,0),(1,1))
      self.Layout()
      print 'constraint button done'
      #Save Constraints button
      self.SaveConstraintsBtn = wx.Button(self,-1,label="Save Constraints",size=(120,25))
      self.SaveConstraintsBtn.SetForegroundColour("#000000")
      self.SaveConstraintsBtn.SetFont(wx.Font(10,wx.DEFAULT,wx.ITALIC,wx.BOLD))
      self.SaveConstraintsBtn.Bind(wx.EVT_BUTTON,self.saveConstraints)
      self.SaveConstraintsBtn.SetToolTipString("Save current constraints to a file for use in later simulations")
      self.Sizer.Add(self.SaveConstraintsBtn,(0,1),(1,1))
      self.Layout()
      print 'save button'
      #Load Constraints file button
      self.LoadConstraintsBtn = wx.Button(self,-1,label="Load Constraints",size=(120,25))
      self.LoadConstraintsBtn.SetForegroundColour("#000000")
      self.LoadConstraintsBtn.SetFont(wx.Font(10,wx.DEFAULT,wx.ITALIC,wx.BOLD))
      self.LoadConstraintsBtn.Bind(wx.EVT_BUTTON,self.loadConstraints)
      self.LoadConstraintsBtn.SetToolTipString("Load constraints from a file")
      self.Sizer.Add(self.LoadConstraintsBtn,(0,2),(1,1))
      self.Layout()
      print 'load button'
      self.Cancelables = []

      #Remove Constraint Button

      #Clear Constraints Button


      #Constraints Grid
      print 'starting constraints grid'
      self.constraintsGrid = wx.grid.Grid(self,size=(1000,500))
      self.constraintsGrid.CreateGrid(0,6)
      self.constraintsGrid.SetLabelFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
      self.Sizer.Add(self.constraintsGrid,(4,0),(1,3))
      self.Layout()
      print 'grid created'

      #Scrolling
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

    def getConstrainableRes(self):
      '''Gets a list of all residues currenty being
      minimized and able to be restrained'''
      logInfo('Geting constrainable residues')

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
      self.constraintTypeText=wx.StaticText(self,-1,'Constraint Type',size=(100,25))
      self.constraintTypeText.SetFont(wx.Font(10,wx.DEFAULT,wx.ITALIC,wx.BOLD))
      self.Sizer.Add(self.constraintTypeText,(1,0),(1,1))
      self.Cancelables.append(self.constraintTypeText)
      self.Layout()
      constraintTypes = ['Choose Type','AtomPair','Angle','Dihedral','Coordinate']
      self.constraintTypeMenu = wx.ComboBox(self,size=(125,25),choices=constraintTypes,style=wx.CB_READONLY)
      self.constraintTypeMenu.Bind(wx.EVT_COMBOBOX,self.setConstraintType)
      self.Sizer.Add(self.constraintTypeMenu,(1,1),(1,1))
      self.Cancelables.append(self.constraintTypeMenu)
      self.Layout()
      #Cancel Button
      self.CancelBtn = wx.Button(self,-1,label='Cancel')
      self.CancelBtn.Bind(wx.EVT_BUTTON,self.cancel)
      self.Sizer.Add(self.CancelBtn,(1,2),(1,1))
      self.Cancelables.append(self.CancelBtn)
      self.Layout()

    def cancel(self,event):
      logInfo('Cancel Button Pressed!')
      print('Cancel Button Pressed!')
      for item in self.Cancelables:
        item.Show(False)
        item.Destroy()
      self.Cancelables = []
      self.Layout()

    def setConstraintType(self,event):
      logInfo('Constraint type selected!')
      if self.constraintTypeMenu.GetStringSelection() == 'Choose Type':
        event.skip()
      try:
        pdbs = ['Choose PDB']
        print pdbs
        for [indx, r, seqpos, poseindx, chainoffset, minType] in self.minPanel.minmap:
          print "poseindx",poseindx
          if len(pdbs) == 0:
            print 'appending',poseindx
            pdbs.append(poseindx)
          isThere = False
          for i in range(0,len(pdbs)):
            if poseindx == pdbs[i]:
              print("%i = %i"%(poseindx,pdbs[i]))
              isThere = True
              break
          if not isThere:
            print 'appending',poseindx
            pdbs.append(poseindx)
        print pdbs
        for i in range(1,len(pdbs)):
          pdbs[i] = str(self.seqWin.poses[pdbs[i]].get_id())
          print pdbs[i]
        print pdbs
        self.PdbText = wx.StaticText(self,-1,'PDB',size=(50,25))
        self.PdbText.SetFont(wx.Font(10,wx.DEFAULT,wx.ITALIC,wx.BOLD))
        self.Cancelables.append(self.PdbText)
        self.Sizer.Add(self.PdbText,(2,0),(1,1))
        self.Layout()
        print 'pdbtext'
        self.PdbMenu = wx.ComboBox(self,size=(125,25),choices=pdbs,style=wx.CB_READONLY)
        self.PdbMenu.Bind(wx.EVT_COMBOBOX,self.setConstraintPDB)
        self.Cancelables.append(self.PdbMenu)
        self.Sizer.Add(self.PdbMenu,(2,1),(1,1))
        self.Layout()
        print 'pdbmenu'
      except Exception as e:
        print(e.message)

    def setConstraintPDB(self,event):
      logInfo("constraint PDB set!")

    def gridClick(self,event):
      logInfo('Constraints Grid Clicked!')

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
      self.ConstraintBtn = wx.Button(self,-1,label="Add Constraint")
      self.ConstraintBtn.SetForegroundColour("#000000")
      self.ConstraintBtn.SetFont(wx.Font(10,wx.DEFAULT,wx.ITALIC,wx.BOLD))
      self.ConstraintBtn.Bind(wx.EVT_BUTTON,self.addConstraint)
      self.ConstraintBtn.SetToolTipString("Add a new constraint to the simulation")
      self.Sizer.Add(self.ConstraintBtn,(0,0),(1,1))
      #self.Layout()
      print 'constraint button done'
      #Save Constraints button
      self.SaveConstraintsBtn = wx.Button(self,-1,label="Save Constraints")
      self.SaveConstraintsBtn.SetForegroundColour("#000000")
      self.SaveConstraintsBtn.SetFont(wx.Font(10,wx.DEFAULT,wx.ITALIC,wx.BOLD))
      self.SaveConstraintsBtn.Bind(wx.EVT_BUTTON,self.saveConstraints)
      self.SaveConstraintsBtn.SetToolTipString("Save current constraints to a file for use in later simulations")
      self.Sizer.Add(self.SaveConstraintsBtn,(0,1),(1,1))
      #self.Layout()
      print 'save button'
      #Load Constraints file button
      self.LoadConstraintsBtn = wx.Button(self,-1,label="Load Constraints")
      self.LoadConstraintsBtn.SetForegroundColour("#000000")
      self.LoadConstraintsBtn.SetFont(wx.Font(10,wx.DEFAULT,wx.ITALIC,wx.BOLD))
      self.LoadConstraintsBtn.Bind(wx.EVT_BUTTON,self.loadConstraints)
      self.LoadConstraintsBtn.SetToolTipString("Load constraints from a file")
      self.Sizer.Add(self.LoadConstraintsBtn,(0,2),(1,1))
      #self.Layout()
      print 'load button'
      self.Cancelables = []
      self.CurrentConstraint = {}

      #Remove Constraint Button

      #Clear Constraints Button


      #Constraints Grid
      print 'starting constraints grid'
      self.constraintsGrid = wx.grid.Grid(self)
      self.constraintsGrid.CreateGrid(0,6)
      self.constraintsGrid.SetLabelFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD))
      self.Sizer.Add(self.constraintsGrid,(6,0),(10,6))
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
#==============================================================================
#       self.constraintTypeText=wx.StaticText(self,-1,'Constraint Type')
#       self.constraintTypeText.SetFont(wx.Font(10,wx.DEFAULT,wx.ITALIC,wx.BOLD))
#       self.Sizer.Add(self.constraintTypeText,(1,0),(1,1))
#       self.Cancelables.append(self.constraintTypeText)
#       #self.Layout()
#==============================================================================
      constraintTypes = ['Constraint Type','AtomPair','Angle','Dihedral','Coordinate']
      self.constraintTypeMenu = wx.ComboBox(self,choices=constraintTypes,style=wx.CB_READONLY)
      self.constraintTypeMenu.Bind(wx.EVT_COMBOBOX,self.setConstraintType)
      self.Sizer.Add(self.constraintTypeMenu,(1,0),(1,1))
      self.Cancelables.append(self.constraintTypeMenu)
      #self.Layout()
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
      self.CurrentConstraint = {}
      self.Layout()

    def setConstraintType(self,event):
      logInfo('Constraint type selected!')
      if self.constraintTypeMenu.GetStringSelection() == 'Choose Type':
        event.skip()
      try:
        self.CurrentConstraint['Constraint Type']=self.constraintTypeMenu.GetStringSelection()
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
#==============================================================================
#         self.PdbText = wx.StaticText(self,-1,'PDB')
#         self.PdbText.SetFont(wx.Font(10,wx.DEFAULT,wx.ITALIC,wx.BOLD))
#         self.Cancelables.append(self.PdbText)
#         self.Sizer.Add(self.PdbText,(2,0),(1,1))
#         #self.Layout()
#         print 'pdbtext'
#==============================================================================
        self.PdbMenu = wx.ComboBox(self,choices=pdbs,style=wx.CB_READONLY)
        self.PdbMenu.Bind(wx.EVT_COMBOBOX,self.setConstraintPDB)
        self.Cancelables.append(self.PdbMenu)
        self.Sizer.Add(self.PdbMenu,(1,1),(1,1))
        self.Layout()
        print 'pdbmenu'
      except Exception as e:
        print(e.message)

    def setConstraintPDB(self,event):
      logInfo("constraint PDB set!")
      self.CurrentConstraint['PDB']=self.PdbMenu.GetStringSelection()
      items = ['Select residue 1']
      for [indx, r, seqpos, poseindx, chainoffset, minType] in self.minPanel.minmap:
        print indx, r, seqpos, poseindx, chainoffset
        chain = self.seqWin.IDs[r][len(self.seqWin.IDs[r])-1]
        if chain == '_':
          chain = ' '
        print self.seqWin.poses[poseindx][0]
        self.CurrentConstraint['poseindx']=poseindx
        print chain
        chain_structure = self.seqWin.poses[poseindx][0][chain]
        print chain_structure.get_id()
        residue = chain_structure[int(seqpos)].resname
        print residue
        #Only considering standard amino acids for the moment
        if residue in "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR":
          residue += str(seqpos)
          residue += chain
          print residue
          items.append(residue)
      print items
      self.Residue1Menu = wx.ComboBox(self,choices=items,style=wx.CB_READONLY)
      self.Residue1Menu.Bind(wx.EVT_COMBOBOX,self.setAtom1Items)
      self.Cancelables.append(self.Residue1Menu)
#==============================================================================
#       self.Res1Text = wx.StaticText(self,-1,'Res 1')
#       self.Res1Text.SetFont(wx.Font(10,wx.DEFAULT,wx.ITALIC,wx.BOLD))
#       self.Cancelables.append(self.Res1Text)
#       self.Sizer.Add(self.Res1Text,(3,0),(1,1))
#==============================================================================
      self.Sizer.Add(self.Residue1Menu,(2,0),(1,1))
      self.Layout()
#==============================================================================
#       self.Atom1Text=wx.StaticText(self,-1,'Res 1 Atom',style=wx.ALIGN_RIGHT)
#       self.Atom1Text.SetFont(wx.Font(10,wx.DEFAULT,wx.ITALIC,wx.BOLD))
#       self.Cancelables.append(self.Atom1Text)
#==============================================================================
      atoms = ['Select Atom 1']
      self.Atom1Menu = wx.ComboBox(self,choices=atoms,style=wx.CB_READONLY)
      self.Atom1Menu.Bind(wx.EVT_COMBOBOX,self.ConstraintSpecifics)
      self.Cancelables.append(self.Atom1Menu)
#      self.Sizer.Add(self.Atom1Text,(3,2),(1,1))
      self.Sizer.Add(self.Atom1Menu,(2,1),(1,1))
      self.Layout()

    def setAtom1Items(self,event):
      residue = self.Residue1Menu.GetStringSelection()
      self.CurrentConstraint['Atom1_ResNum']=residue[3:len(residue)]
      atoms = self.getAtoms(residue)
      print atoms
      self.Atom1Menu.AppendItems(atoms)

    def getAtoms(self,residue):
      results = []
      poseindx = self.CurrentConstraint['poseindx']
      chain = residue[len(residue)-1]
      seqpos = int(residue[3:len(residue)-1])
      for atom in self.seqWin.poses[poseindx][0][chain][seqpos]:
        print atom
        results.append(atom.get_fullname())
      return results


    def ConstraintSpecifics(self,event):
      logInfo('Constraint Specifics!')
      print 'Constraint Specifics'
      #check which constraint method was selected and behave accordingly
      #AtomPair
      #Angle
      #Dihedral
      #Coordinate



    def gridClick(self,event):
      logInfo('Constraints Grid Clicked!')

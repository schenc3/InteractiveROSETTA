#imports
import wx
import wx.grid
import wx.lib.scrolledpanel
import os
import os.path
import time
import platform
import math
import multiprocessing
import Bio.PDB
import webbrowser
import datetime
import gzip
import numpy
from threading import Thread
from tools import *

def cleanPDB(pdbFile):
    '''In order to dock DNA, we must make sure that the pdb file meets a few criteria:
    1. All DNA chains must come AFTER all protein chains
    2. DNA residues must be identified as either ADE GUA THY CYT or A G T C,
    NOT DA DG DT DC <-- may be taken care of by BioPython
    This function will read through the PDB file and fix any of these problems,
    outputting a properly formatted file with the suffix _clean.'''
    # setup stuff
    chains = {} #tracks the chains in the pdb
    Dto3 = {' DA':'ADE',' DT':'THY',' DC':'CYT',' DG':'GUA'}
    DNA = [] #List of DNA chains
    read = open(pdbFile,'r')
    for line in read:
        #read the ATOM lines of the PDB file
        if line[0:6] == 'ATOM  ':
            chain = line[21]
            residue = line[17:20]
            #If it's a bad formatted DNA residue
            if residue in [' DA',' DT',' DC',' DG']:
                residue = Dto3[residue] #Fix it!
            #Add chain to DNA if it contains DNA
            if residue in ['  A','  T','  C','  G','ADE','THY','CYT','GUA']:
                if chain not in DNA:
                    DNA.append(chain)
            #Add line to chains[chain] and create chains[chain] if necessary
            try:
                chains[chain].append('%s%s%s'%(line[0:17],residue,line[20:]))
            except KeyError:
                chains[chain] = []
                chains[chain].append('%s%s%s'%(line[0:17],residue,line[20:]))

    read.close()
    goToSandbox()
    #We have to have fun with slicing to get rid of the .pdb since we don't want
    #example.pdb_clean.pdb
    output = open("%s_clean.pdb"%(pdbFile[:len(pdbFile)-4]))
    #write out all none-DNA chains
    for chain in chains:
        if chain not in DNA:
            for line in chains[chain]:
                output.write(line)
    #write out all the DNA chains
    for chain in DNA:
        for line in chains[chain]:
            output.write(line)
    output.close()

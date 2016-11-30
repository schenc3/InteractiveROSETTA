# This is for spawned PyRosetta processes that do the Rosetta intensive computations
# Why is a new process being spawned -> because Python's Global Interpreter Lock only let's one thread
# perform computations at any time
# Every time a Rosetta job needs to be run, it has to be run in a thread to prevent the GUI from locking up
# The GUI locks up because the job takes a long time to finish and GUI events cannot be handled while it's
# waiting for Rosetta to finish
# Problem is that Rosetta is so computationally heavy that, since the GIL is only allowing one thread to run
# at a time, it hogs all the CPU and the GUI hangs anyway, even in a thread
# So Rosetta has to be run in a separate process to prevent the GUI from hanging

# Only UNIX systems, this is not a big deal because parent processes can fork into children that address the
# same memory space as the parent
# On Windows, unfortunately the children have their own memory space, so all the Rosetta stuff has to be loaded
# in memory in the child process' memory
# However, a working PyRosetta instance uses up to 2.5GB of memory on Windows!!!!!!!!!! (only 1GB on Linux)
# DO NOT ATTEMPT TO SCORE POSES IN THE MAIN GUI OR DO ANYTHING ROSETTA RELATED EXCEPT LOADING AND ACCESSING
# POSES OR YOU WILL BE USING AT LEAST 5GB OF MEMORY WHEN A ROSETTA CHILD IS SPAWNED WHICH WILL SLOW MOST
# COMPUTERS DOWN TREMENDOUSLY SINCE THEY WILL START HAVING TO USE A LOT OF SWAP SPACE!!!!

import time
import os
import os.path
import sys
import traceback
import datetime
import platform
import glob
import gzip
import math
import subprocess
import commands
from threading import Thread
import socket
import shutil
try:
    # Try to import Rosetta
    from pyrosetta import *
    from rosetta import *
    # Extra imports for KIC
    from rosetta.protocols.loops.loop_mover.perturb import *
    from rosetta.protocols.loops.loop_mover.refine import *
except:
    # If it failed, then try to find Rosetta
    # If this already happened once already, then we should have saved the Rosetta path, so let's try to import from there
    print "Rosetta could not be imported.  Attempting to locate the PyRosetta install.  Please be patient..."
    cfgfile = "rosetta.cfg"
    try:
        f = open(cfgfile.strip(), "r")
        rosettadir = "Not Found"
        rosettadb = "Not Found"
        for aline in f:
            if ("[ROSETTAPATH]" in aline):
                rosettapath = aline.split("\t")[1].strip()
            if ("[ROSETTADB]" in aline):
                rosettadb = aline.split("\t")[1].strip()
        f.close()
        if (rosettapath == "Not Found"):
            raise Exception
        else:
            sys.path.append(rosettapath)
            olddir = os.getcwd()
            os.chdir(rosettapath)
            os.environ["PYROSETTA_DATABASE"] = rosettadb
            # Try to import Rosetta
            from rosetta import *
            from pyrosetta import *
            # Extra imports for KIC
            from rosetta.protocols.loops.loop_mover.perturb import *
            from rosetta.protocols.loops.loop_mover.refine import *
            os.chdir(olddir)
            print "Found Rosetta at " + rosettapath.strip() + "!"
            print "Rosetta Database: " + rosettadb.strip()
    except:
        # The error may have been the Rosetta import, which means the file needs to be closed
        try:
            f.close()
        except:
            pass
        # Okay, we still didn't get it, so let's traverse the filesystem looking for it...
        foundIt = False
        if (platform.system() == "Windows"):
            root = "C:\\"
        else:
            root = "/"
        for dpath, dnames, fnames in os.walk(root):
            try:
                if (platform.system() == "Windows"):
                    try:
                        indx = fnames.index("rosetta.pyd") # 64bit
                    except:
                        indx = fnames.index("rosetta.dll") # 32bit
                else:
                    indx = dnames.index("rosetta")
                    files = glob.glob(dpath + "/rosetta/*libmini*")
                    if (len(files) == 0):
                        raise Exception
            except:
                continue
            foundIt = True
            rosettapath = dpath
            for dname in dnames:
                if ("database" in dname):
                    rosettadb = dpath + "/" + dname
                    break
            break
        if (foundIt):
            sys.path.append(rosettapath)
            olddir = os.getcwd()
            os.chdir(rosettapath)
            os.environ["PYROSETTA_DATABASE"] = rosettadb
            try:
                # Try to import Rosetta
                from rosetta import *
                from pyrosetta import *
                # Extra imports for KIC
                from rosetta.protocols.loops.loop_mover.perturb import *
                from rosetta.protocols.loops.loop_mover.refine import *
                # Now let's save these paths so the next time this gets started we don't have to traverse the filesystem again
                data = []
                try:
                    f = open(cfgfile, "r")
                    for aline in f:
                        if (not("[ROSETTAPATH]" in aline) and not("[ROSETTADB]") in aline):
                            data.append(aline.strip())
                    f.close()
                except:
                    pass
                f = open(cfgfile, "w")
                for aline in data:
                    f.write(aline + "\n")
                f.write("[ROSETTAPATH]\t" + rosettapath.strip() + "\n")
                f.write("[ROSETTADB]\t" + rosettadb.strip() + "\n")
                f.close()
                print "Found Rosetta at " + rosettapath.strip() + "!"
                print "Rosetta Database: " + rosettadb.strip()
                os.chdir(olddir)
            except:
                print "PyRosetta cannot be found on your system!"
                print "Until you install PyRosetta, you may only use InteractiveROSETTA to visualize structures in PyMOL"
                exit()

doKICLocally = True # Set to true because the refined KIC step takes a long time, so just have the user do it for now
hostname = socket.gethostname()

def AA3to1(resn):
    indx3 = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR ".find(resn)
    indx = indx3 / 4
    return "ACDEFGHIKLMNPQRSTVWY"[indx]

def initializeRosetta(inputfile):
    # Grab the params files in the inputfile
    # Clean up from last time
    paramsFiles = glob.glob(hostname + "-*.params")
    for params in paramsFiles:
        os.remove(params)
    # Get all the file data for the params files out
    f = open(inputfile, "r")
    readingData = False
    for aline in f:
        if (aline[0:6] == "PARAMS"):
            paramsfile = aline.split("\t")[1].strip()
            f2 = open(hostname + "-" + paramsfile, "w")
        elif (aline[0:17] == "BEGIN PARAMS DATA"):
            readingData = True
        elif (aline[0:15] == "END PARAMS DATA"):
            readingData = False
            f2.close()
        elif (readingData):
            f2.write(aline.strip() + "\n")
    f.close()
    paramsstr = ""
    paramsFiles = glob.glob(hostname + "-*.fa.params")
    for params in paramsFiles:
        paramsstr = paramsstr + params.strip() + " "
    if (len(paramsstr) > 0):
        paramsstr = "-extra_res_fa " + paramsstr.strip()
    init(extra_options=paramsstr + " -ignore_unrecognized_res")
    #goToSandbox()

captured_stdout = ""
stdout_pipe = None
def drain_pipe():
    global captured_stdout
    while True:
        data = os.read(stdout_pipe[0], 1024)
        if not data:
            break
        captured_stdout += data

def doMinimization(inputfile):
    try:
        f = open(inputfile, "r")
    except:
        raise Exception("ERROR: The file \"" + inputfile + "\" is missing!")
    jobs = []
    minmap = []
    readingData = False
    ID = inputfile.split("-")[2] # 0 is the hostname, 1 is jobtype, and 2 is the ID
    for aline in f:
        if (aline[0:3] == "JOB"):
            [pdbfile, strstart, strend] = aline.split("\t")[1:]
            jobs.append([pdbfile.strip(), int(strstart), int(strend)])
            f2 = open(hostname + "-" + pdbfile, "w")
        elif (aline[0:6] == "MINMAP"):
            [strindx, strr, strseqpos, strp, strco, mtype] = aline.split("\t")[1:]
            minmap.append([int(strindx), int(strr), int(strseqpos), int(strp), int(strco), mtype.strip()])
        elif (aline[0:7] == "MINTYPE"):
            minType = aline.split("\t")[1].strip()
        elif (aline[0:8] == "SCOREFXN"):
            weightsfile = aline.split("\t")[1].strip()
        elif (aline[0:14] == "BEGIN PDB DATA"):
            readingData = True
        elif (aline[0:12] == "END PDB DATA"):
            f2.close()
            readingData = False
        elif (aline[0:19] == "BEGIN SCOREFXN DATA"):
            weightsfile = inputfile[inputfile.find("-")+1:] + "-weights"
            f2 = open(weightsfile, "w")
            readingData = True
        elif (aline[0:17] == "END SCOREFXN DATA"):
            f2.close()
            readingData = False
        elif (readingData):
            f2.write(aline.strip() + "\n")
    f.close()
    # Grab the params files in the input file
    initializeRosetta(inputfile)
    scorefxn = ScoreFunction()
    try:
        scorefxn.add_weights_from_file(weightsfile)
    except:
        raise Exception("ERROR: The scoring function weights could not be initialized!")
    f = open(hostname + "-minimizeoutputtemp", "w")
    for [pdbfile, minmapstart, minmapend] in jobs:
        minpose = pose_from_pdb(hostname + "-" + pdbfile)
        mm = MoveMap()
        mm.set_bb(False)
        mm.set_chi(False)
        for [indx, r, seqpos, p, co, mtype] in minmap[minmapstart:minmapend]:
            if (mtype == "BB" or mtype == "BB+Chi"):
                mm.set_bb(indx+1, True)
            if (mtype == "Chi" or mtype == "BB+Chi"):
                mm.set_chi(indx+1, True)
        minmover = MinMover(mm, scorefxn, "dfpmin", 0.01, True)
        if (minType == "Cartesian"):
            minmover.cartesian(True)
        try:
            minmover.apply(minpose)
        except:
            raise Exception("ERROR: The Rosetta minimizer failed!")
        outputpdb = pdbfile.split(".pdb")[0] + "_M.pdb"
        minpose.dump_pdb(hostname + "-" + outputpdb)
        f.write("OUTPUT\t" + outputpdb + "\n")
        f2 = open(hostname + "-" + outputpdb, "r")
        f.write("BEGIN PDB DATA\n")
        for aline in f2:
            f.write(aline.strip() + "\n")
        f.write("END PDB DATA\n")
        f2.close()
        os.remove(hostname + "-" + outputpdb)
        nonzero_scoretypes = scorefxn.get_nonzero_weighted_scoretypes()
        f.write("ENERGY\ttotal_score")
        for scoretype in nonzero_scoretypes:
            f.write("\t" + str(scoretype))
        f.write("\n")
        for res in range(1, minpose.n_residue()+1):
            f.write("ENERGY\t" + str(minpose.energies().residue_total_energy(res)))
            emap = minpose.energies().residue_total_energies(res)
            for scoretype in nonzero_scoretypes:
                f.write("\t" + str(emap.get(scoretype)))
            f.write("\n")
        # Remove the PDB file to keep the server clean
        try:
            os.remove(hostname + "-" + pdbfile)
        except:
            pass
    f.close()
    # So the main GUI doesn't attempt to read the file before the daemon finishes writing its contents
    os.rename(hostname + "-minimizeoutputtemp", "results/minimizeoutput-" + ID)
    os.remove(weightsfile)

def doFixbb(inputfile):
    try:
        f = open(inputfile, "r")
    except:
        raise Exception("ERROR: The file \"" + inputfile + "\" is missing!")
    readingData = False
    ID = inputfile.split("-")[2] # 0 is the hostname, 1 is jobtype, and 2 is the ID
    # Get the pdbfile, resfile, and scorefxn from the input file
    for aline in f:
        if (aline[0:7] == "PDBFILE"):
            pdbfile = aline.split("\t")[1].strip()
            f2 = open(hostname + "-" + pdbfile, "w")
        elif (aline[0:7] == "RESFILE"):
            resfile = aline.split("\t")[1].strip()
        elif (aline[0:8] == "SCOREFXN"):
            weightsfile = aline.split("\t")[1].strip()
        elif (aline[0:14] == "BEGIN PDB DATA"):
            readingData = True
        elif (aline[0:12] == "END PDB DATA"):
            f2.close()
            readingData = False
        elif (aline[0:19] == "BEGIN SCOREFXN DATA"):
            weightsfile = inputfile[inputfile.find("-")+1:] + "-weights"
            f2 = open(weightsfile, "w")
            readingData = True
        elif (aline[0:17] == "END SCOREFXN DATA"):
            f2.close()
            readingData = False
        elif (aline[0:18] == "BEGIN RESFILE DATA"):
            resfile = hostname + "-fixbb.resfile"
            f2 = open(resfile, "w")
            readingData = True
        elif (aline[0:16] == "END RESFILE DATA"):
            f2.close()
            readingData = False
        elif (readingData):
            f2.write(aline.strip() + "\n")
    f.close()
    initializeRosetta(inputfile)
    # Initialize scoring function
    scorefxn = ScoreFunction()
    try:
        print weightsfile
        scorefxn.add_weights_from_file(weightsfile)
    except:
        raise Exception("ERROR: The scoring function weights could not be initialized!")
    f = open(hostname + "-designoutputtemp", "w")
    # Perform fixed backbone design
    pose = pose_from_pdb(hostname + "-" + pdbfile)
    design_pack = TaskFactory.create_packer_task(pose)
    parse_resfile(pose, design_pack, resfile)
    pack_mover = PackRotamersMover(scorefxn, design_pack)
    try:
        pack_mover.apply(pose)
    except:
        raise Exception("ERROR: The Rosetta packer failed!")
    # Now I am going to change the B-factors of the Nbb residues to either 0 for undesigned
    # residues or 100 for designed residues.  That way the user can easily see which residues
    # were designed by looking for red sequence colorings in the sequence viewer
    # Get the designed positions
    varipos = []
    try:
        f2 = open(resfile, "r")
    except:
        raise Exception("ERROR: The file " + resfile + " is missing!")
    for aline in f2:
        if (aline.find("PIKAA") >= 0 or aline.find("NOTAA") >= 0):
            # Not counting NATRO or NATAA as "designed" even though under NATAA the rotamer can change
            varipos.append([int(aline.split()[0]), aline.split()[1].strip()]) # [seqpos, chainID]
    f2.close()
    # Now iterate down the residues and change B-factors
    info = pose.pdb_info()
    for ires in range(1, pose.n_residue()+1):
        try:
            seqpos = int(info.number(ires))
            chain = info.chain(ires)
            if (len(chain.strip()) == 0):
                chain = "_"
            if ([seqpos, chain] in varipos):
                # Designed
                info.temperature(ires, 1, 100.0) # Atom indx 1 is the bb N
            else:
                info.temperature(ires, 1, 0.0)
        except:
            # Maybe this is an NCAA that doesn't have an atom indx of 1?  Don't crash if so
            pass
    # Dump the output
    outputpdb = pdbfile.split(".pdb")[0] + "_D.pdb"
    pose.dump_pdb(hostname + "-" + outputpdb)
    # Now write the output information for the main GUI
    f.write("OUTPUT\t" + outputpdb + "\n")
    f2 = open(hostname + "-" + outputpdb, "r")
    f.write("BEGIN PDB DATA\n")
    for aline in f2:
        f.write(aline.strip() + "\n")
    f.write("END PDB DATA\n")
    f2.close()
    os.remove(hostname + "-" + outputpdb)
    nonzero_scoretypes = scorefxn.get_nonzero_weighted_scoretypes()
    f.write("ENERGY\ttotal_score")
    for scoretype in nonzero_scoretypes:
        f.write("\t" + str(scoretype))
    f.write("\n")
    for res in range(1, pose.n_residue()+1):
        f.write("ENERGY\t" + str(pose.energies().residue_total_energy(res)))
        emap = pose.energies().residue_total_energies(res)
        for scoretype in nonzero_scoretypes:
            f.write("\t" + str(emap.get(scoretype)))
        f.write("\n")
    f.close()
    # Remove the PDB file to keep the server clean
    try:
        os.remove(hostname + "-" + pdbfile)
        os.remove(hostname + "-fixbb.resfile")
    except:
        pass
    # So the main GUI doesn't attempt to read the file before the daemon finishes writing its contents
    os.rename(hostname + "-designoutputtemp", "results/designoutput-" + ID)
    os.remove(weightsfile)

def doScore(inputfile):
    try:
        f = open(inputfile, "r")
    except:
        raise Exception("ERROR: The file \"" + inputfile + "\" is missing!")
    readingData = False
    ID = inputfile.split("-")[2] # 0 is the hostname, 1 is jobtype, and 2 is the ID
    # Get the pdbfile, resfile, and scorefxn from the input file
    for aline in f:
        if (aline[0:7] == "PDBFILE"):
            pdbfile = aline.split("\t")[1].strip()
            f2 = open(hostname + "-" + pdbfile, "w")
        elif (aline[0:8] == "SCOREFXN"):
            weightsfile = aline.split("\t")[1].strip()
        elif (aline[0:14] == "BEGIN PDB DATA"):
            readingData = True
        elif (aline[0:12] == "END PDB DATA"):
            f2.close()
            readingData = False
        elif (aline[0:19] == "BEGIN SCOREFXN DATA"):
            weightsfile = inputfile[inputfile.find("-")+1:] + "-weights"
            f2 = open(weightsfile, "w")
            readingData = True
        elif (aline[0:17] == "END SCOREFXN DATA"):
            f2.close()
            readingData = False
        elif (readingData):
            f2.write(aline.strip() + "\n")
    f.close()
    initializeRosetta(inputfile)
    # Initialize scoring function
    scorefxn = ScoreFunction()
    try:
        scorefxn.add_weights_from_file(weightsfile)
    except:
        raise Exception("ERROR: The scoring function weights could not be initialized!")
    f = open(hostname + "-scoreoutputtemp", "w")
    pose = pose_from_pdb(hostname + "-" + pdbfile)
    # Calculate energy
    total_E = scorefxn(pose)
    # Dump the output
    outputpdb = pdbfile.split(".pdb")[0] + "_S.pdb"
    pose.dump_pdb(hostname + "-" + outputpdb)
    # Now write the output information for the main GUI
    f.write("OUTPUT\t" + outputpdb + "\n")
    f2 = open(hostname + "-" + outputpdb, "r")
    f.write("BEGIN PDB DATA\n")
    for aline in f2:
        f.write(aline.strip() + "\n")
    f.write("END PDB DATA\n")
    f2.close()
    os.remove(hostname + "-" + outputpdb)
    f.write("TOTAL_E\t" + str(total_E) + "\n")
    nonzero_scoretypes = scorefxn.get_nonzero_weighted_scoretypes()
    f.write("ENERGY\ttotal_score")
    for scoretype in nonzero_scoretypes:
        f.write("\t" + str(scoretype))
    f.write("\n")
    info = pose.pdb_info()
    for res in range(1, pose.n_residue()+1):
        # Skip HETATMs/NCAAs
        if (not(pose.residue(res).name1() in "ACDEFGHIKLMNPQRSTVWY")):
            continue
        f.write("ENERGY\t" + str(pose.energies().residue_total_energy(res)))
        emap = pose.energies().residue_total_energies(res)
        for scoretype in nonzero_scoretypes:
            f.write("\t" + str(emap.get(scoretype)))
        f.write("\n")
        chain = info.chain(res)
        if (chain == " " or chain == ""):
            chain = "_"
        f.write("ID\t" + chain + ":" + pose.residue(res).name1() + str(info.number(res)) + "\n")
    f.close()
    # Remove the PDB file to keep the server clean
    try:
        os.remove(hostname + "-" + pdbfile)
    except:
        pass
    # So the main GUI doesn't attempt to read the file before the daemon finishes writing its contents
    os.rename(hostname + "-scoreoutputtemp", "results/scoreoutput-" + ID)
    os.remove(weightsfile)

def doRotamerSearch(inputfile):
    try:
        f = open(inputfile, "r")
    except:
        raise Exception("ERROR: The file \"" + inputfile + "\" is missing!")
    readingData = False
    ID = inputfile.split("-")[2] # 0 is the hostname, 1 is jobtype, and 2 is the ID
    # Get the pdbfile, resfile, and scorefxn from the input file
    for aline in f:
        if (aline[0:7] == "PDBFILE"):
            pdbfile = aline.split("\t")[1].strip()
            f2 = open(hostname + "-" + pdbfile, "w")
        elif (aline[0:8] == "SCOREFXN"):
            weightsfile = aline.split("\t")[1].strip()
        elif (aline[0:7] == "RESTYPE"):
            restype = aline.split("\t")[1].strip()
            resone = AA3to1(restype)
        elif (aline[0:6] == "SEQPOS"):
            data = aline.split("\t")[1].strip()
            chain = data[0]
            seqpos = int(data[3:])
        elif (aline[0:14] == "BEGIN PDB DATA"):
            readingData = True
        elif (aline[0:12] == "END PDB DATA"):
            f2.close()
            readingData = False
        elif (aline[0:19] == "BEGIN SCOREFXN DATA"):
            weightsfile = inputfile[inputfile.find("-")+1:] + "-weights"
            f2 = open(weightsfile, "w")
            readingData = True
        elif (aline[0:17] == "END SCOREFXN DATA"):
            f2.close()
            readingData = False
        elif (readingData):
            f2.write(aline.strip() + "\n")
    f.close()
    initializeRosetta(inputfile)
    # Initialize scoring function
    scorefxn = ScoreFunction()
    try:
        scorefxn.add_weights_from_file(weightsfile)
    except:
        raise Exception("ERROR: The scoring function weights could not be initialized!")
    f = open(hostname + "-rotameroutputtemp", "w")
    try:
        rsd_factory = pose_from_pdb("data/residues.pdb")
    except:
        raise Exception("ERROR: The file \"data/residues.pdb\" is missing!")
    pose = pose_from_pdb(hostname + "-" + pdbfile)
    info = pose.pdb_info()
    # Find the actual Rosetta residue index
    for indx in range(1, pose.n_residue()+1):
        ichain = info.chain(indx)
        iseqpos = int(info.number(indx))
        if (ichain == " " or ichain == ""):
            ichain = "_"
        if (ichain == chain and iseqpos == seqpos):
            break
    phi = pose.phi(indx)
    psi = pose.psi(indx)
    # Read the rotamers
    if (restype != "ALA" and restype != "GLY"):
        if (platform.system() == "Windows"):
            libfile = os.getenv("PYROSETTA_DATABASE") + "\\rotamer\\ExtendedOpt1-5\\" + restype.lower() + ".bbdep.rotamers.lib.gz"
        else:
            libfile = os.getenv("PYROSETTA_DATABASE") + "/rotamer/ExtendedOpt1-5/" + restype.lower() + ".bbdep.rotamers.lib.gz"
        chivals = []
        Elist = []
        rotanames = []
        rotlib = gzip.open(libfile)
        for aline in rotlib:
            if (restype == aline[0:3]):
                lphi = int(aline.split()[1])
                lpsi = int(aline.split()[2])
                chis = aline.split()[9:13]
                if (math.fabs(lphi-phi) < 10 and math.fabs(lpsi-psi) < 10):
                    chivals.append(chis)
                    Elist.append(0.0)
                    rotaname = resone + ":" + aline.split()[4] + aline.split()[5] + aline.split()[6] + aline.split()[7]
                    i = 1
                    while (True):
                        # The r1r2r3r4 code is not unique, so we need an extra counter at the end
                        # to make all the names unique
                        if (not((rotaname + ":" + str(i)) in rotanames)):
                            rotaname = rotaname + ":" + str(i)
                            break
                        i = i + 1
                    rotanames.append(rotaname)
    else:
        chivals = [["-999", "-999", "-999", "-999"]]
        Elist = [0.0]
        rotanames = [restype[0] + ":0000:0"]
    # Mutate to restype (can't use mutate_residue because of the memory problem on Windows)
    resindx = "ACDEFGHIKLMNPQRSTVWY".find(resone) + 2
    try:
        res_mutate = Residue(rsd_factory.residue(resindx))
        res_mutate.place(pose.residue(indx), pose.conformation(), True)
        pose.replace_residue(indx, res_mutate.clone(), True)
    except:
        raise Exception("ERROR: The new residue could not be mutated onto the structure.")
    # Search all the chi values and score them
    k = 0
    nonzero_scoretypes = scorefxn.get_nonzero_weighted_scoretypes()
    for chis in chivals:
        if (chis[0] != "-999"):
            for i in range(0, pose.residue(indx).nchi()):
                pose.residue(indx).set_chi(i+1, float(chis[i]))
        pose.energies().clear_energies() # Necessary otherwise scores are not recalculated
        scorefxn(pose)
        emap = pose.energies().residue_total_energies(indx)
        Elist[k] = [emap.get(total_score)]
        for scoretype in nonzero_scoretypes:
            Elist[k].append(emap.get(scoretype))
        k = k + 1
        # The following line must be here
        # Sometimes after doing multiple set_chis the structure gets really messed up (don't know if I am
        # doing something wrong or if it is a Rosetta bug, it only happens in special occasions) so we have
        # to revert back to the original structure so errors are not propagated
        pose.replace_residue(indx, res_mutate.clone(), True)
    # Now sort according to increasing score
    for i in range(0, len(chivals)-1):
        lowest = i
        for j in range(i+1, len(chivals)):
            if (Elist[j][0] < Elist[lowest][0]):
                lowest = j
        temp = Elist[i]
        Elist[i] = Elist[lowest]
        Elist[lowest] = temp
        temp = chivals[i]
        chivals[i] = chivals[lowest]
        chivals[lowest] = temp
        temp = rotanames[i]
        rotanames[i] = rotanames[lowest]
        rotanames[lowest] = temp
    # Dump the output
    outputpdb = pdbfile.split(".pdb")[0] + "_R.pdb"
    pose.dump_pdb(hostname + "-" + outputpdb)
    # Now write the output information for the main GUI
    f.write("OUTPUT\t" + outputpdb + "\n")
    f2 = open(hostname + "-" + outputpdb, "r")
    f.write("BEGIN PDB DATA\n")
    for aline in f2:
        f.write(aline.strip() + "\n")
    f.write("END PDB DATA\n")
    f2.close()
    os.remove(hostname + "-" + outputpdb)
    f.write("INDEX\t" + str(indx) + "\t")
    info = pose.pdb_info()
    if (len(info.chain(indx).strip()) == 0):
        f.write("_" + str(info.number(indx)) + "\n")
    else:
        f.write(info.chain(indx) + str(info.number(indx)) + "\n")
    for ichi in range(1, pose.residue(indx).nchi()+1):
        chiatoms = pose.residue(indx).chi_atoms()[ichi]
        f.write("CHIATOMS\t" + pose.residue(indx).atom_name(chiatoms[1]).strip() + " " + pose.residue(indx).atom_name(chiatoms[2]).strip() + " " + pose.residue(indx).atom_name(chiatoms[3]).strip() + " " + pose.residue(indx).atom_name(chiatoms[4]).strip() + "\n")
    f.write("ENERGY\ttotal_score")
    for scoretype in nonzero_scoretypes:
        f.write("\t" + str(scoretype))
    f.write("\n")
    for i in range(0, len(Elist)):
        f.write("NAME\t" + rotanames[i] + "\n")
        f.write("ENERGY\t" + str(Elist[i][0]))
        for j in range(1, len(Elist[i])):
            f.write("\t" + str(Elist[i][j]))
        f.write("\n")
        f.write("CHI\t" + chivals[i][0] + "\t" + chivals[i][1] + "\t" + chivals[i][2] + "\t" + chivals[i][3] + "\n")
    f.close()
    # Remove the PDB file to keep the server clean
    try:
        os.remove(hostname + "-" + pdbfile)
    except:
        pass
    # So the main GUI doesn't attempt to read the file before the daemon finishes writing its contents
    os.rename(hostname + "-rotameroutputtemp", "results/rotameroutput-" + ID)
    os.remove(weightsfile)

def doKIC(inputfile):
    # Make a new results folder and unpack the data there
    try:
        os.mkdir("results/" + ID)
    except:
        shutil.rmtree("results/" + ID, ignore_errors=True)
        os.mkdir("results/" + ID)
    os.chdir("results/" + ID)
    # Move the input file here
    os.rename("../../" + inputfile, inputfile.split("-")[1])
    # Submit it to the queue
    sub = subprocess.Popen(["python", "../../rosetta_submit.py", "kic", ID])
    # The Python script will handle everything from here
    # Get back to where we were
    os.chdir("../..")

def old_doKIC(inputfile, stage="Coarse"):
    try:
        f = open(inputfile, "r")
    except:
        raise Exception("ERROR: The file \"" + inputfile + "\" is missing!")
    readingData = False
    ID = inputfile.split("-")[2] # 0 is the hostname, 1 is jobtype, and 2 is the ID
    # Get the pdbfile, resfile, and scorefxn from the input file
    for aline in f:
        if (aline[0:7] == "PDBFILE"):
            # But for the fine grained step the pose comes from repacked.pdb
            if (stage == "Coarse"):
                pdbfile = aline.split("\t")[1].strip()
                f2 = open(hostname + "-" + pdbfile, "w")
            else:
                pdbfile = "repacked.pdb"
        elif (aline[0:8] == "SCOREFXN"):
            weightsfile = aline.split("\t")[1].strip()
        elif (aline[0:7] == "REMODEL"):
            loopType = aline[7:].strip()
        elif (aline[0:8] == "SEQUENCE"):
            sequence = aline.split("\t")[1].strip()
        elif (aline[0:9] == "LOOPBEGIN"):
            loopBegin = int(aline.split("\t")[1])
        elif (aline[0:7] == "LOOPEND"):
            loopEnd = int(aline.split("\t")[1])
        elif (aline[0:14] == "BEGIN PDB DATA" and stage == "Coarse"):
            readingData = True
        elif (aline[0:12] == "END PDB DATA" and stage == "Coarse"):
            f2.close()
            readingData = False
        elif (aline[0:19] == "BEGIN SCOREFXN DATA"):
            weightsfile = inputfile[inputfile.find("-")+1:] + "-weights"
            f2 = open(weightsfile, "w")
            readingData = True
        elif (aline[0:17] == "END SCOREFXN DATA"):
            f2.close()
            readingData = False
        elif (readingData):
            f2.write(aline.strip() + "\n")
    f.close()
    initializeRosetta(inputfile)
    # Initialize scoring function
    scorefxn = ScoreFunction()
    try:
        scorefxn.add_weights_from_file(weightsfile)
    except:
        raise Exception("ERROR: The scoring function weights could not be initialized!")
    pose = pose_from_pdb(hostname + "-" + pdbfile)
    if (loopType == "DE NOVO"):
        if (stage == "Coarse"):
            # Since this is a new sequence being added, we first have to delete all the residues
            # between the beginning and ending points
            for ires in range(loopEnd-1, loopBegin, -1):
                pose.delete_polymer_residue(ires)
            # Now we have to add the sequence using our nifty little "rsd_factory" pose
            # The residues will have coordinates in weird places but it doesn't matter because
            # KIC fixes that and puts them in the right place; they don't need to start out anywhere
            # near being right
            try:
                rsd_factory = pose_from_pdb("data/residues.pdb")
            except:
                raise Exception("ERROR: The file \"data/residues.pdb\" is missing!")
            offset = 0
            for AA in sequence.strip():
                indx = "ACDEFGHIKLMNPQRSTVWY".find(AA) + 2
                pose.append_polymer_residue_after_seqpos(Residue(rsd_factory.residue(indx)), loopBegin+offset, True)
                offset = offset + 1
        # Now maybe the sequence is longer than what was originally the length of the sequence
        # between start and end, so we need to recalculate the loop end
        loopEnd = loopBegin + len(sequence.strip()) + 1
    cutpoint = loopEnd
    loop = Loop(loopBegin, loopEnd, cutpoint, 0, 1)
    loops = Loops()
    loops.add_loop(loop)
    add_single_cutpoint_variant(pose, loop)
    set_single_loop_fold_tree(pose, loop)
    if (stage == "Coarse"):
        # Low res KIC
        sw = SwitchResidueTypeSetMover("centroid")
        try:
            sw.apply(pose)
        except:
            raise Exception("ERROR: The PDB could not be converted to centroid mode!")
        kic_perturb = LoopMover_Perturb_KIC(loops)
        kic_perturb.set_max_kic_build_attempts(1020)
        try:
            kic_perturb.apply(pose)
        except:
            raise Exception("ERROR: The coarse KIC perturber failed!")
        # High res KIC
        sw = SwitchResidueTypeSetMover("fa_standard")
        try:
            sw.apply(pose)
        except:
            raise Exception("ERROR: The PDB could not be converted back to fullatom mode from centroid mode!")
        # Dump it for the repacking daemon
        pose.dump_pdb(hostname + "-torepack.pdb")
    else:
        kic_refine = LoopMover_Refine_KIC(loops)
        try:
            kic_refine.apply(pose)
        except:
            raise Exception("ERROR: The KIC refiner failed!")
        outputpdb = pdbfile.split(".pdb")[0] + "_K.pdb"
        pose.dump_pdb(hostname + "-" + outputpdb)
        # Now score the pose so we have the energy information in the main GUI
        scorefxn(pose)
        f = open(hostname + "-kicoutputtemp", "w")
        f.write("OUTPUT\t" + outputpdb + "\n")
        f2 = open(hostname + "-" + outputpdb, "r")
        f.write("BEGIN PDB DATA\n")
        for aline in f2:
            f.write(aline.strip() + "\n")
        f.write("END PDB DATA\n")
        f2.close()
        os.remove(hostname + "-" + outputpdb)
        f.write("LOOPBEGIN\t" + str(loopBegin) + "\n")
        f.write("LOOPEND\t" + str(loopEnd) + "\n")
        nonzero_scoretypes = scorefxn.get_nonzero_weighted_scoretypes()
        f.write("ENERGY\ttotal_score")
        for scoretype in nonzero_scoretypes:
            f.write("\t" + str(scoretype))
        f.write("\n")
        for res in range(1, pose.n_residue()+1):
            f.write("ENERGY\t" + str(pose.energies().residue_total_energy(res)))
            emap = pose.energies().residue_total_energies(res)
            for scoretype in nonzero_scoretypes:
                f.write("\t" + str(emap.get(scoretype)))
            f.write("\n")
        f.close()
        # Remove the PDB file to keep the server clean
        try:
            os.remove(hostname + "-" + pdbfile)
        except:
            pass
        # So the main GUI doesn't attempt to read the file before the daemon finishes writing its contents
        os.rename(hostname + "-kicoutputtemp", "results/kicoutput-" + ID)
        os.remove(weightsfile)

def doRepack(scorefxninput=""):
    try:
        weightsfile = scorefxninput[scorefxninput.find("-")+1:] + "-weights"
        # Initialize scoring function
        scorefxn = ScoreFunction()
        scorefxn.add_weights_from_file(weightsfile)
    except:
        # Default to Talaris2013
        scorefxn = create_score_function("talaris2013")
    initializeRosetta(scorefxninput)
    # Repack
    pose = pose_from_pdb(hostname + "-repackme.pdb")
    packtask = standard_packer_task(pose)
    packtask.restrict_to_repacking()
    packmover = PackRotamersMover(scorefxn, packtask)
    try:
        packmover.apply(pose)
    except:
        raise Exception("ERROR: The Rosetta packer failed!")
    os.remove(hostname + "-repackme.pdb")
    pose.dump_pdb(hostname + "-repacked.pdb")

def doMSD(inputfile):
    # Unpack the files out of this long file
    f = open(inputfile, "r")
    ID = inputfile.split("-")[2]
    # Make a new results folder and unpack the data there
    try:
        os.mkdir("results/" + ID)
    except:
        shutil.rmtree("results/" + ID, ignore_errors=True)
        os.mkdir("results/" + ID)
    os.chdir("results/" + ID)
    readingData = False
    for aline in f:
        if (len(aline.strip()) == 0):
            continue
        elif (aline.startswith("BEGIN PDB") or aline.startswith("BEGIN CORRESPONDENCE") or aline.startswith("BEGIN 2RESFILE") or aline.startswith("BEGIN ENTITY") or aline.startswith("BEGIN STATES") or aline.startswith("BEGIN FITNESS") or aline.startswith("BEGIN FLAGS")):
            f2 = open(aline.split()[len(aline.split())-1].strip(), "w")
            readingData = True
        elif (aline.startswith("END PDB") or aline.startswith("END CORRESPONDENCE") or aline.startswith("END 2RESFILE") or aline.startswith("END ENTITY") or aline.startswith("END STATES") or aline.startswith("END FITNESS") or aline.startswith("END FLAGS") or aline.startswith("END PARAMS DATA")):
            f2.close()
            readingData = False
        elif (aline.startswith("PARAMS")):
            f2 = open(aline.split()[len(aline.split())-1].strip(), "w")
        elif (aline.startswith("BEGIN PARAMS DATA")):
            readingData = True
        elif (readingData):
            f2.write(aline.strip() + "\n")
    f.close()
    # Submit it to the queue
    sub = subprocess.Popen(["python", "../../rosetta_submit.py", "msd", ID])
    # The Python script will handle everything from here
    # Get back to where we were
    os.chdir("../..")

def doAntibody(inputfile):
    # Make a new results folder and unpack the data there
    try:
        os.mkdir("results/" + ID)
    except:
        shutil.rmtree("results/" + ID, ignore_errors=True)
        os.mkdir("results/" + ID)
    os.chdir("results/" + ID)
    # Move the input file here
    os.rename("../../" + inputfile, inputfile.split("-")[1])
    # Submit it to the queue
    sub = subprocess.Popen(["python", "../../rosetta_submit.py", "antibody", ID])
    # The Python script will handle everything from here
    # Get back to where we were
    os.chdir("../..")

def doDock(inputfile):
    # Make a new results folder and unpack the data there
    try:
        os.mkdir("results/" + ID)
    except:
        shutil.rmtree("results/" + ID, ignore_errors=True)
        os.mkdir("results/" + ID)
    os.chdir("results/" + ID)
    # Move the input file here
    os.rename("../../" + inputfile, inputfile.split("-")[1])
    # Submit it to the queue
    sub = subprocess.Popen(["python", "../../rosetta_submit.py", "dock", ID])
    # The Python script will handle everything from here
    # Get back to where we were
    os.chdir("../..")

def doPMutScan(inputfile):
    # Make a new results folder and unpack the data there
    try:
        os.mkdir("results/" + ID)
    except:
        shutil.rmtree("results/" + ID, ignore_errors=True)
        os.mkdir("results/" + ID)
    os.chdir("results/" + ID)
    # Move the input file here
    os.rename("../../" + inputfile, inputfile.split("-")[1])
    # Submit it to the queue
    sub = subprocess.Popen(["python", "../../rosetta_submit.py", "pmutscan", ID])
    # The Python script will handle everything from here
    # Get back to where we were
    os.chdir("../..")

def doBackrub(inputfile):
    # Make a new results folder and unpack the data there
    try:
        os.mkdir("results/" + ID)
    except:
        shutil.rmtree("results/" + ID, ignore_errors=True)
        os.mkdir("results/" + ID)
    os.chdir("results/" + ID)
    # Move the input file here
    os.rename("../../" + inputfile, inputfile.split("-")[1])
    # Submit it to the queue
    sub = subprocess.Popen(["python", "../../rosetta_submit.py", "backrub", ID])
    # The Python script will handle everything from here
    # Get back to where we were
    os.chdir("../..")

def doFlexPep(inputfile):
    # Make a new results folder and unpack the data there
    try:
        os.mkdir("results/" + ID)
    except:
        shutil.rmtree("results/" + ID, ignore_errors=True)
        os.mkdir("results/" + ID)
    os.chdir("results/" + ID)
    # Move the input file here
    os.rename("../../" + inputfile, inputfile.split("-")[1])
    # Submit it to the queue
    sub = subprocess.Popen(["python", "../../rosetta_submit.py", "flexpep", ID])
    # The Python script will handle everything from here
    # Get back to where we were
    os.chdir("../..")

def doThread(inputfile):
    # Make a new results folder and unpack the data there
    try:
        os.mkdir("results/" + ID)
    except:
        shutil.rmtree("results/" + ID, ignore_errors=True)
        os.mkdir("results/" + ID)
    os.chdir("results/" + ID)
    # Move the input file here
    os.rename("../../" + inputfile, inputfile.split("-")[1])
    # Submit it to the queue
    sub = subprocess.Popen(["python", "../../rosetta_submit.py", "thread", ID])
    # The Python script will handle everything from here
    # Get back to where we were
    os.chdir("../..")

def doCustom(inputfile, protocol):
    # Make a new results folder and unpack the data there
    try:
        os.mkdir("results/" + ID)
    except:
        shutil.rmtree("results/" + ID, ignore_errors=True)
        os.mkdir("results/" + ID)
    os.chdir("results/" + ID)
    # Move the input file here
    os.rename("../../" + inputfile, inputfile.split("-")[1])
    # Submit it to the queue

    sub = subprocess.Popen(["python", "../../rosetta_submit.py", protocol, ID])
    # The Python script will handle everything from here
    # Get back to where we were
    os.chdir("../..")

def writeError(msg, inputfile="None-"):
    # Open a file and write out the error message so the main GUI can tell the user what happened
    # The main GUI needs to check to see if an errreport gets generated and recover from the error
    ID = inputfile.split("-")[2]
    f = open("results/errreport-" + ID, "w")
    f.write(msg + "\n\n")
    f.write(traceback.format_exc() + "\n")
    f.close()
    try:
        os.remove(inputfile[inputfile.find("-")+1:] + "-weights")
    except:
        pass

begintime = datetime.datetime.now()
while (True):
    inputfiles = glob.glob("jobfiles/*")
    for inputfile in inputfiles:
        ID = inputfile.split("-")[len(inputfile.split("-"))-1]
        if (inputfile.startswith("jobfiles/" + hostname + "-minimizeinput")):
            if ("-" in hostname):
                # If the hostname has - in it, it will screw up the splitting later on, so fix it here
                os.rename(inputfile, "jobfiles/" + hostname.replace("-", "") + "-minimizeinput-" + ID)
                inputfile = "jobfiles/" + hostname.replace("-", "") + "-minimizeinput-" + ID
            print "Daemon starting minimization job..."
            try:
                doMinimization(inputfile)
                print "Daemon completed minimization job"
            except Exception as e:
                print "The daemon crashed while performing the minimization job!"
                writeError(e.message, inputfile)
            os.remove(inputfile)
        elif (inputfile.startswith("jobfiles/" + hostname + "-designinput")):
            if ("-" in hostname):
                # If the hostname has - in it, it will screw up the splitting later on, so fix it here
                os.rename(inputfile, "jobfiles/" + hostname.replace("-", "") + "-designinput-" + ID)
                inputfile = "jobfiles/" + hostname.replace("-", "") + "-designinput-" + ID
            print "Daemon starting fixbb design job..."
            try:
                doFixbb(inputfile)
                print "Daemon completed fixbb design job"
            except Exception as e:
                print "The daemon crashed while performing the fixbb job!"
                writeError(e.message, inputfile)
            os.remove(inputfile)
        elif (inputfile.startswith("jobfiles/" + hostname + "-scoreinput")):
            if ("-" in hostname):
                # If the hostname has - in it, it will screw up the splitting later on, so fix it here
                os.rename(inputfile, "jobfiles/" + hostname.replace("-", "") + "-scoreinput-" + ID)
                inputfile = "jobfiles/" + hostname.replace("-", "") + "-scoreinput-" + ID
            print "Daemon starting scoring job..."
            try:
                doScore(inputfile)
                print "Daemon completed scoring job"
            except Exception as e:
                print "The daemon crashed while performing the scoring job!"
                writeError(e.message, inputfile)
            os.remove(inputfile)
        elif (inputfile.startswith("jobfiles/" + hostname + "-rotamerinput")):
            if ("-" in hostname):
                # If the hostname has - in it, it will screw up the splitting later on, so fix it here
                os.rename(inputfile, "jobfiles/" + hostname.replace("-", "") + "-rotamerinput-" + ID)
                inputfile = "jobfiles/" + hostname.replace("-", "") + "-rotamerinput-" + ID
            print "Daemon starting rotamer searching job..."
            try:
                doRotamerSearch(inputfile)
                print "Daemon completed rotamer searching job"
            except Exception as e:
                print "The daemon crashed while performing the rotamer search job!"
                writeError(e.message, inputfile)
            os.remove(inputfile)
        elif (inputfile.startswith("jobfiles/" + hostname + "-coarsekicinput")):
            if ("-" in hostname):
                # If the hostname has - in it, it will screw up the splitting later on, so fix it here
                os.rename(inputfile, "jobfiles/" + hostname.replace("-", "") + "-coarsekicinput-" + ID)
                inputfile = "jobfiles/" + hostname.replace("-", "") + "-coarsekicinput-" + ID
            print "Daemon starting KIC job..."
            doKIC(inputfile)
            #if ("-" in hostname):
                # If the hostname has - in it, it will screw up the splitting later on, so fix it here
                #os.rename(inputfile, "jobfiles/" + hostname.replace("-", "") + "-coarsekicinput-" + ID)
                #inputfile = "jobfiles/" + hostname.replace("-", "") + "-coarsekicinput-" + ID
            #print "Daemon starting coarse KIC loop modeling job..."
            # This is code to pipe stdout to a variable called "captured_stdout" so I can
            # parse the standard output of the coarse KIC perturber and see if it failed to
            # build the loop.  Then I know the sequence was too short
            #stdout_fileno = sys.stdout.fileno()
            #stdout_save = os.dup(stdout_fileno)
            #stdout_pipe = os.pipe()
            #os.dup2(stdout_pipe[1], stdout_fileno)
            #os.close(stdout_pipe[1])
            #captured_stdout = ""
            #t = Thread(target=drain_pipe)
            #t.start()
            #crashed = False
            #try:
                # This function call has the potential to run indefinitely if it cannot find
                # a way to bridge the two endpoint residues (in a de novo loop model)
                # The main GUI is timing the daemon's response though and will kill it and display
                # an error if it doesn't finish before a timeout (usually 3 min)
                #doKIC(inputfile, "Coarse")
                #os.rename(hostname + "-torepack.pdb", hostname + "-repackme.pdb") # So the GUI sees it
                #print "Daemon completed coarse KIC loop modeling job"
            #except Exception as e:
                #print "The daemon crashed while performing the coarse KIC loop modeling job!"
                #os.remove(inputfile)
                #writeError(e.message, inputfile)
                #crashed = True
                # The daemon needs to be killed by the main GUI because the coarse KIC thread is
                # still running and will do so indefinitely
            # Clean up the thread that was reading the standard output
            #os.close(stdout_fileno)
            #t.join()
            # Clean up the pipe and restore the original stdout
            #os.close(stdout_pipe[0])
            #os.dup2(stdout_save, stdout_fileno)
            #os.close(stdout_save)
            #outputline = captured_stdout.split("\n")[len(captured_stdout.split("\n"))-10]
            #try:
                #if (not(crashed) and outputline.find("Attempting loop building") >= 0 and int(outputline.split()[len(outputline.split())-2]) >= 100):
                    #raise Exception("ERROR: The loop sequence is too short and cannot bridge the endpoint residues!")
            #except Exception as e:
                #print "The loop sequence was too short!"
                #os.remove(inputfile)
                #os.remove(hostname + "-repackme.pdb")
                #writeError(e.message, inputfile)
                #continue
            #if (crashed):
                #continue
            #print "Daemon starting rotamer repacking job..."
            #try:
                #doRepack(inputfile)
                #print "Daemon completed rotamer repacking job"
            #except Exception as e:
                #print "The daemon crashed while performing the rotamer repacking job!"
                #writeError(e.message, inputfile)
                #os.remove(inputfile)
                #os.remove("jobfiles/" + hostname + "-repackmetemp.pdb")
                #continue
            #if (doKICLocally):
                # Dump the repacked PDB file to an outputfile that the client can read
                #f = open(hostname + "-coarsekicoutputtemp", "w")
                #f2 = open(hostname + "-repacked.pdb", "r")
                #f.write("OUTPUT\trepackedtemp.pdb\n")
                #f.write("BEGIN PDB DATA\n")
                #for aline in f2:
                #    f.write(aline.strip() + "\n")
                #f.write("END PDB DATA\n")
                #f2.close()
                #f.close()
                #os.rename(hostname + "-coarsekicoutputtemp", "results/coarsekicoutput-" + inputfile.split("-")[2])
                #os.remove(hostname + "-repacked.pdb")
                #os.remove(inputfile)
            #else:
                #print "Daemon starting fine KIC loop modeling job..."
                #try:
                #    doKIC(inputfile, "Fine")
                #    print "Daemon completed fine KIC loop modeling job"
                #except Exception as e:
                #    print "The daemon crashed while performing the fine KIC loop modeling job!"
                #    writeError(e.message, inputfile)
                #os.remove(inputfile)
        elif (inputfile.startswith("jobfiles/" + hostname + "-msdinput")):
            if ("-" in hostname):
                # If the hostname has - in it, it will screw up the splitting later on, so fix it here
                os.rename(inputfile, "jobfiles/" + hostname.replace("-", "") + "-msdinput-" + ID)
                inputfile = "jobfiles/" + hostname.replace("-", "") + "-msdinput-" + ID
            print "Daemon starting multi-state design job..."
            doMSD(inputfile)
            os.remove(inputfile)
        elif (inputfile.startswith("jobfiles/" + hostname + "-antibodyinput")):
            if ("-" in hostname):
                # If the hostname has - in it, it will screw up the splitting later on, so fix it here
                os.rename(inputfile, "jobfiles/" + hostname.replace("-", "") + "-antibodyinput-" + ID)
                inputfile = "jobfiles/" + hostname.replace("-", "") + "-antibodyinput-" + ID
            print "Daemon starting antibody modeling job..."
            doAntibody(inputfile)
        elif (inputfile.startswith("jobfiles/" + hostname + "-coarsedockinput")):
            if ("-" in hostname):
                # If the hostname has - in it, it will screw up the splitting later on, so fix it here
                os.rename(inputfile, "jobfiles/" + hostname.replace("-", "") + "-coarsedockinput-" + ID)
                inputfile = "jobfiles/" + hostname.replace("-", "") + "-coarsedockinput-" + ID
            print "Daemon starting docking job..."
            doDock(inputfile)
        elif (inputfile.startswith("jobfiles/" + hostname + "-scaninput")):
            if ("-" in hostname):
                # If the hostname has - in it, it will screw up the splitting later on, so fix it here
                os.rename(inputfile, "jobfiles/" + hostname.replace("-", "") + "-scaninput-" + ID)
                inputfile = "jobfiles/" + hostname.replace("-", "") + "-scaninput-" + ID
            print "Daemon starting point mutant scanning job..."
            doPMutScan(inputfile)
        elif (inputfile.startswith("jobfiles/" + hostname + "-backrubinput")):
            if ("-" in hostname):
                # If the hostname has - in it, it will screw up the splitting later on, so fix it here
                os.rename(inputfile, "jobfiles/" + hostname.replace("-", "") + "-backrubinput-" + ID)
                inputfile = "jobfiles/" + hostname.replace("-", "") + "-backrubinput-" + ID
            print "Daemon starting backrubbing job..."
            doBackrub(inputfile)
        elif (inputfile.startswith("jobfiles/" + hostname + "-flexpepinput")):
            if ("-" in hostname):
                # If the hostname has - in it, it will screw up the splitting later on, so fix it here
                os.rename(inputfile, "jobfiles/" + hostname.replace("-", "") + "-flexpepinput-" + ID)
                inputfile = "jobfiles/" + hostname.replace("-", "") + "-flexpepinput-" + ID
            print "Daemon starting flexible peptide docking job..."
            doFlexPep(inputfile)
        elif (inputfile.startswith("jobfiles/" + hostname + "-threadinput")):
            if ("-" in hostname):
                # If the hostname has - in it, it will screw up the splitting later on, so fix it here
                os.rename(inputfile, "jobfiles/" + hostname.replace("-", "") + "-threadinput-" + ID)
                inputfile = "jobfiles/" + hostname.replace("-", "") + "-threadinput-" + ID
            print "Daemon starting comparative modeling job..."
            doThread(inputfile)
        elif (inputfile.startswith("jobfiles/" + hostname + "-killinput")):
            fin = open(inputfile, "r")
            filedata = fin.readlines()
            fin.close()
            commandline = "python killjob.py " + filedata[1].strip()
            res, output = commands.getstatusoutput(commandline)
            try:
                os.remove(inputfile)
            except:
                pass
        else:
            if ("-" in hostname):
                inputtype = inputfile.split("-")[len(inputfile.split("-"))-2]
                # If the hostname has - in it, it will screw up the splitting later on, so fix it here
                os.rename(inputfile, "jobfiles/" + hostname.replace("-", "") + "-" + inputtype + "-" + ID)
                inputfile = "jobfiles/" + hostname.replace("-", "") + "-" + inputtype + "-" + ID
            print "Daemon starting custom job..."
            doCustom(inputfile, inputtype.split("input")[0])
    time.sleep(1)
    elapsed = datetime.datetime.now() - begintime
    if (elapsed.seconds >= 60*60*24):
        # Restart the daemon to prevent memory from becoming stale
        break

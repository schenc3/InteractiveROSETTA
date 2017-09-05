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
import platform
import subprocess
import psutil
import glob
import gzip
import math
from tools import *
from wx import Timer
from cStringIO import StringIO
try:
    # Try to import Rosetta
    from rosetta import *
    # Extra imports for KIC
    from rosetta.protocols.loops.loop_mover.perturb import *
    from rosetta.protocols.loops.loop_mover.refine import *
    # Extra import for docking
    import rosetta.protocols.rigid as rigid_moves
    # Extra imports for threading
    import rosetta.protocols.evaluation
    import rosetta.protocols.backrub
    from rosetta.protocols.comparative_modeling import *
    from rosetta.protocols.jd2 import *
    from rosetta.core.scoring.constraints import *
    from rosetta.protocols.pmut_scan import *
except:
    # If it failed, then try to find Rosetta
    # If this already happened once already, then we should have saved the Rosetta path, so let's try to import from there
    print "Rosetta could not be imported.  Attempting to locate the PyRosetta install.  Please be patient..."
    if (platform.system() == "Windows"):
        cfgfile = os.path.expanduser("~") + "/InteractiveROSETTA/seqwindow.cfg"
    else:
        cfgfile = os.path.expanduser("~") + "/.InteractiveROSETTA/seqwindow.cfg"
    try:
        # The main GUI should have already found PyRosetta, so let's look in the config file first
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
            # Nope, wasn't in there
            raise Exception
        else:
            print "Trying to import from " + rosettapath + "..."
            sys.path.append(rosettapath)
            olddir = os.getcwd()
            os.chdir(rosettapath)
            os.environ["PYROSETTA_DATABASE"] = rosettadb
            # Try to import Rosetta
        try:
               from rosetta import *
        except:
            # Extra imports for KIC
            from rosetta.protocols.loops.loop_mover.perturb import *
            from rosetta.protocols.loops.loop_mover.refine import *
            # Extra import for docking
            import rosetta.protocols.rigid as rigid_moves
            # Extra imports for threading
            import rosetta.protocols.evaluation
            import rosetta.protocols.backrub
            from rosetta.protocols.comparative_modeling import *
            from rosetta.protocols.jd2 import *
            os.chdir(olddir)
            # If we got this far, then we successfully imported PyRosetta
            print "Found Rosetta at " + rosettapath.strip() + "!"
            print "Rosetta Database: " + rosettadb.strip()
    except:
        # The error may have been the Rosetta import, which means the file needs to be closed
        try:
            f.close()
        except:
            pass
        raise Exception()
        # Okay, we still didn't get it, so let's traverse the filesystem looking for it...
        foundIt = False
        # Remember that root is C: on Windows, not /
        if (platform.system() == "Windows"):
            root = "C:\\"
        else:
            root = "/"
        for dpath, dnames, fnames in os.walk(root):
            try:
                if (platform.system() == "Windows"):
                    # On Windows, we expect the PyRosetta folder to have either rosetta.pyd or rosetta.dll
                    # in it
                    try:
                        indx = fnames.index("rosetta.pyd") # 64bit
                    except:
                        indx = fnames.index("rosetta.dll") # 32bit
                else:
                    # On Mac/Linux, there should be a folder called "rosetta" in the PyRosetta folder that
                    # contains a shared object file called either libmini.dylib or libmini.so
                    indx = dnames.index("rosetta")
                    files = glob.glob(dpath + "/rosetta/*libmini*")
                    if (len(files) == 0):
                        raise Exception
            except:
                continue
            # If we got here, then we found a potential import
            foundIt = True
            rosettapath = dpath
            for dname in dnames:
                # Figure out what the database is called
                # It will either be "database" or "rosetta_database"
                if ("database" in dname):
                    rosettadb = dpath + "/" + dname
                    break
            break
        if (foundIt):
            # Let's try to import it again
            sys.path.append(rosettapath)
            olddir = os.getcwd()
            os.chdir(rosettapath)
            os.environ["PYROSETTA_DATABASE"] = rosettadb
            try:
                # Try to import Rosetta
                from rosetta import *
                # Extra imports for KIC
                from rosetta.protocols.loops.loop_mover.perturb import *
                from rosetta.protocols.loops.loop_mover.refine import *
                # Extra import for docking
                import rosetta.protocols.rigid as rigid_moves
                # Extra imports for threading
                import rosetta.protocols.evaluation
                import rosetta.protocols.backrub
                from rosetta.protocols.comparative_modeling import *
                from rosetta.protocols.jd2 import *
                # Now let's save these paths so the next time this gets started we don't have to traverse the filesystem again
                data = []
                try:
                    f = open(cfgfile, "r")
                    for aline in f:
                        if (not("[ROSETTAPATH]" in aline) and not("[ROSETTADB]" in aline)):
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
                # Still didn't work, the user will have to figure out what's wrong
                print "PyRosetta cannot be found on your system!"
                print "Until you install PyRosetta, you may only use InteractiveROSETTA to visualize structures in PyMOL"
                exit()
from rotation import *
from threading import Thread
from Bio.Align.Applications import MuscleCommandline
from tools import goToSandbox
from tools import AA3to1

# Function for initializing Rosetta with the parameters files in the params folder in the sandbox
def initializeRosetta(addOn="", extraMutes=""):
    # Grab the params files in the user's personal directory
    os.chdir("params")
    faparamsstr = ""
    faParamsFiles = glob.glob("*.fa.params")
    for params in faParamsFiles:
        faparamsstr = faparamsstr + "params/" + params.strip() + " "
    if (len(faParamsFiles) > 0):
        faparamsstr = "-extra_res_fa " + faparamsstr.strip()
    paramsstr = faparamsstr
    try:
        init(extra_options=paramsstr + " -mute core.kinematics.AtomTree " + extraMutes + " -ignore_unrecognized_res -ignore_zero_occupancy false " + addOn)
        os.chdir("..")
        cm = ChemicalManager.get_instance()
        cm.residue_type_set("fa_standard")
    except:
        print "Water detected already"
        # Take water out, PyRosetta found one elsewhere
        newparamsstr = ""
        for params in paramsstr.split():
            if ("HOH" in params):
                continue
            newparamsstr = newparamsstr + params + " "
        if (newparamsstr.strip() == "-extra_res_fa"):
            paramsstr = ""
        else:
            paramsstr = newparamsstr
        init(extra_options=paramsstr + " -mute core.kinematics.AtomTree " + extraMutes + " -ignore_unrecognized_res -ignore_zero_occupancy false " + addOn)

# This class is used to grab standard output from pmutscan, which is the only output of that protocol
class ScanCapturing(list):
    def __enter__(self):
        self._stdout = sys.stdout
        self.f = open("scanprogress", "w", 1)
        sys.stdout = self.f #self._stringio = StringIO()
        return self

    def __exit__(self, *args):
        self.f.close()
        sys.stdout = self._stdout

captured_stdout = ""
stdout_pipe = None
def drain_pipe():
    global captured_stdout
    while True:
        data = os.read(stdout_pipe[0], 1024)
        if not data:
            break
        print "mydata: " + data.strip()
        captured_stdout += data

# Function for performing energy minimization
def doMinimization():
    # Grab the params files in the user's personal directory
    initializeRosetta()
    import rosetta.protocols.simple_moves
    try:
        f = open("minimizeinput", "r")
    except:
        raise Exception("ERROR: The file \"minimizeinput\" is missing!")
    jobs = []
    minmap = []
    # Read the input data
    # We need:
    #    JOB: contains the pdbfile, and the range in the minmap that corresponds to that pdbfile
    #    MINMAP: A list of all the positions that will be minimized for all jobs and whether BB, Chi, or both will be minimized
    #    MINTYPE: Either minimize in Cartesian or torsion space
    #    SCOREFXN: The file containing the score function weights
    for aline in f:
        if (aline[0:3] == "JOB"):
            [pdbfile, strstart, strend,constraintfile] = aline.split("\t")[1:]
            jobs.append([pdbfile.strip(), int(strstart), int(strend),constraintfile])
        elif (aline[0:6] == "MINMAP"):
            [strindx, strr, strseqpos, strp, strco, mtype,strr_indx] = aline.split("\t")[1:]
            minmap.append([int(strindx), int(strr), int(strseqpos), int(strp), int(strco), mtype.strip(),int(strr_indx.strip())])
        elif (aline[0:7] == "MINTYPE"):
            minType = aline.split("\t")[1].strip()
        elif (aline[0:8] == "SCOREFXN"):
            weightsfile = aline.split("\t")[1].strip()
    f.close()
    # Initialize the scoring function
    scorefxn = ScoreFunction()
    try:
        scorefxn.add_weights_from_file(weightsfile)
    except:
        raise Exception("ERROR: The scoring function weights could not be initialized!")
    f = open("minimizeoutputtemp", "w")
    pyobs = PyMOL_Observer()
    for [pdbfile, minmapstart, minmapend,constraintFile] in jobs:
        # Get a Rosetta pose of the file
        minpose = pose_from_pdb(pdbfile)
        minpose.pdb_info().name("protocol_view")
        pyobs.add_observer(minpose)
        pyobs.pymol.update_energy(True)
        areCST = False
        #Configure Constraints
        if constraintFile.strip() != '':
          try:
                areCST = True
                # print "THERE ARE CONSTRAINTS"
                goToSandbox()
                cstfile = open(constraintFile.strip(),'r')
                # for line in cstfile:
                  # print line
                cstfile.close()
                cstMover = rosetta.protocols.simple_moves.ConstraintSetMover()
                # print "cst file is %s!"%(constraintFile.strip())
                cstMover.constraint_file(constraintFile.strip())
                for cst in [atom_pair_constraint, angle_constraint, dihedral_constraint, coordinate_constraint, constant_constraint]:
                  scorefxn.set_weight(cst,1.0)
          except:
            # print 'cst failed.'
            areCST = False
#          cstMover.apply(minpose)
        # Create the minmap
        mm = MoveMap()
        mm.set_bb(False)
        mm.set_chi(False)
        for [indx, r, seqpos, p, co, mtype,r_indx] in minmap[minmapstart:minmapend]:
            if (mtype == "BB" or mtype == "BB+Chi"):
                #mm.set_bb(indx+1, True)
                mm.set_bb(r_indx,True)
            if (mtype == "Chi" or mtype == "BB+Chi"):
                #mm.set_chi(indx+1, True)
                mm.set_chi(r_indx,True)
        # Minimize using dfpmin always
        minmover = MinMover(mm, scorefxn, "dfpmin", 0.01, True)
#        if constraintFile.strip() != '':
#          minmover.set_cst_file(constraintFile.strip())
        if (minType == "Cartesian"):
            minmover.cartesian(True)
            scorefxn.set_weight(cart_bonded, 0.5)
            scorefxn.set_weight(pro_close, 0)
        try:
            if areCST:
              # print 'areCST!'
              try:
                cstMover.apply(minpose)
              except:
                # print 'pass!'
                pass
            minmover.apply(minpose)
        except:
            raise Exception("ERROR: The Rosetta minimizer failed!")
        # Dump the pdbfile
        outputpdb = pdbfile.split(".pdb")[0] + "_M.pdb"
        minpose.dump_pdb(outputpdb)
        # Update the output file telling the GUI the name of the minimized file and a list of all the energies
        # for display by color in PyMOL
        f.write("OUTPUT\t" + outputpdb + "\n")
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
    f.close()
    # So the main GUI doesn't attempt to read the file before the daemon finishes writing its contents
    renamed = False
    while not renamed:
      try:
        os.rename("minimizeoutputtemp", "minimizeoutput")
        renamed = True
      except:
        renamed = False

# Fixed backbone design protocol
def doFixbb():
    initializeRosetta()
    try:
        f = open("designinput", "r")
    except:
        raise Exception("ERROR: The file \"designinput\" is missing!")
    # Get the pdbfile, resfile, and scorefxn from the input file
    for aline in f:
        if (aline[0:7] == "PDBFILE"):
            pdbfile = aline.split("\t")[1].strip()
        elif (aline[0:7] == "RESFILE"):
            resfile = aline.split("\t")[1].strip()
        elif (aline[0:8] == "SCOREFXN"):
            weightsfile = aline.split("\t")[1].strip()
    f.close()
    # Initialize scoring function
    scorefxn = ScoreFunction()
    try:
        scorefxn.add_weights_from_file(weightsfile)
    except:
        raise Exception("ERROR: The scoring function weights could not be initialized!")
    f = open("designoutputtemp", "w")
    # Perform fixed backbone design
    pose = pose_from_pdb(pdbfile)
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
    pose.dump_pdb(outputpdb)
    # Now write the output information for the main GUI
    f.write("OUTPUT\t" + outputpdb + "\n")
    nonzero_scoretypes = scorefxn.get_nonzero_weighted_scoretypes()
    # Write out the energies of all residues so PyMOL can color by energy
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
    # So the main GUI doesn't attempt to read the file before the daemon finishes writing its contents
    os.rename("designoutputtemp", "designoutput")

# Point mutation scanning protocol
def doPMutScan():
    # We want to mute basic and core so only the protocols standard output is shown
    # We need to grab and parse this standard output and muting the other clears up the clutter
    initializeRosetta(extraMutes="basic core")
    try:
        f = open("scaninput", "r")
    except:
        raise Exception("ERROR: The file \"scaninput\" is missing!")
    # Looking for:
    #    PDBFILE: The structure to analyze
    #    LIST: A list of all single or double point mutants to consider
    #    SCOREFXN: The weights for the scoring function
    #    MUTANTTYPE: Whether we should do single or double mutants
    for aline in f:
        if (aline[0:7] == "PDBFILE"):
            pdbfile = aline.split("\t")[1].strip()
        elif (aline[0:4] == "LIST"):
            listfile = aline.split("\t")[1].strip()
        elif (aline[0:8] == "SCOREFXN"):
            weightsfile = aline.split("\t")[1].strip()
        elif (aline.startswith("MUTANTTYPE")):
            mutanttype = aline.split("\t")[1].strip()
    f.close()
    # Initialize scoring function
    scorefxn = ScoreFunction()
    try:
        scorefxn.add_weights_from_file(weightsfile)
    except:
        raise Exception("ERROR: The scoring function weights could not be initialized!")
    # Perform the scan using the specified mutants
    # First repack the structure otherwise we'll get a lot of hits simply because Rosetta liked repacking
    # the residues into its own energy minimum
    # We want real mutants so the structure has to be stabilized by Rosetta first
    pose = pose_from_pdb(pdbfile)
    design_pack = TaskFactory.create_packer_task(pose)
    design_pack.restrict_to_repacking()
    pack_mover = PackRotamersMover(scorefxn, design_pack)
    pack_mover.apply(pose)
    pose.dump_pdb(pdbfile)
    # The point mutant scanner wants the PDB files in a vector1<string> for some reason
    vec = utility.vector1_string()
    vec.append(pdbfile)
    if (mutanttype == "DOUBLE"):
        double = True
    else:
        double = False
    try:
        PMS = PointMutScanDriver(vec, double, listfile, False)
        PMS.go()
        # Before this function was called, a context manager was started to grab standard output from this
        # function, which is what we want
        # It is being saved in scanprogress
    except:
        pass
    # Append the model name to these data so when the report is saved the GUI can remember what model
    # these data came from
    f = open("scanprogress", "a")
    f.write("#MODEL\t" + pdbfile.strip() + "\n")
    f.close()
    # So the main GUI doesn't attempt to read the file before the daemon finishes writing its contents
    os.rename("scanprogress", "scanoutput")

# Rosetta relax function for generating mild ensembles
def doRelax():
    initializeRosetta()
    try:
        f = open("relaxinput", "r")
    except:
        raise Exception("ERROR: The file \"relaxinput\" is missing!")
    # Looking for:
    #    PDBFILE: The structure for which an ensemble will be generated
    #    NMODELS: The number of models to generate in the ensemble
    #    SCOREFXN: The weights for the scoring function
    for aline in f:
        if (aline[0:7] == "PDBFILE"):
            pdbfile = aline.split("\t")[1].strip()
        elif (aline[0:7] == "NMODELS"):
            nmodels = int(aline.split("\t")[1])
        elif (aline[0:8] == "SCOREFXN"):
            weightsfile = aline.split("\t")[1].strip()
    f.close()
    # Initialize scoring function
    scorefxn = ScoreFunction()
    try:
        scorefxn.add_weights_from_file(weightsfile)
    except:
        raise Exception("ERROR: The scoring function weights could not be initialized!")
    # Perform FastRelax for each model
    inputpose = pose_from_pdb(pdbfile)
    FR = FastRelax(scorefxn)
    # Write the results in an ensb file (basically a large gzipped PDB file with multiple models)
    archive = gzip.open("resultstemp.ensb", "wb")
    for i in range(0, nmodels):
        # Update the progress for the GUI
        f = open("relaxprogress", "w")
        f.write(str(i) + "\n")
        f.close()
        filename = "relax_" + str(i+1) + ".pdb"
        # Get a copy of the original pose
        pose = Pose(inputpose)
        FR.apply(pose)
        # Now save the coordinates in the archive
        pose.dump_pdb("temp.pdb")
        fin = open("temp.pdb", "r")
        archive.write("BEGIN PDB " + filename.strip() + "\n")
        for aline in fin:
            archive.write(aline)
        archive.write("END PDB " + filename.strip() + "\n")
        fin.close()
    archive.close()
    # So the main GUI doesn't attempt to read the file before the daemon finishes writing its contents
    os.rename("resultstemp.ensb", "results.ensb")

# Rosetta backrub protocol for generating more extensive ensembles
def doBackrub():
    initializeRosetta()
    try:
        f = open("backrubinput", "r")
    except:
        raise Exception("ERROR: The file \"relaxinput\" is missing!")
    # Looking for:
    #    PDBFILE: The structure for which an ensemble will be generated
    #    NMODELS: The number of models to generate for the ensemble
    #    MAXRMSD: The RMSD limit for the ensemble gradient
    #    SCOREFXN: The weights for the scoring function
    for aline in f:
        if (aline[0:7] == "PDBFILE"):
            pdbfile = aline.split("\t")[1].strip()
        elif (aline[0:7] == "NMODELS"):
            ntemplates = int(aline.split("\t")[1])
        elif (aline[0:7] == "MAXRMSD"):
            backrubRMSD = int(aline.split("\t")[1])
        elif (aline[0:8] == "SCOREFXN"):
            weightsfile = aline.split("\t")[1].strip()
    f.close()
    # Initialize scoring function
    scorefxn = ScoreFunction()
    try:
        scorefxn.add_weights_from_file(weightsfile)
    except:
        raise Exception("ERROR: The scoring function weights could not be initialized!")
    # Save the ensemble in an ensb file
    archive = gzip.open("resultstemp.ensb", "wb")
    for teamrank in range(0, ntemplates):
        # Update the progress for the GUI
        f = open("backrubprogress", "w")
        f.write(str(teamrank) + "\n")
        f.close()
        # Calculate the target RMSD value for this ensemble member (it is an even gradient between 0
        # and backrubRMSD)
        if (ntemplates == 1):
            myRMSD = float(backrubRMSD)
        else:
            myRMSD = float(teamrank) * (float(backrubRMSD) / float(ntemplates - 1))
        # This member will be generated with a constant seed, so we can easily get the same ensemble
        # given the same input parameters
        RMSDseed = "%7i" % (myRMSD * 1000) # Truncate the decimal
        RMSDseed = RMSDseed.strip()
        addOn = " -constant_seed -seed_offset 0 -use_bicubic_interpolation "
        initializeRosetta(addOn=addOn)
        pose = pose_from_pdb(pdbfile)
        # Allow backrub to potentially pivot all canonical residues
        BM = protocols.backrub.BackrubMover()
        pivot_res = utility.vector1_ulong()
        # Set all canonical AAs to pivot residues
        three = "ALACYSASPGLUPHEGLYHISILELYSLEUMETASNPROGLNARGSERTHRVALTRPTYR"
        for i in range(0, pose.n_residue()):
            if (three.find(pose.residue(i+1).name()) >= 0):
                pivot_res.append(i+1)
        BM.set_pivot_residues(pivot_res)
        # Now minimize the structure to get to an energy minimum in the selected scoring function
        # Use this structure as the reference for RMSD calculations
        mm = MoveMap()
        mm.set_bb(True)
        mm.set_chi(True)
        for i in range(1, pose.n_residue() + 1):
            if ("ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".find(pose.residue(i).name3()) < 0):
                mm.set_bb(i, False)
                mm.set_chi(i, False)
        tol = 0.01
        scorefxn_min = create_score_function("talaris2013")
        minmover = MinMover(mm, scorefxn_min, 'lbfgs_armijo_nonmonotone', tol, True)
        minmover.apply(pose)
        basepose = Pose(pose)
        if (teamrank > 0):
            # If this is template 1, use the inputted PDB and perform no backrub
            attempts = 0
            orig_pose = Pose(basepose)
            addOn = " -constant_seed -seed_offset 0 "
            initializeRosetta(addOn=addOn)
            minmover = MinMover(mm, scorefxn_min, 'lbfgs_armijo_nonmonotone', tol, True)
            addOn = " -constant_seed -seed_offset " + RMSDseed + " -use_bicubic_interpolation "
            initializeRosetta(addOn=addOn)
            # myRMSD is the target RMSD for this member
            # myRange is how close to the target value it needs to get to quit the following loops
            myRMSD = float(teamrank) * (float(backrubRMSD) / float(ntemplates - 1))
            myRange = (float(backrubRMSD) / float(ntemplates - 1)) / 2.0
            if (myRange < 0.025):
                myRange = 0.025 # If you're getting stuck at the backrub diversity step, then this value probably needs to be raised
            backrubTm = 10.0 ** myRMSD # Initial starting value designed to get the Tm close to where it should be to generate the desired RMSD quickly
            Metropolis = rosetta.MonteCarlo(pose, scorefxn, backrubTm)
            myLoop = 10
            lastRMSD = -1
            currRMSD = -1
            prev_pose = Pose(orig_pose)
            while (True):
                # Do some backrub moves to generate diversity
                for i in range(0, myLoop):
                    BM.apply(pose)
                    # If the move is rejected, pose is automatically reverted by this function
                    Metropolis.boltzmann(pose)
                # Minimize to fix clashes
                minmover = MinMover(mm, scorefxn_min, 'lbfgs_armijo_nonmonotone', tol, True)
                minmover.apply(pose)
                attempts = attempts + 1
                # If we've been at this backrub temperature for a while, increase it and do more backrub moves
                if (attempts >= 100):
                    backrubTm = backrubTm * 2.0
                    myLoop = myLoop + 10
                    Metropolis = rosetta.MonteCarlo(pose, scorefxn, backrubTm)
                    attempts = 0
                    lastRMSD = -1
                # If the temperature has gotten too hot or too cold, give up and return the current structure
                if (backrubTm >= 1000000.0):
                    print str(teamrank) + " - Too hot", tol
                    print teamrank, core.scoring.bb_rmsd(pose, orig_pose)
                    break
                if ( backrubTm < 0.0001):
                    print str(teamrank) + " - Too cold", tol
                    pose = Pose(prev_pose)
                    print teamrank, core.scoring.bb_rmsd(pose, orig_pose)
                    break
                # Calculate the current RMSD to the base pose
                calpha_superimpose_pose(pose, orig_pose)
                if (lastRMSD < 0):
                    lastRMSD = core.scoring.bb_rmsd(pose, orig_pose)
                else:
                    currRMSD = core.scoring.bb_rmsd(pose, orig_pose)
                    # If the RMSD did not really change from the last iteration, increase the temperature
                    # so we can get more diversity faster
                    if (currRMSD - lastRMSD == 0 or myRMSD / (currRMSD - lastRMSD) > 10):
                        backrubTm = backrubTm * 2
                        myLoop = myLoop + 10
                        Metropolis = rosetta.MonteCarlo(pose, scorefxn, backrubTm)
                        attempts = 0
                        lastRMSD = -1
                    else:
                        lastRMSD = currRMSD
                if (core.scoring.bb_rmsd(pose, orig_pose) >= myRMSD):
                    # If the actual RMSD is more than halfway towards either of the adjacent templates, then start over and decrease the Tm by a factor of 10
                    if (math.fabs(core.scoring.bb_rmsd(pose, orig_pose) - myRMSD) < myRange):
                        break
                    else:
                        backrubTm = backrubTm / 10.0
                        myLoop = int(float(myLoop) / 2.0)
                        # Revert
                        pose = Pose(prev_pose)
                        Metropolis = rosetta.MonteCarlo(pose, scorefxn, backrubTm)
                        lastRMSD = -1
                else:
                    if (math.fabs(core.scoring.bb_rmsd(pose, orig_pose) - myRMSD) < myRange):
                        break
                    prev_pose = Pose(pose)
        if (teamrank > 0):
            # Now do two Cartesian minimizations to fix the backrub distortions
            # Extensive backrub moves distort ideal bond geometry
            # First, minimize all odd residues, then minimize all even residues, to preserve the RMSD
            # as much as possible
            mm = MoveMap()
            mm.set_bb(True)
            mm.set_chi(True)
            for i in range(1, pose.n_residue() + 1, 2):
                if ("ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".find(pose.residue(i).name3()) < 0):
                    mm.set_bb(i, False)
                    mm.set_chi(i, False)
            minmover = MinMover(mm, create_score_function("talaris2013_cart"), 'lbfgs_armijo_nonmonotone', 0.01, True)
            minmover.cartesian(True)
            minmover.apply(pose)
            mm.set_bb(True)
            mm.set_chi(True)
            for i in range(2, pose.n_residue() + 1, 2):
                if ("ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".find(pose.residue(i).name3()) < 0):
                    mm.set_bb(i, False)
                    mm.set_chi(i, False)
            minmover = MinMover(mm, create_score_function("talaris2013_cart"), 'lbfgs_armijo_nonmonotone', 0.01, True)
            minmover.cartesian(True)
            minmover.apply(pose)
        if (teamrank == 0):
            print "backrub: Backbone 1 Constructed - RMSD(A): 0.0"
        else:
            calpha_superimpose_pose(pose, orig_pose)
            currRMSD = core.scoring.bb_rmsd(pose, orig_pose)
            print "backrub: Backbone " + str(teamrank+1) + " Constructed - RMSD(A): " + str(currRMSD)
        filename = pdbfile.split(".pdb")[0] + "_" + str(teamrank+1) + ".pdb"
        pose.dump_pdb("temp.pdb")
        # Write the coordinates to the archive
        fin = open("temp.pdb", "r")
        archive.write("BEGIN PDB " + filename.strip() + "\n")
        for aline in fin:
            archive.write(aline)
        archive.write("END PDB " + filename.strip() + "\n")
        fin.close()
    archive.close()
    # So the main GUI doesn't attempt to read the file before the daemon finishes writing its contents
    os.rename("resultstemp.ensb", "results.ensb")

def doScore():
    initializeRosetta()
    try:
        f = open("scoreinput", "r")
    except:
        raise Exception("ERROR: The file \"scoreinput\" is missing!")
    # Get the pdbfile, resfile, and scorefxn from the input file
    for aline in f:
        if (aline[0:7] == "PDBFILE"):
            pdbfile = aline.split("\t")[1].strip()
        elif (aline[0:8] == "SCOREFXN"):
            weightsfile = aline.split("\t")[1].strip()
    f.close()
    # Initialize scoring function
    scorefxn = ScoreFunction()
    try:
        scorefxn.add_weights_from_file(weightsfile)
    except:
        raise Exception("ERROR: The scoring function weights could not be initialized!")
    f = open("scoreoutputtemp", "w")
    pose = pose_from_pdb(pdbfile)
    # Calculate energy
    total_E = scorefxn(pose)
    # Dump the output
    outputpdb = pdbfile.split(".pdb")[0] + "_S.pdb"
    # Save atomic charge information in the segment identifier fields (columns 73-76 in PDB)
    # PyMOL reads that as the "segl" field
    info = pose.pdb_info()
    charges = {}
    for ires in range(1, pose.n_residue()+1):
        for iatom in range(1, len(pose.residue(ires).atoms())+1):
            pc = pose.residue(ires).atomic_charge(iatom)
            pc = round(pc * 4.0)
            if (pc < -4):
                pc = -4
            elif (pc > 4):
                pc = 4
            charges[str(info.number(ires)) + pose.residue(ires).atom_name(iatom)] = int(pc)
    pose.dump_pdb(outputpdb)
    # Now manually insert the charges
    fin = open(outputpdb, "r")
    data = []
    for aline in fin:
        data.append(aline)
    fin.close()
    fout = open(outputpdb, "w")
    for aline in data:
        if (aline.startswith("ATOM") or aline.startswith("HETATM")):
            pc = charges[aline[22:26].strip() + aline[12:16]]
            if (pc < 0):
                aline = aline[0:78] + "%1i-" % (pc * -1) + aline[80:]
            else:
                aline = aline[0:78] + "%1i+" % pc + aline[80:]
        fout.write(aline)
    fout.close()
    # Now write the output information for the main GUI
    f.write("OUTPUT\t" + outputpdb + "\n")
    f.write("TOTAL_E\t" + str(total_E) + "\n")
    nonzero_scoretypes = scorefxn.get_nonzero_weighted_scoretypes()
    f.write("ENERGY\ttotal_score")
    for scoretype in nonzero_scoretypes:
        f.write("\t" + str(scoretype))
    f.write("\n")
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
    # So the main GUI doesn't attempt to read the file before the daemon finishes writing its contents
    os.rename("scoreoutputtemp", "scoreoutput")

def doRotamerSearch():
    initializeRosetta()
    try:
        f = open("rotamerinput", "r")
    except:
        raise Exception("ERROR: The file \"rotamerinput\" is missing!")
    # Get the pdbfile, resfile, and scorefxn from the input file
    for aline in f:
        if (aline[0:7] == "PDBFILE"):
            pdbfile = aline.split("\t")[1].strip()
        elif (aline[0:8] == "SCOREFXN"):
            weightsfile = aline.split("\t")[1].strip()
        elif (aline[0:7] == "RESTYPE"):
            restype = aline.split("\t")[1].strip()
            resone = AA3to1(restype)
        elif (aline[0:6] == "SEQPOS"):
            data = aline.split("\t")[1].strip()
            chain = data[0]
            seqpos = int(data[3:])
    f.close()
    # Initialize scoring function
    scorefxn = ScoreFunction()
    try:
        scorefxn.add_weights_from_file(weightsfile)
    except:
        raise Exception("ERROR: The scoring function weights could not be initialized!")
    f = open("rotameroutputtemp", "w")
    try:
        rsd_factory = pose_from_pdb("data/residues.pdb")
    except:
        raise Exception("ERROR: The file \"data/residues.pdb\" is missing!")
    pose = pose_from_pdb(pdbfile)
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
    fout = open("progress", "w")
    fout.write("0/10\n")
    fout.write("Reading rotamer libraries...\n")
    fout.close()
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
                if (restype == "PRO"):
                    chis[2] = "-999.0"
                    chis[3] = "-999.0"
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
        #mutate_residue(pose, indx, resone)
    except:
        raise Exception("ERROR: The new residue could not be mutated onto the structure.")
    # Search all the chi values and score them
    k = 0
    nonzero_scoretypes = scorefxn.get_nonzero_weighted_scoretypes()
    atompos = []
    fout = open("progress", "w")
    fout.write("0/" + str(len(chivals)) + "\n")
    fout.write("Scoring rotamers...\n")
    fout.close()
    count = 0
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
        atomstr = ""
        for i in range(1, pose.residue(indx).natoms()+1):
            atomstr = atomstr + pose.residue(indx).atom_name(i).strip() + " " + str(pose.residue(indx).atom(i).xyz()[0]) + " " + str(pose.residue(indx).atom(i).xyz()[1]) + " " + str(pose.residue(indx).atom(i).xyz()[2]) + "\t"
        atomstr = atomstr.strip()
        atompos.append(atomstr)
        k = k + 1
        # The following line must be here
        # Sometimes after doing multiple set_chis the structure gets really messed up (don't know if I am
        # doing something wrong or if it is a Rosetta bug, it only happens in special occasions) so we have
        # to revert back to the original structure so errors are not propagated
        pose.replace_residue(indx, res_mutate.clone(), True)
        count += 1
        if (count % 10 == 0):
            fout = open("progress", "w")
            fout.write(str(count) + "/" + str(len(chivals)) + "\n")
            fout.write("Scoring rotamers...\n")
            fout.close()
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
        temp = atompos[i]
        atompos[i] = atompos[lowest]
        atompos[lowest] = temp
    # Dump the output
    outputpdb = pdbfile.split(".pdb")[0] + "_R.pdb"
    pose.dump_pdb(outputpdb)
    # Now write the output information for the main GUI
    f.write("OUTPUT\t" + outputpdb + "\n")
    f.write("INDEX\t" + str(indx) + "\t")
    info = pose.pdb_info()
    if (len(info.chain(indx).strip()) == 0):
        f.write("_" + str(info.number(indx)) + "\n")
    else:
        f.write(info.chain(indx) + str(info.number(indx)) + "\n")
    if (restype == "PRO"):
        nchi = 2
    else:
        nchi = pose.residue(indx).nchi()
    #for ichi in range(1, nchi+1):
        #chiatoms = pose.residue(indx).chi_atoms()[ichi]
        #f.write("CHIATOMS\t" + pose.residue(indx).atom_name(chiatoms[1]).strip() + " " + pose.residue(indx).atom_name(chiatoms[2]).strip() + " " + pose.residue(indx).atom_name(chiatoms[3]).strip() + " " + pose.residue(indx).atom_name(chiatoms[4]).strip() + "\n")
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
        #f.write("CHI\t" + chivals[i][0] + "\t" + chivals[i][1] + "\t" + chivals[i][2] + "\t" + chivals[i][3] + "\n")
        f.write("ROTAMER\t" + atompos[i] + "\n")
    f.close()
    # So the main GUI doesn't attempt to read the file before the daemon finishes writing its contents
    os.rename("rotameroutputtemp", "rotameroutput")

def removeNCAAs(pose):
    # Rosetta throws an error if you try to convert a pose that has HETATMS to centroid mode, so we're doing
    # to have to delete all of those HETATMs first
    for ires in range(pose.n_residue(), 0, -1):
        if (not(pose.residue(ires).name3() in "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR")):
            if (ires > 1):
                pose.delete_polymer_residue(ires)
            else:
                pose.replace_residue(1, pose.residue(2), False)
                pose.delete_polymer_residue(2)
    return pose

def doKIC(stage="Coarse"):
    #try:
    initializeRosetta()
    if (stage == "Coarse"):
        try:
            f = open("coarsekicinput", "r")
        except:
            raise Exception("ERROR: The file \"coarsekicinput\" is missing!")
    else:
        try:
            f = open("finekicinput", "r") # Same file as before
        except:
            raise Exception("ERROR: The file \"finekicinput\" is missing!")
    # Get the pdbfile, resfile, and scorefxn from the input file
    pivotOffset = 0
    loopdata = []
    for aline in f:
        if (aline[0:7] == "PDBFILE"):
            # But for the fine grained step the pose comes from repacked.pdb
            pdbfile = aline.split("\t")[1].strip()
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
        elif (aline[0:5] == "PIVOT"):
            pivotOffset = int(aline.split("\t")[1])
        elif (aline[0:4] == "LOOP"):
            # Save the information in a list that we will iterate through later
            loopdata.append(aline.strip().split("\t")[1:])
        elif (aline[0:7] == "NSTRUCT"):
            nstruct = int(aline.split("\t")[1])
        elif (aline[0:7] == "PERTURB"):
            perturbType = aline.split("\t")[1].strip()
        elif (aline[0:9] == "OUTPUTDIR"):
            outputdir = aline.split("\t")[1].strip()
    f.close()
    # Initialize scoring function
    scorefxn = ScoreFunction()
    try:
        scorefxn.add_weights_from_file(weightsfile)
    except:
        raise Exception("ERROR: The scoring function weights could not be initialized!")
    if (stage == "Coarse"):
        pose = pose_from_pdb(pdbfile)
        # We have to take the HETATMs out otherwise it will crash the centroid mover
        # We could try to make centroid params files, but it seems that they sometimes cause Rosetta to seg fault, which we cannot have
        # Later we'll put them back in
        HETATMs = []
        HETATM_indx = []
        info = pose.pdb_info()
        for i in range(pose.n_residue(), 0, -1):
            if (not(pose.residue(i).name3() in "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR ")):
                if (pose.residue(i).is_polymer()):
                    HETATMs.append([i, Residue(pose.residue(i))])
                else:
                    HETATM_indx.append((info.chain(i), info.number(i)))
                # Taking these residues out might change where the loop anchors are
                for j in range(0, len(loopdata)):
                    if (i < int(loopdata[j][2])):
                        loopdata[j][2] = int(loopdata[j][2]) -1
                    if (i < int(loopdata[j][3])):
                        loopdata[j][3] = int(loopdata[j][3]) -1
                    if (i < int(loopdata[j][4])):
                        loopdata[j][4] = int(loopdata[j][4]) -1
                if (i == 1):
                    pose.replace_residue(1, pose.residue(2), False)
                    pose.delete_polymer_residue(2)
                else:
                    pose.delete_polymer_residue(i)
        HETATM_lines = []
        f = open(pdbfile, "r")
        for aline in f:
            if (aline[0:4] == "ATOM" or aline[0:6] == "HETATM"):
                chain = aline[21].strip()
                seqpos = int(aline[22:26])
                if ((chain, seqpos) in HETATM_indx):
                    HETATM_lines.append(aline.strip())
        f.close()
    loops = []
    for [loopType, sequence, begin, pivot, end] in loopdata:
        loopBegin = int(begin)
        pivot = int(pivot)
        loopEnd = int(end)
        if (loopType == "DE NOVO"):
            if (stage == "Coarse"):
                # Since this is a new sequence being added, we first have to delete all the residues
                # between the beginning and ending points
                oldlen = 0
                for ires in range(loopEnd-1, loopBegin, -1):
                    pose.delete_polymer_residue(ires)
                    oldlen += 1
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
                # Now we have to update the begin and end points of the other loops if necessary
                for j in range(0, len(loopdata)):
                    if (loopBegin < int(loopdata[j][2])):
                        loopdata[j][2] = int(loopdata[j][2]) + len(sequence) - oldlen
                    if (loopBegin < int(loopdata[j][3])):
                        loopdata[j][3] = int(loopdata[j][3]) + len(sequence) - oldlen
                    if (loopBegin < int(loopdata[j][4])):
                        loopdata[j][4] = int(loopdata[j][4]) + len(sequence) - oldlen
            # Now maybe the sequence is longer than what was originally the length of the sequence
            # between start and end, so we need to recalculate the loop end
            loopEnd = loopBegin + len(sequence.strip()) + 1
        #elif (loopType == "REFINE"):
            # Here we only want to do high resolution modeling
        #    stage = "Fine"
        #    nstruct = 1
        if (loopType == "DE NOVO" and stage == "Coarse"):
            # This has to be hard-coded, because the loop is not actually there until coarse modeling happens so there's no pivot point
            # other than the loop anchor residues
            cutpoint = loopEnd
        else:
            cutpoint = pivot
        if (loopType == "REFINE"):
            loop = Loop(loopBegin, loopEnd, cutpoint, 0, 0)
        else:
            loop = Loop(loopBegin, loopEnd, cutpoint, 0, 1)
        loops.append(loop)
    if (stage == "Coarse"):
        # Low res KIC
        for decoy in range(0, nstruct):
            sw = SwitchResidueTypeSetMover("centroid")
            try:
                sw.apply(pose)
            except:
                raise Exception("ERROR: The PDB could not be converted to centroid mode!")
            for loop in loops:
                loops_set = Loops()
                loops_set.add_loop(loop)
                add_single_cutpoint_variant(pose, loop)
                set_single_loop_fold_tree(pose, loop)
                kic_perturb = LoopMover_Perturb_KIC(loops_set)
                kic_perturb.set_max_kic_build_attempts(120)
                try:
                    kic_perturb.apply(pose)
                except:
                    raise Exception("ERROR: The coarse KIC perturber failed!")
            if (perturbType == "Perturb Only, Centroid"):
                outputpdb = pdbfile.split(".pdb")[0] + "_K.pdb"
                pose.dump_pdb(outputpdb)
                # Now score the pose so we have the energy information in the main GUI
                scorefxn = create_score_function("cen_std")
                scorefxn(pose)
                f = open("kicoutputtemp", "w")
                f.write("OUTPUT\t" + outputpdb + "\n")
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
                # So the main GUI doesn't attempt to read the file before the daemon finishes writing its contents
                os.rename("kicoutputtemp", "kicoutput")
            else:
                sw = SwitchResidueTypeSetMover("fa_standard")
                try:
                    sw.apply(pose)
                except:
                    raise Exception("ERROR: The PDB could not be converted back to fullatom mode from centroid mode!")
                # Now put the HETATMs back on
                if (len(HETATMs) > 0):
                    ires = 1
                    nres = pose.n_residue()
                    [addres, res] = HETATMs.pop()
                    while (ires <= nres):
                        if (ires == addres):
                            if (res.is_polymer()):
                                pose.prepend_polymer_residue_before_seqpos(res, addres, False)
                            else:
                                pose.append_residue_by_jump(res, addres-1)
                            nres = nres + 1
                            if (len(HETATMs) > 0):
                                [addres, res] = HETATMs.pop()
                            else:
                                break
                        ires = ires + 1
                # Dump it for the repacking daemon
                pose.dump_pdb("torepack_" + str(decoy) + ".pdb")
                # Now append the non-polymer lines
                data = []
                f = open("torepack_" + str(decoy) + ".pdb", "r")
                for aline in f:
                    if (aline[0:4] == "ATOM" or aline[0:6] == "HETATM"):
                        data.append(aline.strip())
                f.close()
                f = open("torepack_" + str(decoy) + ".pdb", "w")
                for aline in data:
                    f.write(aline + "\n")
                for aline in HETATM_lines:
                    f.write(aline + "\n")
                f.close()
        # Now we have to rewrite the input file for refined mode, since if we added or deleted
        # sequence the loop beginning and ending points may have changed
        fin = open("coarsekicinput", "r")
        filedata = []
        for aline in fin:
            if (aline[0:5] == "LOOP\t"):
                continue
            filedata.append(aline)
        fin.close()
        fout = open("coarsekicinput", "w")
        for aline in filedata:
            fout.write(aline)
        for [loopType, sequence, begin, pivot, end] in loopdata:
            fout.write("LOOP\t" + loopType + "\t" + sequence.strip() + "\t" + str(begin) + "\t" + str(pivot) + "\t" + str(end) + "\n")
        fout.close()
    else:
        for decoy in range(0, nstruct):
            if (loopType == "REFINE"):
                pose = pose_from_pdb(pdbfile)
            else:
                pose = pose_from_pdb("repacked_" + str(decoy) + ".pdb")
            for loop in loops:
                loops_set = Loops()
                loops_set.add_loop(loop)
                add_single_cutpoint_variant(pose, loop)
                set_single_loop_fold_tree(pose, loop)
                kic_refine = LoopMover_Refine_KIC(loops_set)
                try:
                    kic_refine.apply(pose)
                except:
                    raise Exception("ERROR: The KIC refiner failed!")
            if (decoy == 0):
                outputpdb = pdbfile.split(".pdb")[0] + "_K.pdb"
                pose.dump_pdb(outputpdb)
                # Now score the pose so we have the energy information in the main GUI
                scorefxn(pose)
                f = open("kicoutputtemp", "w")
                f.write("OUTPUT\t" + outputpdb + "\n")
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
            # Dump all of the files to the indicated directory
            if (nstruct > 1):
                pose.dump_pdb(outputdir + "/" + pdbfile.split(".pdb")[0] + ("_KIC_%4.4i.pdb" % (decoy+1)))
        # So the main GUI doesn't attempt to read the file before the daemon finishes writing its contents
        os.rename("kicoutputtemp", "kicoutput")


def doINDEL(scriptdir):
    # Query database - lookup executable and database files need to be in sandbox
    num_results = 0

    # Find the right binary
    if (platform.system() == "Windows"):
        INDELprogram = scriptdir + "\\bin\\iRosetta_Lookup_Win64.exe"
    elif (platform.system() == "Darwin"):
        INDELprogram = scriptdir + "/bin/iRosetta_Lookup_osx.exe"
    else:
        INDELprogram = scriptdir + "/bin/iRosetta_Lookup_Linux.exe"
    logInfo(INDELprogram)
    print INDELprogram

    runINDEL = INDELprogram

    # Grab the parameters PyRosetta needs
    f = open("INDELinput", "r")
    for line in f:
        if (line[0:8] == "SCAFFOLD"):
            scaffold_pdb = line.split('\t')[1]
        if (line[0:7] == "ANCHORS"):
            anchors = line.split('\t')
            anchors.pop(0)
        if (line[0:8] == "SYMMETRY"):
            symmetry = int(line.split('\t')[1])
        if (line[0:11] == "MAX_RESULTS"):
            max_results = int(line.split('\t')[1])
        if (line[0:8] == "SCOREFXN"):
            selected_scorefxn = line.split('\t')[1]
    f.close()


    # Find where the database files are, append to the input file that iRosetta Lookup uses
    if platform.system() == "Windows":
        dataDir = scriptdir + "\\data\\"
    else:
        dataDir = scriptdir + "/data/"
        print "dataDir:",dataDir

    f = open("INDELinput", "a")

    f.write("PDB_DATA\t"  + "%spdblist.dat"%(dataDir)  + "\n")
    print "PDB_DATA"
    f.write("LOOP_DATA\t" + "%slooplist.dat"%(dataDir) + "\n")
    print "LOOP_DATA"
    f.write("GRID_DATA\t" + "%sgrid.dat"%(dataDir)     + "\n")
    print "GRID_DATA"
    f.write("LOG\t" + "indel_log.txt" + "\n")
    f.flush()
    f.close()
    import commands; print commands.getstatusoutput("cp -v INDELinput INDELinput.log")


    # Run the lookup
    try:
        num_results = int( subprocess.check_output([runINDEL, "INDELinput"]) )
    except:
       writeError("ERROR: The database query failed! Check input parameters and try again")
       return

    print "Finished lookup"


    # If we got results, try to build models with them
    if (num_results > 0):
        try:
            # Get every component of Rosetta started that we can, try to avoid
            # initializing for every loop we try to insert
            print "Found " + str(num_results) + " loops \n"

            initializeRosetta()
            import rosetta.protocols.grafting as graft

            # Set some parameters
            # Just use default score function for now
            #try:
            #    scorefxn = create_score_function(selected_scorefxn)
            #    print "Set score function to " + selected_scorefxn
            #except:
            #    scorefxn = create_score_function("talaris2013")
            #    print "Set score function to default"
            #Shouldn't this use the score function that the user selected???
            scorefxn = create_score_function("talaris2013")

        except:
            writeError("ERROR: Couldn't initialize Rosetta!")
            return

        # Figure out where we're grafting
        start_residue = int(anchors[0])
        stop_residue = int(anchors[1])

        i = 1
        models  = []
        scores  = []
        lengths = []
        goToSandbox()



        # Try to make models from all the loops, or up to the max number specified
        while i < min((num_results + 1), max_results+1):

            if i < 97:
                f = open("progress", 'w')
                f.write("%d/100"%(i))
                f.close()


            print "\n ==================================================== \n"
            print "\t \t \t Attempting to insert loop " + str(i)
            print "\n ==================================================== \n"

            try:
                # Make a copy of the scaffold, load the loop
                loopfile = "loopout_" + str(i) + ".pdb"
                loop = pose_from_pdb(loopfile)
                temp_pose = pose_from_pdb(scaffold_pdb.strip())
                chain_length = temp_pose.total_residue() / symmetry
                graft_dist = stop_residue - start_residue + 1

                # Here's where the actual grafting happens. Graft and add to model list
                # Don't need to repack, AnchoredGraftMover takes care of it
                for j in range(0, symmetry):
                    startgraft = start_residue + (j * chain_length) + (j * (loop.total_residue() - graft_dist ))
                    endgraft   = stop_residue +  (j * chain_length) + (j * (loop.total_residue() - graft_dist ))
                    graftmover = graft.AnchoredGraftMover(startgraft, endgraft , loop, 1, 1)
                    graftmover.set_cycles(50)
                    #graftmover.set_piece(loop, 1, 1)
                    graftmover.apply(temp_pose)
                    graftmover.final_repack()
                    print "Grafted at %d and %d"%(startgraft, endgraft)

                #graftmover = graft.AnchoredGraftMover(start_residue + chain_length, stop_residue + chain_length, loop, 1, 1)
                #graftmover.set_cycles(500)
                #graftmover.set_piece(loop, 1, 1)
                #graftmover.apply(temp_pose)



                #Add the model and its score to their respective lists, as well as the length of each insertion
                temp_score = scorefxn(temp_pose)
                scores.append(temp_score)
                models.append(temp_pose)
                lengths.append(loop.total_residue())
            except:
                writeError("Something went wrong (could be Pyrosetta).")
                return
            i += 1

        # Sort models by energy and write them out so we can look at them
        # Write out a file with info about each loop
        #Console output
        f = open("progress", 'w')
        f.write("100/100")
        f.close()

        i = 1
        f = open("INDELoutputtemp", "w")
        for score, model, length in sorted(zip(scores, models, lengths)):
            outname = "INDELmodel_" + str(i) + ".pdb"
            model.dump_pdb(outname)
            f.write(outname + "\t" + str(score) + "\t" + str(length) + "\n")
            i += 1
        f.close()
        os.rename("INDELoutputtemp", "INDELoutput")

        print "Chain length: %d"%(chain_length)



    else:
        writeError("No results found. Try relaxing search parameters")
        return





    print "Tried to find loops for anchor residues N-term:" + anchors[0] + "\t C-term: " + anchors[1]
    print str(num_results) + " loops found."
    print str(len(models)) + " models built."








def doRepack(scorefxninput="", pdbfile="repackme.pdb", lastStage=False):
    initializeRosetta()
    try:
        f = open(scorefxninput, "r")
        # Get the pdbfile, resfile, and scorefxn from the input file
        for aline in f:
            if (aline[0:8] == "SCOREFXN"):
                weightsfile = aline.split("\t")[1].strip()
            elif (aline[0:9] == "LOOPBEGIN"):
                loopBegin = int(aline.split("\t")[1])
            elif (aline[0:7] == "LOOPEND"):
                loopEnd = int(aline.split("\t")[1])
            elif (aline[0:8] == "SEQUENCE"):
                sequence = aline.split("\t")[1].strip()
        f.close()
        # Initialize scoring function
        scorefxn = ScoreFunction()
        scorefxn.add_weights_from_file(weightsfile)
    except:
        # Default to Talaris2013
        scorefxn = create_score_function("talaris2013")
    # Repack
    pose = pose_from_pdb(pdbfile)
    packtask = standard_packer_task(pose)
    packtask.restrict_to_repacking()
    packmover = PackRotamersMover(scorefxn, packtask)
    try:
        packmover.apply(pose)
    except:
        raise Exception("ERROR: The Rosetta packer failed!")
    os.remove(pdbfile)
    if (lastStage):
        loopEnd = loopBegin + len(sequence.strip()) + 1
        outputpdb = pdbfile.split(".pdb")[0] + "_K.pdb"
        pose.dump_pdb(outputpdb)
        # Now score the pose so we have the energy information in the main GUI
        scorefxn(pose)
        f = open("kicoutputtemp", "w")
        f.write("OUTPUT\t" + outputpdb + "\n")
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
        # So the main GUI doesn't attempt to read the file before the daemon finishes writing its contents
        os.rename("kicoutputtemp", "kicoutput")
    else:
        try:
            outputpdbfile = "repacked" + pdbfile.split(".pdb")[0].split("repackme")[1] + "temp.pdb"
        except:
            outputpdbfile = "repackedtemp.pdb"
        pose.dump_pdb(outputpdbfile)

def dockSuperimpose(pose, pose2, root_res1, root_res2):
    # Function for superimposing two poses for docking purposes
    # It is needed because in the low-resolution stage the polymer HETATMs get dropped
    # If we superimpose the pre-docked fullatom structure with the post coarse-docked fullatom structure
    # we can simply take the polymer resiudes out of the original structure and put them into the docked model
    # Why don't I just create centroid params files for these polymers?
    # Well I tried that for the GFP CRO and it was seg faulting because Rosetta couldn't find dof for certain dihedrals
    # and I couldn't figure out how to fix it, then I was worried this problem might show up more often and don't want
    # the program to be seg faulting for any reason.  This method avoids the possibility of seg faulting
    # How does this work?
    # We only need three equivalent atoms from the two poses to make it happen (root_res specifies the residue and N, CA, and C are the three atoms)
    # First, translate all the atoms such that the first atom is at the origin
    translateToOrigin(pose, root_res1)
    translateToOrigin(pose2, root_res2)
    # Now we calculate a rotation matrix that will rotate atom 2 onto the x-axis, and then rotate all atoms in the poses by this matrix
    unit_x = numpy.array([1.0, 0.0, 0.0]) # X unit vector
    p1 = numpy.array([pose.residue(root_res1).atom(2).xyz()[0], pose.residue(root_res1).atom(2).xyz()[1], pose.residue(root_res1).atom(2).xyz()[2]])
    p2 = numpy.array([pose2.residue(root_res2).atom(2).xyz()[0], pose2.residue(root_res2).atom(2).xyz()[1], pose2.residue(root_res2).atom(2).xyz()[2]])
    p1_v = getUnitVector(p1)
    p2_v = getUnitVector(p2)
    R1 = getRotationMatrix(p1_v, unit_x)
    R2 = getRotationMatrix(p2_v, unit_x)
    if (R1 != "No Rotation"):
        rotatePose(pose, R1)
    if (R2 != "No Rotation"):
        rotatePose(pose2, R2)
    # Now we have to rotate the 3rd atom around the x-axis such that it is in the xy-plane and then we should have superimposed the proteins
    unit_y = numpy.array([0.0, 1.0, 0.0])
    p1 = numpy.array([0.0, pose.residue(root_res1).atom(3).xyz()[1], pose.residue(root_res1).atom(3).xyz()[2]])
    p2 = numpy.array([0.0, pose2.residue(root_res2).atom(3).xyz()[1], pose2.residue(root_res2).atom(3).xyz()[2]])
    p1_v = getUnitVector(p1)
    p2_v = getUnitVector(p2)
    R1 = getRotationMatrix(p1_v, unit_y)
    R2 = getRotationMatrix(p2_v, unit_y)
    if (R1 != "No Rotation"):
        rotatePose(pose, R1)
    if (R2 != "No Rotation"):
        rotatePose(pose2, R2)

def doDock(stage="Coarse"):
    #try:
    if (os.path.isfile("constraints.cst")):
        useConstraints = True
    else:
        useConstraints = False
    if (useConstraints):
        initializeRosetta("-constraints:cst_file constraints.cst")
    else:
        initializeRosetta()
    if (stage == "Coarse"):
        try:
            f = open("coarsedockinput", "r")
        except:
            raise Exception("ERROR: The file \"coarsedockinput\" is missing!")
    else:
        try:
            f = open("finedockinput", "r") # Same file as before
        except:
            raise Exception("ERROR: The file \"finedockinput\" is missing!")
    # Get the pdbfile, resfile, and scorefxn from the input file
    for aline in f:
        if (aline[0:7] == "PDBFILE"):
            # But for the fine grained step the pose comes from repacked.pdb
            if (stage == "Coarse"):
                pdbfile = aline.split("\t")[1].strip()
            else:
                pdbfile = "repackeddock.pdb"
        elif (aline[0:8] == "SCOREFXN"):
            weightsfile = aline.split("\t")[1].strip()
        elif (aline[0:10] == "JUMPCONFIG"):
            jumpconfig = aline.split("\t")[1].strip()
        elif (aline[0:6] == "ORIENT"):
            orient = aline.split("\t")[1].strip()
        elif (aline[0:12] == "COARSEDECOYS"):
            ncoarse = int(aline.split("\t")[1].strip())
        elif (aline[0:13] == "REFINEDDECOYS"):
            nrefined = int(aline.split("\t")[1].strip())
    f.close()
    # Initialize scoring function
    scorefxn = ScoreFunction()
    try:
        scorefxn.add_weights_from_file(weightsfile)
    except:
        raise Exception("ERROR: The scoring function weights could not be initialized!")
    if (useConstraints):
        scorefxn.set_weight(atom_pair_constraint, 1)
    if (stage == "Coarse"):
        pose = pose_from_pdb(pdbfile)
        # Low res docking
        best_decoys = []
        for i in range(0, nrefined):
            best_decoys.append([None, 9999999999.0])
        # Now let's see if there are polymer HETATMs that will be dropped when converting to centroid mode
        # (I am currently not concerned with losing ligands since I think that docking ligands is a separate Rosetta protocol)
        staticpolymers = False
        staticHETATMs = []
        movingpolymers = False
        movingHETATMs = []
        staticchains = jumpconfig.split("_")[0]
        movingchains = jumpconfig.split("_")[1]
        info = pose.pdb_info()
        PDBnumbering = []
        for ires in range(1, pose.n_residue()+1):
            if (not(pose.residue(ires).name3() in "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR") and pose.residue(ires).is_polymer()):
                if (info.chain(ires) in staticchains):
                    staticpolymers = True
                    staticHETATMs.append(ires)
                else:
                    movingpolymers = True
                    movingHETATMs.append(ires)
                PDBnumbering.append(info.number(ires))
            elif (pose.residue(ires).name3() in "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR"):
                PDBnumbering.append(info.number(ires))
        if (staticpolymers):
            posevec = pose.split_by_chain()
            staticpose = Pose(posevec[1])
            for i in range(2, len(staticchains)+1):
                staticpose.append_pose_by_jump(posevec[i], staticpose.n_residue())
            static_root_res = []
            offset = 0
            for i in range(1, staticpose.n_residue()):
                if (not(staticpose.residue(i).name3() in "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR")):
                    offset = offset + 1
                else:
                    static_root_res = [i, i-offset]
                    break
        if (movingpolymers):
            posevec = pose.split_by_chain()
            movingpose = Pose(posevec[jumpconfig.index("_")+1])
            for i in range(jumpconfig.index("_")+2, len(movingchains)+1):
                movingpose.append_pose_by_jump(posevec[i], movingpose.n_residue())
            moving_root_res = []
            offset = 0
            for i in range(1, movingpose.n_residue()):
                if (not(movingpose.residue(i).name3() in "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR")):
                    offset = offset + 1
                else:
                    moving_root_res = [i, i-offset]
                    break
        pose = removeNCAAs(pose)
        orig_pose = Pose(pose)
        scorefxn_low = create_score_function("interchain_cen") # This seems to be necessary for lowres docking, so it's hard-coded
        if (useConstraints):
            scorefxn_low.set_weight(atom_pair_constraint, 1) # To evaluate constraints (NOTE: atom_pair_constraint applies to site_constraints also)
            add_constraints_from_cmdline_to_pose(pose)
        sw = SwitchResidueTypeSetMover("centroid")
        for decoy in range(0, ncoarse):
            pose = Pose(orig_pose)
            try:
                sw.apply(pose)
            except:
                raise Exception("ERROR: The PDB could not be converted to centroid mode!")
            setup_foldtree(pose, jumpconfig, Vector1([1]))
            randomize1 = rigid_moves.RigidBodyRandomizeMover(pose, 1, rigid_moves.partner_upstream)
            randomize2 = rigid_moves.RigidBodyRandomizeMover(pose, 1, rigid_moves.partner_downstream)
            slide = DockingSlideIntoContact(1)
            try:
                # Randomize orientations as per the user's instructions
                if (orient == "Global" or orient == "Fix Mov"):
                    randomize1.apply(pose)
                if (orient == "Global" or orient == "Fix Stat"):
                    randomize2.apply(pose)
                slide.apply(pose)
            except:
                raise Exception("ERROR: The receptor and ligand chains were not able to be brought into contact!")
            dock_lowres = DockingLowRes(scorefxn_low, 1)
            try:
                dock_lowres.apply(pose)
            except:
                raise Exception("ERROR: The coarse docker failed!")
            # Is it better?
            E = scorefxn_low(pose)
            for i in range(0, nrefined):
                if (E < best_decoys[i][1]):
                    for j in range(nrefined-1, i, -1):
                        best_decoys[j] = best_decoys[j-1]
                    best_decoys[i] = [Pose(pose), E]
            # Drop information for the progress bar
            f = open("dock_progress", "w")
            f.write(str(decoy+1) + "\n")
            f.close()
        # Convert back to fa for all the best decoys
        sw = SwitchResidueTypeSetMover("fa_standard")
        for decoy in range(0, min(ncoarse, nrefined)):
            try:
                sw.apply(best_decoys[decoy][0])
            except:
                raise Exception("ERROR: The PDB could not be converted back to fullatom mode from centroid mode!")
            # Cast a magic spell to summon the HETATMs
            if (staticpolymers):
                dockSuperimpose(best_decoys[decoy][0], staticpose, static_root_res[0], static_root_res[1])
                for hetres in staticHETATMs:
                    best_decoys[decoy][0].append_polymer_residue_after_seqpos(Residue(staticpose.residue(hetres)), hetres-1, False)
            if (movingpolymers):
                dockSuperimpose(best_decoys[decoy][0], movingpose, staticpose.n_residue() + moving_root_res[0], moving_root_res[1])
                for hetres in movingHETATMs:
                    best_decoys[decoy][0].append_polymer_residue_after_seqpos(Residue(movingpose.residue(hetres-staticpose.n_residue())), hetres-1, False)
            if (staticpolymers or movingpolymers):
                # Rosetta is annoying, adding those HETATM residues back in renumbered all the PDB numbers, which is going
                # to screw up high resolution docking if there are constraints because the residue numbers are all different
                # Change them back to what they were
                info = best_decoys[decoy][0].pdb_info()
                for ires in range(1, best_decoys[decoy][0].n_residue()+1):
                    info.number(ires, PDBnumbering[ires-1])
            # Dump it for the repacking daemon
            best_decoys[decoy][0].dump_pdb("torepackdock_" + str(decoy) + ".pdb")
            # For some reason Rosetta is refusing to update the PDBnumbering even though the data in pdb_info is correct
            # So I have to change the file manually
            data = []
            f = open("torepackdock_" + str(decoy) + ".pdb", "r")
            ires = -1
            lastres = "0000"
            for aline in f:
                if (aline[0:4] == "ATOM" or aline[0:6] == "HETATM"):
                    if (lastres != aline[22:26]):
                        lastres = aline[22:26]
                        ires = ires + 1
                    data.append(aline[0:22] + "%4i" % PDBnumbering[ires] + aline[26:])
                else:
                    data.append(aline)
            f.close()
            f = open("torepackdock_" + str(decoy) + ".pdb", "w")
            for aline in data:
                f.write(aline)
            f.close()
    else:
        f = open("dockoutputtemp", "w")
        for decoy in range(0, nrefined):
            pdbfile = "repackeddock_" + str(decoy) + ".pdb"
            pose = pose_from_pdb(str(pdbfile))
            scorefxn_dock = create_score_function_ws_patch("docking", "docking_min")
            setup_foldtree(pose, jumpconfig, Vector1([1]))
            if (useConstraints):
                add_constraints_from_cmdline_to_pose(pose)
                scorefxn_dock.set_weight(atom_pair_constraint, 1)
            dock_hires = DockMCMProtocol()
            dock_hires.set_scorefxn(scorefxn_dock)
            dock_hires.set_scorefxn_pack(scorefxn)
            dock_hires.set_partners(jumpconfig)
            try:
                dock_hires.apply(pose)
            except:
                raise Exception("ERROR: The docking refiner failed!")
            outputpdb = pdbfile.split(".pdb")[0] + "_D.pdb"
            pose.dump_pdb(outputpdb)
            # Now score the pose so we have the energy information in the main GUI
            scorefxn(pose)
            f.write("OUTPUT\t" + outputpdb + "\n")
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
        # So the main GUI doesn't attempt to read the file before the daemon finishes writing its contents
        os.rename("dockoutputtemp", "dockoutput")

def doThread(scriptdir):
    try:
        f = open("threadinput", "r")
    except:
        raise Exception("ERROR: The file \"threadinput\" is missing!")
    # Get the pdbfile, resfile, and scorefxn from the input file
    pdbfiles = []
    for aline in f:
        if (aline[0:7] == "PDBFILE"):
            pdbfile = aline.split("\t")[1].strip()
            pdbfiles.append(pdbfile)
        elif (aline[0:4] == "FRAG"):
            fragfile = aline.split("\t")[1].strip()
        elif (aline[0:8] == "SCOREFXN"):
            weightsfile = aline.split("\t")[1].strip()
    f.close()
    # Read the contents of the FRAG file (Note: you need to package up these files using the InteractiveROSETTA code
    # Packaging them up manually will not work
    f = gzip.open(fragfile, "rb")
    # These are "extra_options" for PyRosetta initialization
    frag_files = "-loops:frag_files "
    frag_sizes = "-loops:frag_sizes "
    fasta = "-in:file:fasta "
    psipred = "-in:file:psipred_ss2 "
    for aline in f:
        if (aline[0:5] == "BEGIN"):
            filename = aline.split("\t")[1].strip()
            f2 = open(aline.split("\t")[1].strip(), "w")
            if (filename.endswith(".200_v1_3")):
                frag_files = frag_files + filename + " "
                fragsize = str(int(filename.split("_")[1]))
                frag_sizes = frag_sizes + fragsize + " "
            elif (filename.endswith(".fasta")):
                fastafile = filename
                fasta = fasta + filename + " "
            elif (filename.endswith(".psipred_ss2")):
                psipred = psipred + filename + " "
        elif (aline[0:3] == "END"):
            f2.close()
        elif (len(aline.strip()) > 0):
            f2.write(aline.strip() + "\n")
    f.close()
    # Now we have to do a multiple sequence alignment between the sequence that will be modeled against the
    # sequences of the known templates
    predseq = ""
    f2 = open("muscle.fasta", "w")
    f = open(fastafile)
    # Easy to get the sequence to predict: just read it from the FASTA file
    f2.write(">Pred\n")
    for aline in f:
        if (not(">" in aline)):
            predseq = predseq + aline.strip()
            f2.write(aline.strip() + "\n")
    f.close()
    # Now read the sequences from the PDB files
    for template in pdbfiles:
        f2.write("\n>" + template + "\n")
        f = open(template, "r")
        curr_res = "0000"
        AA3 = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR "
        AA1 = "ACDEFGHIKLMNPQRSTVWY"
        for aline in f:
            if (aline[0:4] == "ATOM"):
                if (aline[22:26] != curr_res):
                    curr_res = aline[22:26]
                    indx = AA3.find(aline[17:20])
                    f2.write(AA1[indx/4])
        f.close()
        f2.write("\n")
    f2.close()
    # Use MUSCLE to do the alignments
    if (platform.system() == "Windows"):
        muscleprogram = scriptdir + "\\bin\\muscle.exe"
    else:
        muscleprogram = scriptdir + "/bin/muscle_unix"
    muscle_cline = MuscleCommandline(muscleprogram, input="muscle.fasta", out="aligned.fasta")
    muscle_cline()
    # Convert this aligned FASTA file to the "general" format used in Rosetta
    f = open("aligned.fasta", "r")
    f2 = open("align.aln", "w")
    f2.write("score 123.456")
    for aline in f:
        if (aline[0] == ">"):
            template = aline[1:].strip()
            f2.write("\n" + template + " 1 ")
        else:
            f2.write(aline.strip())
    f.close()
    f2.close()
    # Initialize scoring function
    initstr = frag_files+frag_sizes+fasta+psipred + "-in:file:template_pdb "
    for template in pdbfiles:
        initstr = initstr + template + " "
    initstr = initstr + "-in:file:alignment align.aln -cm:aln_format general -overwrite"
    if (not(os.path.exists("temp"))):
        os.makedirs("temp")
    initializeRosetta(addOn=initstr)
    scorefxn = ScoreFunction()
    try:
        scorefxn.add_weights_from_file(weightsfile)
    except:
        raise Exception("ERROR: The scoring function weights could not be initialized!")
    f = open("threadoutputtemp", "w")
    # Perform protein threading
    LRT = LoopRelaxThreadingMover()
    LRT.setup()
    #os.chdir("temp")
    JobDistributor.get_instance().go(LRT)
    #os.chdir("..")
    # Get the PDB
    # If the only PDB file in "temp" is the abinitio structure, then we can use glob to get it without knowing
    # the exact filename
    pdbfiles = glob.glob("temp/*.pdb")
    outputpdb = "thread_T.pdb"
    pose = pose_from_pdb(pdbfiles[0])
    scorefxn(pose)
    pose.dump_pdb(outputpdb)
    for pdbfile in pdbfiles:
        os.remove(pdbfile)
    # Now write the output information for the main GUI
    f.write("OUTPUT\t" + outputpdb + "\n")
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
    # So the main GUI doesn't attempt to read the file before the daemon finishes writing its contents
    os.rename("threadoutputtemp", "threadoutput")

def doCustom(jobtype):
    # Runs a custom protocol
    # Let's get the module first (it should be at modules/<jobtype>/daemon/__init__.py
    try:
        sys.path.append(os.getcwd() + "/modules/" + jobtype)
        module = __import__("job")
        sys.path.pop()
    except:
        raise Exception("Could not import modules/" + jobtype + "/job!\nDoes modules/" + jobtype + "/job/__init__.py exist?")
    module.runJob(jobtype)

def writeError(msg):
    # Open a file and write out the error message so the main GUI can tell the user what happened
    # The main GUI needs to check to see if an errreport gets generated and recover from the error
    f = open("errreport", "w")
    f.write(msg + "\n\n")
    f.write(traceback.format_exc() + "\n")
    f.close()

def daemonLoop():
    scriptdir = os.getcwd()
    goToSandbox()
    sys.path.append(os.getcwd() + "/modules")
    stillrunning = True
    while (stillrunning):
        if (os.path.isfile("minimizeinput")):
            print "Daemon starting minimization job..."
            try:
                doMinimization()
                print "Daemon completed minimization job"
            except Exception as e:
                print "The daemon crashed while performing the minimization job!"
                writeError(e.message)
            os.remove("minimizeinput")
        elif (os.path.isfile("designinput")):
            print "Daemon starting fixbb design job..."
            try:
                doFixbb()
                print "Daemon completed fixbb design job"
            except Exception as e:
                print "The daemon crashed while performing the fixbb job!"
                writeError(e.message)
            os.remove("designinput")
        elif (os.path.isfile("scaninput")):
            print "Daemon starting point mutation scanning job..."
            try:
                with ScanCapturing() as output:
                    doPMutScan()
                print "Daemon completed point mutation scanning job"
            except Exception as e:
                print "The daemon crashed while performing the scanning job!"
                writeError(e.message)
            os.remove("scaninput")
        elif (os.path.isfile("relaxinput")):
            print "Daemon starting relaxing job..."
            try:
                doRelax()
                print "Daemon completed relaxing job"
            except Exception as e:
                print "The daemon crashed while performing the relaxing job!"
                writeError(e.message)
            os.remove("relaxinput")
        elif (os.path.isfile("backrubinput")):
            print "Daemon starting backrubbing job..."
            try:
                doBackrub()
                print "Daemon completed backrubbing job"
            except Exception as e:
                print "The daemon crashed while performing the backrubbing job!"
                writeError(e.message)
            os.remove("backrubinput")
        elif (os.path.isfile("scoreinput")):
            print "Daemon starting scoring job..."
            try:
                doScore()
                print "Daemon completed scoring job"
            except Exception as e:
                print "The daemon crashed while performing the scoring job!"
                writeError(e.message)
            os.remove("scoreinput")
        elif (os.path.isfile("rotamerinput")):
            print "Daemon starting rotamer searching job..."
            try:
                doRotamerSearch()
                print "Daemon completed rotamer searching job"
            except Exception as e:
                print "The daemon crashed while performing the rotamer search job!"
                writeError(e.message)
            os.remove("rotamerinput")
        elif (os.path.isfile("coarsekicinput")):
            print "Daemon starting coarse KIC loop modeling job..."
            try:
                # This function call has the potential to run indefinitely if it cannot find
                # a way to bridge the two endpoint residues (in a de novo loop model)
                # The main GUI is timing the daemon's response though and will kill it and display
                # an error if it doesn't finish before a timeout (usually 3 min)
                doKIC("Coarse")
                f = open("coarsekicinput", "r")
                for aline in f:
                    if (aline[0:7] == "PERTURB"):
                        perturbType = aline.split("\t")[1].strip()
                    elif (aline[0:7] == "REMODEL"):
                        loopType = aline[7:].strip()
                    elif (aline[0:7] == "NSTRUCT"):
                        nstruct = int(aline.split("\t")[1].strip())
                f.close()
                if (perturbType == "Perturb Only, Centroid"): # or loopType == "REFINE"):
                    os.remove("coarsekicinput")
                else:
                    os.rename("coarsekicinput", "finekicinputtemp")
                    for decoy in range(nstruct-1, -1, -1):
                        os.rename("torepack_" + str(decoy) + ".pdb", "repackmetemp_" + str(decoy) + ".pdb") # So the GUI sees it
                print "Daemon completed coarse KIC loop modeling job"
            except Exception as e:
                print "The daemon crashed while performing the coarse KIC loop modeling job!"
                os.remove("coarsekicinput")
                writeError(e.message)
        elif (os.path.isfile("INDELinput")):
            print "Daemon starting INDEL loop modeling job..."
            try:
                doINDEL(scriptdir)
            except Exception as e:
                print "The daemon crashed while performing the INDEL loop modeling job!"
                writeError(e.message)
                raise(e)
            try:
                os.remove("INDELinput")
            except Exception as e:
                print "The daemon crashed while performing the INDEL loop modeling job!"
                writeError(e.message)

        elif (os.path.isfile("repackme_0.pdb")):
            print "Daemon starting rotamer repacking job..."
            try:
                f = open("finekicinputtemp", "r")
                for aline in f:
                    if (aline[0:7] == "PERTURB"):
                        perturbType = aline.split("\t")[1].strip()
                    elif (aline[0:7] == "NSTRUCT"):
                        nstruct = int(aline.split("\t")[1].strip())
                f.close()
                if (perturbType == "Perturb Only, Fullatom"):
                    for decoy in range(0, nstruct):
                        doRepack("finekicinputtemp", pdbfile="repackme_" + str(decoy) + ".pdb", lastStage=True)
                    os.remove("finekicinputtemp")
                else:
                    for decoy in range(nstruct-1, -1, -1):
                        doRepack("finekicinputtemp", pdbfile="repackme_" + str(decoy) + ".pdb")
                        os.rename("repacked_" + str(decoy) + "temp.pdb", "repacked_" + str(decoy) + ".pdb") # So the GUI sees it
                print "Daemon completed rotamer repacking job"
            except Exception as e:
                print "The daemon crashed while performing the rotamer repacking job!"
                writeError(e.message)
        elif (os.path.isfile("finekicinput")):
            print "Daemon starting fine KIC loop modeling job..."
            try:
                with ScanCapturing() as output:
                    doKIC("Fine")
                print "Daemon completed fine KIC loop modeling job"
            except Exception as e:
                print "The daemon crashed while performing the rotamer search job!"
                writeError(e.message)
            os.remove("finekicinput")
        elif (os.path.isfile("coarsedockinput")):
            print "Daemon starting coarse docking job..."
            try:
                doDock("Coarse")
                f = open("coarsedockinput", "r")
                for aline in f:
                    if (aline[0:13] == "REFINEDDECOYS"):
                        nrefined = int(aline.split("\t")[1].strip())
                        break
                f.close()
                os.rename("coarsedockinput", "finedockinputtemp")
                for i in range(0, nrefined):
                    os.rename("torepackdock_" + str(i) + ".pdb", "repackmedocktemp_" + str(i) + ".pdb") # So the GUI sees it
                print "Daemon completed coarse docking job"
            except Exception as e:
                print "The daemon crashed while performing the coarse docking job!"
                os.remove("coarsedockinput")
                writeError(e.message)
                # The daemon needs to be killed by the main GUI because the coarse KIC thread is
                # still running and will do so indefinitely
        elif (os.path.isfile("repackmedock_0.pdb")):
            print "Daemon starting rotamer repacking job..."
            try:
                f = open("finedockinputtemp", "r")
                for aline in f:
                    if (aline[0:13] == "REFINEDDECOYS"):
                        nrefined = int(aline.split("\t")[1].strip())
                f.close()
                for i in range(0, nrefined):
                    doRepack("finedockinputtemp", "repackmedock_" + str(i) + ".pdb")
                print "Daemon completed rotamer repacking job"
                for i in range(0, nrefined):
                    os.rename("repackeddock_" + str(i) + "temp.pdb", "repackeddock_" + str(i) + ".pdb") # So the GUI sees it
            except Exception as e:
                print "The daemon crashed while performing the rotamer repacking job!"
                writeError(e.message)
        elif (os.path.isfile("finedockinput")):
            print "Daemon starting refined docking job..."
            try:
                doDock("Fine")
                print "Daemon completed refined docking job"
            except Exception as e:
                print "The daemon crashed while performing the rotamer search job!"
                writeError(e.message)
            os.remove("finedockinput")
        elif (os.path.isfile("threadinput")):
            print "Daemon starting threading job..."
            try:
                doThread(scriptdir)
                print "Daemon completed threading job"
            except Exception as e:
                print "The daemon crashed while performing the threading job!"
                writeError(e.message)
            os.remove("threadinput")
        else:
            # CUSTOM MODULE
            for inputfile in glob.glob("*input"):
                jobtype = inputfile.split("input")[0]
                print "Daemon starting " + jobtype + " job..."
                try:
                    doCustom(jobtype)
                    print "Daemon completed " + jobtype + " job"
                except Exception as e:
                    print "The daemon crashed while performing the " + jobtype + " job!"
                    writeError(e.message)
                os.remove(inputfile)
        time.sleep(1)
        # See if the main GUI is still running and have the daemon terminate if the main GUI was closed
        # Since this is a separate process, the user can exit out of the main GUI but this daemon
        # will still be running until they close the Python window, which they may forget to do
        # So the daemon should check to see if the main GUI is still running and if not it should
        # exit this loop
        stillrunning = False
        count = 0
        for proc in psutil.process_iter():
            try:
                if (platform.system() == "Windows"):
                    if (len(proc.cmdline()) >= 2 and proc.cmdline()[0].find("python") >= 0 and proc.cmdline()[1].find("InteractiveROSETTA.py") >= 0):
                        stillrunning = True
                        break
                else:
                    # On Unix systems you just have to make sure two instances of python are running
                    # because there isn't any forking information in the daemon's instance of python
                    if (len(proc.cmdline()) >= 2 and proc.cmdline()[0].find("python") >= 0 and proc.cmdline()[1].find("InteractiveROSETTA") >= 0):
                        count = count + 1
            except:
                # In Windows it will crash if you try to read process information for the Administrator
                # Doesn't matter though since InteractiveROSETTA is run by a non-Administrator
                # But we need to catch these errors since we don't know which processes are admin ones
                pass
        if (platform.system() != "Windows" and count == 2):
            stillrunning = True

if (__name__ == "__main__"):
    daemonLoop()

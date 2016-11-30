import sys
import math
import glob
from mpi4py import MPI
try:
    # Try to import Rosetta
    from rosetta import *
    import rosetta.protocols.backrub
except:
    # If it failed, then try to find Rosetta
    # If this already happened once already, then we should have saved the Rosetta path, so let's try to import from there
    print "Rosetta could not be imported.  Attempting to locate the PyRosetta install.  Please be patient..."
    if ("/" in sys.argv[0]):
        scriptdir = sys.argv[0][0:sys.argv[0].rfind("/")]
    else:
        scriptdir = "."
    cfgfile = scriptdir + "/rosetta.cfg"
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
            # Extra imports for KIC
            from rosetta.protocols.loops.loop_mover.perturb import *
            from rosetta.protocols.loops.loop_mover.refine import *
            os.chdir(olddir)
            print "Found Rosetta at " + rosettapath.strip() + "!"
            print "Rosetta Database: " + rosettadb.strip()
    except:
        print "Where, oh where is PyRosetta?"
        exit()

noderank = MPI.COMM_WORLD.Get_rank()
nodesize = MPI.COMM_WORLD.Get_size()
minimize = False
if (len(sys.argv) != 4):
    print "Usage: python backrub.py pdbfile ntemplates backrubRMSD"
    exit()
individual = True
pdbfile = sys.argv[1].strip()
ntemplates = int(sys.argv[2])
backrubRMSD = float(sys.argv[3])
init(extra_options="-ignore_unrecognized_res -mute all")
scorefxn = create_score_function("talaris2013")
# Get the params files
paramsfiles = glob.glob("*.fa.params")
paramsstr = " -extra_res_fa "
if (len(paramsfiles) == 0):
    paramsstr = ""
else:
    for filename in paramsfiles:
        paramsstr += filename + " "
# Load the given PDB file and then use a PyRosetta Backrub mover to generate the template
MPI.COMM_WORLD.Barrier()
if (noderank == 0):
    print "backrub: Starting ensemble generation"
for teamrank in range(noderank, ntemplates, nodesize):
    if (ntemplates == 1):
        myRMSD = float(backrubRMSD)
    else:
        myRMSD = float(teamrank) * (float(backrubRMSD) / float(ntemplates - 1))
    RMSDseed = "%7i" % (myRMSD * 1000) # Truncate the decimal
    RMSDseed = RMSDseed.strip()
    init(extra_options="-ignore_unrecognized_res " + paramsstr + " -constant_seed -seed_offset 0 -use_bicubic_interpolation -mute all")
    pose = pose_from_pdb(pdbfile)
    BM = protocols.backrub.BackrubMover()
    pivot_res = utility.vector1_Size()
    # Set all canonical AAs to pivot residues
    three = "ALACYSASPGLUPHEGLYHISILELYSLEUMETASNPROGLNARGSERTHRVALTRPTYR"
    for i in range(0, pose.n_residue()):
        if (three.find(pose.residue(i+1).name()) >= 0):
            pivot_res.append(i+1)
    BM.set_pivot_residues(pivot_res)
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
        init(extra_options="-ignore_unrecognized_res " + paramsstr + " -constant_seed -seed_offset 0 -mute all")
        minmover = MinMover(mm, scorefxn_min, 'lbfgs_armijo_nonmonotone', tol, True)
        init(extra_options="-ignore_unrecognized_res " + paramsstr + " -constant_seed -seed_offset " + RMSDseed + " -use_bicubic_interpolation -mute all")
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
            for i in range(0, myLoop):
                BM.apply(pose)
                # If the move is rejected, pose is automatically reverted by this function
                Metropolis.boltzmann(pose)
            minmover = MinMover(mm, scorefxn_min, 'lbfgs_armijo_nonmonotone', tol, True)
            minmover.apply(pose)
            attempts = attempts + 1
            if (attempts >= 100):
                backrubTm = backrubTm * 2.0
                myLoop = myLoop + 10
                Metropolis = rosetta.MonteCarlo(pose, scorefxn, backrubTm)
                attempts = 0
                lastRMSD = -1
            if (backrubTm >= 1000000.0):
                print str(teamrank) + " - Too hot", tol
                print teamrank, core.scoring.bb_rmsd(pose, orig_pose)
                break
            if ( backrubTm < 0.0001):
                print str(teamrank) + " - Too cold", tol
                pose = Pose(prev_pose)
                print teamrank, core.scoring.bb_rmsd(pose, orig_pose)
                break
            calpha_superimpose_pose(pose, orig_pose)
            if (lastRMSD < 0):
                lastRMSD = core.scoring.bb_rmsd(pose, orig_pose)
            else:
                currRMSD = core.scoring.bb_rmsd(pose, orig_pose)
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
    pose.dump_pdb(pdbfile.split(".pdb")[0] + "_" + str(teamrank+1) + ".pdb")
MPI.COMM_WORLD.Barrier()
MPI.Finalize()

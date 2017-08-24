import numpy
import os
import sys
import platform
try:
    # Try to import Rosetta
    from pyrosetta import *
    from pyrosetta.rosetta import *
    # Extra imports for KIC
    from pyrosetta.rosetta.protocols.loops.loop_mover.perturb import *
    from pyrosetta.rosetta.protocols.loops.loop_mover.refine import *
    # Extra import for docking
    import pyrosetta.rosetta.protocols.rigid as rigid_moves
    # Extra imports for threading
    import pyrosetta.rosetta.protocols.evaluation
    from pyrosetta.rosetta.protocols.comparative_modeling import *
    from pyrosetta.rosetta.protocols.jd2 import *
    from pyrosetta.rosetta.core.scoring.constraints import *
except:
    #TODO rewrite to prompt user to properly install PyRosetta 4
    # If it failed, then try to find Rosetta
    # If this already happened once already, then we should have saved the Rosetta path, so let's try to import from there
    print "Rosetta could not be imported.  Attempting to locate the PyRosetta install.  Please be patient..."
    if (platform.system() == "Windows"):
        cfgfile = os.path.expanduser("~") + "/InteractiveROSETTA/seqwindow.cfg"
    else:
        cfgfile = os.path.expanduser("~") + "/.InteractiveROSETTA/seqwindow.cfg"
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
            os.environ["PYROSETTA_DATABASE"] = rosettadb
            # Try to import Rosetta
            from rosetta import *
            # Extra imports for KIC
            from rosetta.protocols.loops.loop_mover.perturb import *
            from rosetta.protocols.loops.loop_mover.refine import *
            # Extra import for docking
            import rosetta.protocols.rigid as rigid_moves
            # Extra imports for threading
            import rosetta.protocols.evaluation
            from rosetta.protocols.comparative_modeling import *
            from rosetta.protocols.jd2 import *
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
                from rosetta.protocols.comparative_modeling import *
                from rosetta.protocols.jd2 import *
                # Now let's save these paths so the next time this gets started we don't have to traverse the filesystem again
                data = []
                f = open(cfgfile, "r")
                for aline in f:
                    if (not("[ROSETTAPATH]" in aline) and not("[ROSETTADB]") in aline):
                        data.append(aline.strip())
                f.close()
                f = open(cfgfile, "w")
                for aline in data:
                    f.write(aline + "\n")
                f.write("[ROSETTAPATH]\t" + rosettapath.strip() + "\n")
                f.write("[ROSETTADB]\t" + rosettadb.strip() + "\n")
                f.close()
                print "Found Rosetta at " + rosettapath.strip() + "!"
                print "Rosetta Database: " + rosettadb.strip()
            except:
                print "PyRosetta cannot be found on your system!"
                print "Until you install PyRosetta, you may only use InteractiveROSETTA to visualize structures in PyMOL"
                exit()

def translateToOrigin(pose, root_res):
    transv = numpy.array([pose.residue(root_res).atom(1).xyz()[0], pose.residue(root_res).atom(1).xyz()[1], pose.residue(root_res).atom(1).xyz()[2]])
    for r in range(1, pose.n_residue()+1):
        for a in range(1, len(pose.residue(r).atoms())+1):
            newx = pose.residue(r).atom(a).xyz()[0] - transv[0]
            newy = pose.residue(r).atom(a).xyz()[1] - transv[1]
            newz = pose.residue(r).atom(a).xyz()[2] - transv[2]
            pose.residue(r).atom(a).xyz(numeric.xyzVector_double(newx, newy, newz))

def getUnitVector(v):
    norm = numpy.linalg.norm(v)
    return v / norm

def getRotationMatrix(move_v, static_v):
    v = numpy.cross(move_v, static_v)
    s = numpy.linalg.norm(v)
    c = numpy.dot(move_v, static_v)
    vx = numpy.zeros(9).reshape(3, 3)
    vx[0][1] = -1.0 * v[2]
    vx[0][2] = v[1]
    vx[1][0] = v[2]
    vx[1][2] = -1.0 * v[0]
    vx[2][0] = -1.0 * v[1]
    vx[2][1] = v[0]
    I = numpy.zeros(9).reshape(3, 3)
    I[0][0] = 1.0
    I[1][1] = 1.0
    I[2][2] = 1.0

    if (s**2 == 0):
        return "No Rotation" # Because the cross product of move_v and static_v is 0, in which case no rotation is necessary
    R = I + vx + (vx.dot(vx))*((1.0-c)/(s**2))
    return R

def rotatePose(pose, R):
    for r in range(1, pose.n_residue()+1):
        for a in range(1, len(pose.residue(r).atoms())+1):
            v = numpy.array([pose.residue(r).atom(a).xyz()[0], pose.residue(r).atom(a).xyz()[1], pose.residue(r).atom(a).xyz()[2]])
            newv = R.dot(v)
            pose.residue(r).atom(a).xyz(numeric.xyzVector_double(newv[0], newv[1], newv[2]))
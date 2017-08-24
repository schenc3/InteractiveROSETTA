import subprocess
import sys
import os

# This script is an attempt to use PyMOL to do the RosettaAntibody superimposition
# You have to get permission to use ProFit and the PyRosetta script they provided does not seem to work
if (len(sys.argv) != 2):
        print "Usage: python superimpose_pymol.py profit.in"
        exit()
profitin = sys.argv[1].strip()
prefix = profitin[0:profitin.rfind("/")]

# Read the ProFit input data
atoms = ["ca"]
zones = []
f = open(profitin, "r")
for aline in f:
        if ("reference" in aline.lower()):
                refpdb = aline.split()[1].strip()
        elif ("mobile" in aline.lower()):
                mobpdb = aline.split()[1].strip()
        elif ("write" in aline.lower()):
                outpdb = aline.split()[1].strip()
        elif ("atoms" in aline.lower()):
                atoms = aline.split()[1].split(",")
        elif ("zone" in aline.lower()):
                if (":" in aline):
                        zones.append(aline.split()[1].split(":"))
                else:
                        zones.append([aline.split()[1]])
f.close()
f = open(prefix + "/pymol_superimp.py", "w")
f.write("import pymol\n")
f.write("pymol.cmd.load(" + refpdb + ", \"ref\")\n")
f.write("pymol.cmd.load(" + mobpdb + ", \"mob\")\n")
refstr = "\""
mobstr = "\""
# Convert the zone data into PyMOL selection strings
for zone in zones:
        if (len(zone) == 1):
                bounds = zone[0].split("-")
                if (not(bounds[0][0] in "0123456789")):
                        chain = "chain " + bounds[0][0]
                        bounds[0] = bounds[0][1:]
                        bounds[1] = bounds[1][1:]
                else:
                        chain = ""
                if (len(chain) > 0):
                        refstr = refstr + "(model ref and " + chain + " and resi " + bounds[0] + "-" + bounds[1]
                        mobstr = mobstr + "(model mob and " + chain + " and resi " + bounds[0] + "-" + bounds[1]
                else:
                        refstr = refstr + "(model ref and resi " + bounds[0] + "-" + bounds[1]
                        mobstr = mobstr + "(model mob and resi " + bounds[0] + "-" + bounds[1]
                if (atoms[0] != "*"):
                        for atom in atoms:
                                refstr = refstr + " and name " + atom
                                mobstr = mobstr + " and name " + atom
                refstr = refstr + ") or "
                mobstr = mobstr + ") or "
        else:
                pass
refstr = refstr[0:len(refstr)-4] + "\"" # Trim off the trailing "or"
mobstr = mobstr[0:len(mobstr)-4] + "\"" # Trim off the trailing "or"
f.write("pymol.cmd.align(" + mobstr + ", " + refstr + ")\n")
f.write("pymol.cmd.save(" + outpdb + ", \"mob\")\n")
f.close()
p = subprocess.Popen("pymol -qc " + prefix + "/pymol_superimp.py", shell=True, stdout=subprocess.PIPE)
p.communicate()
os.remove(prefix + "/pymol_superimp.py")
